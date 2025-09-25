clc, clear, close all

load('parameters.mat');
%load('vantage_params.mat');

field_init(-1);

f0 = 5.2083e6;                   % Transducer center frequency [Hz]
fs = 20.833e6;
c = opt_struct.patient.c;        % 1487.5 [m/s]
element_width = 2.07e-4;         % Width of element
kerf = 2.3e-5;                   % Kerf [m]
pitch = element_width+kerf;      % Pitch of transducer

N_elements = 192;
active_elements_tx = 64;
tx_order = opt_struct.B_mode.xmt_placement;

N_frames = 16;
N_vs = 17;
folder_name = 'seq_0005';

t1 = sarus_rx.record.start_time; % receive start time
padDelay = round(t1 * fs);       % amount of padding

%% rf_data grouping

% Loop over the 272 transmissions
for frame_idx = 1:N_frames
    for vs_idx = 1:N_vs

        % Get the correct TX index from xmt_placement
        tx_index = tx_order(vs_idx, frame_idx);

        % create the file name
        vs_filename = sprintf('elem_data_em%04d.mat', tx_index);
        full_path = fullfile(folder_name, vs_filename);

        % Load the file [time sample x 192]
        data = load(full_path);

        %Concatenate the zero-padding matrix [padDelay x 192] at the beginning of the signal
        v_padded = cat(1, zeros(padDelay, N_elements),  data.samples);

        % Store samples
        rf_data(:, :, frame_idx, vs_idx) = v_padded;
    end

end

save(['exp_RFdata_',folder_name,'.mat'], 'rf_data','-v7.3');

%% SA beamforming settings and rf_data preprocessing
clear rf_data
filename = sprintf('exp_RFdata_%s.mat',folder_name);
load(filename);

vs_x = (-64:8:64)*pitch;           % VS lateral position based on the spacing (= pitch)
no_lines = N_vs;                   % num transmission = num vs
vs_z = -active_elements_tx*pitch;    % VS depth based on F# = -1;
VS = [vs_x.', repmat(vs_z, no_lines, 1)]; % Add depth vs_z to create the virtual sources coords

% SA Beamforming parameters
lateral_range = 44.16e-3; % ±22.08 mm
depth_min = opt_struct.patient.start_depth; % start depth: 20 mm
depth_max = opt_struct.patient.end_depth;   % end depth: 50 mm

pixel_size = 0.1e-3;
x_grid = -lateral_range/2:pixel_size:lateral_range/2;  % image lateral dimension
z_grid = depth_min:pixel_size:depth_max;               % image depth dimension
[X_grid, Z_grid] = meshgrid(x_grid, z_grid);

pixelMap = cat(3, X_grid, Z_grid);
apod = hamming(N_elements);

% Define the elements position: 192 elements separated by kerf (pitched array)
arrayPos = (-(N_elements-1)/2:(N_elements-1)/2)*pitch;
arrayPos = [arrayPos.', zeros(N_elements,1)]; % attach 192 zeros to add the z-axis coords

% Precompute constant variable for all pixels within ImageFormation()
si = size(pixelMap);
t = (1:length(rf_data))/fs;
[X, T] = meshgrid(1:length(arrayPos), t);   % X is the array positions, T the number of samples
Xq = repmat(1:length(arrayPos),[si(2), 1]); % Precompute X2 for interpolation

% Initialize LRIs matrix as: [depth range x lateral range x transmissions x frames]
LRIs = zeros(si(1), si(2), no_lines, N_frames);

% Apply Hilbert transform
rf_data_hil = hilbert(rf_data);

%% Beamforming process

for frame = 1:size(rf_data,3)

    disp(['Frame ', num2str(frame),':']);

    % Select rf_data for the current frame and remove the frame dimension as we
    % process each frame separately
    rf_data_frame = squeeze(rf_data_hil(:,:,frame,:));

    % Loop over transmissions as each of them corresponds to one different VS
    parfor trans = 1:no_lines

        disp(['Generating LRI ', num2str(trans), '/', num2str(no_lines)]);

        % Get RF data for current transmission
        rf_data_vs = rf_data_frame(:, :, trans);

        % Generate the LRIs for that vs
        LRIs(:, :, trans, frame) = imageFormation(rf_data_vs, pixelMap, VS(trans,:), arrayPos, c, pitch, X,T,Xq, x_grid);
    end

end

% save the file to avoid repeating this lenghty process
save(['exp_LRIs_',folder_name,'.mat'], 'LRIs','-v7.3');

%% Load beaformed LRIs
clear LRIs
filename = sprintf('exp_LRIs_%s.mat',folder_name);
load(filename);

PRF = opt_struct.patient.fprf ; % 5000 Hz
Tprf = 1/PRF;                   % [sec]
theta = 15/180*pi;              % angle lateral to flow

% Define vessel bounds
vessel_radius = opt_struct.patient.radius; %6 mm
vessel_depth = opt_struct.patient.vessel_depth; % 35 mm

% Define vessel boundaries location
z_lower = vessel_depth - vessel_radius+1.5e-3;
z_upper = vessel_depth + vessel_radius-0.5e-3;

% Logical mask for depth points inside vessel
inside_vessel = (z_grid >= z_lower) & (z_grid <= z_upper);

peak_velocity = -0.275;      % -0.4912 for seq_0005 [m/s]

%% HRI formation

% Create HRIs matrix [depth x lateral x frames] by summing along the VS dimension (3rd)
HRIs = squeeze(sum(LRIs,3));           % 16 frames
%HRIs = squeeze(sum(LRIs(:,:,:,1:2), 3));  % 2 frames

% Stationary echo cancellation
%HRIs = HRIs - mean(HRIs,3);

HRIs_rot = imrotate(HRIs, -15, 'bicubic', 'crop');
max_ref = max(sum(abs(LRIs), 3), [], 'all');
HRIs_dB = 20 * log10(abs(HRIs_rot)/max_ref);

% HRI frames visualization
figure;
while true
    for i = 1:N_frames
        imagesc(x_grid*1e3, z_grid*1e3 ,HRIs_dB(:,:,i));
        colormap gray;
        clim([-60, 0]);
        cb = colorbar;
        xlabel('Lateral position [mm]');
        ylabel('Depth [mm]');
        ylabel(cb, 'dB');
        axis ij equal tight
        title(sprintf('Frame %d', i));

        % Add vertical lines at roi locations
        hold on;
        x_lines = -5:0.5:5; %mm
        for x = x_lines
            line([x x], [z_lower*1e3 z_upper*1e3], 'Color', 'r', 'LineWidth', 1);
        end
        drawnow;
    end

    % Pause and wait until the figure window is closed
    if ~ishandle(gcf)
        break; % exit the loop when the figure is closed
    end
end

%% Autocorrelation-based velocity estimator

% Define the range of interest in meters
roi_range = (-5:0.5:5)*1e-3;

% Initialize velocity profile arrays
velocity_flow = zeros(length(z_grid), length(roi_range));
expected_profile = zeros(length(z_grid), length(roi_range));

% Loop over roi_range positions
for xi = 1:length(roi_range)

    % Extract lateral position value
    x_val = roi_range(xi);

    % Find the index in the grid of that value
    [~, x_idx] = min(abs(x_grid - x_val));

    % Extract the signal at the current lateral location and only inside
    % the vessel
    signals = squeeze(HRIs_rot(inside_vessel, x_idx, :));

    % Compute lag-one autocorrelation between signal of frames
    Rxx = sum(conj(signals(:,1:end-1)).*signals(:,2:end),2);
    phase_diff = angle(Rxx);

    % velocity values along the beam direction
    velocity_beam = (phase_diff * c) / (4 * pi * f0 * Tprf);

    % NOTE: measured velocity = real velocity along vessel * cos(90-theta)
    % thus, real velocity along vessel = measured velocity/cos(90-theta)

    % velocity values along the vessel direction
    velocity_flow(inside_vessel, xi) = velocity_beam/sin(theta);

    % FILTER (OPTIONAL)
    %velocity_flow(inside_vessel) = smoothdata(velocity_flow(inside_vessel), 'gaussian', 5); % window size = 5

    % Generate expected velocity profile inside the vessel
    z_rel = z_grid(inside_vessel) - vessel_depth-0.5e-3;
    expected_profile(inside_vessel,xi) = peak_velocity * (1 - (z_rel /(vessel_radius-1e-3)).^2);
end

% Average profiles across lateral points (depth x n° of lateral points)
avg_estimated_profile = mean(velocity_flow, 2);
avg_expected_profile = mean(expected_profile, 2);

%%
figure;

% ESTIMATED PROFILE
h1 = plot(avg_estimated_profile, z_grid* 1e3, 'r', 'LineWidth', 1.2);
hold on

% EXPECTED PROFILES
h2 = plot(avg_expected_profile, z_grid* 1e3, 'b', 'LineWidth', 1);

% EXPECTED PEAKS
h3 = plot(peak_velocity*ones(length(z_grid),1), z_grid* 1e3, 'b--', 'LineWidth', 1);

set(gca, 'YDir', 'reverse');
xlabel('Velocity [m/s]');
ylabel('Depth [mm]');
title(['Estimated vs Expected profiles (x = ± 5 mm): ' folder_name], 'Interpreter', 'none');
legend([h1, h2, h3], ...
    {'Estimated profile', ...
    'Expected profile', ...
    'Expected peak'},...
    'Location', 'northwest');

grid on;

%%
Nframes=size(HRIs,3);
time = (0:4) * Tprf;
velocity_trend = all_profiles(156, :);

plot(time, velocity_trend);
    
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend show;
title('Pulsatile Velocity Trend at Vessel Center');
grid on;
