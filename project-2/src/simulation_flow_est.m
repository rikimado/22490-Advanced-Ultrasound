clc, clear, close all

field_init(-1);

f0 = 5e6;                        % Transducer center frequency [Hz]
M = 5;                           % Number of cycles in emitted pulse
fs = 100e6;                      % Sampling frequency [Hz]
c = 1540;                        % Speed of sound [m/s]
element_height = 5e-3;           % Height of element [m]
element_width = 0.2e-3;          % Width of element
kerf = 0.03e-3;                  % Kerf [m]
pitch = element_width+kerf;      % Pitch of transducer
elevation_focus = [0, 0, 30e-3]; % Fixed focal point [m]
N_elements = 192;

active_elements_tx = 64;

% Set the sampling frequency
set_sampling(fs);

% Set the impulse response and excitation of the emit aperture
t_excitation = 0:1/fs:2/f0;
excitation = sin(2*pi*f0*t_excitation);
t_impulse = 0:1/fs:M/f0;
impulse_response = sin(2*pi*f0*t_impulse);

plotTimeDomain(t_excitation, excitation,t_impulse, impulse_response,fs, f0);
impulse_response = impulse_response.*hann(numel(impulse_response))';
%% Create transducer array

% Build aperture for emission
emit_aperture = xdc_linear_array(N_elements, element_width, element_height, kerf, 1, 1, elevation_focus);

% Set how the transmit aperture responds to an excitation
xdc_impulse(emit_aperture, impulse_response);

% Build the receive aperture
receive_aperture = xdc_linear_array(N_elements, element_width, element_height, kerf, 1, 1, elevation_focus);

% Set how the receive aperture responds to incoming signals
xdc_impulse(receive_aperture, impulse_response);

% Set the system excitation
xdc_excitation(emit_aperture, excitation);

% Set a hanning apodization on the active elements during transmission
apo = hann(active_elements_tx)';

%% Scanning parameters
Nframes = 2; % Number of frame to simulate

vs_spacing = pitch;
vs_x = (-64:8:64)*vs_spacing;             % VS lateral position based on the spacing
no_lines = numel(vs_x);                   % num transmission = num vs
vs_z = -active_elements_tx*vs_spacing;    % VS depth based on F# = -1;
VS = [vs_x.', repmat(vs_z, no_lines, 1)]; % Add depth vs_z to create the virtual sources coords

%% Create a digital phantom flow
vessel_length = 15/1000;    % [m] length of the vessel
vessel_height = 8/1000;     % [m] diameter (height) of the vessel
peak_velocity = 0.3;        % [m/s]
PRF = 5000;                 % Pulse repetition frequency [Hz]
Tprf = 1/PRF;               % [sec]
num_scatterers = 15000;
depth_offset = 30/1000;
velocity_profile = 'plug';  % plug or parabolic
theta = 30/180*pi;     % angle transducer to flow

[x,y,z,amp,velocity] = make_flow_phantom(vessel_length, vessel_height, peak_velocity, num_scatterers, velocity_profile);

%% Simulate RF acquisition of the phantom flow
clear rf_data times

% preallocation
times = zeros(Nframes,no_lines);

for i = 1:no_lines

    currentVS_idx = mod(i-1, no_lines) + 1;

    disp(['VS ', num2str(currentVS_idx), '/', num2str(no_lines),': ', num2str(vs_x(currentVS_idx))]);

    % Center the transmit aperture on the current vs position
    xdc_center_focus(emit_aperture, [vs_x(currentVS_idx) 0 0]);

    % Place the transmit focal point based on F#
    xdc_focus(emit_aperture, 0, [vs_x(currentVS_idx), 0, vs_z]);

    % Center the receive aperture at the center of the transducer?
    % and set the time for which this point becomes valid (0s, meaning since the beginning)
    xdc_center_focus(receive_aperture, [0 0 0]);

    % Set the position of the receive focus and the time in which becomes valid
    xdc_focus(receive_aperture, 0, [0, 0, 5]); % 5m ensure no delay is applied on the receivers

    % Clear the current apo window and activate the next 64 elements
    current_apo_vector = zeros(N_elements,1);
    current_apo_vector( (1:active_elements_tx) + (currentVS_idx-1)*8 ) = apo;

    % Apply the apodization specifing the time in which it becomes valid
    xdc_apodization(emit_aperture, 0, current_apo_vector.');

    for j = 1:Nframes

        disp([num2str(j),'/', num2str(Nframes),' frames being processed...']);

        % Rotate the coordinates (x, z) by θ in the X-Z plane.
        % The angle θ is used such that the rotation occurs from the X toward the Z-axis
        % according to the provided code!

        % Compute the new x-coordinate (xnew) and z-coordinate (znew) after rotation
        % by applying the 2D rotation matrix to the original (x, z) coordinates.
        % The depth_offset is added to the z-coordinate to shift the vessel along the Z-axis.
        xnew = x*cos(theta) + z*sin(theta);
        znew = z*cos(theta) - x*sin(theta) + depth_offset;
        scatterers = [xnew; y; znew]';

        % Calculate the receive response
        [v, t1] = calc_scat_multi(emit_aperture, receive_aperture, scatterers, amp');

        % Amount of padding
        padDelay = round(t1 * fs);
        
        % Concatenate the zero-padding matrix at the beginning of the signal v
        v_padded = cat(1, zeros(padDelay, N_elements), v);

        % Store padded signal v_padded and times:
        % - 3rd dimension (j): frame index (2 or 16)
        % - 4th dimension (i): transmission index (17)
        rf_data(1:size(v_padded,1),:, j, i) = v_padded;
        times(j, i) = t1;

        % Move the scatterers and align them to lie within the correct
        % range, to simulate their movement and acquire
        x = x + velocity * Tprf;
        outside_range = (x > vessel_length/2);
        x = x-vessel_length*outside_range;

    end
end

% save the file to avoid repeating this lenghty process
save([velocity_profile,num2str(Nframes),'_RFdata.mat'], 'rf_data','-v7.3');

%% SA beamforming settings and rf_data preprocessing
clear rf_data
filename = sprintf('%s%d_RFdata.mat', velocity_profile, Nframes);
load(filename);

lateral_range = 20e-3;  % ±10 mm
depth_min = 20e-3;      % start depth
depth_max = 40e-3;      % end depth

pixel_size = 0.1e-3;

% image lateral dimension
x_grid = -lateral_range/2:pixel_size:lateral_range/2;

% image depth dimension
z_grid = depth_min:pixel_size:depth_max;             
[X_grid, Z_grid] = meshgrid(x_grid, z_grid);

pixelMap = cat(3, X_grid, Z_grid);
apod = hamming(N_elements);

% Define the elements position: 192 elements separated by kerf (pitched array)
arrayPos = (-(N_elements-1)/2:(N_elements-1)/2)*pitch;

% attach 192 zeros to add the z-axis coords
arrayPos = [arrayPos.', zeros(N_elements,1)];

% Precompute constant variable for all pixels within ImageFormation()
si = size(pixelMap);
t = (1:length(rf_data))/fs;

% X is the array positions, T the number of samples
[X, T] = meshgrid(1:length(arrayPos), t);

% Precompute Xq for interpolation
Xq = repmat(1:length(arrayPos),[si(2), 1]);

% Initialize LRIs matrix as: [depth range x lateral range x transmissions x frames]
LRIs = zeros(si(1), si(2), no_lines, Nframes);

% Apply Hilbert transform
rf_data_hil = hilbert(rf_data*1e23);

%% SA rf_data beamforming process

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
save([velocity_profile, num2str(Nframes),'_LRIs.mat'], 'LRIs', '-v7.3');

%% Load beaformed LRIs
clear LRIs

% Load the LRIs
filename = sprintf('%s%d_LRIs.mat', velocity_profile, Nframes);
load(filename);

% Define vessel bounds
vessel_radius = vessel_height / 2; 
z_center = depth_offset+0.7e-3; % Vessel center at z = 30.7 mm

% Define vessel boundaries location
z_lower = z_center - vessel_radius;
z_upper = z_center + vessel_radius;

% Logical mask for depth points inside vessel
inside_vessel = (z_grid >= z_lower) & (z_grid <= z_upper);

%% HRI formation

% Create HRIs matrix
HRIs = squeeze(sum(LRIs,3));
HRIs_rot = imrotate(HRIs, -30, 'bicubic', 'crop');

max_ref = max(sum(abs(LRIs), 3), [], 'all');
HRIs_dB = 20 * log10(abs(HRIs_rot)/max_ref);

% HRI visualization
figure;
while true
    for i = 1:Nframes
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
        x_lines = 0;
        %x_lines = -5:0.5:5; %mm
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
roi_range = 0;
%roi_range = (-5:0.5:5)*1e-3;

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

    % velocity values along the cross section
    velocity_flow(inside_vessel, xi) = velocity_beam/sin(theta);

    % FILTER (OPTIONAL)
    %velocity_flow(inside_vessel) = smoothdata(velocity_flow(inside_vessel), 'gaussian', 5); % window size = 5

    % Generate expected velocity profile inside the vessel
    if strcmp(velocity_profile, 'plug')

        % Compute the mean value along the cross section
        mean_velocity_flow = mean(velocity_flow(inside_vessel),1);

        % Generate the expected profile (0.3 m/s)
        expected_profile(inside_vessel,xi) = peak_velocity;

    elseif strcmp(velocity_profile, 'parabolic')
        z_rel = z_grid(inside_vessel) - z_center-0.5e-3;
        expected_profile(inside_vessel,xi) = peak_velocity * (1 - (z_rel /(vessel_radius-1e-3)).^2);
    end
end

% Average profiles across lateral points (depth x n° of lateral points)
avg_estimated_profile = mean(velocity_flow, 2);
avg_expected_profile = mean(expected_profile, 2);
%%
figure;

% ESTIMATED PROFILE
h1 = plot(avg_estimated_profile, z_grid* 1e3, 'r', 'LineWidth', 1.2);
hold on

% ESIMATED MEAN VELOCITY (ONLY FOR PLUG)
h2 = plot(mean_velocity_flow*ones(length(z_grid),1), z_grid* 1e3,'r--', 'LineWidth', 1.2);

% EXPECTED PROFILES
h3 = plot(avg_expected_profile, z_grid* 1e3, 'b', 'LineWidth', 1);

% EXPECTED MEAN(plug)/PEAK(parabolic)
h4 = plot(peak_velocity*ones(length(z_grid),1), z_grid* 1e3, 'b--', 'LineWidth', 1);

set(gca, 'YDir', 'reverse');
xlabel('Velocity [m/s]');
ylabel('Depth [mm]');

title(['Estimated vs Expected profiles (x = 0 mm): ' velocity_profile]);  % or ±5 mm
legend([h1, h2, h3, h4], ...
    {'Estimated velocity profile', ...
    'Estimated mean velocity',...
    'Expected velocity profile', ...
    'Expected mean velocity'},...
    'Location', 'northwest');

grid on;

%% Functions

% Delay-and-Sum Beamforming (DAS process)
function LRI = imageFormation(rf_data, pixelMap, VS_focus, arrayPos, c, pitch, X,T, Xq, x_grid)

LRI = zeros(size(pixelMap,1),size(pixelMap,2));

for i = 1:size(pixelMap,1) % loop over depth values (rows of pixelMap)

    % Initialize a NaN matrix for storing computed delays for the pixels
    % in row i. Each pixel can cover at most 192 element (transducer physical
    % boundary), thus the matrix is [lateral length x transducer elements]
    delay_line = nan(size(pixelMap,2), length(arrayPos));

    % The apodization has to work on the aperture created for each pixel.
    % For the same reason as before, can be at most the same size of the
    % delay matrix.
    apod_line = zeros(size(pixelMap,2), length(arrayPos));

    % loop over the pixels of the selected row along the lateral direction
    for j = 1:size(pixelMap,2)

        % extract pixel at the j-th position
        pixel = [pixelMap(i,j,1), pixelMap(i,j,2)];

        % Compute delays for that pixel
        [pixel_delays, pixel_apod] = delayCal(pixel, [i,j], VS_focus, arrayPos, c, x_grid, pitch);

        % Stores delays and apodization window size for that pixel
        delay_line(j,:) = pixel_delays';
        apod_line(j,:) = pixel_apod';
    end

    % Interpolate with interp2: Columns = elements, rows = time
    value = interp2(X, T, rf_data, Xq, delay_line, 'cubic',0);
    value_weighted = value .* apod_line;

    LRI(i,:) = sum(value_weighted, 2, 'omitnan');
end

end
