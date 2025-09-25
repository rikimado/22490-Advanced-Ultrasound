clc, clear, close all

field_init(-1);

% -------- ATTENTION -----------------------------

% speed of sound
c = 1540;   % in tissue
%c = 1480;  % in water (only for scanner points)

% sampling rate [Hz]
fs = 100e6;    % for simulation
%fs = 62.5e6;  % for scanner

% ------------------------------------------------

% Transducer parameters
f0 = 7e6;                      % center frequency (Hz)
lambda = c/f0;                 % wavelength [m]
element_height = 5e-3;         % element heigth in y direction
element_width = 0.2e-3;        % element width in x direction
kerf = 0.03e-3;                % gap between elements in x direction
pitch = element_width + kerf;  % element width + kerf
elevation_focus = [0,0,30e-3];
num_elements = 192;

active_elements_tx = 64;

% Set the sampling frequency
set_sampling(fs);

% Set the impulse response and excitation of the emit aperture
t_excitation = 0:1/fs:1/f0;
excitation = sin(2*pi*f0*t_excitation);
t_impulse = 0:1/fs:2/f0;
impulse_response = sin(2*pi*f0*t_impulse);
plotTimeDomain(t_excitation, excitation,t_impulse, impulse_response,fs, f0);
impulse_response = impulse_response.*hann(numel(impulse_response))';

%% Create transducer array

% Build aperture for emission
emit_aperture = xdc_linear_array(num_elements, element_width, element_height, kerf, 1, 1, elevation_focus);

% Set how the transmit aperture responds to an excitation
xdc_impulse(emit_aperture, impulse_response);

% Build the receive aperture
recv_aperture = xdc_linear_array(num_elements, element_width, element_height, kerf, 1, 1, elevation_focus);

% Set how the receive aperture responds to incoming signals
xdc_impulse(recv_aperture, impulse_response);

% Set the excitation
xdc_excitation(emit_aperture, excitation);

%% Scanning parameters
num_transmissions = 129; % this also represents the n° of vs, as at each transmission we change vs

% VS settings:
vs_spacing = pitch;  % spacing between each vs
vs_span = (num_transmissions-1)*vs_spacing;
x = -vs_span/2;      % Position of the first vs
vs_x = linspace(x, -x, num_transmissions);
vs_z = -active_elements_tx * vs_spacing;   % vs depth based on F# = -1 = focal_depth / transmitter aperture

% Time in which the focuses are valid
t0 = 0;

% Define the apodization window to isolate the 64 active elements
apo_vector = hamming(active_elements_tx);

%% Choose and define simulation target

% Simulate points
phantom_positions = [0 0 10e-3; 0 0 20e-3; 0 0 30e-3; 0 0 40e-3; 0 0 50e-3];%
phantom_amplitudes = ones(size(phantom_positions, 1), 1)*1e6;

% Simulate cyst
%[phantom_positions, phantom_amplitudes] = cystFormation(10000);

%% Simulate RF acquisition for each transmission

clear rf_data times
times = zeros(1,num_transmissions); % preallocation

for i=1:num_transmissions

    % Center the transmit aperture on the current vs position
    xdc_center_focus(emit_aperture, [x 0 0]);

    % Place the transmit focal point based on F#
    xdc_focus(emit_aperture, t0, [x ,0, vs_z]);

    % Center the receiver aperture on the current vs position 
    xdc_center_focus(recv_aperture, [x 0 0]);

    % Set the position of the receive focus and the time in which becomes valid
    xdc_focus(recv_aperture, t0, [x, 0, 5]); % 5m ensure no delay is applied on the receivers

    % Clear the current apo window and activate the next 64 elements
    currentApod = zeros(num_elements,1);
    currentApod((1:active_elements_tx) +i- 1) = apo_vector;

    % Apply the apodization specifing the time in which it becomes valid
    xdc_apodization(emit_aperture, t0, currentApod.');

    % NOT REQUIRED AS ALL THE 192 ELEMENTS ARE ACTIVE DURING RECEPTION
    %xdc_apodization(recv_aperture, t0, currentApod.');

    % Calculate the received response
    [v, t1] = calc_scat_multi(emit_aperture, recv_aperture, phantom_positions, phantom_amplitudes);

    % Zero padding matrix
    padDelay = round(t1 * fs);

    % Concatenate the zero-padding matrix at the beginning of the signal v
    v_padded = cat(1, zeros(padDelay, num_elements), v);

    % Store padded signal and times:
    rf_data(1:length(v_padded),:,i) = v_padded;
    times(i) = t1;

    % Shift the vs
    x = x + vs_spacing;
end

%% Run only for SIMULATED CYST
%  load simulation_cyst_RFdata.mat
%  Add white noise
%  SNR_dB = 20;
%  rf_data = rf_data + randn(size(rf_data)) * (max(abs(rf_data(:))) / (10^(SNR_dB/20)));

%% LOAD scanner rf_data
% rf_data = load('pointsScanner_RFdata.mat').raw; % For scanner points
% rf_data = load('cystScanner_RFdata.mat').raw;   % For scanner cyst

% rf_data = rf_data(90:end,:,:);                  % For both scanner data

%% SA beamforming settings and rf_data preprocessing

% For simulated points
x_range = 20e-3;        % Lateral range (±22 mm)
z_min = 0e-3;           % Start image depth
z_max = 55e-3;          % End image depth

% For scanner points
  % x_range = 30e-3;          % Lateral range (±15 mm)
  % z_min = 0;
  % z_max = 65e-3;

% For simulated cyst
 % x_range = 20e-3;         % Lateral range (±10 mm)
 % z_min = 20e-3;
 % z_max = 40e-3;

% For scanner cyst
 % x_range = 50e-3;         % Lateral range (±25 mm)
 % z_min = 0;
 % z_max = 55e-3;

pixel_size = 0.1e-3;
x_grid = -x_range/2:pixel_size:x_range/2;  % lateral dimension
z_grid = z_min:pixel_size:z_max;           % depth dimension
[X_grid, Z_grid] = meshgrid(x_grid, z_grid);

pixelMap = cat(3, X_grid, Z_grid);
apod = hamming(num_elements); % could be changed

% Define the elements position: 192 elements separated by kerf (pitched array)
arrayPos = (-(num_elements-1)/2:(num_elements-1)/2)*pitch;
arrayPos = [arrayPos.', zeros(num_elements,1)]; % attach 192 zeros to define the z-axis coords

%VS_x = x_start + (0:(num_transmissions-1)) * vs_spacing; % VS positions along the transducer
VS = [vs_x.', repmat(vs_z, num_transmissions, 1)];  % VS positions adding also the depth z

% Precompute constant variable for all pixels within ImageFormation()
si = size(pixelMap);
t = (1:length(rf_data))/fs;
[X, T] = meshgrid(1:length(arrayPos), t);   % X is the array positions, T the number of samples
Xq = repmat(1:length(arrayPos),[si(2), 1]); % Precompute X2 for interpolation

% Initialize LRI
LRIs = zeros(si(1), si(2), num_transmissions);

% Apply Hilbert transform
rf_data_hil = hilbert(rf_data);

%% SA rf_data beamforming process

% Loop over each transmission as each of them corresponds to one different VS
parfor numTrans = 1:num_transmissions
    
    disp(['Processing LRI ', num2str(numTrans), '/', num2str(num_transmissions)]);
    
    % Get RF data for current transmission
    rf_data_VS = rf_data_hil(:, :, numTrans);  
    
    LRIs(:, :, numTrans) = imageFormation(rf_data_VS, pixelMap, VS(numTrans,:), arrayPos, c, pitch, X,T,Xq, x_grid);
end

% save the output to avoid repeating this lenghty process
save('simulation_points_LRIs.mat', 'LRIs');