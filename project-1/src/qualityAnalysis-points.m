clc, clear, close all

f0 = 7e6;

%load simulation_points_LRIs.mat;
load scanner_points_LRIs.mat;

% For simulated points
% x_range = 20e-3;        % Lateral range (±10 mm)
% z_min = 0e-3;              % Start depth
% z_max = 55e-3;          % End depth
% c = 1540;
% target_depths = [10.3, 20.3, 30.3, 40.3, 50.3]/1000;

% For scanner points
x_range = 30e-3;          % Lateral range (±15 mm)
z_min = 0;
z_max = 65e-3;
target_depths = [27.8, 52.4] / 1000;
c = 1480;

pixel_size = 0.1e-3;
x_grid = -x_range/2:pixel_size:x_range/2;  % lateral dimension
z_grid = z_min:pixel_size:z_max;  % depth dimension
[X_grid, Z_grid] = meshgrid(x_grid, z_grid);

lambda = c/f0;
max_ref = max(sum(abs(LRIs), 3), [], 'all');

%% HRI visualzation
HRI = sum(LRIs,3);
HRI_dB = 20 * log10(abs(HRI)/max_ref);

figure;
imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
clim([-60 0]);
c = colorbar;
c.Label.String = 'dB';
colormap gray ;
xlabel('Lateral Position [mm]');
ylabel('Depth [mm]');
title('High Resolution Image');
axis ij equal tight;

%% FWHM calculation (only for simulated points)

% Convert target depths to pixel indices
[~, depth_indices] = min(abs(z_grid - target_depths'), [], 2);

% Plot HRI with lateral lines
figure;
imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
colormap gray;
clim([-60, 0]);
c = colorbar;
c.Label.String = 'dB';
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
title('HRI');
axis ij equal tight;
hold on;

% Upsample laterla profile before FWHM calculation
scaling_factor = 10;
x_fine = linspace(x_grid(1), x_grid(end), length(x_grid)*scaling_factor);

lateral_lines = zeros(length(depth_indices), length(x_fine));
fwhm_values = zeros(length(depth_indices), 1);

half_max = 0.5;

for i = 1:length(depth_indices)

    % Extract the lateral line at the specified depth
    lateral_line = abs(HRI(depth_indices(i), :));

    % Interpolate the lateral line to the finer x_grid
    lateral_line_fine = interp1(x_grid, lateral_line, x_fine, 'linear', 'extrap');

    % Normalize the line between 0 and 1
    line_norm = lateral_line_fine - min(lateral_line_fine);
    line_norm = line_norm / max(line_norm);

    % update lateral_line with its norm vlaue for plotting
    lateral_lines(i, :) = line_norm;

    % Find and extract peak location
    [peak_val, peak_idx] = max(line_norm);
    peak_pos = x_fine(peak_idx);

    % Search for last leftmost idx crossing the half-peak value
    left_idx = find(line_norm(1:peak_idx) <= half_max*peak_val, 1, 'last');

    if isempty(left_idx)
        left_pos = x_fine(1); % Edge case
    else

        left_x = x_fine(left_idx:left_idx+1);
        left_y = line_norm(left_idx:left_idx+1);
        left_pos = interp1(left_y, left_x, half_max*peak_val, 'linear', 'extrap');
    end

    % Search for right half-max crossing (values above peak)
    right_idx = find(line_norm(peak_idx+1:end) <= half_max*peak_val, 1, 'first')+peak_idx;

    if isempty(right_idx)
        right_pos = x_fine(end); % Edge case
    else
        right_x = x_fine(right_idx-1:right_idx);
        right_y = line_norm(right_idx-1:right_idx);
        right_pos = interp1(right_y, right_x, half_max*peak_val, 'linear', 'extrap');
    end

    % Calculate FWHM in m
    fwhm_values(i) = right_pos - left_pos;

    % Plot the lateral line on the HRI
    plot(x_fine * 1e3, ones(size(x_fine)) * z_grid(depth_indices(i)) * 1e3, 'r-', 'LineWidth', 0.3);

end

figure;
hold on
grid on
plot_handles = gobjects(1, length(depth_indices));
for i = 1:length(depth_indices)
    
    % Plot the lateral line intensity profile
    plot_handles(i) = plot(x_fine * 1e3, lateral_lines(i, :), 'Linewidth', 1.5, 'DisplayName', ...
        [num2str((target_depths(i)-0.3e-3)*1000,'%.0f') ' mm,',' FWHM: ' num2str(fwhm_values(i) * 1000, '%.2f') ' mm']);

    x_min = x_grid(91) * 1e3;  %96
    x_max = x_grid(111) * 1e3;  %106
    xlim([x_min, x_max]);

    y_max = max(lateral_lines(i, :));
    ylim([0, y_max*1.1]);
end

% Add labels and title
xlabel('Lateral position [mm]');
ylabel('Intensity (normalized)');
title('Intensity Profiles');

% Create the legend after the loop using the handles
legend(plot_handles, 'Location', 'northeast');
hold off; % Release the hold to allow further plotting if needed

% Plot of lateral resolution vs depth
figure;
plot(target_depths * 1000, fwhm_values*1000, '-o');
title('FWHM vs depth');
xlabel('Depth [mm]');
ylabel('Lateral resolution [mm]');
grid on;

%% ROI DEFINITION

% general_roi = [x1 (left), x2 (right); z1 (top), z2 (bottom)]

% SIMULATION POINTS
%total_roi = [-5,5;5,15]*1e-3; % meters

% SCANNER POINTS
 total_roi = [-5.5,4.5;23,33]*1e-3; % meters

% Visualize first ROI
HRI = sum(LRIs,3);
HRI_dB = 20 * log10(abs(HRI)/ max_ref);

figure;
imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
clim([-60 0]);
c = colorbar;
c.Label.String = 'dB';
colormap gray ;
xlabel('Lateral Position [mm]');
ylabel('Depth [mm]');
axis ij equal tight;
hold on;

r1 = rectangle('Position', [total_roi(1,1), total_roi(2,1), total_roi(1,2)-total_roi(1,1), total_roi(2,2)-total_roi(2,1)]*1000, 'EdgeColor', 'r', 'LineWidth', 1);
plot(NaN, NaN, 'r', 'LineWidth', 2);

title('HRI');
legend({'ROI'});

%% PSF PLOT

% Define contour levels at 6 dB intervals
min_val = -60;
max_val = 0;
levels = min_val:6:max_val;

% Create contour plot
figure;
contour(x_grid*1000, z_grid*1000, HRI_dB, levels);
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
axis ij equal tight;
title('Point PSF');
c = colorbar;
c.Label.String = 'dB';
hold on;

r = 2.5 * lambda;

% Generate circle points
theta = linspace(0, 2*pi, 100); % 100 points for smoothness
x_circle = (total_roi(1,1)+total_roi(1,2))/2 + r * cos(theta);
z_circle = (total_roi(2,1)+total_roi(2,2))/2 + r * sin(theta);

% Convert to mm and plot
plot(x_circle * 1000, z_circle * 1000, 'r', 'LineWidth', 1);
grid on

ylim([z_grid(251), z_grid(351)]*1000);

hold off;

%% CR CALCULATION FOR BOTH POINTS CASES

summation_rate = [10,8,6,4,2,1];

% For subplots
if mod(length(summation_rate),2) == 0
    cols = length(summation_rate)/2;
    rows = length(summation_rate)/cols;
else
    cols = (length(summation_rate)+1)/2;
    rows = (length(summation_rate)+1)/cols;

end

CR_values = zeros(length(summation_rate), length(target_depths));

figure;
for s = 1:length(summation_rate)

    step = summation_rate(s);
    indices = 1:step:size(LRIs,3);
    HRI = sum(LRIs(:,:,indices), 3);

    for i = 1:length(target_depths)

        CR = crCalc(x_grid, z_min, total_roi, pixel_size, lambda, HRI);
        CR_values(s,i) = CR;

        % Shift the ROI to next point target before
        % the new calculation:

        % for simulation
        %total_roi(2,:)=total_roi(2,:)+10e-3;

        % for scanner
        total_roi(1,:)=total_roi(1,:)+1e-3;
        total_roi(2,:)=total_roi(2,:)+25e-3;
    end

    % Restore the ROI to the initial location before
        % evaluating a new summation strategy

       % total_roi = [-5,5;5,15]*1e-3; %mm  for simulation
        total_roi = [-5.5,4.5;22.5,32.5]*1e-3; %mm for scanner

    % Show how the contrast changes
    HRI_dB = 20 * log10(abs(HRI)/max_ref);
    subplot(rows, cols, s); % Select the subplot for the current summation type
    imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
    colormap gray;
    colorbar;
    clim([-60 0]);
    title(sprintf('%d LRIs [1/%d]', length(indices), step));
    xlabel('Lateral Position [mm]');

    if s==1 || s==4
        ylabel('Depth [mm]');
    else
    end
    
    axis ij equal tight;
    hold on;
end

%% CR TREND PLOTS

% Create a n° of labels based on the n° of summation type
summation_labels = cell(1, length(summation_rate));

for s = 1:length(summation_rate)

    if summation_rate(s) == 1
        summation_labels{s} = 'All frames';
    else
        summation_labels{s} = ['1/', num2str(summation_rate(s)), ' frames'];
    end
end

figure;
colors = lines(length(target_depths));

for i = 1:length(target_depths)

    %scatter(1:length(summation_rate), CR_values(:,i), 50, colors(i,:), 'filled');
    plot(1:length(summation_rate), CR_values(:,i), 'Color', colors(i,:), 'LineWidth', 1);
    hold on;
end

legend_labels = arrayfun(@(d) sprintf('%.0f mm', d), target_depths*1000, 'UniformOutput', false);
legend(legend_labels, 'Location', 'best');

xticks(1:length(summation_rate));
xticklabels(summation_labels);
xtickangle(45);
xlabel('Summation Strategy');
ylim([min(CR_values(:))-1, max(CR_values(:))+1])
ylabel('CR (dB)');
title('CR over different summation Strategies');
grid on;

%% ITERATIVE LRIs VISUALIZATION

fig = figure;
ax = axes('Parent', fig);

% Initialize progressive summation with the first LRI
HRI = LRIs(:,:,1);
HRI_dB = 20 * log10(abs(HRI) / max_ref);

hImg = imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
colormap('gray');
clim([-120, 0]);
colorbar;
xlabel('Lateral position [mm]');
ylabel('Depth [mm]');
title('HRI');
axis ij equal tight;

% Create slider for scrolling through summation steps
slider = uicontrol('Style', 'slider', ...
    'Min', 1, ...
    'Max', size(LRIs,3), ...
    'Value', 1, ...
    'SliderStep', [1/(size(LRIs,3)-1) , 0.1], ...
    'Position', [400 45 20 330], ...
    'BackgroundColor', [0.8 0.8 0.8], ... % Change slider color
    'Callback', @updateImage);

% Create text label for slider
slider_label = uicontrol('Style', 'text', ...
    'Position', [361.5 375 100 20], ...
    'String', sprintf('LRIs: %d', 1),...
    'FontWeight', 'bold');

% Store required data in `guidata`
data.LRIs = LRIs;
data.x_grid = x_grid;
data.z_grid = z_grid;
data.ax = ax;
data.hImg = hImg;
data.slider_label = slider_label;

guidata(fig, data); % Store data in figure

%% Functions

% Callback function to update the image based on slider position
function updateImage(src, ~)
data = guidata(src);
numToSum = round(src.Value);  % Get slider value as an integer
HRI = sum(data.LRIs(:,:,1:numToSum), 3); % Sum first 'numToSum' LRIs
max_ref = max(sum(abs(data.LRIs), 3), [], 'all');
HRI_dB = 20 * log10(abs(HRI) / max_ref);
set(data.hImg, 'CData', HRI_dB); % Update image
set(data.slider_label, 'String', sprintf('LRIs: %d', numToSum)); % Update label
end
