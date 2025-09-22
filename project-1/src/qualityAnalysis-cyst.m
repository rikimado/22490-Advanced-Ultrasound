clc, clear, close all

load cystSimulation_LRIs.mat;
%load cystScanner_LRIs.mat;

% For simulated cyst
  x_range = 20e-3;       % Lateral range (±10 mm)
  z_min = 20e-3;
  z_max = 40e-3;

% For scanner cyst
  % x_range = 50e-3;        % Lateral range (±25 mm)
  % z_min = 0;
  % z_max = 55e-3;

pixel_size = 0.1e-3;
x_grid = -x_range/2:pixel_size:x_range/2;  % lateral dimension
z_grid = z_min:pixel_size:z_max;  % depth dimension
[X_grid, Z_grid] = meshgrid(x_grid, z_grid);
%%
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

%% ROI DEFINITION

% general_roi = [x1 (left), x2 (right); z1 (top), z2 (bottom)]

% SIMULATION CYST
 cyst_roi = [0,2;29.5,31.5]; %mm
 tissue_roi = [3,5;29.5,31.5];

% SCANNER CYST 1
 % cyst_roi = [-1,1;15.5,17.5]; %mm
 %  tissue_roi = [2.5,4.5;15.5,17.5];

% SCANNER CYST 2
  % cyst_roi_2 = [0,2;46,48]; %mm
  % tissue_roi_2 = [3.5,5.5;46,48];

% Visualize ROIs
HRI = sum(LRIs,3);
max_ref = max(abs(HRI(:)));
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

% For simulated cyst and scanner cyst 1
r1 = rectangle('Position', [cyst_roi(1,1), cyst_roi(2,1), cyst_roi(1,2)-cyst_roi(1,1), cyst_roi(2,2)-cyst_roi(2,1)], 'EdgeColor', 'r', 'LineWidth', 1.5); % Cyst ROI
r2 = rectangle('Position', [tissue_roi(1,1), tissue_roi(2,1), tissue_roi(1,2)-tissue_roi(1,1), tissue_roi(2,2)-tissue_roi(2,1)], 'EdgeColor', 'b', 'LineWidth', 1.5); % Tissue ROI
plot(NaN, NaN, 'r', 'LineWidth', 2); % red for cyst
plot(NaN, NaN, 'b', 'LineWidth', 2); % green for tissue

% Add also below for scanner cyst 2
 % r3 = rectangle('Position', [cyst_roi_2(1,1), cyst_roi_2(2,1), cyst_roi_2(1,2)-cyst_roi_2(1,1), cyst_roi_2(2,2)-cyst_roi_2(2,1)], 'EdgeColor', 'r', 'LineWidth', 1); % Cyst ROI
 % r4 = rectangle('Position', [tissue_roi_2(1,1), tissue_roi_2(2,1), tissue_roi_2(1,2)-tissue_roi_2(1,1), tissue_roi_2(2,2)-tissue_roi_2(2,1)], 'EdgeColor', 'g', 'LineWidth', 1); % Tissue ROI
 % plot(NaN, NaN, 'r', 'LineWidth', 2); % red line for cyst
 % plot(NaN, NaN, 'g', 'LineWidth', 2); % green line for tissue

title('HRI');
legend({'Cyst ROI', 'Tissue ROI'});

%% CNRs CALCULATION

summation_rate = [10,8,6,4,2,1];

% For subplots
if mod(length(summation_rate),2) == 0 % even
    cols = length(summation_rate)/2;
    rows = length(summation_rate)/cols;
else
    cols = (length(summation_rate)+1)/2;
    rows = (length(summation_rate)+1)/cols;

end

CNR_values = zeros(length(summation_rate), 1);  % for simulated cyst
%CNR_values = zeros(length(summation_rate), 2);  % for scanner cyst

figure;
for s = 1:length(summation_rate)

    step = summation_rate(s);
    indices = 1:step:size(LRIs,3);
    HRI = sum(LRIs(:,:,indices), 3);
    max_ref = max(abs(HRI(:)));

    % CNR for simulation and scanner cyst
    CNR = cnrCalc(x_grid,z_min, pixel_size, HRI/max_ref, cyst_roi,tissue_roi);
    CNR_values(s,1) = CNR;

    % Add also below for scanner cyst 2
     % CNR_2 = cnrCalc(x_grid,z_min, pixel_size, HRI/max_ref, cyst_roi_2,tissue_roi_2);
     % CNR_values(s,2) = CNR_2;

    % Show how the contrast changes
    HRI_dB = 20 * log10(abs(HRI) / max_ref);
    subplot(rows, cols, s); % Select the correct subplot for the current summation type
    imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
    colormap gray;
    colorbar;
    clim([-60 0]); % Adjust dynamic range
    title(sprintf('%d LRIs [1/%d]', length(indices), step));
    xlabel('Lateral Position [mm]');
    hold on;
    
    if s==1 || s==4
        ylabel('Depth [mm]');
    else
    end

    axis ij equal tight;
end


%% CNR TREND PLOTS

% Create a n° of labels based on the n° of summation type
summation_labels = cell(1, length(summation_rate));

for s = 1:length(summation_rate)
    if summation_rate(s) == 1
        summation_labels{s} = 'All frames';  % Label for "All Frames"
    else
        summation_labels{s} = ['1/', num2str(summation_rate(s)), ' frames'];  % Other cases
    end
end

figure;
hold on;

%for simulated and scanner cyst 1
%scatter(1:length(summation_rate), CNR_values(:,1),80, 'r', 'filled');
plot(1:length(summation_rate), CNR_values(:,1), 'r');

% Add also below for scanner cyst 2
 %scatter(1:length(summation_rate), CNR_values(:,2),80, 'b', 'filled');
 plot(1:length(summation_rate), CNR_values(:,2), 'b');
 legend({'Upper cyst', 'Bottom cyst'}, 'Location', 'best');

xticks(1:length(summation_rate));
xticklabels(summation_labels);
xtickangle(45);

xlabel('Summation Strategy');
ylabel('CNR (dB)');
title('CNR over different summation Strategies');
grid on;


%% ITERATIVE LRIs SUMMATION

fig = figure;
ax = axes('Parent', fig);

% Initialize progressive summation with the first LRI
HRI = LRIs(:,:,1);
HRI_dB = 20 * log10(abs(HRI) / max_ref);

hImg = imagesc(x_grid * 1e3, z_grid * 1e3, HRI_dB);
colormap('gray');
clim([-60, 0]);
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
