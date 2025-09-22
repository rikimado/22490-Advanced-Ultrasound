function CR = crCalc(x_grid, z_min, total_roi, pixel_size, lambda, HRI)

% Define first (x1,z1) and last (x2,z2) location of lateral and depth edge of ROI as: 
% general_roi = [x1, x2; z1, z2] 

x_start = round((total_roi(1,1) - x_grid(1)) / pixel_size)+1;  % Start index based on total_roi(1,1)
x_end   = round((total_roi(1,2) - x_grid(1)) / pixel_size)+1;  % End index based on total_roi(1,2)

% Define total_ROI pixel indices
%total_ROI_x = (center_idx + total_roi(1,1)/pixel_size):(center_idx + total_roi(1,2)/pixel_size);
total_ROI_x = x_start:x_end;
total_ROI_z = 1+(total_roi(2,1) - z_min)/pixel_size:1+(total_roi(2,2) - z_min)/pixel_size;

% Define circle radius in pixel
point_radius = 2.5*(lambda/pixel_size);

% Extract signal within the ROI
total_ROI = abs(HRI(round(total_ROI_z), round(total_ROI_x)));

% Create meshgrid of ROI positions
[X, Z] = meshgrid((1:length(total_ROI_x)), (1:length(total_ROI_z)));

% Compute distance of each ROI pixel from the center
distances = sqrt((X - length(total_ROI_x)/2).^2 + (Z - length(total_ROI_z)/2).^2);

% Mask pixels that are outside the circular point ROI
out_mask = distances > point_radius;

% Extract signal values outside the circle
out_ROI = total_ROI(out_mask);
%out_ROI = out_ROI(out_ROI ~= 0);

% Compute energies
E_out = sum(out_ROI(:).^2);
E_total = sum(total_ROI(:).^2);

CR = sqrt(E_out/E_total);
CR = 20 * log10(CR);
end