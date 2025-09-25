% Delay-and-Sum Beamforming (DAS process) Function
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


function [delay, apod_vec] = delayCal(pixel, idx, vs, arrayPos, c, x_grid, pitch)

% Computes delay vector relative to the aperture elements activated for a
% specific pixel, setting NaNs for inactive elements of the transducer

% Estimate central aperture element based on current lateral index j
center = round(idx(2)/length(x_grid)*length(arrayPos));

% Determine number of active elements based on depth
active_elements = round(pixel(2)/pitch);

% Compute aperture limits

% START = CENTER - ACTIVE ELEMENTS/2
start_idx = center - round((active_elements/2));
if start_idx < 1
    start_idx = 1;
end

% END = CENTER + ACTIVE ELEMENTS/2
end_idx = center + round((active_elements/2));
if end_idx > length(arrayPos)
    end_idx = length(arrayPos);
end

aperture_idx = start_idx:end_idx;
active_array = arrayPos(aperture_idx, :);

tx_path = sqrt(sum((vs - pixel).^2)) + vs(2);
rx_path = sqrt(sum((active_array - pixel).^2, 2));

% Initialize delay vector with NaNs for inactive elements
delay = nan(length(arrayPos),1);

% Compute delay for active elements only
delay(aperture_idx) = (tx_path + rx_path) / c;

% Apodization vector: Hamming for active elements, 0 elsewhere
apod_vec = zeros(length(arrayPos),1);
apod_vec(aperture_idx) = hamming(length(aperture_idx));
end