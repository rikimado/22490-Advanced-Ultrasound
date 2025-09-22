function cnr_value = cnrCalc(x_grid,z_min, pixel_size, HRI, roi_cyst, roi_bgnd)

% general_roi definition = [x1, x2; z1, z2]

% center of roi_cyst region at half of x_grid
center_idx = round(length(x_grid)/2);

% X-range for cyst region
cyst_ROI_x = center_idx+(roi_cyst(1,1)*1e-3)/pixel_size:center_idx+(roi_cyst(1,2)*1e-3)/pixel_size;

% Z-range for cyst region
cyst_ROI_z = 1+(roi_cyst(2,1)*1e-3 -z_min)/pixel_size:1+(roi_cyst(2,2)*1e-3-z_min)/pixel_size;

% X-range for tissue region
tissue_ROI_x = center_idx+(roi_bgnd(1,1)*1e-3)/pixel_size:center_idx+(roi_bgnd(1,2)*1e-3)/pixel_size;

% z-range for tissue region
tissue_ROI_z = (1+(roi_bgnd(2,1)*1e-3 -z_min)/pixel_size):(1+(roi_bgnd(2,2)*1e-3 -z_min)/pixel_size);

cyst_ROI = HRI(round(cyst_ROI_z), cyst_ROI_x); 
meanCyst = mean(abs(cyst_ROI(:)));
stdCyst = std(abs(cyst_ROI(:)));

tissue_ROI = HRI(round(tissue_ROI_z), tissue_ROI_x);
meanTissue = mean(abs(tissue_ROI(:)));
stdTissue = std(abs(tissue_ROI(:)));

%stdTissue = std(abs(tissue_ROI(:)));

% Take the abs of the difference to ensure to capture the contrast 
% regardless of whether the cyst is brighter or darker than the 
% surrounding tissue.
cnr_value = abs(meanCyst - meanTissue)/sqrt(stdCyst^2+stdTissue^2);
%cnr_value = abs(meanCyst - meanTissue)/stdTissue;
%cnr_value = abs(meanCyst - meanTissue);
cnr_value = 20*log10(cnr_value);
%fprintf('CNR = %.2f\n', cnr_value);
end
