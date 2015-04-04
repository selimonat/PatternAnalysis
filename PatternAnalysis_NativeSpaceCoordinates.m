function [xyz,atlas]=PatternAnalysis_NativeSpaceCoordinates(subject,roi,threshold)
%[xyz,atlas]=PatternAnalysis_NativeSpaceCoordinates(subject,roi,threshold);
%
%   Returns XYZ coordinates for roi ROI in the native space binarized at a
%   threshold level THRESHOLD (default = 50%). Threshold defines the
%   prcntile. You can directly call spm_get_data with xyz values.
%
%   See also: PatternAnalysis_GetNativeData

if nargin < 3
    threshold = 50;
end
h          = Atlas2NativeSpace(subject);
atlas      = h(roi);
mask       = spm_read_vols(atlas);
% probability values of the mask is now a little bit distorded wrt to the
% original mni mask. that's why I use percentile metric to threshold
% voxels.
threshold  = prctile(mask(mask ~= 0),threshold);

[X Y Z]    = ind2sub(atlas.dim,find(mask > threshold));
xyz        = [X Y Z ones(length(X),1)]';