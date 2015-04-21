function [volume_handle,roi,xyz]=pa_GetAtlas(atlas,threshold,varargin)
%[volume_handle,roi,xyz]=pa_GetAtlas(atlas,threshold,varargin)
%
%   Returns VOLUME_HANDLE for the all the ROIs present in the ATLAS folder.
%   It basically loads the 4D uncompressed nifti image, read by spm_vol. It
%   adds to this structure the field .roi_names. These are just the volume
%   handles therefore threshold has no effect.
%
%   Optionally the third output argument may contain the binary mask
%   thresholded at THRESHOLD. This will contain .d and .name fields.
%   VARARGIN can contain the index for a roi, if you know which ROI to load
%   exactly, this could speed up loading times.
%
%   XYZ will return the indices of the selected voxels as 4xN matrix
%   suitable for spm_get_data.
%
%   Dependency: spm_vol, spm_read_vols
%
%
%   See also: GetOxfordRoi LoadMergedAtlas, LoadREMAtlas, FeargenROIPrepare
%


%these two are the products of the FeargenROIPrepare.m, which is located at
%the atlas directory.
if ismac
    atlas_path    = ['/Users/onat/Documents/fsldata/atlases/' atlas '/'];
else
    atlas_path    = ['/home/schenk/Documents/MATLAB/RSA/' atlas '/'];
end
volume_handle = spm_vol([atlas_path 'AtlasCleaned.nii']);
name          = textread([atlas_path 'LabelsCleaned.txt'],'%s');
%
for n = 1:length(volume_handle)
    volume_handle(n).roi_name = name{n};
end

if nargout > 1%mask required
    
    if nargin == 1%all ROIs are required.
        roi_index = 1:length(volume_handle);
    else%a specific ROI is required
        roi_index = varargin{1};
    end
    
    roi      = [];
    if roi_index <= length(volume_handle)
        roi.d    = spm_read_vols(volume_handle(roi_index));
        if ~isempty(threshold)%thresholding required.
            roi.d    = roi.d > threshold;
        end
        roi.name = name(roi_index);
    else
        fprintf('ROI index: %03d too big\n',roi_index);
        return
    end
    %get the voxel indices
    if nargout > 2
        [X Y Z]    = ind2sub(volume_handle(1).dim,find(roi.d));
        xyz         = [X Y Z ones(sum(roi.d(:)),1)]';
    end
end
fclose('all');
