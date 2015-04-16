function [volume_files]=pa_GetBOLDFiles(subject_id,phase,pattern)
%[volume_files]=pa_GetBOLDFiles(subject_id,phase,pattern)
%
%   Will gather all the NII files from SUBJECT recorded during PHASE. Will
%   ignore the HR files. The files that are returned matches the
%   PATTERN.
%
%   This function is meaningful after dicom conversion has been done.
%   VOLUME_FILES will be a cell array of size 1.
%
%   Dependency: pa_GetRoot, spm_select
%


volume_files=[];
if nargin == 3
    %get all subject path
    data_path        = sprintf('%ssub%03d/phase%02d/mrt/nii/',pa_GetRoot,subject_id,phase);   
    %
    if ~isempty(data_path)        
        volume_files  = spm_select('FPList',data_path,pattern);
        %this is necessary for SPM
        volume_files  = cellstr(volume_files);
    else
        fprintf('No nii folder is found for subject %d and phase %d\n',subject_id,phase);
    end
else
    fprintf('Need 3 arguments...\n');    
end
