function [Y]=pa_GetY(subject,phase,nroi,threshold,pattern)
%[Y]=pa_GetY(subject,phase,nroi,threshold,pattern)
%
%   Will return the raw data for SUBJECT recorded at phase number PHASE from
%   the ROI which is THRESHOLDed from .nii files defined with the PATTERN.
%   On the first call the data will be cached to
%   pa_defaults('save_path'), next calls will simply read this data.
%
%   Dependency: pa_GetBOLDFiles, pa_GetAtlas, pa_defaults, spm_get_data
%
%   Example Usage:
%   for ns = subject_list;Y = pa_GetY(ns,2,1,75,pattern);end

%% save all the time courses

% cleaned pattern necessary for saving
cleaned_pattern = pattern(isstrprop(pattern,'alphanum'));
save_path       = sprintf('%sdata_matrices/s%02d_p%02d_r%03d_t%03d_%s.mat',pa_defaults('save_path'),subject,phase,nroi,threshold,cleaned_pattern);
%
Y               = [];%init
GM              = 100;%global target mean
if exist(save_path) == 0
    fprintf('File doesn''t exist so need to save the data first...\n');
    %% load files with a given pattern
    volume_files  = pa_GetBOLDFiles(subject,phase,pattern);
    M             = spm_get_space(volume_files{1});
    %% iterate through ROIs and save the time-series one by one...
    for nroi2 = nroi
        %xyz is in voxel space of the atlas and beta images
        [~,d,xyz]   = pa_GetAtlas(threshold,nroi2);        
        %
        Y           = spm_get_data(volume_files,xyz);                
        %
        save_path2  = sprintf('%sdata_matrices/s%02d_p%02d_r%03d_t%03d_%s.mat',pa_defaults('save_path'),subject,phase,nroi2,threshold,cleaned_pattern);
        %% save the time x voxel matrices
        save(save_path2,'Y');
    end    
    %call it self to load from the cache this time, this way global
    %correction will be applied
    [Y]=pa_GetY(subject,phase,nroi,threshold,pattern);
else
    fprintf('%s: File will be loaded from cache...\n',mfilename);
    load(save_path);%will spawn Y
    % it might be possible that at this threshold level Y has no valid
    % voxels.
    if isempty(Y) == 1
        Y = [];
        return;
    end    
    %method 1 (spm):
    % global normalization: multiply with a factor that makes the global
    % mean equal to 100 (session specific grand mean scaling in spmish)
    % this could be realized differently so that the mean of the ROI is
    % taken as reference rather than spm_global's mean.
    g      = pa_GetGlobal(subject,phase,pattern);
    factor = GM./mean(g);
    Y      = Y*factor;
    
% %     %method 2 (roi specific):
% %     Y = Y/mean(Y(:))*GM;
end
