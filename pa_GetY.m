function [Y]=PatternAnalysis_GetY(subject,phase,nroi,threshold,pattern)
% [Y]=PatternAnalysis_GetY(subject,phase,nroi,threshold,pattern)
%
%
% Will return the data for SUBJECT, PHASE, NROI, THRESHOLD from files
% defined with the PATTERN.
%
%for ns = subject_list;Y = pa_GetY(ns,2,1,75,pattern);Y = pa_GetY(ns,4,1,75,pattern);end

%% save all the time courses
% threshold       = 75;
% pattern         = '^swRealigned.*nii$';
cleaned_pattern = pattern(isstrprop(pattern,'alphanum'));
save_path       = sprintf('/Users/onat/Documents/feargen_midlevel/data_matrices/s%02d_p%02d_r%03d_t%03d_%s.mat',subject,phase,nroi,threshold,cleaned_pattern);
Y               = [];
GM              = 100;%global mean
if exist(save_path) == 0
    fprintf('File doesn''t exist so need to save the data first...\n');
    %% load files with a given pattern
    volume_files  = GetNiiFiles(subject,phase,pattern);
    M             = spm_get_space(volume_files{1});
    %% iterate through ROIs and save the time-series one by one...
    for nroi2 = 1:96
        %xyz is in voxel space of the atlas and beta images
        [~,d,xyz]   = LoadCommonAtlas(threshold,nroi2);
        %
        Y           = spm_get_data(volume_files,xyz);                
        % remove the all zero lines.
        i           = sum(abs(Y))==0;
        fprintf('Roi: %03d, removed %04d voxels coz zero...\n',nroi2,sum(i));
        Y           = Y(:,~i);
        save_path2  = sprintf('/Users/onat/Documents/feargen_midlevel/data_matrices/s%02d_p%02d_r%03d_t%03d_%s.mat',subject,phase,nroi2,threshold,cleaned_pattern);
        %% save the time x voxel matrices
        save(save_path2,'Y');
    end    
    %call it self to load from the cache this time.
    [Y]=pa_GetY(subject,phase,nroi,threshold,pattern);
else
    fprintf('%s: File will be loaded from cache...\n',mfilename);
    load(save_path);%will spawn Y
    
% %     %method 1 (spm):
% %     % global normalization: multiply with a factor that makes the global
% %     % mean equal to 100 (session specific grand mean scaling in spmish)
% %     % this could be realized differently so that the mean of the ROI is
% %     % taken as reference rather than spm_global's mean.
% %     g      = PatternAnalysis_GetGlobal(subject,phase,pattern);
% %     factor = GM./mean(g);
% %     Y      = Y*factor;
% %     
    %method 2 (roi specific):
    Y = Y/mean(Y(:))*GM;
end
