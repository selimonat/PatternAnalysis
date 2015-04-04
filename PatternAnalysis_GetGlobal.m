function [g]=PatternAnalysis_GetGlobal(subject,phase,pattern)
% [g]=PatternAnalysis_GetGlobal(subject,phase,pattern)
%
% Will return the global activity for SUBJECT, PHASE, from files
% defined with the PATTERN.
%


%% save all the time courses
cleaned_pattern = pattern(isstrprop(pattern,'alphanum'));
save_path       = sprintf('/Users/onat/Documents/feargen_midlevel/globals/s%02d_p%02d_%s.mat',subject,phase,cleaned_pattern);
g               = [];
if exist(save_path) == 0
    fprintf('File doesn''t exist so need to save the data first...\n');
    %% load files with a given pattern
    volume_files  = GetNiiFiles(subject,phase,pattern);
    volume_files  = strvcat(volume_files);    
    %% iterate through ROIs and save the time-series one by one...
    g             = spm_global(spm_vol(volume_files));    
    %% save the time x voxel matrices
    save(save_path,'g');         
else
    fprintf('%s: File will be loaded from cache...\n',mfilename);
    load(save_path);
end
