function [tc,name]=pa_GetCSF(subject,phase,pattern)
%[tc,name]=pa_GetCSF(subject,phase,pattern)
%
%   Returns two separate time-series representing the average time-course of
%   the CSF voxels (threshold = 0). Make only sense with MNI normalized
%   volumes. PATTERN is the filename pattern for the files to 
%   extract the time-courses from. 
%
%   Time-courses are z-score transformed.
%
%   Dependency: pa_defaults, pa_GetY

%% 
tc        = [];
threshold = 0;
for ventricle  = pa_defaults('ventricle');%these are ventricle indices of the oxford roi.    
    %oxford indices of the volumes
    Y   = pa_GetY(subject,phase,ventricle,threshold,pattern);
    Y   = mean(Y,2);%average across subjects
    tc  = [tc zscore(Y)];%concat
end
%
name = repmat({'CSF'},1,2);
