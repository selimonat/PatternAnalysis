function [D,xyz]=PatternAnalysis_GetNativeData(subject,phase,roi,threshold)
% [D,xyz]=PatternAnalysis_GetNativeData(subject,phase,roi,threshold)
% 
% D is a time x space matrix.
%
% See also: PatternAnalysis_NativeSpaceCoordinates


if nargin == 3
        threshold = 50;
end
filename    = sprintf('%ssub%03d/phase%02d/midlevel/%s_%03d_%03d.mat',GetProjectRootFeargen2,subject,phase,mfilename,roi,threshold)
D           = [];
xyz         = [];
%
if exist(filename) == 0        
    %
    fprintf('%s: Getting native data for:\n Subject: %03d, Phase: %02d, ROI: %03d, Threshold: %03d.\n',mfilename,subject,phase,roi,threshold);
    %
    [xyz,atlas]  = PatternAnalysis_NativeSpaceCoordinates(subject,roi,threshold);
    p            = sprintf('%ssub%03d/phase%02d/mrt/nii/',GetProjectRootFeargen2,subject,phase);%path
    files        = spm_select('FPList',p,'^Realigned.*nii$');
    sprintf('Here are few files:\n');
    display(files(1:10,:));
    fprintf('...\n');
    D            = spm_get_data(files,xyz);
    %remove the full zero time-series
    i            = sum(D == zeros(size(D))) == size(D,1);
    D(:,i)       = [];
    fprintf('%03d voxels were found to be outside of the brain in the native space.\n',sum(i));
    xyz(:,i)     = [];%remove also the coordinates.
    %add outlier detection maybe
    
    save(filename,'D','xyz');
else
    fprintf('Loading file:\n%s\nfrom cache...\n',filename);
    load(filename);
end