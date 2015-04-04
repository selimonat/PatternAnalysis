function [X] = PatternAnalysis_GetDesignMatrix(subject,phase);
% [X] = PatternAnalysis_GetDesignMatrix(subject,phase);
% Will return the Design matrix, X for subject.
%
% TR is fixed to 2.02


%% load the onsets
load(GetSubjectData(subject,phase,'stimulation'));
csp               = p.stim.cs_plus;
dm                = GetDesignMatrix(subject,phase);
dm                = [circshift(dm(:,1:8),[0 4-csp]) dm(:,9:10)];
[onsets conds]    = find(dm);
%this should always give us the correct number of files.
k                 = size(dm,1);
%
% discard always the first trial
onsets(1)         = [];
conds(1)          = [];
%%
fMRI_T                = spm_get_defaults('stats.fmri.t');
fMRI_T0               = spm_get_defaults('stats.fmri.t0');
xBF.T                 = fMRI_T;
xBF.T0                = fMRI_T0;
xBF.dt                = 2.02/xBF.T;
xBF.UNITS             = 'scans';
xBF.Volterra          = 1;
xBF.name              = 'hrf';
xBF                   = spm_get_bf(xBF);

%%
for i = 1:10;%one regressor for each condition
    Sess.U(i).dt        = xBF.dt;%- time bin (seconds)
    
    Sess.U(i).ons       = onsets( conds == i );%- onsets    (in SPM.xBF.UNITS)
    Sess.U(i).name      = {sprintf('%02d',i)};%- cell of names for each input or cause
    
    %no parametric modulation here
    Sess.U(i).dur    =  repmat(0,length(Sess.U(i).ons),1);%- durations (in SPM.xBF.UNITS)
    Sess.U(i).P.name =  'none';
    Sess.U(i).P.P    =  'none';
    Sess.U(i).P.h    =  0;%- order of polynomial expansion
    Sess.U(i).P.i    =  1;%- sub-indices of u pertaining to P
end
SPM.xBF                 = xBF;
SPM.nscan               = k;
SPM.Sess                = Sess;
SPM.Sess.U              = spm_get_ons(SPM,1);
U                       = SPM.Sess.U; %for convinience
%%
% Convolve stimulus functions with basis functions
[X,Xn,Fc]               = spm_Volterra(SPM.Sess.U,SPM.xBF.bf,SPM.xBF.Volterra);
% Resample regressors at acquisition times (32 bin offset)
X                       = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
