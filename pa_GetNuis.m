function [N]=pa_GetNuis(subject,phase,pattern)
%[N]=pa_GetNuis(subject,phase,pattern)
%
% this function DEPENDS ON THE ATLAS USED as far as extraction of ventricle
% time-courses are concerned.
%
% 

%% Get the motion nuissance variables.
[Nuis]                = GetNuissanceVariables(subject,phase,pattern,{'mc_diff_square'});
tNuis                 = length(Nuis);
%% Get also the ventricle time-courses.
ventricle             = 36;%necessary for constructing a nuissance regressor
threshold             = 75;
Nuis(tNuis + 1).val   = zscore(mean(pa_GetY(subject,phase,ventricle,threshold,pattern),2));
Nuis(tNuis + 2).val   = zscore(mean(pa_GetY(subject,phase,mod(ventricle-1,46)+1+46,threshold,pattern),2));
Nuis(tNuis + 1).name  = 'RH_csf';
Nuis(tNuis + 2).name  = 'LH_csf';
% final output
N                     = [Nuis(:).val];