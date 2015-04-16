function [folder]=pa_GetSubjectPath(subject,phase)
%[folder]=pa_GetSubjectPath(subject,phase);
%
%   It returns the path to the data of SUBJECT recorded at PHASE.
%
%   See also: pa_GetSubjectPath
%
%   Dependency: pa_GetRoot
%
%   See also: pa_GetSubjectDataPath


folder = sprintf('%ssub%03d/phase%02d/',pa_GetRoot,subject,phase);
