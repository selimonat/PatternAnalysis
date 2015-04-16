function [ProjectRoot]=pa_GetRoot
%[ProjectRoot]=pa_GetRoot
%
%   Returns the project path.


if ismac
    ProjectRoot = '/Volumes/feargen2/feargen2/data/';
elseif isunix
    ProjectRoot = '/projects/feargen2/data/';
end
