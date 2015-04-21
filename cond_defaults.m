function [default]=cond_defaults(varargin)
%[default]=pa_defaults(varargin)
%
%Returns either all default values or simply the requested one in
%VARARGIN. If there is no match [] is returned.

%
default.atlas         = 'CommonAtlas';%the name of the folder where atlas is located.
default.troi          = 96;
default.threshold     = 50; %threshold for binarizing probability maps
%
default.fov           = [121 145 121];
%
if ismac 
    default.project_path  = '/Volumes/feargen2/cond/';
else
    default.project_path  = '/projects/cond/';
end
%subject groups
default.group_context     = [24 25 27 28 29 31 32 34 37 39 42 43 45 46 47 49 53 54 56 58 59 60 63 65];
default.group_nocontext   = [7 9 10 11 15 16 17 20 30 33 35 36 38 40 41 44 48 50 52 55 57 62 64 66];

default.metrics           = {'mean' 'cov' 'correlation' 'euclidean' 'seuclidean'};
%
%if requested simply return the field
if nargin == 1    
    if isfield(default,varargin{1})
        default = default.(varargin{1});
    else
        default = [];
        fprintf('no such field present\n');
    end    
end
