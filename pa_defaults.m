function [default]=pa_defaults(varargin)
%[default]=pa_defaults(varargin)
%
%Returns either all default values or simply the requested one in
%VARARGIN. If there is no match [] is returned.

%data path
default.save_path = '/Users/onat/Documents/pa_midlevel/';

%atlas business
default.atlas     = 'FeargenMerged';%the name of the folder where atlas is located.
default.troi      = 92;
default.ventricle = [36 36+default.troi/2];
default.threshold = 50;

%subjects
default.subject_list = [2 3 4 5 6 8 9 11 12 13 14 16 18 19 20 21 22 25 26 27 28 29 30 31 32 33 34 35 36];
default.gs_index     = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
default.gs           = default.subject_list(default.gs_index);%tuned subjects
default.bs           = setdiff(default.subject_list,default.gs);%untuned subjects

%fmri
default.HParam       = 128;
default.RT           = 2.02;

%stimuli
default.design       = @pa_GetDesignMatrix;%will return design, onsets
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
