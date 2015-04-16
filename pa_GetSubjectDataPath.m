function [file,folder,existence]=pa_GetSubjectDataPath(subject,phase,datatype)
%[file,folder,existence]=pa_GetSubjectDataPath(subject,phase,datatype)
%
%   Returns the path to various data files defined in DATATYPE for SUBJECT and
%   PHASE. DATATYPE is a string and can take the following values: stimulation,
%   triads, scr, eye. EXISTENCE is true if the file exists. FOLDER is the
%   path to the folder which is subject and phase specific.
%
%   Assumes that your data is conform with '%ssub%03d/phase%02d/%s/'
%
%   No dependency: pa_GetRoot
%
%   See also: pa_GetSubjectPath


folder    = sprintf('%ssub%03d/phase%02d/%s/',pa_GetRoot,subject,phase,datatype);
file      = sprintf('%sdata.mat',folder);
existence = exist(file,'file');


    
