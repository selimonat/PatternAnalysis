function [nuisance,names]=pa_GetNuissance(subject,phase,pattern,what)
%[nuisance,names]=pa_GetNuissance(subject,phase,pattern,what)
%
%   Returns time-series usable as regressors of nuissance. These can be
%   subject's head-motion (as detected by spm realignement), their absolute
%   values, their derivatives, retroicor heartbeat and respiration
%   regressors and time-course of CSF voxels. 
%
%   WHAT is a string of cells that defines what is to be returned.
%   Possibilities are mc, mc_abs, mc_mc2_diff_diff2, mc_diff_square,
%   mc_square, mc_diff, hb, resp, csf. If WHAT is omitted, then
%   'mc_diff_square' 'hb' 'resp' 'csf' are returned.
%
%   NAMES contains human-readable descriptions.
%
%   Use hb and resp with your own risk.
%
%   See also: pa_GetMotion
%   
%   Dependency: pa_GetCSF, pa_GetMotion,
%   pa_GetRetroicorHeartbeat, pa_GetRetroicorRespiration

%default behavior
if nargin == 3
    what = {'mc_diff_square' 'hb' 'resp' 'csf'};%WHAT to compute
end
%
c = 0;
for n = 1:length(what)
    filename = sprintf('%ssub%03d/phase%02d/midlevel/%s_%s.mat',pa_GetRoot,subject,phase,mfilename,what{n});
    if exist(filename) == 0 || 1
        fprintf('%s: Getting %s Nuis.Var.\n',mfilename,what{n})        
        %get the data
        if strcmp(what{n}(1:2),'mc')            
            [dummy, names] = pa_GetMotion(subject,phase,what{n});
        elseif strcmp(what{n},'hb')
            [dummy,names]  = pa_GetRetroicorHeartbeat(subject,phase);
        elseif strcmp(what{n},'resp')
            [dummy,names]  = pa_GetRetroicorRespiration(subject,phase);
        elseif strcmp(what{n},'csf')
            [dummy,names]  = pa_GetCSF(subject,phase,pattern);
        end
        %store in a way SPM likes.
        for nreg = 1:size(dummy,2)
            c = c+1;
            nuisance(c).val  = dummy(:,nreg);
            nuisance(c).name = sprintf('%s_%02d',names{nreg},nreg);
        end
        %save it
        dummy = [];
        save(filename,'nuisance','names');
    else
        fprintf('%s:\nFile cached: %s\nLoading from cache...\n',mfilename,filename);
        load(filename);
    end
end
