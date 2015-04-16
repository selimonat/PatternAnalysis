function [mc,name]=pa_GetMotion(subject,phase,what)
%[mc,name]=pa_GetMotion(subject,phase,what)
%
%   Returns motion parameters computed by the SPMs
%   realignement routine. 
%
%   Based on these it further derives new motion
%   features, as indicated in WHAT. These can be the derivatives, squares, squares of
%   derivatives, or absolute values, etc (best read the code).
%
%   For example if WHAT is mc_diff_square, then basic motion parameters,
%   their derivatives and its square are returned.
%
%   The basic motion parameters but NOT the derived ones are z-scored. All
%   time-courses are demeaned.
%
%   Dependency: pa_GetRoot, spm_select

%%
subject_path        = sprintf('%ssub%03d/phase%02d/mrt/nii',pa_GetRoot,subject,phase);
file                = cellstr(spm_select('FPList',subject_path,'^rp_'));
mc                  = load(file{1});
mc                  = zscore(mc);
name                = repmat({'mc'},1,6);

%%
if nargin > 2
    if     strcmp(what,'mc_diff')
        mc                  = demean([mc [diff(mc) ;repmat(0,1,size(mc,2))]]);
        name                = [name repmat({'mc_diff'},1,6)];
    elseif strcmp(what,'mc_square')
        mc                  = demean([mc mc.^2]);    
        name                = [name repmat({'mc_square'},1,6)];
    elseif strcmp(what,'mc_diff_square')%the one that I used for most of the analyses
        mc                  = demean([mc [diff(mc) ;repmat(0,1,size(mc,2))] mc.^2]);
        name                = [name repmat({'mc_diff'},1,6) repmat({'mc_square'},1,6)];
    elseif strcmp(what,'mc_mc2_diff_diff2')%new
        mc                  = demean([mc mc.^2 [diff(mc) ;repmat(0,1,size(mc,2))] [diff(mc).^2; repmat(0,1,size(mc,2))]]);
        name                = [name repmat({'mc_square'},1,6) repmat({'mc_diff'},1,6) repmat({'mc_diff_square'},1,6) ]
    elseif strcmp(what,'mc_abs')
        mc                  = demean([mc abs(mc)]);
        name                = [name 'mc_abs'];
    end
end
