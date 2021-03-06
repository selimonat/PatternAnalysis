function [beta]=pa_GetBetas(X,Y,N)
%[beta]=pa_GetBetas(X,Y,N)
%
% Returns the weights for the design matrix X that fits best the observed
% BOLD time-series in Y. Y and X are high-pass filtered as in spm. Can be
% update to contain the residuals.
%
% Dependency: spm_sp, spm_filter
%



%% Get the High-pass filter by default
K  = struct('HParam', pa_defaults('HParam') , 'row',    1:size(X,1) , 'RT',     pa_defaults('RT') );
Y  = spm_filter(K,Y);
%%
% N  = [N K.X0];%combine covariates
% N  =  N - repmat(mean(N),size(N,1),1);%and demean them

%% first residualize with respect to Nuissance parameters.
% Yr   = Y - N*(N\Y);
% beta = (X\Yr)';
DM        = [X N ones(size(X,1),1)];
DM        = spm_filter(K,DM);
DM        = spm_sp('Set',DM);
DM        = spm_sp('x-',DM);% projector;
beta      = DM*Y;
%ignore the N coefficients.
beta      = beta(pa_defaults('conditions'),:);