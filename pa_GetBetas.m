function [beta]=pa_GetBetas(X,Y,N,HParam)
%[beta]=pa_GetBetas(X,Y,N,HParam)
%
% Return beta values best fitting the design matrix X to observed BOLD time
% courses in Y, after residualizing with the Nuissance matrix N, which will
% be extended with spm-like high-pass filter matrix with parameter HPARAM.
%
% important: TR is fixed in the function body.



%% plot the correlation matrix of X
% %             figure(1);
% %             subplot(1,2,1);
% %             imagesc(corrcoef(X).^2,[0 1]);colorbar;
% %             SaveFigure(sprintf('%ssub%03d/phase%02d/figures/%s.png',GetProjectRootFeargen2,subject,phase,'DesignCorrelation'));

%% Get the High-pass filter by default
K = struct('HParam', HParam , 'row',    1:size(X,1) , 'RT',     2.02);
K = spm_filter(K);
%%
N = [M K.X0];%combine covariates
N  = N - repmat(mean(N),size(N,1),1);%and demean them
%% plot the nuissance correlation
%             figure(1)
%             subplot(1,2,2);
%             imagesc(corrcoef(M).^2,[0 1]);
%             colorbar;
%             SaveFigure(sprintf('%ssub%03d/phase%02d/figures/%s.png',GetProjectRootFeargen2,subject,phase,'NuissanceCorrelation'));
%% first residualize with respect to Nuissance parameters.
Yr   = Y - N*(N\Y);
beta = (X\Yr)';