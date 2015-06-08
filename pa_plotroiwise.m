%function pa_plonroiwise(simmat)
% % addpath('\\samba\schenk\Documents\MATLAB\RSA\globalfunctions\');
% % results_path='\\samba\schenk\Documents\MATLAB\RSA\';
% % load(sprintf('%ssimmatpain1regr.mat',results_path));
% % load(sprintf('%sbetaspain1regrneu2.mat',results_path));
%%
simmat=sim.raw.cov; %[condition,condition,subject,roi,phase]
% load(sprintf('%sbetaspain1regrneu1.mat',results_path));
nsub=size((sim.mc.cov),3);
%% atlas business
betamat=betas.mc.cov;
rsqmat=rsq;
if ismac
    atlas=cond_defaults('atlas');
else
    atlas=cond_defaults('atlas');
end

nroi      = size(simmat,4);
ncond     = size(simmat,1);
nreg = size(betamat,1)-2;
% try
%     [row,col] = GetSubplotNumber(nroi/2);
%     row=row+1;
% catch
row=5;
col=4;
%     col=9;
% end
[d u]     = GetColorMapLimits(simmat(:),3);
fontsize  = 7;

%% model
V1=[0.4 0.8 0.6 0.6]; % temp
V2=[1 0 1 0]; %placcontext
V3=[0 1 0 1]; %kontrollcontext
V4=[1 1 -1 -1]; %session

O1=pa_outerproduct(V1);
O2=pa_outerproduct(V2);
O3=pa_outerproduct(V3);
O4=pa_outerproduct(V4);

% O1        = [0 -1 0 0;-1 0 0 0;0 0 0 1;0 0 1 0;];
% O2        = [0 0 1 0;0 0 0 0;1 0 0 0;0 0 0 0];
% O3        = [0 0 0 0;0 0 0 1;0 0 0 0;0 1 0 0];
% O4        = [0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0];%sess

D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
K= ones(4,4); % konstant
nreg=6; %5 %6

%%%%%%%%%%%%%%%%%%%% figure
figure;
%% plot model

h            = subplot(row,col,1);
subplotChangeSize(h,.009,.009);
imagesc(O1);
title('models');
axis image
axis off
h            = subplot(row,col,2);
subplotChangeSize(h,.009,.009);
imagesc(O2);
axis image
axis off
h            = subplot(row,col,3);
subplotChangeSize(h,.009,.009);
imagesc(O3);
axis image
axis off
h            = subplot(row,col,4);
subplotChangeSize(h,.009,.009);
imagesc(O4);
axis image
axis off

for croi = 1:(nroi/2)
    %% plot data
    
    [~,roi_name] = pa_GetAtlas(atlas,[],croi);
    roi_name     = regexprep(roi_name.name,'_[L,R]_','_');
    
    matlg1 = simmat(:,:,:,croi,1);matlg1 = reshape(zscore(matlg1(:)),size(matlg1));
    matlg2 = simmat(:,:,:,croi,2);matlg2 = reshape(zscore(matlg2(:)),size(matlg2));
    matrg1 = simmat(:,:,:,croi+nroi/2,1);matrg1 = reshape(zscore(matrg1(:)),size(matrg1));
    matrg2 = simmat(:,:,:,croi+nroi/2,2);matrg2 = reshape(zscore(matrg2(:)),size(matrg2));
    
    h            = subplot(row,col,5);
    subplotChangeSize(h,.009,.009);
    imagesc(mean( matlg1,3),[d u]);
    colorbar
    title(roi_name{1},'interpreter','none','fontsize',fontsize);
    axis image
    axis off
    h            = subplot(row,col,6);
    subplotChangeSize(h,.009,.009);
    imagesc(mean( matrg1,3),[d u]);
    colorbar
    title('grp 1 R');
    axis image
    axis off
    h            = subplot(row,col,7);
    subplotChangeSize(h,.009,.009);
    imagesc(mean( matlg2,3),[d u]);
    colorbar
    title('grp 2 L');
    axis image
    axis off
    h            = subplot(row,col,8);
    subplotChangeSize(h,.009,.009);
    imagesc(mean( matrg2,3),[d u]);
    colorbar
    title('grp 2 R');
    axis image
    axis off
    
    %% plot betas
    %betamat(:,croi,gr).b(1)

        Bp1=[];
        Bp2=[];
        Bp3=[];
        Bp4=[];
        Bp5=[];
        Bp6=[];
        Bp7=[];
        Bp8=[];
        Bk1=[];
        Bk2=[];
        Bk3=[];
        Bk4=[];
        Bk5=[];
        Bk6=[];
        Bk7=[];
        Bk8=[];        
        for csub = 1:nsub
            Bp1=cat(1,Bp1,betamat(csub,croi,1).b(1));
            Bp2=cat(1,Bp2,betamat(csub,croi,1).b(2));
            Bp3=cat(1,Bp3,betamat(csub,croi,1).b(3));
            Bp4=cat(1,Bp4,betamat(csub,croi,1).b(4));
            Bp5=cat(1,Bp5,betamat(csub,croi+nroi/2,1).b(1));
            Bp6=cat(1,Bp6,betamat(csub,croi+nroi/2,1).b(2));
            Bp7=cat(1,Bp7,betamat(csub,croi+nroi/2,1).b(3));
            Bp8=cat(1,Bp8,betamat(csub,croi+nroi/2,1).b(4));
            Bk1=cat(1,Bk1,betamat(csub,croi,2).b(1));
            Bk2=cat(1,Bk2,betamat(csub,croi,2).b(2));
            Bk3=cat(1,Bk3,betamat(csub,croi,2).b(3));
            Bk4=cat(1,Bk4,betamat(csub,croi,2).b(4));
            Bk5=cat(1,Bk5,betamat(csub,croi+nroi/2,2).b(1));
            Bk6=cat(1,Bk6,betamat(csub,croi+nroi/2,2).b(2));
            Bk7=cat(1,Bk7,betamat(csub,croi+nroi/2,2).b(3));
            Bk8=cat(1,Bk8,betamat(csub,croi+nroi/2,2).b(4));
        end
        
        h            = subplot(row,col,9);
        bar([mean(Bp1) mean(Bp2) mean( Bp3) mean( Bp4)]);
        title('betas');
        
        h            = subplot(row,col,10);
        bar([ mean(Bp5) mean( Bp6) mean( Bp7) mean( Bp8)]);
        
        h            = subplot(row,col,11);
        bar([mean(Bk1) mean(Bk2) mean( Bk3) mean( Bk4)]);
        title('betas');
        
        h            = subplot(row,col,12);
        bar([ mean(Bk5) mean( Bk6) mean( Bk7) mean( Bk8)]);
    
    %% test for significance
    
    [h1T1L,p1T1L]=ttest(Bp1);
    [h1T2L,p1T2L]=ttest(Bk1);
    [h1TBL,p1TBL]=ttest2(Bp1,Bk1);
    [h2T1L,p2T1L]=ttest(Bp2);
    [h2T2L,p2T2L]=ttest(Bk2);
    [h2TBL,p2TBL]=ttest2(Bp2,Bk2);
    [h3T1L,p3T1L]=ttest(Bp3);
    [h3T2L,p3T2L]=ttest(Bk3);
    [h3TBL,p3TBL]=ttest2(Bp3,Bk3);
    [h4T1L,p4T1L]=ttest(Bp4);
    [h4T2L,p4T2L]=ttest(Bk4);
    [h4TBL,p4TBL]=ttest2(Bp4,Bk4);
    
    [h1T1R,p1T1R]=ttest(Bp5);
    [h1T2R,p1T2R]=ttest(Bk5);
    [h1TBR,p1TBR]=ttest2(Bp5,Bk5);
    [h2T1R,p2T1R]=ttest(Bp6);
    [h2T2R,p2T2R]=ttest(Bk6);
    [h2TBR,p2TBR]=ttest2(Bp6,Bk6);
    [h3T1R,p3T1R]=ttest(Bp7);
    [h3T2R,p3T2R]=ttest(Bk7);
    [h3TBR,p3TBR]=ttest2(Bp7,Bk7);
    [h4T1R,p4T1R]=ttest(Bp8);
    [h4T2R,p4T2R]=ttest(Bk8);
    [h4TBR,p4TBR]=ttest2(Bp8,Bk8);
    
    %% plot ttests
    
    h            = subplot(row,col,13);
    bar([p1T1L p2T1L p3T1L p4T1L]);
    line([0 6],[0.05 0.05]);
    ylim([0 0.1]);
    xlim([0 5]);
    title('ttest grp1 L');
    h            = subplot(row,col,14);
    bar([p1T1R p2T1R p3T1R p4T1R]);
    line([0 6],[0.05 0.05]);
    ylim([0 0.1]);
    xlim([0 5]);
    title('ttest grp1 R');
    h            = subplot(row,col,15);
    bar([p1T2L p2T2L p3T2L p4T2L]);
    line([0 6],[0.05 0.05]);
    ylim([0 0.1]);
    xlim([0 5]);
    title('ttest grp2 L');
    h            = subplot(row,col,16);
    bar([p1T2R p2T2R p3T2R p4T2R]);
    line([0 6],[0.05 0.05]);
    ylim([0 0.1]);
    xlim([0 5]);
    title('ttest grp2 R');
    
    %% plot r2 and group comparison
    
    h            = subplot(row,col,17);
    R1=mean(rsq(:,croi,1));
    R2=mean(rsq(:,croi,1));
    R3=mean(rsq(:,croi+nroi/2,2));
    R4=mean(rsq(:,croi+nroi/2,2));
    bar([R1 R2 R3 R4]);
    title('R2 for each model');
    ylim([0.5 1]);
    
    h            = subplot(row,col,18);
    bar([p1TBL p2TBL p3TBL p4TBL]);
    line([0 6],[0.05 0.05]);
    ylim([0 0.1]);
    xlim([0 5]);
    title('ttest grp compare L');
    
    h            = subplot(row,col,19);
    bar([p1TBR p2TBR p3TBR p4TBR]);
    line([0 6],[0.05 0.05]);
    ylim([0 0.1]);
    xlim([0 5]);
    title('ttest grp compare R');
    
    % continue to next roi
    input('Press return key to continue');
    
end
