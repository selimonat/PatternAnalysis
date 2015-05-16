%% Go 
clear all
%results_path='/projects/cond/midlevel/pa_spmbeta2cov_45_158ef17';
% results_path='/projects/cond/midlevel/pa_spmbeta2cov_50_18b3013';
% samba_path='/home/schenk/Documents/MATLAB/RSA';
results_path='\\samba\schenk\Documents\MATLAB\RSA';

% load simmat
% load(sprintf('%s/simmatcue.mat',results_path));
% load(sprintf('%s/simmatpain.mat',results_path));
% Bigfield=1;
% load(sprintf('%s/simmatcue1regr.mat',results_path));
% load(sprintf('%s/simmatpain1regr.mat',results_path));
load(sprintf('%s/simmatpain1regr.mat',results_path));
Bigfield=0;
%% define models
%collapse
lut      = kron(reshape([1:16],[4 4])',[ones(9)]);
nsub=size((sim.mc.cov),3);
% nsub=size((sim.raw.cov),3);
%X matrix with diagonals;
% X1        = [1 0 .5 .5;0 1 .5 .5;.5 .5 1 1;0.5 .5 1 1];
% X2        = [1 0 1 0;0 1 0 1   ; 1 0 1 0; 0 1 0 1];
% Sess      = [1 1 0 0; 1 1 0 0  ; 0 0 1 1 ; 0 0 1 1 ];%sess
% D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
%%X matrix without diagonals;
% X1        = [0 0 .5 .5;0 0 .5 .5;.5 .5 0 1;0.5 .5 1 0];
% X2        = [0 0 1 0;0 0 0 1   ; 1 0 0 0; 0 1 0 0];
% Sess      = [0 1 0 0; 1 0 0 0  ; 0 0 0 1 ; 0 0 1 0 ];%sess
% D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
V1=[40 80 60 60]; % temp
V2=[1 0 1 0]; %placcontext
V3=[0 1 0 1]; %kontrollcontext
V4=[1 1 -1 -1]; %session

O1=pa_outerproduct(V1);
O2=pa_outerproduct(V2);
O3=pa_outerproduct(V3);
O4=pa_outerproduct(V4);

D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
K= ones(4,4); % konstant
nreg=6; %5 %6

%% fit models
% betas.mc.cov = nan(nreg,24,124,2);
% betas.raw.cov = nan(nreg,24,124,2);
bfit=zeros(1,16);
bres=zeros(1,16);
for gr  = 1:2;
    for roi = 1:124;
        for csub = 1:nsub
           
            %[cond,cond,sub,group,roi,laterality,metric]
            Ymc  = mean(sim.mc.cov(:,:,csub,roi,gr),3);  
            Yraw  = mean(sim.raw.cov(:,:,csub,roi,gr),3);
            if Bigfield ==1
            Ycollapsmc = reshape(accumarray(lut(:),Ymc(:)),[4 4]);
            Ycollapsraw = reshape(accumarray(lut(:),Yraw(:)),[4 4]);
            else
                Ycollapsmc=Ymc;
                Ycollapsraw=Yraw;
            end
            %
%             DM = [X1(:) X2(:) Sess(:) D(:) K(:)];%design matrix
DM = [O1(:) O2(:) O3(:) O4(:) D(:) K(:)];%design matrix
% DM = [O4(:) D(:) K(:)];%design matrix
            [betas.mc.cov(csub,roi,gr).b,betas.mc.cov(csub,roi,gr).dev,betas.mc.cov(csub,roi,gr).stats] = glmfit(DM,Ycollapsmc(:));
            [betas.raw.cov(csub,roi,gr).b,betas.mc.cov(csub,roi,gr).dev,betas.raw.cov(csub,roi,gr).stats] = glmfit(DM,Ycollapsraw(:));
            
            % the plotroiwise function still uses this data, but can be removed
            % and .stats can be used once glmfit works
            bfit(:)=DM*betas.mc.cov(csub,roi,gr).b;
            bres(:)=Ycollapsmc(:)-bfit(:);
            SSresid = sum(sum(bres.^2));
            SStotal = (length(Ycollapsmc(:))-1) * var(Ycollapsmc(:));
            rsq(csub,roi,gr) = 1 - SSresid/SStotal;
        end
        disp(sim.roiname{roi});
        disp(mean(rsq(:,roi,gr)));
    end
end
disp('residuals');
disp(mean(rsq(:,[1:108 110:124],:)));
betas.roiname=sim.roiname;
% save(sprintf('%s/betascue.mat',results_path),'betas');
% save(sprintf('%s/betascue.mat',samba_path),'betas');
% save(sprintf('%s/betaspain.mat',results_path),'betas');
% save(sprintf('%s/betaspain.mat',samba_path),'betas');
% save(sprintf('%s/betascue1regr.mat',results_path),'betas');
% save(sprintf('%s/betascue1regr.mat',samba_path),'betas');
save(sprintf('%s/betaspain1regrneu1.mat',results_path),'betas', 'rsq');
% save(sprintf('%s/betaspain1regrneu1.mat',samba_path),'betas', 'rsq');

