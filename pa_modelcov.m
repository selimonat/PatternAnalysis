%% Go 
clear all
results_path='/projects/cond/midlevel/pa_spmbeta2cov_45_158ef17';
% load simmat
load(sprintf('%s/simmatcue.mat',results_path));
% load(sprintf('%s/simmatpain.mat',results_path));
%% define models
%collapse
lut      = kron(reshape([1:16],[4 4])',[ones(9)]);
nsub=size((sim.mc.cov),3);
%X matrix with diagonals;
% X1        = [1 0 .5 .5;0 1 .5 .5;.5 .5 1 1;0.5 .5 1 1];
% X2        = [1 0 1 0;0 1 0 1   ; 1 0 1 0; 0 1 0 1];
% Sess      = [1 1 0 0; 1 1 0 0  ; 0 0 1 1 ; 0 0 1 1 ];%sess
% D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
%%X matrix without diagonals;
X1        = [0 0 .5 .5;0 0 .5 .5;.5 .5 0 1;0.5 .5 1 0];
X2        = [0 0 1 0;0 0 0 1   ; 1 0 0 0; 0 1 0 0];
Sess      = [0 1 0 0; 1 0 0 0  ; 0 0 0 1 ; 0 0 1 0 ];%sess
D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
K= ones(4,4); % konstant
%% fit models
betas.mc.cov = nan(5,24,96,2);
for gr  = 1:2;
    for roi = 1:96;
        for csub = 1:nsub
           
            %[cond,cond,sub,group,roi,laterality,metric]
            Y  = mean(sim.mc.cov(:,:,csub,roi,gr),3);     
            Ycollaps = reshape(accumarray(lut(:),Y(:)),[4 4]);
            %%
            DM                 = [X1(:) X2(:) Sess(:) D(:) K(:)];%design matrix
            betas.mc.cov(:,csub,roi,gr) = DM\Ycollaps(:);
        end
    end
end
betas.roiname=sim.roiname;
save(sprintf('%s/betascue.mat',results_path),'betas');
% save(sprintf('%s/betaspain.mat',results_path),'betas');