%% Go 
gr = 2;
Y  = mean(sim.d.mc(:,:,:,gr,roi,lr,2),3);
%% define models
%collapse
lut      = kron(reshape([1:16],[4 4])',[ones(9)]);
Ycollaps = reshape(accumarray(lut(:),Y(:)),[4 4]);
%X matrix with diagonals;
X1        = [1 0 .5 .5;0 1 .5 .5;.5 .5 1 1;0.5 .5 1 1];
X2        = [1 0 1 0;0 1 0 1   ; 1 0 1 0; 0 1 0 1];
Sess      = [1 1 0 0; 1 1 0 0  ; 0 0 1 1 ; 0 0 1 1 ];%sess
D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
%%X matrix without diagonals;
X1        = [0 0 .5 .5;0 0 .5 .5;.5 .5 0 1;0.5 .5 1 0];
X2        = [0 0 1 0;0 0 0 1   ; 1 0 0 0; 0 1 0 0];
Sess      = [0 1 0 0; 1 0 0 0  ; 0 0 0 1 ; 0 0 1 0 ];%sess
D         = [1 0 0 0; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];%self similarity
%% fit models
betas = nan(5,2,48,2);
for gr  = 1:2;
    for roi = 1:48;
        for lr  = 1:2;            
            %[cond,cond,sub,group,roi,laterality,metric]
            Y        = mean(sim.d.mc(:,:,:,gr,roi,lr,2),3);            
            %%
            DM                 = [X1(:) X2(:) Sess(:) D(:) ones(length(D(:)),1)];%design matrix
            betas(:,gr,roi,lr) = DM\Ycollaps(:);
        end
    end
end