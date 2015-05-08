% load beta file

for gr  = 1:2;
    for roi = 1:96;
        for csub = 1:nsub
           
            
            if 
            disp(strcat('group :',gr));
            disp(betas.roiname{roi});
            disp('is significantly different from zero');
            
            %[cond,cond,sub,group,roi,laterality,metric]
            Y  = mean(sim.mc.cov(:,:,csub,roi,gr),3);     
            Ycollaps = reshape(accumarray(lut(:),Y(:)),[4 4]);
            %%
            DM                 = [X1(:) X2(:) Sess(:) D(:) K(:)];%design matrix
            betas(:,csub,roi,gr) = DM\Ycollaps(:);
        end
    end
end