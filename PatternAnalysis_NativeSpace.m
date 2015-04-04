clear all
load /Volumes/feargen2/feargen2/data/midlevel/selectedsubjects.mat
threshold = 50;
Dc        = [];
%
sub_c  = 0;
for nsubject = subject_list;
    sub_c = sub_c + 1;
    %% load the beta images
    c = 0;
    for nc = [1:10 47:56]
        c = c+1;
        beta(c) = spm_vol(sprintf('/Volumes/feargen2/feargen2/data/sub%03d/phase04/mrt/spm/firstlevel/chrf_derivs_00/Realigned/beta_%04d.img',nsubject,nc));
    end
    beta_vol = spm_read_vols(beta);
    s        = size(beta_vol);
    beta_vol = reshape(beta_vol,prod(s(1:3)),s(4));
    %% load the time courses
    load(sprintf('/Volumes/feargen2/feargen2/data/sub%03d/phase04/mrt/spm/firstlevel/chrf_derivs_00/Realigned/SPM.mat',nsubject))
    YY        = spm_read_vols(spm_vol(SPM.xY.VY));
    s         = size(YY);
    YY        = reshape(YY,prod(s(1:3)),s(4));
    %%
    h        = Atlas2NativeSpace(nsubject);
    %%    
    for nroi      = 41;
        [nsubject nroi]               
        ROIs     = spm_read_vols(h(nroi));
        ROIs     = ROIs > threshold;        
        %%
        D        = beta_vol(ROIs(:),:);
        %%
        Y        = YY(ROIs(:),:)';
        KWY      = spm_filter(SPM.xX.K,SPM.xX.W*Y);
        res      = spm_sp('r',SPM.xX.xKXs,KWY);
        Dpw      = cov(res)^-.5*D;
        %%
        Dc(:,:,sub_c,  nroi)  = nancov(D);
        Dpwc(:,:,sub_c,nroi)  = nancov(Dpw);
        %         figure(1);
        %         imagesc(squeeze(mean(Dc(:,:,:,nroi),3)));
        %         drawnow;
    end
end
save('~/Desktop/Native.mat','Dc','Dpwc');
% % %%
% % for n = 1:45
% %     subplot(8,6,n)
% %     imagesc([squeeze(mean(Dc(:,:,:,n),3)) squeeze(mean(Dc(:,:,:,n+45),3))] )
% %     [~,roi]=LoadFeargenAtlas(50,n);
% %     title(roi.name,'interpreter','none')
% % end


