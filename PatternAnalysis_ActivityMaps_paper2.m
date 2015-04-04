function PatternAnalysis_ActivityMaps(SPMidcell,conditions,subjects,rois,cmap_limit,method,zsc,id,thr)
%PatternAnalysis_ActivityMaps(SPMidcell,thr,subjects,rois,cmap_limit,method,zsc,id)
%
%   Will plot the evoked activity maps. ZSC is a flag for zscoring
%   indivudal patterns.
%
%   As input an SPM structure is entered, therefore the relevant analysis
%   has to be carried out already. From the SPM file the first-level data
%   will be conviniently extracted.
%
%   ID is a string that is appendend to the end of the filename.
%
%CMAP is about 3.5, the colormap boundary, number of std from mean.
%
%TYPE is the common pattern correction method.
%
%Example:
%gs           = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
%rois         = [[2 10 15 16 17 23 26 29 30 31 35 36 41 42 43];
%PatternAnalysis_ActivityMaps({'anova_cov_',4,'chrf_derivs_00_fullcovariate','swRealigned'},[1:10;11:20]',1:29,rois,2.65,'',1,'all');
%PatternAnalysis_ActivityMaps({'anova_cov_',4,'chrf_derivs_00_fullcovariate','swRealigned'},[1:10;11:20]',gs,rois,2.65,'',1,'gs');
%
%   to run this anew you must have already the
%   PatternAnalysis_SPM2ExplainedVariance and
%   PatternAnalysis_BetaFiles2Matrix functions ran.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_folder  = sprintf('~/Pictures/results_rsa/%s/%s/phase%02d/%s/%s/',mfilename,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4});
if exist(save_folder) == 0;mkdir(save_folder);end
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
col        = size(conditions,1);
row        = size(conditions,2)*3;
tpatt      = size(conditions,2);
% close all;
%Load the data, DM is a cell with left right hemispheres.
%PatternAnalysis_BetaFiles2Matrix('load',50,'anova_cov_',4,'chrf_derivs_00_fullcovariate','swRealigned','',1,1:29)
[DM1,~,~,W]      = PatternAnalysis_BetaFiles2Matrix('load',thr,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4},method,rois(1),1:29);
[DM2,~,~,W]      = PatternAnalysis_BetaFiles2Matrix('load',thr,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4},method,rois(2),1:29);
for lr = 1:2
    ffigure;
    clf
    [roi1, mask1]         = LoadFeargenAtlas(thr,rois(1)+(45*(2-lr)));
    [roi2, mask2]         = LoadFeargenAtlas(thr,rois(2)+(45*(2-lr)));
    mask_i1            = find(mask1.d);
    mask_i2            = find(mask2.d);
    
    %Load the data and insert it in a fake brain...
    brain             = zeros(91*109*91,size(DM1{lr},2),length(subjects));
    for ns = subjects
        for nc = 1:size(DM{lr},2)
            %                 brain(mask_i,nc,ns) = diag(1./sqrt(W{lr}(:,1,ns)))*DM{lr}(:,nc,ns);
            brain(mask_i1,nc,ns) = DM1{lr}(:,nc,ns);
            brain(mask_i2,nc,ns) = DM2{lr}(:,nc,ns);
        end
    end
    %%
    phase = 0;
    for conds = conditions
        phase         = phase +1;
        mean_brain    = reshape(mean(brain(:,conds ,:),3),[91 109 91 col]);
        amap{phase}   = Analysis_FMRI_ActivityMaps(mask,mean_brain,cmap_limit);
    end
    %%
    fields_alpha = {'x_alpha' 'y_alpha' 'z_alpha' };
    fields_cmap  = {'xcmap' 'ycmap' 'zcmap' };
    fields_data  = {'x' 'y' 'z'};
    for n = 1:col;
        for nrow = 0:row-1
            
            nfields = mod(nrow,3)+1;%will be 1,2,3 for each next row
            npatt   = floor(nrow/3)+1;%will be 1, 1, 1, 2, 2, 2
            N       = n+col*nrow;%will be 1,6,11,.. if col is 5
            %
            subplot(row,col,N);
            imagesc(amap{npatt}.(fields_data{nfields})(:,:,n), amap{npatt}.(fields_cmap{nfields})(n,:) );
            
            alpha(double(amap{npatt}.(fields_alpha{nfields})));
            Bellish;
            drawnow;
        end
    end
end 
    %
    supertitle( sprintf('%s #vx: %d, %s',mask.name{1}, sum(mask.d(:)),method)    ,1,'interpreter','none','fontsize',20);
    %         SaveFigure(sprintf('%s%s_Maps_zs%d_%s_%03d_%s.png',save_folder,mask.name{1},zsc,method,thr,id),'-r120','-nocrop');
end


    function [roi]=Analysis_FMRI_ActivityMaps(mask,data_mat,sigma_cmap)
        
        %needs [Y,X,Z,Cond] volume information in DATA_MAT
        
        tcond = size(data_mat,4);
        %
        x = [min(find(sum(squeeze(sum(mask.d,3)),2) ~= 0)) max(find(sum(squeeze(sum(mask.d,3)),2) ~= 0))];
        y = [min(find(sum(squeeze(sum(mask.d,3))) ~= 0)) max(find(sum(squeeze(sum(mask.d,3))) ~= 0))];
        z = [min(find(sum(squeeze(sum(mask.d))) ~= 0)) max(find(sum(squeeze(sum(mask.d))) ~= 0))];
        %mask the data
        data_mat            = data_mat.*repmat(mask.d,[1 1 1 tcond]);
        
        %vectorize it so that we can get the contribution of all voxels to the
        %colormap computation
        
        %for the first COL conditions
        s     = size(data_mat(:,:,:,1:col));
        dummy = reshape(data_mat(:,:,:,1:col),prod(s(1:3)),s(4));
        dummy = mean(dummy(mask.d(:),:),2);
        [roi.cmap(1:col,1) roi.cmap(1:col,2)]  = GetColorMapLimits(dummy(:),sigma_cmap);
        
        
        for n = 1:tcond
            current          = data_mat(x(1):x(2),y(1):y(2),z(1):z(2),n);
            if zsc
                SSS     = size(current);
                current = reshape(nanzscore(current(:)),SSS);
            end
            current_mask     = mask.d(x(1):x(2),y(1):y(2),z(1):z(2));
            
            % get the data
            roi.x(:,:,n)     = squeeze(nanmean(current,1));
            roi.y(:,:,n)     = squeeze(nanmean(current,2));
            roi.z(:,:,n)     = squeeze(nanmean(current,3));
            
            % get the alpha masks
            roi.x_alpha      = squeeze(mean(double(current_mask),1) ~=  0);
            roi.y_alpha      = squeeze(mean(double(current_mask),2) ~=  0);
            roi.z_alpha      = squeeze(mean(double(current_mask),3) ~=  0);
            %
        end
        
        [roi.xcmap(1:col-2,1) roi.xcmap(1:col-2,2)]   = GetColorMapLimits(Vectorize(roi.x(:,:,1:col-2)),sigma_cmap);
        [roi.ycmap(1:col-2,1) roi.ycmap(1:col-2,2)]   = GetColorMapLimits(Vectorize(roi.y(:,:,1:col-2)),sigma_cmap);
        [roi.zcmap(1:col-2,1) roi.zcmap(1:col-2,2)]   = GetColorMapLimits(Vectorize(roi.z(:,:,1:col-2)),sigma_cmap);
        
        [roi.xcmap(col-1:col,1) roi.xcmap(col-1:col,2)]   = GetColorMapLimits(Vectorize(roi.x(:,:,col-1:col)),sigma_cmap);
        [roi.ycmap(col-1:col,1) roi.ycmap(col-1:col,2)]   = GetColorMapLimits(Vectorize(roi.y(:,:,col-1:col)),sigma_cmap);
        [roi.zcmap(col-1:col,1) roi.zcmap(col-1:col,2)]   = GetColorMapLimits(Vectorize(roi.z(:,:,col-1:col)),sigma_cmap);
    end

    function Bellish
        
        %         axis off;
        box off
        if N <= col
            title(sprintf('C:%02d',n));
        end
        thincolorbar('horizontal');
        color = repmat((npatt/tpatt)*.5+.5,1,3);
        set(gca,'color',color,'xtick',[],'ytick',[],'linewidth',1,'xcolor',color,'ycolor',color);
        axis image;
    end
end
