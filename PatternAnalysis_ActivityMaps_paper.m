function [amap]=PatternAnalysis_ActivityMaps_paper(SPMidcell,conditions,subjects,rois,cmap_limit,method,zsc,id,thr)
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
%
%PatternAnalysis_ActivityMaps_paper({'anova_cov_',4,'faces/chrf_derivs_00','swRealigned'},[1:10;11:20]',gs,33,3,'',1,'all',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_folder  = sprintf('/Users/onat/Pictures/results_rsa/%s/%s/phase%02d/%s/%s/',mfilename,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4});
if exist(save_folder) == 0;mkdir(save_folder);end
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
colnames   = {'-3' '-2' '-1' 'CS+' '+1' '+2' '+3' 'CS-' '' 'UCS'};
col        = size(conditions,1);
row        = size(conditions,2)*3;
tpatt      = size(conditions,2);
% close all;

for nroi = rois
    ffigure;
    %Load the data, DM is a cell with left right hemispheres.
             %PatternAnalysis_BetaFiles2Matrix('load',50,'anova_cov_',4,'chrf_derivs_00_fullcovariate','swRealigned','',1,1:29)
    [DM,~,~,W]      = PatternAnalysis_BetaFiles2Matrix('load',thr,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4},method,nroi,1:29);
    for lr = [2]
        clf
         [roi, mask]      = LoadFeargenAtlas(thr,nroi+(45*(2-lr)));
         mask_i           = find(mask.d);        
        %Load the data and insert it in a fake brain...        
        brain             = zeros(91*109*91,size(DM{lr},2),length(subjects));
        for ns = subjects
            for nc = 1:size(DM{lr},2)
                
%                 brain(mask_i,nc,ns) = diag(1./sqrt(W{lr}(:,1,ns)))*DM{lr}(:,nc,ns);
                  brain(mask_i,nc,ns) = DM{lr}(:,nc,ns);
                                    
            end
        end
        %%        
        phase = 0;
        for conds = conditions
            phase         = phase +1;
            mean_brain    = reshape(mean(brain(:,conds ,:),3),[91 109 91 col]);
%             %mark the voxels that are peak values.
            mean_brain(24,73,35,:) = NaN;
            mean_brain(29,76,39,:) = NaN;
            %                        
            amap{phase}         = Analysis_FMRI_ActivityMaps(mask,mean_brain,cmap_limit);
        end        
        %        
        amap{2}.x(:,:,end)  =  amap{1}.x(:,:,end);
        amap{2}.y(:,:,end)  =  amap{1}.y(:,:,end);
        amap{2}.z(:,:,end)  =  amap{1}.z(:,:,end);
        %
        amap{1}.x(:,:,end)  =  mean(amap{1}.x(:,:,1:8),3);
        amap{1}.y(:,:,end)  =  mean(amap{1}.y(:,:,1:8),3);
        amap{1}.z(:,:,end)  =  mean(amap{1}.z(:,:,1:8),3);
        %%
        %%get the colormaps limits, the above computations separates the
        %%phases which I don't want
        conds         = [Vectorize(conditions(1:8,:)) ;9];
        col           = length(conds);
        mean_brain    = reshape(mean(brain(:,conds,:),3),[91 109 91 col]);
        amap2         = Analysis_FMRI_ActivityMaps(mask,mean_brain,cmap_limit);
        amap2.xcmap(10,:)=amap2.xcmap(end,:);
        amap2.ycmap(10,:)=amap2.ycmap(end,:);
        amap2.zcmap(10,:)=amap2.zcmap(end,:);
        
        fields_alpha  = {'x_alpha' 'y_alpha' 'z_alpha' };
        fields_cmap   = {'xcmap' 'ycmap' 'zcmap' };
        fields_data   = {'x' 'y' 'z'};
        col           = size(conditions,1);
        for n = 1:col;
            for nrow = 0:row-1
                if n ~= (col-1)%skip one column
                    
                    nfields = mod(nrow,3)+1;%will be 1,2,3 for each next row
                    npatt   = floor(nrow/3)+1;%will be 1, 1, 1, 2, 2, 2
                    N       = n+col*nrow;%will be 1,6,11,.. if col is 5
                    %
                    sph(nrow+1,n) = subplot(row,col,N);
                    
                    %if ~((n == 10) && (nrow > 2));%exclude the UCS of the test phase
                        if ismember(nrow,[0 3 1 4])
                            imagesc(rot90(fliplr(amap{npatt}.(fields_data{nfields})(:,:,n)),-1), amap2.(fields_cmap{nfields})(n,:));
                            alpha(rot90(fliplr(double(amap{npatt}.(fields_alpha{nfields}))),-1));
                        else
                            imagesc(rot90((amap{npatt}.(fields_data{nfields})(:,:,n)),2), amap2.(fields_cmap{nfields})(n,:));
                            alpha(rot90((double(amap{npatt}.(fields_alpha{nfields}))),2));
                        end
                    %end                                        
                    Bellish;
                    drawnow;
                end
            end
        end        
        %
        supertitle( sprintf('%s #vx: %d, %s',mask.name{1}, sum(mask.d(:)),method)    ,1,'interpreter','none','fontsize',20);
        SaveFigure(sprintf('%s%s_Maps_zs%d_%s_%03d_%s.png',save_folder,mask.name{1},zsc,method,thr,id));
    end
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
% % %         s     = size(data_mat(:,:,:,1:col));
% % %         dummy = reshape(data_mat(:,:,:,1:col),prod(s(1:3)),s(4));
% % %         dummy = mean(dummy(mask.d(:),:),2);
% % %         [roi.cmap(1:col,1) roi.cmap(1:col,2)]  = GetColorMapLimits(dummy(:),sigma_cmap);

        data = [];
        for n = 1:tcond
            current          = data_mat(x(1):x(2),y(1):y(2),z(1):z(2),n);
            if zsc
                SSS     = size(current);
                current = reshape(nanzscore(current(:)),SSS);
            end
            current_mask     = mask.d(x(1):x(2),y(1):y(2),z(1):z(2));
            %Nanize outside mask voxels
%             current(~current_mask)= NaN;
            % get the data
            roi.x(:,:,n)     = squeeze(mean(current,1));
            roi.y(:,:,n)     = squeeze(mean(current,2));
            roi.z(:,:,n)     = squeeze(mean(current,3));
            
            % get the alpha masks
            roi.x_alpha      = squeeze(mean(double(current_mask),1) ~=  0); 
            roi.y_alpha      = squeeze(mean(double(current_mask),2) ~=  0);
            roi.z_alpha      = squeeze(mean(double(current_mask),3) ~=  0);            
            %            
        end
        
        dummy = [reshape(roi.x,size(roi.x,1)*size(roi.x,2),size(roi.x,3))];
        data  = dummy(roi.x_alpha,:);            
        dummy = reshape(roi.y,size(roi.y,1)*size(roi.y,2),size(roi.y,3));
        data  = [data ; dummy(roi.y_alpha,:)];        
        dummy = reshape(roi.z,size(roi.z,1)*size(roi.z,2),size(roi.z,3));
        data  = [data ; dummy(roi.z_alpha,:)];
        
        
        
% %         [roi.xcmap(1:tcond-1,1) roi.xcmap(1:tcond-1,2)]   = GetColorMapLimits(Vectorize(data(:,1:tcond-1)),sigma_cmap);
% %         [roi.ycmap(1:tcond-1,1) roi.ycmap(1:tcond-1,2)]   = GetColorMapLimits(Vectorize(data(:,1:tcond-1)),sigma_cmap);
% %         [roi.zcmap(1:tcond-1,1) roi.zcmap(1:tcond-1,2)]   = GetColorMapLimits(Vectorize(data(:,1:tcond-1)),sigma_cmap);
% %         
% %         [roi.xcmap(tcond,1) roi.xcmap(tcond,2)]   = GetColorMapLimits(Vectorize(data(:,tcond)),sigma_cmap);
% %         [roi.ycmap(tcond,1) roi.ycmap(tcond,2)]   = GetColorMapLimits(Vectorize(data(:,tcond)),sigma_cmap);
% %         [roi.zcmap(tcond,1) roi.zcmap(tcond,2)]   = GetColorMapLimits(Vectorize(data(:,tcond)),sigma_cmap);
        %the colormap adapted by hand to be perceptually best
        cmap = [-1.5 2];%operculum
        cmap = [-.5 2];
        roi.xcmap(1:tcond-1,:) = repmat(cmap,tcond-1,1);
        roi.ycmap(1:tcond-1,:) = repmat(cmap,tcond-1,1);
        roi.zcmap(1:tcond-1,:) = repmat(cmap,tcond-1,1);
        
        cmap = [-1 3];%operculum
        cmap = [-.5 2];
        roi.xcmap(tcond,:) = repmat(cmap,1,1);
        roi.ycmap(tcond,:) = repmat(cmap,1,1);
        roi.zcmap(tcond,:) = repmat(cmap,1,1);
        %
        roi.x_alpha      = Scale(log10(1+squeeze(sum(current_mask,1))));%squeeze(mean(double(current_mask),1) ~=  0); 
        roi.y_alpha      = Scale(log10(1+squeeze(sum(current_mask,2))));%squeeze(mean(double(current_mask),2) ~=  0);
        roi.z_alpha      = Scale(log10(1+squeeze(sum(current_mask,3))));%squeeze(mean(double(current_mask),3) ~=  0);            
        
    end

    function Bellish
            
%         axis off;
        box off
        if N <= col
%             title(colnames{n},'interpreter','none','FontSize',26);
        end        
        color = repmat((npatt/tpatt)*.5+.5,1,3);
%         set(gca,'color',color,'xtick',[],'ytick',[],'linewidth',1,'xcolor',color,'ycolor',color);
        axis image;
        axis off;
        subplotChangeSize(gca,.016,.016);
        %
        if ismember(n,[8 10]);%colorbar only for condition 8 and UCS
            if nrow == 0%colorbar on the first row only
            ap = get(gca','position');
%             thincolorbar('horizontal');
            if n == 8%put the colorbar etiher to the left or right
                h = colorbar('position',[ap(1)+ap(3)+0.01 ap(2) 0.01 ap(4) ]);
                colormap(flipud(lbmap(256,'RedBlue')));
            set(h,'box','off','TickLength',[0.01 0.01],'FontSize',20,'ycolor','k','FontName','Arial','ytick',[ceil(10*amap2.xcmap(n,1))/10 floor(10*amap2.xcmap(n,2))/10] );
%             else
%                 h = colorbar('position',[ap(1)-0.025 ap(2) 0.01 ap(4) ]);
            end
            
            
%             set(gca,'position',ap);            
            end
        end
    end

end
