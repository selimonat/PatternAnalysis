function pa_PlotSimMat(pattern,metric)
%pa_PlotSimMat(pattern,metric)
%
%   Plots similarity matrices for each ROI in a panel where left/right
%   hemispheres and baseline/test phases are shown together.


C         = pa_GetSimMat(pattern,metric);
C         = C.(metric{1});%[condition,condition,subject,roi,phase]
%
troi      = pa_defaults('troi');
[row,col] = GetSubplotNumber(troi/2);
%

%
[d u]     = GetColorMapLimits(C(:),3);
fontsize  = 8;
ffigure;
for nroi = 1:troi/2
    %
    [~,roi_name] = pa_GetAtlas('FeargenMerged',[],nroi);
    roi_name     = regexprep(roi_name.name,'_[L,R]_','_');
    %
    h         = subplot(row,col,nroi);
    subplotChangeSize(h,.009,.009);
    %
    imagesc(mean(cat(2,C(:,:,:,nroi),C(:,:,:,nroi+troi/2)),3));
    axis image    
    axis off
    hold on    
    plot([16 16]+.5,[0.5 16.5],'k')                    
    title(roi_name{1},'interpreter','none','fontsize',fontsize);
    hcb     = thincolorbar('vertical');
    set(hcb,'Fontsize',fontsize);    
end
% SaveFigure(sprintf('%s.png',fname),'-nocrop')