function pa_PlotSimMat(simmat)
%pa_PlotSimMat(simmat)
%
%   Plots similarity matrices for each ROI in a panel where left/right
%   hemispheres and baseline/test phases are shown together. SIMMAT is
%   [condition,condition,subject,roi]. The plotted values do not represent
%   real values as L and R hemispheres are separately scaled for plotting
%   purposes.
%
%
%   Dependency: pa_defaults, 


C         = simmat;%[condition,condition,subject,roi,phase]
%
troi      = size(simmat,4);
[row,col] = GetSubplotNumber(troi/2);
%
tcond     = size(C,1);
%
[d u]     = GetColorMapLimits(C(:),3);
fontsize  = 8;
ffigure;
for nroi = 1:troi/2
    %
    [~,roi_name] = pa_GetAtlas(pa_defaults('atlas'),[],nroi);
    roi_name     = regexprep(roi_name.name,'_[L,R]_','_');
    %
    h            = subplot(row,col,nroi);
    subplotChangeSize(h,.009,.009);
    %
    mat1 = C(:,:,:,nroi);mat1 = reshape(zscore(mat1(:)),size(mat1));
    mat2 = C(:,:,:,nroi+troi/2);mat2 = reshape(zscore(mat2(:)),size(mat2));
    imagesc(mean(cat(2,mat1,mat2),3));
    axis image    
    axis off
    hold on    
    plot([tcond tcond]+.5,[0 tcond]+.5,'k')                    
    title(roi_name{1},'interpreter','none','fontsize',fontsize);
    hcb     = thincolorbar('vertical');
    set(hcb,'Fontsize',fontsize);    
end
% SaveFigure(sprintf('%s.png',fname),'-nocrop')