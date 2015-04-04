function PatternAnalysis_ROIShow(mcmethod,simmethod,what,matrix,both)
%PatternAnalysis_ROIShow(mcmethod,matrix,simmethod,both,what)
%
%   MCMETHOD : method for correction
%
%   MATRIX : 0 or 1 depending on whether the result to be plotted is a
%   matrix or vector.
%
%   SIMMETHOD: similarity metric.
%
%   BOTH     : plot B and T simultaneously or take the difference
%
%   WHAT     : what part of the covariance matrix you like to plot (see
%   PatternAnalysis_DataMatrix2Similarity.m for more info)

gs         = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
[col row]  = GetSubplotNumber(length(GetNicedROIs));
x          = deg2rad(linspace(-135,180,8));
xx         = [x x+2*pi];
close all;
ffigure;

roi_c = 0;
for nroi = GetNicedROIs
    roi_c = roi_c + 1;
    subplot(col,row,roi_c);
    %
    [DM roiname]   = PatternAnalysis_BetaFiles2Matrix(nroi,gs,mcmethod);
    DM             = DM([2 1]);%make left hemisphere appear on the left side.
    [SM fun]       = PatternAnalysis_DataMatrix2Similarity(DM,simmethod);    
    ts             = size(SM{1},3);
    %
    if ~matrix
        PlotVector(both,what);
    else
        PlotMatrix;
    end    
    %
    drawnow;
end
    function PlotMatrix
        imagesc([mean(SM{1}(1:16,1:16,:),3) mean(SM{2}(1:16,1:16,:),3) ]);
        subplotChangeSize(gca,.009,.009);
        axis off
        hold on
        plot([16 16]+.5,[1 16]-.5,'k')
        hold off
        title(roiname,'interpreter','none');                
        thincolorbar('vertical')
    end


    function PlotVector(both,what)
        %
        if both == 1%plot B and T in one plot
            %
            MB = [fun.B.(what).mean(1) fun.B.(what).mean(2)];
            MT = [fun.T.(what).mean(1) fun.T.(what).mean(2)];
            SB = [fun.B.(what).sem(1) fun.B.(what).sem(2)];
            ST = [fun.T.(what).sem(1) fun.T.(what).sem(2)];
            %
            plot(xx,MB ,'k','linewidth',2);
            hold on;
            plot(xx,MT,'r','linewidth',2);
            errorbar(xx,MB,SB,'ok');
            errorbar(xx,MT,ST,'or');
            
        elseif both == 0%plot only the difference, D.
            
            
            Mdiff = [fun.D.(what).mean(1) fun.D.(what).mean(2)];
            Sdiff = [fun.D.(what).sem(1) fun.D.(what).sem(2)];
            
            plot(xx,Mdiff,'r','linewidth',2);
            hold on;
            errorbar(xx,Mdiff,Sdiff,'or');
            
        end
        
        axis tight;
        plot(repmat(mean(xx(8:9)),1,2),ylim,'k');
        plot(xlim,[0 0],'k');
        plot(xx([4 4]),ylim,'k:');
        plot(xx([12 12]),ylim,'k:')
        box off;
        title(roiname,'interpreter','none');
        set(gca,'xticklabel','');
        hold off;
        subplotChangeSize(gca,.009,.009);
    end
end