%%
SPMid      = {'anova_cov_',4,'faces/chrf_derivs_00','swRealigned'};
gs        = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
bs        = setdiff(1:29,gs);
D         = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',1:29,'swRealigned');
Dgs       = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',gs,'swRealigned');
Dbs       = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',bs,'swRealigned');
%%
tt       = {'Combined ROI' 'Insula', 'Frontal Operculum'};
figure(9);
set(gcf,'position',[  103   229   965   678])
colors   = {'k' 'k'}
face_color = [.1 .1 .1]*6;
edge_color = face_color;[.1 .1 .1]+.9;
lr       = 2;
field_c  = 0;
addon    = [0 .1];
fs       = 30;
clf;
cc = 0;
%
for field = {'B2T' 'UCS2T'};
    %    
    field_c = field_c  + 1;
    pval    = 0.05;    
    nroi_c  = 0;
    for nroi    = [3 1 2]
        %
        nroi_c   = nroi_c+1;
        cc       = cc + 1
        hsub(cc) = subplot(2,3,cc);
        %
        X    = deg2rad(linspace(-135,180,8)')+addon(field_c);
        Xp   = deg2rad(linspace(-135,180,8)')+.1;
        data = (D{lr}.(field{1})(:,:,nroi));
        B2B  = D{lr}.B2B(:,nroi);
        M    = mean(data,2);
        S    = prctile(data,[2.5 97.5],2);
        h    = bar(X,M,1);
        SetFearGenBarColors(h)
%         set(h,'edgecolor',edge_color,'facecolor',face_color)
        
        hold on
        a    = errorbar(X,M,S(:,1)-M,S(:,2)-M,'color',colors{field_c},'marker','o','linestyle','none','linewidth',2);
        errorbar_tick(a,Inf);
        Est = FitGauss(X,M,3);
        if Est.pval >= -log10(pval)            
            plot(Est.xsup,Est.fitup,'linewidth',3,'color',colors{field_c})
        else
            plot(Est.xsup,repmat(mean(M),1,length(Est.xsup)),'-','linewidth',3,'color',colors{field_c})
        end
        %        
        if field_c == 1
            title(tt{nroi_c} ,'interpreter','none','fontsize',fs);
        end
        axis tight
        xlim([X(1)-.5,X(end)+.5])
        ylim([0,1])
        SetTickNumber(gca,3,'y');
%         ylim([-.67 1]);
        
        set(gca,'color','none','xtick',[X(4) X(end)],'xticklabel',{'CS+' 'CS-'},'ytick',[0 1]);
        axis square;
        set(gca,'fontsize',fs,'linewidth',3)
        box off;
        %
        if field_c == 1
            plot(xlim,repmat(mean(B2B),1,2),'--','color','k','linewidth',1.5);
            if cc == 1
                ylabel('r')
            end
        end
        %  
        if cc == 2 | cc == 3 | cc == 5 | cc == 6
            set(gca,'yticklabel',[]);
        end
         if cc <= 6 & cc > 3
             set(gca,'ylim',[-.35 1]);
             if cc == 4
                ylabel('r')
             end
         end        
        subplotChangeSize(hsub(cc),0.1/3,0.1/3);
%         pause;
%         lr   = 2;
%         data = (D{lr}.(field{1})(:,:,nroi));
%         B2B  = D{lr}.B2B(:,nroi);
%         M    = mean(data,2);
%         S    = prctile(data,[2.5 97.5],2);
%         % bar(X,M,.5,'r');
%         hold on
%         a = errorbar(Xp,M,S(:,1)-M,S(:,2)-M,'b','marker','o','linestyle','none');
%         errorbar_tick(a,Inf);
%         Est = FitGauss(X,M,3);
%         if Est.pval >= -log10(pval)
%             plot(Est.xsup,Est.fitup,'linewidth',3,'color','b')
%         else
%             plot(Est.xsup,mean(M),'linewidth',3,'color','b')
%         end
%         %
%         plot([Xp(1) Xp(end)],repmat(mean(B2B),1,2),'-b');
%         %
%         hold off
        
%         axis tight;
    end    
end