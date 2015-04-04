function [r]=PatternAnalysis_EvokedSimilarityAnalysis(SPMidcell,pool,how,gs,subfolder)
%Implemtents single subject OR EVOKED analysis of pattern similarity
%between Test patterns and Average Baseline Pattern.
%
%POOL is a logical about left/right pooling
%
%GS is a list of subjects
%
%FOLDERID is string added to the folder for saving
%
%PARAM contains the fitted parameters
%
%If you use a column vector for subject vector, it will make an evoked
%analysis.
%%
thr = 10;
folder = ...
    sprintf('/Users/onat/Pictures/results_rsa/%s/thr_%02d/%s_%s/',mfilename,thr,how,subfolder);
if exist(folder) ==0;mkdir(folder);end;
%
xbase       = linspace(0,2*pi-2*pi/8,8)-pi/4*3;
%gs         = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
%gs         = 1:29;
if strcmp(how,'evoked')
    gs = gs(:);
else
    gs = gs(:)';
end

%%
roi_counter  = 0;
r            = [];
for nroi = [2 31 33 34]%setdiff(GetNicedROIs('pmod'),44);
    roi_counter = roi_counter + 1;
    %%
    [data,dataw,roiname,resms] = PatternAnalysis_BetaFiles2Matrix('load',thr,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4},'',nroi,1:29);
%     [data,dataw,roiname,resms]= PatternAnalysis_BetaFiles2Matrix('load',thr,SPMidcell{1},SPMidcell{2},SPMidcell{3},SPMidcell{4},'',nroi,1:29);    
    %%    
    %
    ns          = 0;
    x           = [];
    ts          = length(gs);    
    for subjects = gs%
        ns = ns + 1;
        for lr = 1:2
            r{lr}.roi(roi_counter)                 = nroi;
            i = isnan(mean(dataw{lr},3));
            data{lr}(i,:,:) = [];
            dataw{lr}(i,:,:) = [];
            
            if strcmp(how,'induced')
                %
                ODD         = mean(data{lr}(:,10,subjects),3)';
                UCS         = mean(data{lr}(:,9,subjects),3)';
                B           = mean(data{lr}(:,1:8,subjects),3)';
                T           = mean(data{lr}(:,11:18,subjects),3)';
                %weight vector
%                 w           = mean(dataw{lr}(:,1,subjects),3)';
                          w           = ones(size(dataw{lr},1),1)*size(dataw{lr},1)^-1;
                %%
                for nface = 1:8
%                     r{lr}.B2T(nface,ns)     = pdist2(mean(B,1),T(nface,:), @(X,Y) wcorr2(X,Y,w));
                    r{lr}.B2T(nface,ns)     = pdist2(mean(B,1),T(nface,:), @(X,Y) wcorr2(X,Y,w));
                    r{lr}.UCS2T(nface,ns)   = pdist2(UCS,T(nface,:), @(X,Y) wcorr2(X,Y,w));
                end
                %%
                r{lr}.mean(:,ns,:)       = cat(3,mean(B,2) , mean(T,2));
                r{lr}.var(:,ns,:)        = cat(3,var(B,1,2), var(T,1,2));
                
            elseif strcmp(how,'evoked')
                tsub    = length(subjects);
                %weight each pattern by the inverse-sqrt of the resms
%                 for ns = subjects(:)'
%                     for nc = 1:size(data{lr},2)
%                         data{lr}(:,nc,ns) = diag(1./sqrt(resms{lr}(:,1,ns)))*data{lr}(:,nc,ns);
%                     end
%                 end
                %%
                w       = mean(dataw{lr}(:,1,subjects),3)';
%                 w       = ones(size(data{lr},1),1)*size(data{lr},1)^-1;
                for nbs = 1:1000
                    ODD         = mean(zscore(data{lr}(:,10,randsample(subjects,tsub,1))),3)';
                    UCS         = mean(zscore(data{lr}(:,9,randsample(subjects,tsub,1))),3)';
                    B           = mean(zscore(data{lr}(:,1:8,randsample(subjects,tsub,1))),3)';                    
                    T           = mean(zscore(data{lr}(:,11:18,randsample(subjects,tsub,1))),3)';
                    %
                    for nface = 1:8
                        r{lr}.B2T(nface,nbs,roi_counter)     = pdist2(mean(B,1),T(nface,:), @(X,Y) wcorr2(X,Y,w));
                        r{lr}.UCS2T(nface,nbs,roi_counter)   = pdist2(UCS,T(nface,:), @(X,Y) wcorr2(X,Y,w));                              
                    end
                    r{lr}.mean(:,nbs,:)       = cat(3,mean(B,2) , mean(T,2));
                    r{lr}.var(:,nbs,:)        = cat(3,var(B,1,2), var(T,1,2));                    
                end    
                
                %% compute B2B correlation between all possible B pairs                
                C = 0;
                for nface1 = 1:8
                    for nface2 = 1:8
                        if nface1 ~= nface2
                            C = C+1;
                            r{lr}.B2B(C,roi_counter) = pdist2( mean(data{lr}(:,nface1,subjects),3)' , mean(data{lr}(:,nface2,subjects),3)' , @(X,Y) wcorr2(X,Y,w));
                        end
                    end
                end
                %%
            elseif strcmp(how,'evoked_weighted')
                
                keyboard
%                 ODD = diag(1./sqrt(resms{lr}(:,1,ns))*data{lr}(:,10,ns);
                
                
            end
        end
    end
    
    %%
% % %     x  = repmat(xbase(:),1,size(r{lr}.B2T,2));
% % %     %%
% % %     if pool
% % %         dummy       = r;
% % %         r           = [];
% % %         r{1}.B2T    = (dummy{1}.B2T   + dummy{2}.B2T)*.5;
% % %         r{1}.UCS2T  = (dummy{1}.UCS2T + dummy{2}.UCS2T)*.5;
% % %         r{1}.mean   = (dummy{1}.mean  + dummy{2}.mean)*.5;
% % %         r{1}.var    = (dummy{1}.var   + dummy{2}.var)*.5;
% % %     end
% % %     %%
% % %     top = [];bottom=[];
% % %     [~, roi]  = LoadFeargenAtlas(thr,nroi);
% % %     figure;clf;
% % %     set(gcf,'position',[248         280        1195         504]);
% % %     for lr = 1:length(r)        
% % %         %
% % %         %PLOT THE SIMILARITY MEASURES
% % %         f  = {'B2T' 'UCS2T'};%fields to plot in a subplot
% % %         o  = [1 2];%subplot order of plotting
% % %         for sp = 1:length(f);
% % %             data     = r{lr}.(f{sp});
% % %             %
% % %             dummy = subplot(2,4,o(sp)+2*(lr-1));
% % %             top   = [top dummy];
% % %             hold on;
% % %             %
% % % %             fit  = FitGauss(x(:),data(:),2);
% % % %             fitN = FitGauss(x(:),data(:),1);
% % % %             dof  = fit.dof - fitN.dof;
% % % %             pval = -log10(1-chi2cdf(-2*(fit.Likelihood - fitN.Likelihood),dof)+eps);
% % %             %
% % %             bar(xbase , mean(data,2) , .85 ,'facecolor',[.5 .5 .5],'linestyle','none');
% % %             xlim(xbase([1 8])+[-pi/8 pi/8]);box off
% % %             errorbar(xbase(:) ,mean(data,2),std(data,1,2)./sqrt(ts),'ro');
% % % %             plot(fit.xsup(:),fit.fitup(:),'r');
% % % %             title([f{sp} ': ' mat2str(pval,2)]);
% % %             drawnow;
% % %             %store the fit variables to output
% % % %             Param.(f{sp}).(fit.funname).roi{roi_counter}       = roi.name{1};
% % % %             Param.(f{sp}).(fit.funname).Est(roi_counter,:,lr)  = [fit.Est pval];
% % %         end
% % %         
% % %         
% % %         %PLOT THE MEAN AND VARIANCE MEASURES
% % %         f  = {'mean' 'var'};%fields to plot in a subplot
% % %         o  = [5 6];%subplot order of plotting
% % %         for sp = 1:length(f);
% % %             data   = r{lr}.(f{sp})(:,:,2) - r{lr}.(f{sp})(:,:,1);
% % %             dummy  = subplot(2,4,o(sp)+2*(lr-1));
% % %             bottom = [bottom dummy];
% % %             hold on;
% % %             %
% % % %             fit    = FitGauss(x(:),data(:),2);
% % % %             fitN   = FitGauss(x(:),data(:),1);
% % % %             dof    = fit.dof - fitN.dof;
% % % %             pval   = -log10(1-chi2cdf(-2*(fit.Likelihood - fitN.Likelihood),dof)+eps);
% % %             %
% % %             %bar(xbase , mean(data,2) , .85 ,'facecolor',[.5 .5 .5],'linestyle','none');
% % %             xlim(xbase([1 8])+[-pi/8 pi/8]);box off
% % %             bar(xbase , mean(data,2) , .85 ,'facecolor',[.5 .5 .5],'linestyle','none');
% % %             errorbar(xbase(:) ,mean(data,2),std(data,1,2)./sqrt(ts),'ro');
% % % %             plot(fit.xsup(:),fit.fitup(:),'r');
% % % %             title([f{sp} ': ' mat2str(pval,2)]);
% % %             drawnow;
% % %             %store the fit variables to output
% % % %             Param.(f{sp}).(fit.funname).roi{roi_counter}       = roi.name{1};
% % % %             Param.(f{sp}).(fit.funname).Est(roi_counter,:,lr)  = [fit.Est pval];
% % %         end
% % %         hold off
% % %     end
% % %     EqualizeSubPlotYlim(100,top);
% % %     %EqualizeSubPlotYlim(100,bottom);
% % %     supertitle(regexprep( roi.name{1},'_[R,L]_',''),1,'interpreter','none');
% % %     %
% % % %     SaveFigure(sprintf('%s%s.png',folder,roi.name{1}),'-transparent','-nocrop');
end


