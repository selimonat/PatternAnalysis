function [R RV]=PatternCorrelation(f)
%[R RV]=PatternCorrelation(f)


R  = [];
RV = [];
%
for nf = 1:length(f);
    load(f{nf});
    for nroi = 1:64;
        ccc = 0;
        for s = [1 2]
            GM = mean(zscore([roi(nroi).pat{s,:}]),2);
            
            for c1 = 1:8
                ccc = ccc + 1;
                Y = zscore(roi(nroi).pat{s,c1}) - GM;
                %Y = (roi(nroi).pat{s,c1} - repmat(GM,1,size(roi(nroi).pat{s,c1},2)));
                %Y = zscore(zscore(roi(nroi).pat{s,c1} ')');%.^2-repmat(GM,1,size(roi(nroi).pat{s,c1},2));
                for c2 = 1:8
                    [nf nroi s c1 c2]
                    if c2 >= c1;
                        %
                        X = zscore(roi(nroi).pat{s,c2}) - GM;
                        %X = (roi(nroi).pat{s,c2} - repmat(GM,1,size(roi(nroi).pat{s,c2},2)));
                        nbs = 0;
                        r = [];
                        if size(X,2) ~= 1 && size(Y,2) ~= 1
                            while nbs < 250
                                %
                                nbs    = nbs + 1;
                                %if c1 == c2;
                                    x      = randsample(34,34);
                                    r(nbs) = corr2(mean(Y(:,x(1:17)),2), mean(X(:,x(18:end)),2));
                                %else
                                %    x      = randsample(34,17);
                                %    r(nbs) = corr2(mean(Y(:,x),2), mean(X(:,x),2));
                                %end
                            end
                        elseif size(X,2) == 1 && size(Y,2) == 1
                            
                            r = corrcoef(X,Y);
                            r = r(2,1);
                            
                        end
                        %
                        R(c1,c2,s,nroi,nf) = mean(r);
                        R(c2,c1,s,nroi,nf) = mean(r);
                        %
                        RV(c1,c2,s,nroi,nf) = std(r);
                        RV(c2,c1,s,nroi,nf) = std(r);
                        %
                        %
                    end
                end
                %
% %                 figure(2);
% %                 if ccc == 1
% %                     clf;
% %                 end
% %                 subplot(4,4,ccc);
% %                 imagesc(Y);
                
            end
        end
% %         figure(1)
% %         clf
% %         subplot(2,2,1);imagesc(R(:,:,1,nroi,nf));colorbar;
% %         subplot(2,2,2);imagesc(R(:,:,2,nroi,nf));colorbar
% %         subplot(2,2,3);imagesc(RV(:,:,1,nroi,nf));colorbar
% %         subplot(2,2,4);imagesc(RV(:,:,2,nroi,nf));colorbar
% %         supertitle(roi(nroi).name,1,'interpreter','none','fontsize',25)
        %
        %pause;
    end
    
end