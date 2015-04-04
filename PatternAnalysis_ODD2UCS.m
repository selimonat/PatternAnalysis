%%
load /Volumes/feargen2/feargen2/data/midlevel/selectedsubjects.mat
load /Volumes/feargen2/feargen2/data/midlevel/Behavior
%%
dtype    = 1;
pv       = -log10(subject_view_FGProfile(s,2,dtype,subject_list));
[y subs] = sort(pv,'descend');
param    = 3;
subjects = subs(1:14);
tsub     = length(subjects);
%%
r = [];
roi = 0;
rois = [2 23 41 31 17 10 36 26 15 16 30 35 42 43]
troi = length(rois);
for nroi= rois
    roi = roi + 1;
    for lr = 1:2
        [data,dataw,roiname]= PatternAnalysis_BetaFiles2Matrix('load',50,'anova_cov_', 4 , 'faces/chrf_derivs_00_fullcovariate' , 'swRealigned','',nroi,1:29)    
        w                    = mean(dataw{lr}(:,1,subjects),3)';
        for nbs = 1:100
            ODD         = mean(data{lr}(:,18,randsample(subjects,tsub,1)),3)';
            UCS         = mean(data{lr}(:,17,randsample(subjects,tsub,1)),3)';
            
            r{lr}.UCS2ODD(roi,nbs)     = pdist2(UCS,ODD, @(X,Y) wcorr2(X,Y,w));
        end        
    end
end
%%
figure;
[~, roi]  = LoadFeargenAtlas(50,rois);
for lr = 1:2
    subplot(1,2,lr)
    barh(1:troi, mean(r{lr}.UCS2ODD,2));
    set(gca,'yticklabel',roi.name)
    hold on
    herrorbar(mean(r{lr}.UCS2ODD,2),1:troi,std(r{1}.UCS2ODD,1,2),'ro')
    hold off
end