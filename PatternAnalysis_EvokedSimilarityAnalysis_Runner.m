%% CREATING THE ACTIVITY MAPS
thr     = 30;
preprop = 'swRealigned';
phase   = 4;
folder  = 'chrf_derivs_00'; 
secondlevelfolder = 'faces';
ttype   = 'anova_cov_';
%%
[datalr,w_r2lr,roiname] = PatternAnalysis_BetaFiles2Matrix('cache',thr,ttype,phase,sprintf('%s/%s',secondlevelfolder,folder),preprop);
%%
zsc = 0;
PatternAnalysis_ActivityMaps({ttype,phase,sprintf('%s/%s',secondlevelfolder,folder),preprop},[1:10;11:20]',1:29,41,2.65,'',zsc,'swRealigned',thr);
%%
zsc       = 0;
gs        = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
PatternAnalysis_ActivityMaps({ttype,phase,sprintf('%s/%s',secondlevelfolder,folder),preprop},[1:9;11:19]',gs,[18],2.65,'',zsc,'gs_testnew',thr);

%% BOLD FG ANALYSIS
load /Volumes/feargen2/feargen2/data/midlevel/selectedsubjects.mat
load /Volumes/feargen2/feargen2/data/midlevel/Behavior
%%
SPMid      = {'anova_cov_',4,'faces/chrf_derivs_00','swRealigned'};
paramnames = {'amp' 'sigma' 'ttest' '' 'pval'};
%%
PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',1:29,'swRealigned');
PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked_weighted',1:29,'all');
PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',1:29,'AllSubjects');
%% sort by the difference between CS+ and CS- during Conditioning on SCR
dtype    = 1;
param    = 3;
pv       = -log10(subject_view_FGProfile(s,2,dtype,subject_list));
[y subs] = sort(pv,'descend');
%
%evoked
%MSparams(dtype,param).high = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',subs(1:15),sprintf('%s_SCR_+',paramnames{param}));
%MSparams(dtype,param).low  = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',subs(16:end),sprintf('%s_SCR_-',paramnames{param}));
%induced
MSparams(dtype,param).high = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',subs(1:15),sprintf('%s_SCR_+',paramnames{param}));
MSparams(dtype,param).low  = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',subs(16:end),sprintf('%s_SCR_-',paramnames{param}));
%% PLOT GOOD GAUSSIAN SUBJECTS AND BAD GAUSSIAN SUBJECTS
%sort by SCR goodness of the model of the TEST phase.
dtype        = 1;%scr
phase        = 3;%test
param        = 5;%
%
[Est]        = subject_GetFitParam(s,subject_list,dtype,phase);
[y i]        = sort(Est(:,end-1),'descend');%pval,
goodsubjects = i(1:14);
badsubjects  = setdiff(1:29,goodsubjects);%EST is in register with the subject vector 1:29;
% MSparams(dtype,param).high = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',goodsubjects,sprintf('%s_SCR_+',paramnames{param}));
% MSparams(dtype,param).low  = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',badsubjects,sprintf('%s_SCR_-',paramnames{param}));
  MSparams(dtype,param).high = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',goodsubjects,sprintf('%s_SCR_+',paramnames{param}));
  MSparams(dtype,param).low  = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',badsubjects,sprintf('%s_SCR_-',paramnames{param}));

%% sort by Rating sigma of the Conditioning phase
dtype = 5;%rating
phase = 2;%test phase
%these are selected according to the picture in
%/Users/onat/Pictures/Ratings_TestPhase.png/, I do it by hand because there
%are some subjects which are just flat
param = 2;
%this part selects subjects with a nice Gaussian rating profiles AND
%sorts them according to their sigma width
[Est] = subject_GetFitParam(s,subject_list,dtype,phase);
[y i] = sort(Est(:,end-1),'descend');
invalid = i(find(y < -log10(.05)));%the indices of the invalid subjects
%
Est(invalid,:) = NaN;
[y i]          = sort(Est(:,param),'descend');
i(isnan(y))    = [];
tsub           = round(sum(~isnan(i))/2);
subs_high      = i(1:tsub);
subs_low       = i(tsub+1:end);
%
MSparams(dtype,param).high = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',subs_high,sprintf('%s_Cond_Rating_+wide',paramnames{param}));
MSparams(dtype,param).low  = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'evoked',subs_low,sprintf('%s_Cond_Rating_-tight',paramnames{param}));
% %     MSparams(dtype,param).high = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',subs_high,sprintf('_induced_%s_Cond_Rating_GoodAndWide',paramnames{param}));
% %     MSparams(dtype,param).low  = PatternAnalysis_EvokedSimilarityAnalysis(SPMid,0,'induced',subs_low,sprintf('_induced_%s_Cond_Rating_GoodAndTight',paramnames{param}));
