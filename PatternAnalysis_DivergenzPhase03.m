%% BLOCK ANALYSIS
clear all;
nroi     = 6;
subjects = [1 4 6 7 8 9 10 11 12 13 14 16 18 23 25 26];
subjects = 1:28;
tsub     = length(subjects);
DM3      = PatternAnalysis_BetaFiles2Matrix('load',50,'anova_cov_',3,'faces/chrf_derivs_00_tBlock_05_AbsoluteTime','swRealigned','',nroi,subjects);
% DM2      = PatternAnalysis_BetaFiles2Matrix('load',50,'anova_cov_',4,'faces/chrf_derivs_00_fullcovariate','swRealigned','',nroi,subjects);

for lr       = 1:2;    
    p2p      = [];p2n      = [];u2p      = [];
    tbs      = 1000;
    for ncond = 1:5
        [lr ncond]
        for nbs = 1:tbs            
            %
            s1          = randsample(1:tsub,tsub,1);
            s2          = s1;%randsample(1:tsub,tsub,1);
            %from the baseline
            CSPm        = mean(mean(DM3{lr}(:,1:5 ,s1),3)');
            CSNm        = mean(mean(DM3{lr}(:,6:10,s2),3)');
            %from the conditioning phase
            CSP         = mean(DM3{lr}(:,ncond  , s1),3)';
            CSN         = mean(DM3{lr}(:,ncond+5, s2),3)';
            UCS2        = mean(mean(DM3{lr}(:,11:15,s2),3)');
            %
%             r(lr).p2UCS2(nbs,ncond) = corr2(CSP,UCS2);
%             r(lr).n2UCS2(nbs,ncond) = corr2(CSN,UCS2); %CSN and CSP to the UCS of the phase3
            r(lr).p2UCS(nbs,ncond)  = corr2(CSP,UCS2);  %CSP to UCS of baseline
            r(lr).n2UCS(nbs,ncond)  = corr2(CSN,UCS2);  %CSN to UCS of baseline
            r(lr).p2B(nbs,ncond)    = corr2(CSP,CSNm);    %CSP to B
            r(lr).n2B(nbs,ncond)    = corr2(CSN,CSPm);
            
            r(lr).p2p(nbs,ncond)    = corr2( ...
                nanmean(DM3{lr}(:,1,s1),3), ...
                nanmean(DM3{lr}(:,ncond,s2),3) ...
                );
            r(lr).n2n(nbs,ncond)    = corr2( ...
                nanmean(DM3{lr}(:,6,s1),3), ...
                nanmean(DM3{lr}(:,ncond+5,s2),3) ...
                );
            r(lr).p2n(nbs,ncond)    = corr2( ...
                nanmean(DM3{lr}(:,ncond,s1),3), ...
                nanmean(DM3{lr}(:,ncond+5,s2),3) ...
                );
        end
    end
end
%%
figure;
set(gcf,'position',[680 49 527 1027])
for lr = 1:2
    subplot(4,2,lr)
    errorbar(1:5,mean(r(lr).p2p),std(r(lr).p2p),'k.-','linewidth',2);
    hold on
    errorbar(1:5,mean(r(lr).n2n),std(r(lr).n2n),'r.-','linewidth',2)
    legend({'csp2csp' 'csn2csn'});legend boxoff
    set(gca,'xtick',1:5)
    %
    subplot(4,2,lr+2)
    errorbar(1:5,mean(r(lr).p2n),nanstd(r(lr).p2n),'k.-','linewidth',2);    
    legend({'csp2csn'});legend boxoff
    set(gca,'xtick',1:5)
    %
    subplot(4,2,lr+4)
    errorbar(1:5,nanmean(r(lr).p2UCS),nanstd(r(lr).p2UCS),'r.-','linewidth',2);    
    hold on
    errorbar(1:5,nanmean(r(lr).n2UCS),nanstd(r(lr).n2UCS),'k.-','linewidth',2);    
    legend({'csp2UCS' 'csn2UCS'});legend boxoff
    set(gca,'xtick',1:5)
    %
    subplot(4,2,lr+6)
    errorbar(1:5,nanmean(r(lr).p2B),nanstd(r(lr).p2B),'r.-','linewidth',2);    
    hold on
    errorbar(1:5,nanmean(r(lr).n2B),nanstd(r(lr).n2B),'k.-','linewidth',2);    
    legend({'csp2B' 'csn2B'});legend boxoff;
    set(gca,'xtick',1:5)
end

%% SINGLE TRIAL CORRELATION ANALYSIS

CSP = load_nii('/Volumes/feargen2/feargen2/data/midlevel/singletrial/phase02/singletrial_cond_04.nii.gz');
CSP = CSP.img;s = size(CSP);CSP = reshape(CSP,prod(s(1:3)),s(4));i = ~isnan(sum(CSP,2));CSP_B = CSP(i,:);
CSN = load_nii('/Volumes/feargen2/feargen2/data/midlevel/singletrial/phase02/singletrial_cond_08.nii.gz');
CSN = CSN.img;s = size(CSN);CSN = reshape(CSN,prod(s(1:3)),s(4));CSN_B = CSN(i,:);
UCS = load_nii('/Volumes/feargen2/feargen2/data/midlevel/singletrial/phase02/singletrial_cond_09.nii.gz');
UCS = UCS.img;s = size(UCS);UCS = reshape(UCS,prod(s(1:3)),s(4));UCS_B = UCS(i,:);
%
CSP = load_nii('/Volumes/feargen2/feargen2/data/midlevel/singletrial/phase03/singletrial_cond_04.nii.gz');
CSP = CSP.img;s = size(CSP);CSP = reshape(CSP,prod(s(1:3)),s(4));i = ~isnan(sum(CSP,2));CSP_C = CSP(i,:);
CSN = load_nii('/Volumes/feargen2/feargen2/data/midlevel/singletrial/phase03/singletrial_cond_08.nii.gz');
CSN = CSN.img;s = size(CSN);CSN = reshape(CSN,prod(s(1:3)),s(4));CSN_C = CSN(i,:);
UCS = load_nii('/Volumes/feargen2/feargen2/data/midlevel/singletrial/phase03/singletrial_cond_09.nii.gz');
UCS = UCS.img;s = size(UCS);UCS = reshape(UCS,prod(s(1:3)),s(4));UCS_C = UCS(i,:);

%%




%%
C  = corrcoef([CSP_B CSN_B UCS_B CSP_C CSN_C UCS_C ]);
figure;
imagesc(C)
%
figure(2);clf;
% C      = fisherz(C);
subplot(2,2,1);
dummy  = ( mean(C(1:32,33:64)));
dummys = ( std(C(1:32,33:64)));
errorbar(linspace(0,1,length(dummy)),dummy,dummys,'ro-','marker','none');
%
hold on;
dummy  = (mean(C(1:32,65:74)));
dummys = (std(C(1:32,65:74)));
errorbar(linspace(0,1,length(dummy)),dummy,dummys,'ko-','marker','none');
%
dummy  = (mean(C(33:64,65:74)));
dummys = (std(C(33:64,65:74)));
errorbar(linspace(0,1,length(dummy)),dummy,dummys,'mo-','marker','none');
hold off
axis tight
set(gca,'ylim',[-1 1])
legend({'CSPvsCSN' 'UCSvsCSP' 'UCSvsCSN'})

subplot(2,2,4);
dummy  = (mean(C(75:103,104:162)));
dummys = (std(C(75:103,104:162)));
% errorbar(linspace(0,1,length(dummy)),dummy,dummys,'ro-','marker','none');
plot(linspace(0,1,length(dummy)),dummy,'ro-')
hold on;
%
dummy  = (mean(C(75:103,163:end)));
dummys = (std(C(75:103,163:end)));
% errorbar(linspace(0,1,length(dummy)),dummy,dummys,'ko-','marker','none');
plot(linspace(0,1,length(dummy)),dummy,'ko-')
%
dummy  = (mean(C(104:162,163:end)));
dummys = (std(C(104:162,163:end)));
% errorbar(linspace(0,1,length(dummy)),dummy,dummys,'mo-','marker','none');
plot(linspace(0,1,length(dummy)),dummy,'mo-')
hold off
axis tight
set(gca,'ylim',[-1 1])
legend({'CSPvsCSN' 'UCSvsCSP' 'UCSvsCSN'})

subplot(2,2,2);
dummy  = (mean(C(1:32,75:103)));
dummys = (std(C(1:32,75:103)));
errorbar(linspace(0,1,length(dummy)),dummy,dummys,'ro-','marker','none');
%
hold on;
dummy  = (mean(C(33:64,104:162)));
dummys = (std(C(33:64,104:162)));
errorbar(linspace(0,1,length(dummy)),dummy,dummys,'ko-','marker','none');
%
dummy  = (mean(C(65:74,163:end)));
dummys = (std(C(65:74,163:end)));
errorbar(linspace(0,1,length(dummy)),dummy,dummys,'mo-','marker','none');

axis tight
set(gca,'ylim',[-1 1])

legend({'CSPvsCSP' 'CSNvsCSN' 'UCSvsUCS'})


