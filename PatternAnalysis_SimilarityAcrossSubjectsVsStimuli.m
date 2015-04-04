
gs              = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
[data,roiname]=PatternAnalysis_BetaFiles2Matrix(41,gs,'raw');
%

for lr = 1:2
    
    r_acrossstim = [];    
    for nstim = 1:16;%exclude the UCS and ODDball conditions
        mat          = squeeze(data{lr}(:,nstim,:));
        r_acrossstim = [r_acrossstim 1-pdist(mat','correlation')];
    end
    subplot(2,2,1+lr-1)
    hist(r_acrossstim,100)
    title('correlation across subjects');
    %    
    r_acrosssubjects = [];    
    for nsub = 1:size(data{lr},3)
        mat              = squeeze(data{lr}(:,1:16,nsub));
        r_acrosssubjects = [r_acrosssubjects 1-pdist(mat','correlation')];
    end
    subplot(2,2,3+lr-1)
    hist(r_acrosssubjects,100)
    title('correlation across stimuli');
end