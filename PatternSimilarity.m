function [Sim]=PatternSimilarity(subjects)
%% Condition based fmri analysis
% in this analysis trials are not considered, similarity is measured between
% condition specific patterns and similarity values are averaged across
% subjects.
%%
Sim = [];
for subject         = subjects;
    fprintf('Subject: %03d...\n',subject)
    %%
    phase           = 4;
    pattern         = '^swRealigned.*nii$';
    roi             = 19;
    %%
    %first create an SPM array
    SPM             = GetSPMArray(subject,phase,pattern,roi);
    
    
    [Y,xyz]         = PatternAnalysis_GetNativeData(subject,phase,roi);
    
    SPM.xX.xKXs     = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.X));       % KWX
    SPM.xX.xKXs.X   = full(SPM.xX.xKXs.X);
    SPM.xX.pKX      = spm_sp('x-',SPM.xX.xKXs);
    KWY             = spm_filter(SPM.xX.K,Y);
    beta            = SPM.xX.pKX*KWY;
    D               = beta(1:8,:);
    Sim             = cat(3,Sim,cov(D'));
end
    
end