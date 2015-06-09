load('/Users/onat/Dropbox/project_cond/simmatpain1regr.mat')

%% get model
O1 = squareform([1 0 1 1; 0 1 1 1; 1 1 1 2; 1 1 2 1]-1)';
O2 = squareform(ones(4)-eye(4))';
%%
n       = size(sim.mc.cov,1);
ind     = tril(logical(ones(n)));%matrix 2 vector wo redundant entries
tentry  = sum(ind(:));
nroi    = 17;%amy or 58
tsub    = size(sim.mc.cov,3);
M       = [];
%%
model = [];
M = [];
for ngroup = 1%:2;
    for nsub    = 1:tsub;
        nsub
        %
        covmat      = sim.mc.cov(:,:,nsub,nroi,ngroup);        
        [~,corrmat] = cov2corr(covmat);
        dmat        = 1- corrmat;
        dmat([1 6 11 16]) = 0;
        dmat     = squareform(dmat);
        %
        M       = [[dmat(:) O1 repmat(nsub,6,1)]];        
        tbl     = array2table(M,'VariableNames',{'Y' 'Temp' 'Sub'});%'Low' 'High' 'Sess' 'Diag' 'Sub' 'Group'
        model{nsub}       = fitlme(tbl,'Y ~ 1 + Temp');  
        beta(:,nsub)      = fixedEffects(model{nsub});
        fitted            = squareform(    model{nsub}.fitted);
        imagesc(1-[squareform(dmat) fitted])
        drawnow
        pause;
        %%        
    end
end
%%


