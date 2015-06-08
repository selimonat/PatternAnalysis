%% load the data first.
load /Volumes/feargen2/cond/midlevel/pa_spmbeta2cov_57_f63684f/simmatpain.mat
[version, ID]     = GetGit(fileparts(which('cond_defaults')));
results_path      = sprintf('%smidlevel/%s_%s_%s',cond_defaults('project_path'),'pa_fit',version,ID);
if exist(results_path) ==0;mkdir(results_path);end
%% put the data (and its residualized version) into a table
% get necessary variables
S                 = sim.raw.correlation;%correlation is actually 1-correlation.
roi_names =[];
for n = 1:length(sim.roiname)
    roi_names{n}  = sim.roiname{n}{1};
end
tgroup            = size(S,5);
troi              = size(S,4);
tsub              = size(S,3);
trialpercond      = size(S,1)/4;
n_elements        = (size(S,1)^2-size(S,1))/2;%number of independent elements in the similarity matrix.
t_elements        = n_elements*tsub*tgroup*troi;
% model regressors
reg_names         = {'session' 'temperature' 'perception'}; 
reg(:,:,1)        = kron([1 1 0 0; 1 1 0 0; 0 0 1 1; 0 0 1 1 ],ones(trialpercond));
reg(:,:,2)        = kron([1 -1 0 0; -1 1 0 0; 0 0 1 1; 0 0 1 1 ],ones(trialpercond));
reg(:,:,3)        = kron([1 -1 1 -1; -1 1 -1 1; 1 -1 1 -1; -1 1 -1 1 ],ones(trialpercond));
% prepare regressors
M = struct();
for r = 1:size(reg,3);
    M.(reg_names{r}) = squareform(1-reg(:,:,r))';
end
%% Create Y, Yres, Indicator vectors
Y           = nan(t_elements,1);
Yres        = nan(t_elements,1);
I.subject   = nan(t_elements,1);
I.roi       = nan(t_elements,1);
I.group     = nan(t_elements,1);
nmat = 0;
DM   = [ones(n_elements,1) M.session];
for ngroup = 1:tgroup
    for nroi = 1:troi
        fprintf('Collecting the data from group: %d, roi: %d\n',ngroup,nroi);
        for nsub = 1:tsub
            %
            nmat        = nmat + 1;
            dummy_Y     = squareform(S(:,:,nsub,nroi,ngroup))';            
            dummy_I     = repmat([nsub nroi ngroup],n_elements,1);
            dummy_Y_res = dummy_Y - DM*(DM\dummy_Y);
            %
            i               = [1:n_elements]+((nmat-1)*n_elements); 
            Y(i)            = dummy_Y;
            Yres(i)         = dummy_Y_res;
            I.subject(i)    = repmat(nsub  ,n_elements,1);
            I.roi(i)        = repmat(nroi  ,n_elements,1);
            I.group(i)      = repmat(ngroup,n_elements,1);
        end
    end
end
%% beta values
beta.temperature = nan(2,length(unique(I.subject)'),length(unique(I.roi)'),length(unique(I.group)'));
beta.perception  = nan(2,length(unique(I.subject)'),length(unique(I.roi)'),length(unique(I.group)'));
DMt              = [ones(630,1) M.temperature];
DMp              = [ones(630,1) M.perception];
DMtp              = [ones(630,1) M.temperature M.perception];
for ngroup = unique(I.group)'    
    for nroi = unique(I.roi)'    
        fprintf('Collecting the data from group: %d, roi: %d\n',double(ngroup),double(nroi));
        for nsub = unique(I.subject)'
            i = logical(I.group == ngroup) & logical(I.roi == nroi) & logical(I.subject == nsub);
            beta.temperature(:,nsub,nroi,ngroup) = DMt\Yres(i);
            beta.perception(:,nsub,nroi,ngroup)  = DMp\Yres(i);
            beta.joint(:,nsub,nroi,ngroup)       = DMtp\Yres(i);
            %
            r.sess(:,nsub,nroi,ngroup)           = corr(Y(i),M.session,'type','Spearman');
            r.temperature(:,nsub,nroi,ngroup)    = corr(Y(i),M.temperature,'type','Spearman');
            r.perception(:,nsub,nroi,ngroup)     = corr(Y(i),M.perception,'type','Spearman');
        end
    end
end
%%
% Plot
betas       = beta.joint;
troi        = 124;
test_limits = [0.0001;0.0005;0.001;0.005;0.01;0.05];
nbeta       = 3;
for ngroup  = 1:2;    
    %make the test and get best TROI
    [h pval] = ttest(betas(:,:,ngroup,nbeta));
    [~,i]    = sort(pval);
    B        = mean(betas(:,i(1:troi),ngroup,nbeta));
    pval     = pval(i(1:troi));
    %
    subplot(1,2,ngroup);
    n_asterix  = sum(repmat(log(pval),length(test_limits),1) < repmat(log(test_limits),1,troi))
    barh(B,'k');
    for nroi = 1:troi
        text(B(nroi)-B(nroi)*.1,nroi,repmat('*',1,n_asterix(nroi)),'color','r','fontsize',25)
    end
    set(gca,'ytick',1:troi,'yticklabel',roi_names(i(1:troi)),'ygrid','on')
    box off;
    xlabel('beta weight')
    title('Temperature')
    axis tight;
end
EqualizeSubPlotXlim(gcf)
%% now analyze the difference
%compare groups
[h, pval] = ttest2(beta.joint(:,:,1,3),beta.joint(:,:,2,3));
[~,i]     = sort(pval);
barh(mean(beta.joint(:,i,1,3))-mean(beta.joint(:,i,2,3)));
    set(gca,'ytick',1:troi,'yticklabel',roi_names(i(1:troi)),'ygrid','on')


