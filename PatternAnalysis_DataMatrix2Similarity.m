function [SM,fun]=PatternAnalysis_DataMatrix2Similarity(DM,method)
%
%SM: similarity matrix
%DM: data matrix
%
%FUN is used to extract subcomponents from the SM. 


%COMPUTE THE SIMILARITY MATRIX FOR EACH SUBJECT
for lr = 1:2
    SM{lr} = nan(size(DM{lr},2),size(DM{lr},2),size(DM{lr},3));
    for ns = 1:size(DM{lr},3)
        dm = DM{lr}(:,:,ns);
        if strcmp(method,'correlation');            
            sm = squareform(  1- pdist(dm',method));
        elseif any(strcmp(method,{'euclidean' 'seuclidean' 'mahalonobis'}))            
            sm = squareform(pdist(dm',method));
        elseif strcmp(method,'covariance') | strcmp(method,'cov')                        
            sm = cov(dm);
            
            % %     elseif strcmp(method,'robust_covariance')
            % %
            % %         sm = mcdcov(dm');
            % %
            % %     elseif strcmp(method,'symm_dkl');
            % %
            % %
            % %     elseif strcmp(method,'mutual_info');                        
        elseif strcmp(method,'mean') || strcmp(method,'average');            
            sm = diag(mean(dm));            
        elseif strcmp(method,'var') || strcmp(method,'variance');            
            sm = diag(var(dm));            
        else
            SM = [];
            return
        end
        SM{lr}(:,:,ns) = sm;
    end
end

%% DEFINE ALL THE FUNCTIONS THAT EXTRACTS SUBCOMPONENTS FROM THE SIM MATRIX

%CSP vs ALL
fun.B.CSPvsAll.raw     = @(lr) squeeze(SM{lr}(1:8,4,:))';
fun.T.CSPvsAll.raw     = @(lr) squeeze(SM{lr}(9:16,12,:))';
fun.D.CSPvsAll.raw     = @(lr) squeeze(SM{lr}(9:16,12,:))' - squeeze(SM{lr}(1:8,4,:))';

fun.B.CSPvsAll.mean    = @(lr) mean(squeeze(SM{lr}(1:8,4,:))');
fun.T.CSPvsAll.mean    = @(lr) mean(squeeze(SM{lr}(9:16,12,:))');
fun.D.CSPvsAll.mean    = @(lr) mean(squeeze(SM{lr}(9:16,12,:))') - mean(squeeze(SM{lr}(1:8,4,:))');

fun.B.CSPvsAll.sem     = @(lr) sem(squeeze(SM{lr}(1:8,4,:))');
fun.T.CSPvsAll.sem     = @(lr) sem(squeeze(SM{lr}(9:16,12,:))');
fun.D.CSPvsAll.sem     = @(lr) sem(squeeze(SM{lr}(9:16,12,:))') - sem(squeeze(SM{lr}(1:8,4,:))');


%diagonals
fun.B.diag.raw         = @(lr) reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(1:8,1,size(SM{lr},3)),repmat(1:8,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))';
fun.T.diag.raw         = @(lr) reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(9:16,1,size(SM{lr},3)),repmat(9:16,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))';
fun.D.diag.raw         = @(lr) reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(9:16,1,size(SM{lr},3)),repmat(9:16,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))' - reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(1:8,1,size(SM{lr},3)),repmat(1:8,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))';

fun.B.diag.mean        = @(lr) mean(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(1:8,1,size(SM{lr},3)),repmat(1:8,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))');
fun.T.diag.mean        = @(lr) mean(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(9:16,1,size(SM{lr},3)),repmat(9:16,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))');
fun.D.diag.mean        = @(lr) mean(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(9:16,1,size(SM{lr},3)),repmat(9:16,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))') -mean(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(1:8,1,size(SM{lr},3)),repmat(1:8,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))');

fun.B.diag.sem         = @(lr) sem(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(1:8,1,size(SM{lr},3)),repmat(1:8,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))');
fun.T.diag.sem         = @(lr) sem(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(9:16,1,size(SM{lr},3)),repmat(9:16,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))');
fun.D.diag.sem         = @(lr) sem(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(9:16,1,size(SM{lr},3)),repmat(9:16,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))') - sem(reshape(SM{lr}(sub2ind([size(SM{lr})],repmat(1:8,1,size(SM{lr},3)),repmat(1:8,1,size(SM{lr},3)),Vectorize(repmat(1:size(SM{lr},3),8,1))')),8,size(SM{lr},3))');

%cross comparison
fun.D.Cross.raw     = @(lr) squeeze(SM{lr}(1:8,12,:))';
fun.D.Cross.mean    = @(lr) mean(squeeze(SM{lr}(1:8,12,:))');
fun.D.Cross.sem     = @(lr) sem(squeeze(SM{lr}(1:8,12,:))');

%cross comparison integrated
fun.D.CrossAll.raw     = @(lr) squeeze(mean(SM{lr}(1:8,9:16,:)))';
fun.D.CrossAll.mean    = @(lr) mean(squeeze(mean(SM{lr}(1:8,9:16,:)))');
fun.D.CrossAll.sem     = @(lr) sem(squeeze(mean(SM{lr}(1:8,9:16,:)))');





