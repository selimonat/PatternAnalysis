function [data]=PatternAnalysis_CommonPatternCorrection(data,w,method)
%[data]=PatternAnalysis_CommonPatternCorrection(data,w,method);
%
%Will apply different pattern correction methods to the DATA. it assumes
%that DATA is [Voxel, Condition, Subject]. Will apply the common pattern
%correction separately for each single subject. It assumes as well that
%conditions 1:8 are baseline and 9:16 are test phase. W has the same size
%as DATA and represents the weights for each voxel. 
%
%METHOD is a string.


fprintf('Will correct for common pattern according to the selected method: %s.\n', method);
for ns = 1:size(data,3)
    if strcmp(method,'baseline');
        %will remove the baseline average from everything
        CP           = squeeze(mean(data(:,1:8,ns),2));
        data(:,:,ns) = data(:,:,ns) - repmat(CP,1,size(data,2));        
    elseif strcmp(method,'best_baseline');        
        %will remove the best fitting baseline average from everything
        CP           = [squeeze(mean(data(:,1:8,ns),2)) ones(size(data,1),1)];
        betas        = CP\data(:,:,ns);
        data(:,:,ns) = data(:,:,ns) - CP*betas;
        
    elseif strcmp(method,'all');      
        %will remove the average computed across all faces...
        CP           = squeeze(mean(data(:,:,ns),2));
        data(:,:,ns) = data(:,:,ns) - repmat(CP,1,size(data,2));
    elseif strcmp(method,'') | strcmp(method,'raw')
        fprintf('No correction requested.\n');
        return;
        
% % %     elseif strcmp(method,'custom')
% % %         CP           = squeeze(mean(data(:,w,ns),2));
% % %         data(:,:,ns) = data(:,:,ns) - repmat(CP,1,size(data,2));
    else        
        fprintf('No such method recognized\n');
        return
    end
end