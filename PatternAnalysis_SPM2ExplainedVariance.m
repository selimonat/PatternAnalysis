function [explained_variance]=PatternAnalysis_SPM2ExplainedVariance(SPM,type,regressors,coor)
%[explained_variance]=PatternAnalysis_SPM2ExplainedVariance(SPM,type,coor,regressors)
%
%if TYPE is test: will plot the time-series of the voxels at COOR and will
%displays histogram of time-seris for raw data, fitted data and the
%residuals, will also print explained variance. Will compare the SPM values
%of explained variance to my computation. This is just for sanity check
%purpose.
%
%if TYPE is ALL, it will ignore the COOR and directly compute the explained
%variance for all voxels, will also write this as a volume to the SPM
%directory. In this case EXPLAINED_VARIANCE is a vol handle.
%
%If REGRESSORS are provided then the explained variance due to only the
%REGRESSORS is computed.
%
%   This is called automatically from the PatternAnalysis_BetaFiles2Matrix
%   routine, using the full model.
%
%   
%Selim Onat, 2014
total_SS = [];
res_SS   = [];
fit_SS   = [];
fit      = [];
here     = pwd;
SPM.swd  = regexprep(SPM.swd,'/Users/onat/Desktop','/Volumes/feargen2/feargen2/data');
cd(SPM.swd);
%
if ~isempty(regressors)
    fprintf('Will consider a reduced model:\n')
    fname = 'reduced';
else
    fprintf('Model is full model\n');
    regressors = SPM.xX.iC;
    fname = 'full';
end
%
if strcmp(type,'all')%compute a r2 volume
    
    [res_SS]      = GetResMS('all');%returns the residuals for all voxels.
    %
    %returns the betas for all voxels.
    betas_all     = spm_read_vols(spm_vol(SPM.Vbeta));
    s             = size(betas_all);
    betas_all     = reshape(betas_all,[prod(s(1:3)), s(4)])';
    s(end)        = [];
    %
    %returns the valid voxels.
    mask      = spm_read_vols(spm_vol(SPM.VM));
    mask      = reshape(mask,[prod(s(1:3)) 1])';
    voxels    = find(mask);
    %
    %init the final variable
    r2        = nan(1,prod(s));
    for nvox = voxels
        betas     = betas_all(SPM.xX.iC,nvox);        
        [fit_SS]  = GetModelSS('matrixmultip',betas,regressors);%will use betas        
        r2(nvox)  = fit_SS./(res_SS(nvox)+fit_SS);
    end
    %write the r2 to a file
    r2 = reshape(r2,[s(1) s(2) s(3)]);
    %
    explained_variance         = SPM.VResMS;
    explained_variance.pinfo   = [1 0 0]';
    if isempty(regressors)
        explained_variance.fname   = sprintf('%s/r2_%smodel.img',SPM.swd,fname);
    else
        explained_variance.fname   = sprintf('%s/r2_%smodel.img',SPM.swd,fname);
    end
    explained_variance.descrip = sprintf('%s: Explained Variance',mfilename);
    fprintf('Will write the volume to:\n %s\n',explained_variance.fname)
    spm_write_vol(explained_variance,r2);
    
    
    
elseif strcmp(type,'test')
    
    %load the necessary stuff for these GetTotalSS, GetModelSS and GetResMS
    %functions
    %read the data and low pass filter it
    
    %load the betas
    betas    = spm_get_data(SPM.Vbeta(:),coor);
    %
    Y        = spm_get_data(SPM.xY.VY(:),coor);
    KWY      = spm_filter(SPM.xX.K,Y);
    KWY      = KWY - SPM.xX.xKXs.X(:,SPM.xX.iB)*betas(SPM.xX.iB);
    
    
    fprintf('Total Variance directly from data: %g\n',GetTotalSS(KWY));
    fprintf('Total Variance: res_SS + fit_SS: %g\n',GetResMS('projection') + GetModelSS('matrixmultip',betas,regressors));
    fprintf('residual SS from SPM ResMS.img: %f\n',GetResMS('likespm'));
    fprintf('residual SS from me: %f\n',GetResMS('projection'));
    %
    explained_variance = GetModelSS('matrixmultip',betas,regressors)/GetTotalSS(KWY);
    %
    fprintf('Explained Variance: %g\n',explained_variance );
    figure(103);
    %
    subplot(2,3,1:3)
    plot(KWY,'ro');
    hold on
    plot(fit);
    plot(res,'g');
    axis tight;
    legend('raw','fit','res');
    hold off;
    tbin = 100;
    h(1) = subplot(2,3,4);
    hist(KWY,tbin);title('raw data');
    h(2) = subplot(2,3,5);
    hist(fit,tbin);title('fit');
    h(3) = subplot(2,3,6);
    hist(res,tbin);title('residuals');
    EqualizeSubPlotXlim(103,h);
end
cd(here);
    function [res_SS]=GetResMS(how)
        if strcmp(how,'likespm');%residual variance from spm
            res_SS   = spm_get_data(SPM.VResMS,coor)*SPM.VResMS.pinfo(1)^-1;
        elseif strcmp(how,'projection');
            res     = spm_sp('r',SPM.xX.xKXs,KWY);
            res_SS  = res'*res;
        elseif strcmp(how,'all');
            res_SS     = spm_read_vols(SPM.VResMS)*SPM.VResMS.pinfo(1)^-1;
            res_SS     = res_SS(:);
        end
    end

    function [total_SS]=GetTotalSS(KWY)
        
        %correct for the mean shift, this is still real data
        %KWY     = KWY - SPM.xX.xKXs.X(:,SPM.xX.iB)*betas(SPM.xX.iB);
        %total variance from the raw data (corrected for session shifts)
        total_SS = KWY'*KWY;
    end

    function [fit_SS]=GetModelSS(how,betas,RegORData)
        
        if strcmp(how,'matrixmultip');                        
            %RegORData is REGRESSOR indices right now
            fit        = SPM.xX.xKXs.X(:,RegORData)*betas(RegORData);
        elseif strcmp(how,'likespm');            
            %RegORData is the DATA right now
            fit        = spm_sp('op',SPM.xX.xKXs,RegORData);
        end
        fit_SS   = fit'*fit;
    end
end