function [C]=pa_GetSimMat(pattern,metrics)
%[C]=pa_GetSimMat(pattern,what)
%
%   Returns a similarity matrix organized as
%   [condition,condition,subject,roi,phase]. Similarity is computed using
%   different metrics as precised in METRICS. It can be any of the
%   following: {'mean' 'cov' 'correlation' 'euclidean' 'seuclidean'}, and
%   must be a cell array. If nothing is given then all metrics are
%   computed.
%
%   Dependency: pa_defaults, pa_GetY, pa_GetSPMDesignMatrix,
%   pa_GetNuissance, pa_GetBetas

threshold  = pa_defaults('threshold');
if nargin == 1
    metrics = {'mean' 'cov' 'correlation' 'euclidean' 'seuclidean'};
end
%clean pattern string
cleaned_pattern = pattern(isstrprop(pattern,'alphanum'));
tcond           = pa_defaults('tcond');
tphase          = length(pa_defaults('phases'));
for metric = metrics
    %init
    covmat = zeros(tcond*2,tcond*2,length(pa_defaults('gs')),pa_defaults('troi'));
    save_path  = sprintf('%s%s/%s_%s_%02d.mat',pa_defaults('save_path'),mfilename,cleaned_pattern,metric{1},threshold);
    if exist(save_path) == 0;
        
        for roi = 1:pa_defaults('troi')
            ns = 0;
            for subject = pa_defaults('gs')
                fprintf('%s: Getting covmat %02d-%02d.\n',mfilename,subject,roi);
                ns = ns + 1;
                np = 0;
                B  = [];
                for phase = pa_defaults('phases')
                    %%
                    np = np + 1;
                    Y     = pa_GetY(subject,phase,roi,threshold,pattern);
                    X     = pa_GetSPMDesignMatrix(subject,phase);
                    N     = pa_GetNuissance(subject,phase,pattern,{'mc_diff_square' 'csf'});
                    N     = [N(:).val];
                    B     = [B pa_GetBetas(X,Y,N)'];
                end
                
                %% This is the core where different metrics are computed
                if any(strcmp(metric,{'cov' 'mean'}))
                    fh       = str2func(metric{1});
                    cmat     = fh(B);
                    if isvector(cmat)%this is necessary to store the mean vector
                        cmat = diag(cmat);
                    end
                elseif any(strcmp(metric{1},{'correlation' 'euclidean' 'seuclidean'}))
                    cmat     = squareform(pdist(B,metric{1}));
                end
                %% store
                covmat(:,:,ns,roi) = cmat;
            end
        end
        save(save_path,'covmat')
    else
        load(save_path);%will spawn covmat
    end
    C.(metric{1}) = covmat;
end





