function [C]=pa_GetCovMats(pattern)

threshold = pa_defaults('threshold');
C = zeros(10,10,length(pa_defaults('gs')),2,pa_defaults('troi'));
for roi = 1:pa_defaults('troi')
    ns = 0;
    for subject = pa_defaults('gs')
        ns = ns + 1;
        np = 0;
        for phase = [2 4]
            %%
            np = np + 1;
            Y     = pa_GetY(subject,phase,roi,threshold,pattern);
            X     = pa_GetSPMDesignMatrix(subject,phase);
            N     = pa_GetNuissance(subject,phase,pattern,{'mc_diff_square' 'csf'});
            N     = [N(:).val];
            B     = pa_GetBetas(X,Y,N)';
            %cov:
            cmat     = cov(B);
        end
        C(:,:,ns,np,roi) = cmat;
    end
end





