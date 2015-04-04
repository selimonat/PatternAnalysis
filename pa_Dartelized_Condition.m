%%
cd /Users/onat/Documents/Code/Matlab/PatternAnalysis;
save_path         = '/Users/onat/Pictures/PatternAnalysis/'
ventricle         = 38;%necessary for constructing a nuissance regressor
HParam            = 128;%high pass filter parameter
threshold         = 75;
pattern           = '^swRealigned.*nii$';
gs                = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
load /Volumes/feargen2/feargen2/data/midlevel/selectedsubjects.mat%
%get the version of number of the analysis
[~,version_count] = system('git rev-list --count --first-parent HEAD')
[~,version_id]    = system('git describe --always')
for phase     = [2 4];
    for roi       = 1:96;        
        betas     = [];
        for subject = subject_list(gs);
            %%
            % get data
            Y = pa_GetY(subject,phase,roi,threshold,pattern);
            % remove mean
            Y = Y - repmat(mean(Y),size(Y,1),1);
            % get the design matrix
            X = pa_GetDesignMatrix(subject,phase);
            %% plot the correlation matrix of X
% %             figure(1);
% %             subplot(1,2,1);
% %             imagesc(corrcoef(X).^2,[0 1]);colorbar;
% %             SaveFigure(sprintf('%ssub%03d/phase%02d/figures/%s.png',GetProjectRootFeargen2,subject,phase,'DesignCorrelation'));
            %% Get the nuissance variables.
            [Nuis]                = GetNuissanceVariables(subject,phase,pattern,{'mc_diff_square'});            
            tNuis                 = length(Nuis);
            Nuis(tNuis + 1).val   = zscore(mean(pa_GetY(subject,phase,ventricle,75,pattern),2));
            Nuis(tNuis + 2).val   = zscore(mean(pa_GetY(subject,phase,mod(ventricle-1,48)+1+48,75,pattern),2));
            Nuis(tNuis + 1).name  = 'RH_csf';
            Nuis(tNuis + 2).name  = 'LH_csf';
            % Nuissance Matrix
            M                     = [Nuis(:).val];
            %% Get the High-pass filter as well            
            K = struct('HParam', HParam , ...
                'row',    1:size(X,1) , ...
                'RT',     2.02);
            K = spm_filter(K);
            %%
            N = [M K.X0];%combine covariates            
            N  = N - repmat(mean(N),size(N,1),1);%and demean them
            %% plot the nuissance correlation
%             figure(1)
%             subplot(1,2,2);
%             imagesc(corrcoef(M).^2,[0 1]);
%             colorbar;
%             SaveFigure(sprintf('%ssub%03d/phase%02d/figures/%s.png',GetProjectRootFeargen2,subject,phase,'NuissanceCorrelation'));
            %% first residualize with respect to Nuissance parameters.
            Yr = Y - N*(N\Y);
            %% get the betas as (voxel,condition,subject) matrix
            betas = cat(3,betas,(X\Yr)');            
        end
        filename = sprintf('%sver%s_id%s/R%02d_P%02d.png' , save_path , deblank(version_count), deblank(version_id), roi , phase );
        if exist(fileparts(filename)) == 0
            mkdir(fileparts(filename));
        end
        %make the MDS analysis and save the plot
        pa_MDS(betas,filename);
    end
end