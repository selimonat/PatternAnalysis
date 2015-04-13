%%
cd /Users/onat/Documents/Code/Matlab/PatternAnalysis;
save_path         = '/Users/onat/Pictures/PatternAnalysis/';
HParam            = 128;%high pass filter parameter
threshold         = 75;
pattern           = '^swRealigned.*nii$';
gs                = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
load /Volumes/feargen2/feargen2/data/midlevel/selectedsubjects.mat%
%get the version of number of the analysis
[~,version_count] = system('git rev-list --count --first-parent HEAD')
[~,version_id]    = system('git describe --always')
for phase     = [2 4];
    for roi       = 1:90;        
        betas     = [];
        for subject = subject_list(gs);
            %%
            % get data
            Y = pa_GetY(subject,phase,roi,threshold,pattern);
            if isempty(Y) == 1
                % it might be possible that at this threshold level Y has
                % no valid voxels. 
                fprintf('No data for ROI: %02d, subject:%02d, phase: %02d\n',roi,subject,phase);
                break;
            end
            % remove mean
            Y     = Y - repmat(mean(Y),size(Y,1),1);
            % get the design matrix
            X     = pa_GetDesignMatrix(subject,phase);
            beta  = pa_GetBetas(X,Y);            
            %% get the betas as (voxel,condition,subject) matrix
            betas = cat(3,betas,beta);
        end
        filename = sprintf('%sver%s_id%s/R%02d_P%02d.png' , save_path , deblank(version_count), deblank(version_id), roi , phase );
        if exist(fileparts(filename)) == 0
            mkdir(fileparts(filename));
        end        
        %make the MDS analysis and save the plot
        pa_MDS(betas,filename);
    end
end