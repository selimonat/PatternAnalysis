function SPM = PatternAnalysis_ConditionBased_fMRI_design(subject,phase,pattern,roi)
%SPM = PatternAnalysis_ConditionBased_fMRI_design(subject,phase,pattern,roi)
%
% Will return an first-level SPM array ready for manual estimation. This is
% the condition-based version of the
% PatternAnalysis_SingleTrial_fMRI_design function which generates many SPM
% structure arrays for single trial analysis.
%
% See also: PatternAnalysis_SingleTrial_fMRI_design

        fprintf('%s: Collecting the SPM structure arrays.\n', mfilename);
        
        p                 = '/Volumes/feargen2/feargen2/';
        files             = spm_select('FPList',sprintf('%s/data/sub%03d/phase%02d/mrt/nii/',p,subject,phase),pattern);
        k                 = length(files);%for future convienience.
        if isempty(files)
            fprintf('No files found ... \n');
            SPMs = [];
            return
        else
            fprintf('Here are few of the selected files:\n');
            display(files(1:10,:));
            fprintf('...\n');
        end
        load(GetSubjectData(subject,phase,'stimulation'));
        csp               = p.stim.cs_plus;
        dm                = GetDesignMatrix(subject,phase);
        dm                = [circshift(dm(:,1:8),[0 4-csp]) dm(:,9:10)];
        [onsets conds]    = find(dm);
        % only take face conditions and discard always the first trial
        onsets(1)         = [];
        conds(1)          = [];
        fprintf('Running Subject: %02d in Phase %02d, CSP: %02d, total of %03d stimuli\n',subject,phase,csp,length(onsets));
        %% which ROI? accordingly select a mask
        mask              = Atlas2NativeSpace(subject);
        mask              = mask(roi);
        %% get all the files
        xY.P              = files;
        xY.VY             = spm_vol(xY.P);
        xY.RT             = 2.02;
        if ~spm_check_orientations(xY.VY)
            fprintf('These images don''t have the same orientation\n');
            return
        end
        %% compute globals and normalize
        GM                = 100;
        q                 = length(xY.VY);
        g                 = zeros(q,1);
        for i = 1:q
            g(i) = spm_global(xY.VY(i));
        end
        gSF = repmat(GM./mean(g),k,1);
        for i = 1:q
            xY.VY(i).pinfo(1:2,:) = xY.VY(i).pinfo(1:2,:)*gSF(i);
        end
        %-place global variates in global structure
        xGX.rg          = g;
        xGX.GM          = GM;
        xGX.gSF         = gSF;
        xGX.iGXcalc     = 'None';
        xGX.sGXcalc     = 'mean voxel value';
        xGX.sGMsca      = 'session specific';
        %%
        TH          = g.*gSF*0.8;
        xM          = struct('T',    ones(q,1),...
            'TH',   TH,...
            'I',    0,...
            'VM',   mask,...
            'xs',   struct('Masking','analysis threshold'));
        %% intrinsic autocorrelations (Vi)
        xVi.V       = speye(k);
        xVi.form    = 'i.i.d';
        %% High-pass filtering and serial correlations
        % create and set filter struct
        HParam = 128;
        K = struct('HParam', HParam , ...
            'row',    1:k , ...
            'RT',     xY.RT);
        K = spm_filter(K);%will be placed in SPM.xX.K
        %% get the nuissance covariates.
        [Nuis]                = GetNuissanceVariables(subject,phase,pattern,{'mc_diff_square'});
%         tNuis                 = length(Nuis);
%         Nuis(tNuis + 1).val   = zscore(mean(PatternAnalysis_GetNativeData(subject,phase,roi,75),2));
%         Nuis(tNuis + 2).val   = zscore(mean(PatternAnalysis_GetNativeData(subject,phase,mod(roi-1,48)+1+48,75),2));
%         Nuis(tNuis + 1).name  = 'RH_csf';
%         Nuis(tNuis + 2).name  = 'LH_csf';
        
        %% Create a SPM structure for each single trial, time consuming things has to be done before this for loop
        SPM       = [];
        SPM.SPMid =  'SPM8: spm_fmri_spm_ui (v4178)';
        % for t = trials;%current trial
        %% place some useful information about trials
        SPM.external_trial_id = conds;
        SPM.external_onsets   = onsets;
        
        %% first place all what is needed to be placed
        SPM.xY   = xY;
        SPM.xGX  = xGX;
        SPM.xX.K = K;
        SPM.xM   = xM;
        SPM.xVi  = xVi;
        
        %%
        fMRI_T                = spm_get_defaults('stats.fmri.fmri_t');
        fMRI_T0               = spm_get_defaults('stats.fmri.fmri_t0');
        %% xBF structure containing the basis function information
        xBF.T                 = fMRI_T;
        xBF.T0                = fMRI_T0;
        xBF.dt                = SPM.xY.RT/xBF.T;
        xBF.UNITS             = 'scans';
        xBF.Volterra          = 1;
        xBF.name              = 'hrf';
        xBF                   = spm_get_bf(xBF);
        rep                   = 0;%replication of sessions
        %% Session field for 2 trials, no parametric modulator
        for i = 1:10;%two regressors (first one with one trial, and the other one with all the rest)
            Sess.U(i).dt        = xBF.dt;%- time bin (seconds)
            
            Sess.U(i).ons       = onsets( conds == i );%- onsets    (in SPM.xBF.UNITS)
            Sess.U(i).name      = {sprintf('%02d',i)};%- cell of names for each input or cause
            
            %no parametric modulation here
            Sess.U(i).dur    =  repmat(0,length(Sess.U(i).ons),1);%- durations (in SPM.xBF.UNITS)
            Sess.U(i).P.name =  'none';
            Sess.U(i).P.P    =  'none';
            Sess.U(i).P.h    =  0;%- order of polynomial expansion
            Sess.U(i).P.i    =  1;%- sub-indices of u pertaining to P
        end
        % I rely on spm_get_ons to create the stimulus matrix.
        % so we need a first merge to SPM structure here.
        SPM.xBF                 = xBF;
        SPM.nscan               = k;
        SPM.Sess                = Sess;
        SPM.Sess.U              = spm_get_ons(SPM,1);
        U                       = SPM.Sess.U; %for convinience
        %% get the design matrix and resample it
        Xx                      = [];
        Xb                      = [];
        Xname                   = {};
        Bname                   = {};
        % Convolve stimulus functions with basis functions
        [X,Xn,Fc]               = spm_Volterra(U,xBF.bf,xBF.Volterra);
        % Resample regressors at acquisition times (32 bin offset)
        X                       = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
        %% append the nuissance variables to the Session
        C.C                     = zscore([Nuis.val]);
        C.name                  = {Nuis.name};
        % Get nuissance variables
        X                       = [X spm_detrend(C.C)];
        Xn                      = [Xn C.name];
        % Session structure array
        SPM.Sess.C              = C;
        SPM.Sess.row            = (1:k);
        SPM.Sess.col            = 1:size(X,2);
        SPM.Sess.Fc             = Fc;
        % Append names
        for i = 1:length(Xn)
            Xname{end + 1}      = Xn{i};
        end
        %
        %% Confounds: Session effects
        B                       = ones(k,1);
        Bn                      = {'constant'};
        for i = 1:length(Bn)
            Bname{end + 1}      = Bn{i};
        end
        %% append into Xx and Xb
        Xx                      = X;
        Xb                      = B;
        % finished, now merge
        SPM.xX.X                = [Xx Xb];
        SPM.xX.iH               = [];
        SPM.xX.iC               = 1:size(Xx,2);
        SPM.xX.iB               = (1:size(Xb,2)) + size(Xx,2);
        SPM.xX.iG               = [];
        SPM.xX.name             = {Xname{:} Bname{:}};
        %%
        SPM.xsDes = struct(...
            'Basis_functions',      SPM.xBF.name,...
            'Number_of_sessions',   sprintf('%d',1),...
            'Trials_per_session',   sprintf('%-3d',length(onsets)),...
            'Interscan_interval',   sprintf('%0.2f {s}',SPM.xY.RT),...
            'High_pass_Filter',     sprintf('Cutoff: %d {s}',SPM.xX.K(1).HParam),...
            'Global_calculation',   SPM.xGX.sGXcalc,...
            'Grand_mean_scaling',   SPM.xGX.sGMsca,...
            'Global_normalisation', SPM.xGX.iGXcalc);
    end