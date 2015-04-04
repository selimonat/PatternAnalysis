function [datalr,w_r2lr,roiname,w_resmslr]=PatternAnalysis_BetaFiles2Matrix(what2do,thr,ttype,phase,folder,preprop,varargin)
%[datalr,w_r2lr,roiname]=PatternAnalysis_BetaFiles2Matrix(what2do,thr,ttype,phase,folder,preprop,varargin)
%[datalr,w_r2lr,roiname]=PatternAnalysis_BetaFiles2Matrix(what2do,thr,ttype,phase,folder,preprop,method,ROI)
%[datalr,w_r2lr,roiname]=PatternAnalysis_BetaFiles2Matrix(what2do,thr,ttype,phase,folder,preprop,method,ROI,SUBJECTS)
%
%WHAT2DO: cache or load, if CACHE: will cache the data to a handy matlab
%format on the SSD drive. ROIW_* images are also saved, which are the ResMS
%images. if LOAD: then specific ROI, subject and method of correction can
%be requested.
%
%THR: ROI binarization threshold
%
%TYPE, PHASE, FOLDER, PREPROP are necessary to find out where the SPM file
%is located.
%
%TTYPE: Second level Analysis folder, mainly useless for the results here,
%but necessary to load the SPM file
%PHASE, FOLDER and PREPROP: Self evident
%
%METHOD is used to CPC, common pattern correction, this is done on a
%subject by subject basis, as indicated by results that showed that CPs are
%specific to subjects. Options are mean, best, all.
%
%VARARGIN: If LOAD, then the first varargin should contain the ROI index to
%load the data from.
%
%Additionnaly the 2nd VARARGIN can contain the subject indices, to do the
%analysis only on a selected subset. %SUBJECTS is in reference with the 3rd
%dimension of the ROI matrix, and NOT the real index of the subject. Good
%subjects are 
%gs           = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
%
%Example for caching:
%PatternAnalysis_BetaFiles2Matrix('cache',30,'anova_cov_',4,'faces/chrf_derivs_00','swRealigned')


cache_path = sprintf('~/Documents/feargen_midlevel/%s_Thr%02d/%s/phase%02d/%s/%s/',mfilename,thr,ttype,phase,folder,preprop);
%
w_fname   = 'r2_fullmodel.img';
resms     = 'ResMS.img';
%
roiname = [];datalr = [];w_r2lr=[];w_resms = [];
if strcmp(what2do,'cache')
    if exist(cache_path) == 0;mkdir(cache_path);end;
    fprintf('Will load the single subject data and save the individual files ROI by ROI.\n');
    %%
    spm_path = sprintf('/Volumes/feargen2/feargen2/data/spm/secondlevel/%s/phase%02d/%s/%s/SPM.mat',ttype, phase,folder,preprop);
    load(spm_path)
    %Get the number of regressors per session on the original SPM
    a          = load(sprintf('%s/SPM.mat',fileparts(SPM.xY.P{1})));
    RegPerSess = length(a.SPM.Sess(1).col);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GET THE FILEPATHS AND ASSOCIATED INFORMATION
    %run through all the files and for each file gather the subject index,
    %condition index, phase index and the filename.
    for nfile = 1:length(SPM.xY.VY)
        numbers        = regexp( SPM.xY.VY(nfile).fname,'\d\d\d','match');
        i_sub(nfile)   = str2num(numbers{1});
        numbers        = regexp( SPM.xY.VY(nfile).fname ,'\d\d\d\d','match');
        i_cond(nfile)  = mod(str2num(numbers{1})-1,RegPerSess)+1;
        i_phase(nfile) = ceil(str2num(numbers{1})./RegPerSess);
        fname{nfile}   = SPM.xY.VY(nfile).fname;
        %get as well the residuals mean square
        R2{nfile}      = [fileparts(SPM.xY.VY(nfile).fname) '/' w_fname];
        ResMS{nfile}   = [fileparts(SPM.xY.VY(nfile).fname) '/ResMS.hdr'];
    end
    figure;plot(i_cond)
    %%
    [subject_list i]   = unique(i_sub);
    ts                 = length(subject_list);    
    %this will give us one ResMS image per subject        
    R2                 = R2(i);
    ResMS              = ResMS(i);
    %%
    for r = R2;
        if exist(r{1}) == 0
            fprintf('Need to compute the explained variance\n');
            load(regexprep(r{1},'r2_.*','SPM.mat'))
            PatternAnalysis_SPM2ExplainedVariance(SPM,'all',[]);            
        end
    end
        
    
    
    %% LOAD THE DATA IN A MEANINGFUL MATRIX: BRAIN [voxel,condition,subject]    
    sc    = 0;
    for ns = subject_list
        sc               = sc + 1;
        fprintf('Reading whole brain from subject %02d.\n',ns)
        i                = find((i_sub == ns));
        dummy            = spm_read_vols(spm_vol(strvcat(fname{i})));
        %
        %compute the datamatrix
        fprintf('Storing pattern from subject %02d.\n',ns);
        brain(:,:,sc)    = reshape(dummy,[91*109*91 size(dummy,4)]);
        %
        %load as well the explained variance images;
        R2{sc}           = regexprep(R2{sc},'projects','Volumes/feargen2');
        dummy            = spm_read_vols( spm_vol(R2{sc}));
        fprintf('loading ResMS from subject %02d.\n',ns);
        brain_r2(:,:,sc) = reshape(dummy,[91*109*91 size(dummy,4)]);
        
        %load as well the resms images
        ResMS{sc} = regexprep(ResMS{sc},'projects','Volumes/feargen2');
        dummy         = spm_read_vols( spm_vol(ResMS{sc}));
        
        fprintf('laoding ResMS from subject %02d.\n',ns);
        brain_resms(:,:,sc) = reshape(dummy,[91*109*91 size(dummy,4)]);
        
    end
    fprintf('Saving the brain data.\n',ns);
    filename     = sprintf('%sBrain_VOX_COND_SUB.mat',cache_path);
    save(filename,'brain');    
           
    
    %% STORE DATA PER ROI and LATERALITY
    for nroi = [2 31 33 34]
        for lr = 0:1
            
            fprintf('Storing pattern for ROI: %02d, Side: %02d.\n',nroi,lr);
            [~, roi]          = LoadFeargenAtlas(thr,nroi+45*lr);            
            mask              = roi.d;
            mask_i            = find(mask);
            %
            [xyz mat name maskold roiold] = GetOxfordRoi(1 , nroi, 50 );
            %extract the pattern from whole brain volume
            data              = [];
            for ns = 1:ts
                data(:,:,ns)  = brain(mask_i,:,ns);
            end
            %            
            %
            filename          = sprintf('%sROI_%02d_%01d',cache_path,nroi,lr);
            save(filename,'data');
            %% extract as well the R2 data.
            data              = [];
            for ns = 1:ts
                data(:,:,ns)  = brain_r2(mask_i,:,ns);
            end
            
            filename          = sprintf('%sROI_%02d_%01d_W%s',cache_path,nroi,lr,w_fname(1:end-4));
            save(filename,'data');
            %% extract as well the ResMS data.
            data              = [];
            for ns = 1:ts
                data(:,:,ns)  = brain_resms(mask_i,:,ns);
            end
            
            filename          = sprintf('%sROI_%02d_%01d_%s',cache_path,nroi,lr,resms(1:end-4));
            save(filename,'data');
        end
    end
    
elseif strcmp(what2do,'load')
    
    %
    if isempty(varargin)
        fprintf('With LOAD, at least one VARARGIN will be necessary..\n')
        return
    else
        method   = varargin{1};
        roi      = varargin{2};        
        subjects = varargin{3};        
    end
    %
    fprintf('Will load the data for Roi: %02d-%02d from disk.\n',roi)
    %
    for lr = 0:1
        %LOAD THE BOLD DATA
        data              = [];
        filename          = sprintf('%sROI_%02d_%01d.mat',cache_path,roi,lr);
        load(filename);
        %
        [~, atlas]        = LoadFeargenAtlas(thr,roi+45*lr);
        roiname           = atlas.name;
        roiname           = regexprep(roiname{1},'_[L,R]_','_');
        %%LOAD THE explained variance
        filename          = sprintf('%sROI_%02d_%01d_W%s',cache_path,roi,lr,w_fname(1:end-4));
        r2                = load(filename);
        
        %%LOAD THE RESMS
        filename          = sprintf('%sROI_%02d_%01d_%s',cache_path,roi,lr,'ResMS');
        resms             = load(filename);
        
        fprintf('Will only include the requested subjects.\n');
        data              = data(:,:,subjects);
        w_r2              = r2.data(:,:,subjects);        
        w_resms           = resms.data(:,:,subjects);
        
        %
        if nargin > 6
             data         = PatternAnalysis_CommonPatternCorrection(data,w_r2,method);
        end
        %create a cell array to output
        %{2} = Right;
        %{1} = Left;
        datalr{2-lr}      = data;
        w_r2lr{2-lr}      = w_r2;
        w_resmslr{2-lr}     = w_resms;
    end
end