%This script will return a similarity matrix of size
%[condition,condition,subject,group,roi,laterality,metric]
if ismac    
    addpath('/Users/onat/Documents/Code/Matlab/PatternAnalysis')
else
    addpath('/home/schenk/Documents/MATLAB/RSA')
    addpath('/home/schenk/Documents/MATLAB/RSA/project_cond')
    addpath('/home/schenk/Documents/MATLAB/RSA/globalfunctions')
    addpath('/home/schenk/Documents/MATLAB/RSA/PatternAnalysis')
end
%% initial variables
rois              = 1:cond_defaults('troi')/2;%all rois; 1:46
troi              = cond_defaults('troi')/2;%effective roi number 46
roi_constant      = [troi 0]; %[46 0]
vs                = cond_defaults('fov'); %[121 145 121]
project_path      = cond_defaults('project_path');%/where the data is '/projects/cond/'
analysis_path     = fileparts(which('cond_defaults'));%where the analysis scripts are located. '/home/schenk/Documents/MATLAB/RSA'
[Version, ID]     = GetGit(analysis_path);
group_context     = cond_defaults('group_context');%subject groups
group_nocontext   = cond_defaults('group_nocontext');
atlas             = cond_defaults('atlas');
threshold         = cond_defaults('threshold');
metrics           = cond_defaults('metrics');
vis               = 1;%visuals
start             = tic; % start clock % what is the start variable used for?

%% Get the beta file list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this part is still project specific, could be made more general in the
%future.
if ismac
    a = load(sprintf('%s15_04_16s4_f10_Mpain_ffMpains1_dummy/SPM.mat',project_path));
    b = load(sprintf('%s15_04_16s4_f10_Mpain_ffMpains2_dummy/SPM.mat',project_path));
else    
    a = load(sprintf('%s2ndlevel/15_04_21s4_f12_Mcue_ffMcue1_dummy/SPM.mat',project_path));
    b = load(sprintf('%s2ndlevel/15_04_21s4_f12_Mcue_ffMcue2_dummy/SPM.mat',project_path));
end

beta_list = [strvcat(a.SPM.xY.VY(:).fname); strvcat(b.SPM.xY.VY(:).fname)];
%adapt to ganglion
if ismac
    dummy = [];
    for n = 1:size(beta_list ,1)
        %replace path prefix and remove the folder name (because we didn't
        %copy it)
        dummy= [dummy; regexprep(regexprep(beta_list (n,:),'/projects','/Volumes/feargen2'),'15_04_15s4_f10_Mpain/','')];
    end
    beta_list = dummy;
    dummy = [];%garbage
end

%% Get useful information from the BETA_LIST
i_sub = [];i_cond=[];i_group = [];fname=[];
for nfile = 1:length(beta_list) %'numbers' is just a changing loop variable
    numbers        = regexp( cellstr(regexp( beta_list(nfile,:),'sub\d\d','match')) ,'[0-9]*','match');
    i_sub(nfile)   = str2double(numbers{1});
    numbers        = regexp( beta_list(nfile,:),'\d\d\d\d','match');
    i_cond(nfile)  = str2double(numbers{1})-2;
    i_group(nfile) = ismember(i_sub(nfile),group_nocontext);    % Info 0 = context, 1= nocontext: so that in raw context=1 and nocontext =2
end
% %sort everything by groups so that first 0s then 1s i
beta_list = [beta_list(i_group == 0,:); beta_list(i_group == 1,:)];
i_sub  = [i_sub(i_group == 0) i_sub(i_group == 1)];
i_cond = [i_cond(i_group == 0) i_cond(i_group == 1)];
i_group= [i_group(i_group == 0) i_group(i_group == 1)];

%sanity check
if vis
    figure(3);
    subplot(3,1,1);plot(i_sub,'o');subplot(3,1,2);plot(i_cond,'o');subplot(3,1,3);plot(i_group,'o')
end
%% Compute similarity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ts             = length(unique(i_sub))/2;%total subjects per group
tc             = length(unique(i_cond));%total conditions
for nroi       = rois;%run across all ROIs, and fill up the similarity matrix
    for lr = 1:2        
        %Get roiname and voxel indices
        [dummy, mask,xyz]       = pa_GetAtlas(atlas,threshold,nroi+roi_constant(lr));
        tvox                    = sum(mask.d(:));%total number of voxels
        fprintf('ROI %02d (%s): Side: %01d, %s...\n',nroi,mask.name{1},lr,datestr(now,'HH:MM:SS'));
        for groups = {group_context group_nocontext}%run across two groups
            subs = groups{1};%vectorize the cell
            s_c  = 0;%counter init
            for ns = subs  
                fprintf('.')
                s_c = s_c + 1;                
                dummy = beta_list((i_sub == ns),:);%36 conditions for subject NS
                %store the matrix once as raw and then after common pattern
                %correction
                gr                 = unique(i_group(i_sub == ns)) + 1;%+1 to avoid zero indexing
                roi.raw            = spm_get_data(spm_vol(dummy),xyz)';
                %clean from voxels that do not change across conditions
                nz                 = ~any(diff(roi.raw'));
                removed_zeros(s_c) = sum(nz); % nr of voxels with 0 change over conditions
                roi.raw(nz,:)      = [];
                %motion correction
                common_pattern     = mean( roi.raw ,2);%average across conditions (but not across sub or group)
                roi.mc             = roi.raw - repmat(common_pattern,[1 tc 1 1]);
                % at this stage we have a cell array ROI, which stores beta values for
                % all conditions, subjects and groups. Left/Right has to be distributed
                % to different cells because the number of voxels is different, so it
                % can not be stored in a matrix.
                %         roi.raw is:
                %         [voxel,condition,subject,group]
                %         group 1: context
                %         group 2: nocontext
                %% compute the similarity using different metrics
                for dtype = fieldnames(roi)';
                    data     = roi.(dtype{1}); % data of one sub
                    c_metric = 0;
                    for metric = metrics                        
                        smat = [];
                        if any(strcmp(metric{1},{'cov' 'mean'}))
                            fh       = str2func(metric{1});%get the function handle
                            smat     = fh(data); % applies function to data
                            if isvector(smat)%this is necessary to store the mean vector
                                smat = diag(smat); %puts means into a diagonal
                            end
                        elseif any(strcmp(metric{1},{'correlation' 'euclidean' 'seuclidean'}))
                            try %not all machines have toolbox for pdist
                                smat     = squareform(pdist(data',metric{1})); % i do not get this line (and its not working)
                            end
                        end
                        if ~isempty(smat)
                            %store
                            c_metric             = c_metric + 1;
                            sim.metric{c_metric} = metric{1};
                            sim.(dtype{1})(:,:,s_c,gr,nroi,lr,c_metric) = smat;
                        end
                        
                    end
                end             
            end
            %display results on removed voxels for sanity.
            fprintf('%3g (of %3g) voxels were removed.\n',mean(removed_zeros),length(nz));
            if mean(removed_zeros) ~= 0
                fprintf('\tfrom subjects: ');fprintf('%g, ',find(removed_zeros));fprintf('\n');
            end
        end
        fprintf('\n');
    end
    sim.name = mask.name{1};%store the name
end
fprintf('Analysis finished in %.2g minutes.\n',toc/60);



