function [sim]=Analysis_FMRI_Multivariate2(force)

rois      = 1:46;%all rois;
troi      = length(rois);
validc    = [1:8 11:18];%will be used to exclude ucs & oddball conditions
load /Volumes/feargen2/feargen2/data/midlevel/selectedsubjects.mat
gs              = [2 4 7 8 9 11 12 14 16 18 20 21 24 27];
subject_list    = subject_list(gs);
ts              = length(subject_list);%total subject
tc              = 10;%total conditions
vs              = [91 109 91];
%     load /Volumes/feargen2/feargen2/data/spm/secondlevel/chrf_derivs_00/swRealigned/SPM.mat
load /Volumes/feargen2/feargen2/data/spm/secondlevel/anova_cov_/phase04/faces/chrf_derivs_00/swRealigned/SPM.mat
%%
if exist(roi_fname) == 0 | force
    %%
    % function [brain]=GetBrain(SPM,subjectlist,loadcache);
    brainrecompute  = 0;
    i_sub = [];i_cond=[];i_phase = [];fname=[];
    for nfile = 1:length(SPM.xY.VY)
        numbers        = regexp( SPM.xY.VY(nfile).fname,'\d\d\d','match');
        i_sub(nfile)   = str2num(numbers{1});
        numbers        = regexp( SPM.xY.VY(nfile).fname ,'\d\d\d\d','match');
        i_cond(nfile)  = mod(str2num(numbers{1})-1,30)+1;
        i_phase(nfile) = ceil(str2num(numbers{1})./30);
        fname{1,nfile} = SPM.xY.VY(nfile).fname;
    end
    %%
    %create ucs and odd data
    fname   = [fname   regexprep(fname(i_phase == 1 & i_cond == 8),'0008','0009') regexprep(fname(i_phase == 1 & i_cond == 8),'0008','0010') regexprep(fname(i_phase == 2 & i_cond == 8),'0038','0039')   regexprep(fname(i_phase == 2 & i_cond == 8),'0038','0040')]';
    i_sub   = [i_sub   repmat(subject_list,1,4)];
    i_cond  = [i_cond  ones(1,ts)*9 ones(1,ts)*10 ones(1,ts)*9 ones(1,ts)*10];
    i_phase = [i_phase ones(1,ts*2)*1 ones(1,ts*2)*2];
    %% load the whole brain in a matrix [voxel,face,subject,phase]
    brain   = zeros(prod(vs),tc,ts,2);%init the brain
    sub      = 0;
    for ns = subject_list
        sub = sub + 1;
        for ph = [1 2]
            fprintf('Collecting pattern from subject %02d phase %02d.\n',ns,ph)
            i                 = find((i_sub == ns).*(i_phase == ph));
            dummy             = spm_read_vols(spm_vol(strvcat(fname{i})));
            %compute the datamatrix
            brain(:,:,sub,ph) = reshape(dummy,[prod(vs) size(dummy,4)]);
        end
    end
    %I dont want to save BRAIN
    %% function [roi]=brain2ROI(brain);
    fprintf('Extracting ROIs:\n')
    roi_fname = sprintf('%s/midlevel/SimilarityPerROI_common1_16.mat',GetProjectRootFeargen2);
    roi            = [];
    roi_constant   = [troi 0];
    ts             = size(brain,3);%total subject
    tc             = size(brain,2);%total conditions
    for nroi       = rois;
        fprintf('%d ',nroi)
        for lr = 1:2
            %store left and right data
            [~, mask]           = LoadFeargenAtlas(50,nroi+roi_constant(lr));
            tv                  = sum(mask.d(:));
            %store the name
            roi{nroi,lr}.name   = mask.name{1};
            %brain is [voxel,face,subject,phase]
            %raw   is [voxel,face,subject];
            %   phase dimension is collapsed into one single to ease
            %   computation of the similarity matrices
            %   the VX, space dimension is distributed on different ROIs as
            %   cells. The reason we have to cell-split the data is that at
            %   this point the ROIs have different voxel sizes. However
            %   when we compute the similarity metrics all the data will be
            %   stored on a huge matrix.
            raw                 = [brain(mask.d(:),1:8,:,1) brain(mask.d(:),1:8,:,2) brain(mask.d(:),9:10,:,1)  brain(mask.d(:),9:10,:,2)];
            %
            roi{nroi,lr}.d.raw  = raw;
            %get the MC (this is subject specific, we don't like fanciness, we
            %make a for loop
            for ns = 1:ts
                common_pattern              = mean( raw(:,1:16,ns) ,2);
                roi{nroi,lr}.d.mc(:,:,ns)   = raw(:,:,ns) - repmat(common_pattern,[1 tc*2 1]);
                %make a more sophisticated MC, find the best fitting common pattern
                %and subtract it
                beta                        = lscov([common_pattern ones(tv,1)], raw(:,:,ns));
                best_common_pattern         = [common_pattern ones(tv,1)]*beta;
                %
                roi{nroi,lr}.d.bmc(:,:,ns)  = [raw(:,:,ns) - best_common_pattern];
            end
        end
    end
    fprintf('\n');
    save(roi_fname,'roi');
else
    fprintf('Loading ROI from cache\n');
    load(roi_fname);
    ts             = size(roi{1}.d.raw,3);
end

%% compute the similarity metrics
metric          = 0;
fprintf('Extracting ROIs:\n')
for nroi = rois
    fprintf('%02d ',nroi)
    sim.roi_name{nroi}                        = roi{nroi,1}.name;
    for lr = 1:2
        for dtype = fieldnames(roi{nroi,lr}.d)'
            data                              = roi{nroi,lr}.d.(dtype{1});
            data_combined                     = reshape( data , [size(data,1)*2 size(data,2)/2 size(data,3)]);
            for sub = 1:ts
                %FIRST ORDER STUFF, simply the mean
                metric                                              = 5;
                sim.d.(dtype{1})(:,:,sub,nroi,lr,metric)            = diag(mean(data(:,:,sub)));
                sim.metric{metric}                                  = 'mean';
                %I donnow what is this cross
                %sim.cross.(dtype{1})(:,:,sub,nroi,lr,metric)       = cov(data_combined(:,:,sub));
                %sim.metric{metric}                                  = 'cov';
                %COVARIANCE
                metric                                              = 1;
                sim.d.(dtype{1})(:,:,sub,nroi,lr,metric)            = cov(data(:,:,sub));
                %sim.cross.(dtype{1})(:,:,sub,nroi,lr,metric)        = cov(data_combined(:,:,sub));
                sim.metric{metric}                                  = 'cov';
                %CORRELATION
                metric                                              = 2;
                sim.metric{metric}                                  = 'corr';
                sim.d.(dtype{1})(:,:,sub,nroi,lr,metric)            = squareform(pdist(data(:,:,sub)','correlation'));
                %EUCLIDEAN
                metric                                              = 3;
                sim.metric{metric}                                  = 'euclidian';
                sim.d.(dtype{1})(:,:,sub,nroi,lr,metric)            = squareform(pdist(data(:,:,sub)','euclidean'));
                %SCALED EUCLIDEAN
                metric                                              = 4;
                sim.metric{metric}                                  = 's-euclidian';
                sim.d.(dtype{1})(:,:,sub,nroi,lr,metric)            = squareform(pdist(data(:,:,sub)','seuclidean'));
            end
        end
    end
end
fprintf('\n');
%% Extract the interesting stuff from the covariance matrices
for nroi = rois;
    for lr = 1:2
        %FG1: simply the difference of the CS+ column between B and T.
        %the subject-wise difference of two matrices
        %extract the 4th column
        fg = 1;
        data                               = squeeze(sim.d.bmc(4,1:8,:,nroi,lr) - sim.d.bmc(12,1:8,:,nroi,lr))';
        sim.result(fg).bmc.mean(:,nroi,lr) = mean(data);
        sim.result(fg).bmc.std(:,nroi,lr)  = std(data)./sqrt(size(data,1));
        sim.result(fg).bmc.title           = 'CovDiff';
        
        %FG2: The forth column of th Cxy matrix
        fg = 2;
        data                               = squeeze(sim.d.bmc(1:8,12,:,nroi,lr))' ;
        sim.result(fg).bmc.mean(:,nroi,lr) = mean(data);
        sim.result(fg).bmc.std(:,nroi,lr)  = std(data)./sqrt(size(data,1));
        sim.result(fg).bmc.title           = 'C_xy4thcolumn';
        
        %FG = 3; The sum of the C_xy matrix
        fg = 3;
        data                               = squeeze(sum(sim.d.bmc(1:8,9:16,:,nroi,lr)))' ;
        sim.result(fg).bmc.mean(:,nroi,lr) = -mean(data);
        sim.result(fg).bmc.std(:,nroi,lr)  = std(data)./sqrt(size(data,1));
        sim.result(fg).bmc.title           = 'C_xySummed';
    end
end
%% plot all the matrices
fs = 8;%fontsize;
savepath = '~/Pictures/results_rsa/SimilarityMatrices/';
mkdir(sprintf('%s',savepath));
[col row] = GetSubplotNumber(troi);
figure(1)
set(1,'position',get(0,'ScreenSize'));
for colormap_equalize = [0 1]
    for mc_type   = {'mc'};%{'raw' 'mc' 'bmc'};
        for metric    = 1:4;
            figure(1);clf
            [d u]     = GetColorMapLimits(Vectorize(sim.d.(mc_type{1})(1:16,1:16,:,:,:,metric)),2);
            for nroi = rois
                h=subplot(row,col,nroi);
                subplotChangeSize(h,.009,.009);
                if colormap_equalize
                    imagesc([mean(sim.d.(mc_type{1})(1:16,1:16,:,nroi,1,metric),3) mean(sim.d.(mc_type{1})(1:16,1:16,:,nroi,2,metric),3)],[d u])
                else
                    imagesc([mean(sim.d.(mc_type{1})(1:16,1:16,:,nroi,1,metric),3) mean(sim.d.(mc_type{1})(1:16,1:16,:,nroi,2,metric),3)])
                end
                axis off
                axis image;
                hold on
                plot([16 16]+.5,[1 16]-.5,'k')
                hold off
                title(regexprep(sim.roi_name{nroi},'_[L,R]_','_'),'interpreter','none','fontsize',fs);
                h = thincolorbar('vertical');
                set(h,'fontsize',fs)
            end
            supertitle(sprintf('Metric: %s, MC type: %s',sim.metric{metric},mc_type{1}),1,'fontsize',fs)
            SaveFigure(sprintf('%scmeq_%g_%s_%s.png',savepath,colormap_equalize,sim.metric{metric},mc_type{1}),'-r200')
        end
    end
end
close all;
% % % % %% now we plot all the individual components of these matrices
% % % % close all;
% % % % savepath  = '~/Pictures/results_rsa/FGprofiles/';
% % % % mkdir(sprintf('%s',savepath));
% % % % troi      = size(sim.d.bmc,4);
% % % % [col row] = GetSubplotNumber(troi);
% % % % ts = size(sim.d.bmc,3);%total subjects
% % % % for nr = 1:length(sim.result);
% % % %     figure;
% % % %     for mc_type   = {'bmc'};%{'raw' 'mc' 'bmc'};
% % % %         for metric = 1;
% % % %             %[d u]     = GetColorMapLimits(Vectorize(data),2);
% % % %             spc       = 0;%subplot counter...
% % % %             for nroi  = 1:troi
% % % %                 spc = spc + 1;
% % % %                 h   = subplot(row,col,spc);
% % % %                 for lr = 1:2;
% % % %                     fprintf('%02d ',nroi);
% % % %                     bar((1:8)+8*(lr-1),sim.result(nr).(mc_type{1}).mean(:,nroi,lr),1,'k','linestyle','none');
% % % %                     hold on;
% % % %                     errorbar((1:8)+8*(lr-1),sim.result(nr).(mc_type{1}).mean(:,nroi,lr),sim.result(nr).(mc_type{1}).std(:,nroi,lr),'bo');
% % % %                 end
% % % %                 title(regexprep(sim.roi_name{nroi},'_[L,R]_','_'),'interpreter','none');
% % % %                 axis tight;
% % % %                 box off;
% % % %                 set(gca,'xticklabel','');
% % % %                 plot([9 9]-.5,ylim,'r');
% % % %             end
% % % %             set(gcf,'position',get(0,'ScreenSize'));
% % % %             supertitle(sim.result(nr).(mc_type{1}).title,1,'interpreter','none');
% % % %             SaveFigure(sprintf('%scmeq_%g_%s_%s_%s.png',savepath,0,sim.metric{metric},mc_type{1},sim.result(nr).(mc_type{1}).title),'resolution',200);
% % % %             %calibrate the y-axis
% % % %             %set(findall(2,'type','axes'),'ylim',[d u]);
% % % %         end
% % % %     end
% % % % end
% % % % %% now we have all the similarity matrices, we would like to evaluate
% % % % %  how similar are matrices between Baseline and Test Phases for different
% % % % %  regions of interest... This analysis will result in exclusion of some
% % % % %  areas, and quantification of left right differences.
% % % % nroi = 1;
% % % % B    = mean(sim.d.bmc(1:8,1:8,:,nroi,1,1),3);
% % % % T    = mean(sim.d.bmc(1:8,1:8,:,nroi,2,1),3);
% % % %
% % % %
% % % %
% % % %
% % % % %% get common ylim values for each ROI
% % % % mc_type     = 'bmc';
% % % % sim.bmc.ylim = [];
% % % % f = 4;
% % % % for metric = 1:4;
% % % %     for nroi = 1:troi
% % % %         data                   = sim.d.bmc(1:16,1:16,:,:,:,metric);
% % % %         m                      = mean(Vectorize(data(1:8,4,1,:)));
% % % %         s                      = std(Vectorize(data(1:8,4,1,:)));
% % % %         sim.bmc.ylim(metric,:) = [m-f*s  m+f*s];
% % % %         %for the diffs
% % % %         data                   = sim.d.bmc(1:8,1:8,:,:,:,metric) - sim.d.bmc(9:16,9:16,:,:,:,metric);
% % % %         m                      = mean(Vectorize(data(1:8,4,1,:)));
% % % %         s                      = std(Vectorize(data(1:8,4,1,:)));
% % % %         sim.bmc.ylim_diff(metric,:) = [m - f*s m+f*s]
% % % %         %for the off diagonales
% % % %     end
% % % % end
% % % %
% % % %
% % % % %% plot
% % % % savepath = '~/Pictures/results_rsa/GeneralizationGradients/';
% % % % mkdir(sprintf('%s',savepath));
% % % % x           = deg2rad(linspace(-135,180,8));
% % % % row         = 3;
% % % % mc_type     = 'bmc';
% % % % tmetric     = size(sim.d.bmc);tmetric = tmetric(end);
% % % % col         = tmetric*3;
% % % % laterality  = {'L' 'R'};
% % % %
% % % % metric = 1;
% % % % nroi   = 1;
% % % % [suby,subx,w,h] = SubPlotCoor(row,col);
% % % % lrshift=[0 .009 0.009];
% % % % set(gcf,'position',[0 0 1438 595]);
% % % % for nroi = 1:troi
% % % %     clf
% % % %     for metric = 1:4
% % % %         for lr = 1:2
% % % %             subplot('position',[subx(lr+(metric-1)*3)+lrshift(lr) suby(1) w h]);
% % % %             hold on
% % % %             data   = sim.d.(mc_type)(1:16,1:16,:,nroi,lr,metric);
% % % %             m      = mean(data,3);
% % % %             s      = std(data,1,3)./sqrt(ts);
% % % %             %%the first row
% % % %             plot(x,m(1:8,4),'k')
% % % %             errorbar(x,m(1:8,4),s(1:8,4),'ko');
% % % %             plot(x,m(9:end,12),'r')
% % % %             errorbar(x,m(9:end,12),s(1:8,4),'ro');
% % % %             hold off
% % % %             axis tight
% % % %             box off
% % % %             ylim(sim.bmc.ylim(metric,:));
% % % %             title(sprintf('%s',laterality{lr}))
% % % %             if lr == 1;ylabel(sim.metric{metric});end
% % % %             set(gca,'color','none')
% % % %             %%differences
% % % %             subplot('position',[subx(lr+(metric-1)*3)+lrshift(lr) suby(2) w h]);hold on
% % % %             data   = sim.d.(mc_type)(9:16,9:16,:,nroi,lr,metric) - sim.d.(mc_type)(1:8,1:8,:,nroi,lr,metric);
% % % %             m      = mean(data,3);
% % % %             s      = std(data,1,3)./sqrt(ts);
% % % %             ylim(sim.bmc.ylim_diff(metric,:));
% % % %             %
% % % %             plot(x,m(1:8,4),'k');
% % % %             errorbar(x,m(1:8,4),s(1:8,4),'ko');
% % % %             axis tight
% % % %             box off
% % % %             ylim(sim.bmc.ylim_diff(metric,:));
% % % %             set(gca,'color','none')
% % % %             if lr == 1;ylabel(['Diff ' sim.metric{metric}]);end
% % % %
% % % %             %% plot the matrices
% % % %             subplot('position',[subx(lr+(metric-1)*3)+lrshift(lr) suby(3)-0.05 w h]);hold on;cla
% % % %             axis off;axis ij
% % % %             data = [ShiftCorrMat(mean(sim.d.(mc_type)(1:8,1:8,:,nroi,lr,metric),3));ShiftCorrMat(mean(sim.d.(mc_type)(9:16,9:16,:,nroi,lr,metric),3))];
% % % %             imagesc(data);
% % % %             %%
% % % %             if lr == 2
% % % %                 lr = 3;
% % % %                 subplot('position',[subx(lr+(metric-1)*3)+lrshift(lr) suby(3) w.*0.85 h.*0.85]);hold on;
% % % %                 m      = mean(data(1:8,1:8));
% % % %                 s      = std(data(1:8,1:8))./sqrt(ts);
% % % %                 plot(x,m,'k');
% % % %                 errorbar(x,m,s,'ko');
% % % %                 m      = mean(data(9:16,1:8));
% % % %                 s      = std(data(9:16,1:8))./sqrt(ts);
% % % %                 plot(x,m,'r');
% % % %                 errorbar(x,m,s,'ro');
% % % %                 axis tight;box off;set(gca,'color','none')
% % % %                 ylim(sim.bmc.ylim(metric,:));
% % % %             end
% % % %         end
% % % %     end
% % % %     supertitle(regexprep(sim.roi_name{nroi},'_[L,R]_','_'),1,'interpreter','none')
% % % %     SaveFigure(sprintf('%s%s_%s_%s.png',savepath,regexprep(sim.roi_name{nroi},'_[L,R]_','_'),sim.metric{metric},mc_type),'resolution',200)
% % % % end
% % % %


