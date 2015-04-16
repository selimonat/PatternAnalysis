function [dm, onsets]= pa_GetDesignMatrix(subject,phase)
%[DM, onsets,dm_spm]= GetDesignMatrix(subject,phase)
%
%   Returns a binary design matrix based on the pulses recorded at the
%   physiocomputer by the CED. DM has the same size as the number of
%   volumes recorded and conditions. If the DM size is different than the
%   number of recorded volumes, it enters debugging mode. The result is
%   cached.
%   
%   Dependency: pa_GetPhysioData, pa_defaults, pa_GetRoot
%
%   see also: pa_GetSPMDesignMatrix


dm     = [];
onsets = [];
save_path    = sprintf('%ssub%03d/phase%02d/midlevel/%s.mat',pa_GetRoot,subject,phase,mfilename);
%load from cache or recompute it
if exist(save_path) == 0 ;
    TR           = pa_defaults('RT');
    spike2       = pa_GetPhysioData(subject,phase);
    %get the number of presented trials
    load(pa_GetSubjectDataPath(subject,phase,'stimulation'));
    cond_ind     = p.presentation.cond_id;
    tStim        = length(cond_ind);
    tUCS         = sum(p.presentation.cond_id == 9);
    %get the times
    stim_times   = spike2.Ch6.times;
    ucs_times    = spike2.Ch8.times;
    pulse_times  = spike2.Ch4.times;%pulse sanity check is done in Spike2Matlab
    %discard stim times which are before the first pulse
    stim_times(stim_times < pulse_times(1)) = [];
    ucs_times(ucs_times < pulse_times(1))   = [];
    %and take the first 293 or 123, basically discard the rating stimulus onsets
    stim_times   = stim_times(1:tStim);
    %make the occurence of the 7th volume zero time
    if subject == 31 && phase == 2        
        stim_times   = stim_times - pulse_times(29);
        ucs_times    = ucs_times  - pulse_times(29);
    else
        stim_times   = stim_times - pulse_times(7);
        ucs_times    = ucs_times  - pulse_times(7);        
    end
    %now replace CSP onsets with the UCS onset, we would like to have GLM that
    %models the UCS onset not the CSP onset.
    i                                   = [1 ; find(diff(ucs_times) > 1)+1 ];
    %TAke the first tUCS entries. in Subject 19, the last input in UCS_TIMES is
    %wrong, and counted as a UCS, which causes an error. on line 45
    i                                   = i(1:tUCS);
    stim_times(p.presentation.ucs == 1) = ucs_times(i);
    % % see the result
    % % figure;clf;
    % % plot(stim_times_,'bo-')
    % % hold on;
    % % plot(stim_times,'ro-')
    % % hold off
    
    %stim_times to volume index
    volume_index     = floor(stim_times./TR)+1;
    volume_index_dec = (stim_times./TR)+1;
    
    %create a design matrix that has the same size as the number of recorded volumes.
    volume_files = pa_GetBOLDFiles(subject,phase,'^fTRIO.*nii$');
    if isempty(volume_files)
        keyboard;
    end
    dm           = zeros(length(volume_files),10);
    %fill it in with ones where there is an event
    for ntrial = 1:length(cond_ind)
        dm(volume_index(ntrial),cond_ind(ntrial)) = 1;
    end
    %
    %prepare a cell array of onsets for each condition
    for ncond = unique(cond_ind(:)')
        onsets{ncond} = volume_index_dec(cond_ind == ncond);
    end
    save(save_path,'onsets','dm');
else
    load(save_path);
    
end
