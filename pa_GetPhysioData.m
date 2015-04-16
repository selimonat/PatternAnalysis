function [spike2]=pa_GetPhysioData(subject,phase)
%[spike2]=pa_GetPhysioData(subject,phase)
%
%   Reads in Matlab the Spike2 recorded data (exported as matlab format).
%   Field names are made human readable and compatible across different
%   sessions of recordings.
%
%   Dependency: pa_GetSubjectDataPath
%
%   See also: SCR_Spike2Ledalab
%

%%
fprintf('%s: Reading subject %d, phase %d\n',mfilename,subject,phase);
scrfile   = pa_GetSubjectDataPath(subject,phase,'scr');
p         = load(pa_GetSubjectDataPath(subject,phase,'stimulation'));
dummy     = p.p.mrt.dummy_scan;%number of dummy scans

if exist(scrfile) ~= 0
    
    % load the data
    spike2    = load(scrfile);
    spike2    = Spike2_BeautifyFieldnames(spike2);
    
    %% Clean The pulses        
    p           = spike2.Ch4.times(:);
    %p = [2 3 4 5 9 (5:10)*2.02 229 ]
    % figure;
    % plot(diff(p),'o-');
    tVolumes    = length(ListFiles(regexprep(scrfile,'scr/data.mat','mrt/nii/fTRIO*')));
    tPulses     = length(p);
    
    
    %1/ Remove the first/last pulses that might occur before/after recordings starts/stops,
    %After step 1, the first pulse should be REALLY the pulse of the first dummy scan
    d           = diff(p);
    mid         = round(length(p)./2);
    %start from the middle and go both directions until you hit a too late pulse
    pulse_validity= zeros(1,length(p));
    for c = mid:tPulses-1
        if d(c) < 2.025
            pulse_validity(c+1) = 1;%if diff at point c is TRUE than c+1th pulse is valid
        else
            break
        end
    end
    for c = (mid-1):-1:1
        if d(c) < 2.025
            pulse_validity(c+1) = 1;
        else
            break
        end
    end
    %the one before the valid is always valid.
    pulse_validity(find(pulse_validity==1,1)-1)=1;
    
    % %take only the stuff between stop and start
    p                 = spike2.Ch4.times(pulse_validity == 1);
    d                 = diff(p);
    
    %% 2/ Remove pulses that come too quick
    invalid           = d < 2.015;
    revalid           = find( diff(d < 2.015) == -1);
    invalid(end+1)    = 0;
    invalid(revalid)  = 0;
    %
    p(find(invalid==1)+1) = [];
    %cleaned pulses!
    spike2.Ch4.times    = p;
    spike2.Ch4.length   = length(p);
    % figure
    % plot(diff(p),'o-')
    
    After = length(p);
    if (After - tVolumes-dummy) > 1
        fprintf('Problem with pulse detection:\n Number of Volumes do not match number of pulses...\n');
        keyboard
    else
        fprintf('Pulse detection succeeded!\n%d volumes recorded\n %d pulses recorded\n %d final pulses\n',tVolumes+dummy,tPulses,After)
    end
    
    
else
    fprintf('scrfile: %s doesn''t exist\n', scrfile);
    spike2 = [];
end

