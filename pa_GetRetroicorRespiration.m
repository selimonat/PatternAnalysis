function [regressor]=RetroicorRegressor_Respiration(subject,phase)
%[regressor]=RetroicorRegressor_Respiration(subject,phase)
%   
%   Assigns to each pulse a phase value based on peak to peak distances
%   between hearth beats. This is based on Glover GH et al., Magnetic
%   Resonance in Medicine, 2000.
%
%   Dependency: pa_GetBOLD, pa_GetPhysioData
%



%% load the data, downsample FACTOR times, get the pulse times.
factor          = 8;%downsampling factor
rw              = 10;%running average window in samples
addon           = 5;%i seconds, it includes the data plus minus ADDON seconds from the first and last pulse time on 
nbin            = 100;
tOrder          = 3;
dummy_scans     = 6;
%
volume_files    = pa_GetBOLDFiles(subject,phase,'^fTRIO.*nii$');
spike2          = pa_GetPhysioData(subject,phase);
sampling_period = spike2.Ch3.interval*factor;
%
pulse_times     = spike2.Ch4.times;
%
r_ori           = spike2.Ch3.values;
%invert the signal so that peaks are peaks inhalation moments
r_ori           = max(r_ori) - r_ori;
r_ori_time      = [(0:length(r_ori)-1).*sampling_period/factor]';
%
[r_ori K]       = HighPassFilter(r_ori, range(pulse_times)*10, sampling_period);
% figure;imagesc(K.X0);
%downsample
r_downs         = decimate(r_ori,factor);
r_downs_time    = [(0:length(r_downs)-1).*sampling_period]';
%smooth
r_downs         = conv(r_downs,ones(1,rw)./rw,'same');
%
%%
%detect window of interest
i               = (r_downs_time > (pulse_times(1)-addon)) & (r_downs_time < (pulse_times(end)+addon));
r               = r_downs(i);
r_time          = r_downs_time(i);
% normalize the range in the temporal window of interest to [0 1]
r               = (r - min(r));
r               = (r./max(r));
%take the derivative
r_diff          = [diff(r);NaN];
%
y               = hist(r,nbin);
transfer        = [ 0 cumsum(y)]'./sum(y);
%
phases          = pi*transfer(round(r*100)+1).*sign(r_diff);
%%
% for each pulse gather a phase value based on inter-heart-beat position of
% the pulse.
pulse_samples   = [];
pulse_phase     = [];
for n = 1:length(volume_files);
    pulse_i                  = n + dummy_scans;%pulse index
    %find the sample that is closest to this pulse
    [mini i]            = min(abs(pulse_times(pulse_i) - r_time));        
    pulse_samples(n,1)  = i;
    pulse_phase(n,1)    = phases(i);
end
%%
%gather the regressors
regressor = [];
for nOrder = 1:tOrder
    regressor = [regressor cos(nOrder*pulse_phase) sin(nOrder*pulse_phase)];
end
name = repmat({mfilename},1,size(regressor,2));
%%
%plot the stuff
% rz      = zscore(r);
% phasesz = zscore(phases(1:end-1));%last value is not defined
% figure
% plot(r_time, rz, 'k');
% hold on
% plot(r_ori_time, zscore(r_orihp), 'b');
% plot(r_time(1:end-1), phasesz, 'r');
% plot(r_time(pulse_samples), phasesz(pulse_samples),'mo');
% hold off



