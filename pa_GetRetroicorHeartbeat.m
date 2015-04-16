function [regressor,name]=pa_GetRetroicorHeartbeat(subject,phase)
%[regressor,name]=pa_GetRetroicorHeartbeat(subject,phase)
%   
%   Assigns to each pulse a phase value based on peak to peak distances
%   between hearth beats. HIGHPASS is 1 if you like to high-pass filter the
%   data, however this is not really always beneficial, even maybe
%   detrimentral.
%
%   This is based on Glover GH et al., Magnetic Resonance in Medicine, 2000.
%
%   Dependency: pa_GetBOLD, pa_GetPhysioData

%order of the regressors
tOrder          = 3;
%load the data, downsample FACTOR times, get the pulse times.
factor          = 4;%so we are effectively working on 50 Hz.
dummy_scans     = 6;
%
%%get the total number of volume files so we can limit the size of the
%%regressors to the design matrix.
volume_files    = pa_GetBOLDFiles(subject,phase,'^fTRIO.*nii$');
t_volume_files  = length(volume_files);
%
spike2          = pa_GetPhysioData(subject,phase);
sampling_period = spike2.Ch2.interval*factor;
%
h               = spike2.Ch2.values;%pulseoxy meter
h               = decimate(h,factor);
h_time          = (0:length(h)-1).*sampling_period;
%
pulse_times     = spike2.Ch4.times;%scanner pulses

% % % % high pass filter the data to remove the low frequency oscilations. This
% % % % is not always beneficial.
% % % h              = HighPassFilter(h, 20, sampling_period);
%Detect the peaks and their times
peaks           = DetectPeaks(h);%peaks are in samples indices
peak_times      = h_time(peaks);%in register with the heart beat vector

% for each pulse gather a phase value based on inter-heart-beat position of
% the pulse.
pulse_phase = [];
pulse_samples = [];
for n = 1:length(volume_files);
    pulse_i                         = n + dummy_scans;%pulse index
    %find the peak that preceeds the pulse.
    [last_peak_time i]              = max(peak_times(peak_times < pulse_times(pulse_i)));
    pulse_samples(n,1)              = i;
    %in words: the difference between the last peak time and the current
    %pulse normalized by the time between this peak and the next peak.
    pulse_phase(n,1)                = 2*pi*(pulse_times(pulse_i) - last_peak_time)./(peak_times(i+1) - peak_times(i));
end

%gather the regressors
regressor = [];
for nOrder = 1:tOrder
    regressor = [regressor cos(nOrder*pulse_phase) sin(nOrder*pulse_phase)];
end
name = repmat({mfilename},1,size(regressor,2));
%%
% % hz = zscore(h);
% % figure
% % plot(h_time,hz,'r')
% % hold on
% % plot(h_time(peaks),hz(peaks),'sk');
% % plot(h_time(pulse_samples),zscore(pulse_phase),'b.-');
% % hold off







