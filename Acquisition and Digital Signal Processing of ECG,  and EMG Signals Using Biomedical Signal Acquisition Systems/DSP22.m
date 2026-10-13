clc; clear; close all;

%% EXPERIMENT 22: ACQUISITION AND DSP OF ECG, EEG AND EMG SIGNALS
% Loads real signals or generates synthetic examples, filters them,
% displays time/frequency plots, and extracts basic signal features.
% Requires MATLAB Signal Processing Toolbox.

%% 1. INPUT PARAMETERS
fs_ecg = 360;                  % ECG sampling frequency (Hz)
fs_eeg = 256;                  % EEG sampling frequency (Hz)
fs_emg = 1000;                 % EMG sampling frequency (Hz)

band_ecg = [0.5 40];           % ECG band-pass band (Hz)
band_eeg = [1 40];             % EEG band-pass band (Hz)
band_emg = [20 450];           % EMG band-pass band (Hz)

T = 10;                        % Record duration (s)
mains = 50;                    % Power-line frequency in India (Hz)
filtOrder = 4;                 % Butterworth prototype order

useRealData = false;           % true to load real .mat files
ecgFile = 'ecg_data.mat';      % Must contain ecg_raw
eegFile = 'eeg_data.mat';      % Must contain eeg_raw
emgFile = 'emg_data.mat';      % Must contain emg_raw

rng(1);                        % Repeatable synthetic noise

%% 2. ACQUIRE OR LOAD SIGNALS
if useRealData
    S = load(ecgFile);
    assert(isfield(S,'ecg_raw'), 'ECG file must contain ecg_raw.');
    ecg_raw = S.ecg_raw(:).';

    S = load(eegFile);
    assert(isfield(S,'eeg_raw'), 'EEG file must contain eeg_raw.');
    eeg_raw = S.eeg_raw(:).';

    S = load(emgFile);
    assert(isfield(S,'emg_raw'), 'EMG file must contain emg_raw.');
    emg_raw = S.emg_raw(:).';

    assert(numel(ecg_raw) >= T*fs_ecg, 'ECG data is shorter than T seconds.');
    assert(numel(eeg_raw) >= T*fs_eeg, 'EEG data is shorter than T seconds.');
    assert(numel(emg_raw) >= T*fs_emg, 'EMG data is shorter than T seconds.');

    ecg_raw = ecg_raw(1:T*fs_ecg);
    eeg_raw = eeg_raw(1:T*fs_eeg);
    emg_raw = emg_raw(1:T*fs_emg);
else
    [ecg_raw, ecg_clean] = synth_ecg(fs_ecg,T,mains);
    [eeg_raw, eeg_clean] = synth_eeg(fs_eeg,T,mains);
    [emg_raw, emg_clean] = synth_emg(fs_emg,T,mains);
end

t_ecg = (0:numel(ecg_raw)-1)/fs_ecg;
t_eeg = (0:numel(eeg_raw)-1)/fs_eeg;
t_emg = (0:numel(emg_raw)-1)/fs_emg;

%% 3. DIGITAL FILTERING
ecg_f = bp_filter(ecg_raw,band_ecg,fs_ecg,filtOrder);
eeg_f = bp_filter(eeg_raw,band_eeg,fs_eeg,filtOrder);
emg_f = bp_filter(emg_raw,band_emg,fs_emg,filtOrder);

% Remove 50-Hz mains interference from EMG with a notch filter.
% The notch frequency and bandwidth are normalized to Nyquist.
w0 = mains/(fs_emg/2);
bw = w0/35;
[bn,an] = iirnotch(w0,bw);
emg_f = filtfilt(bn,an,emg_f);

%% 4. FREQUENCY SPECTRA
[fE, XE_raw] = fft_spectrum(ecg_raw,fs_ecg);
[~, XE_filt] = fft_spectrum(ecg_f,fs_ecg);

[fG, XG_raw] = fft_spectrum(eeg_raw,fs_eeg);
[~, XG_filt] = fft_spectrum(eeg_f,fs_eeg);

[fM, XM_raw] = fft_spectrum(emg_raw,fs_emg);
[~, XM_filt] = fft_spectrum(emg_f,fs_emg);

%% 5. FEATURE EXTRACTION

% ECG: R-peaks, RR intervals, and heart rate
[pk,loc] = findpeaks(ecg_f, ...
    'MinPeakHeight',0.5*max(ecg_f), ...
    'MinPeakDistance',round(0.3*fs_ecg));

rPeakTimes = (loc-1)/fs_ecg;     % Correct zero-based time axis
RR = diff(rPeakTimes);
HR_all = 60./RR;

if isempty(RR)
    HR_avg = NaN;
else
    HR_avg = 60/mean(RR);
end

% EEG: delta, theta, alpha, and beta band powers using Welch PSD
[Pxx,fw] = pwelch(eeg_f,hamming(2*fs_eeg),fs_eeg, ...
    2*fs_eeg,fs_eeg);

bands = struct('name',{'Delta','Theta','Alpha','Beta'}, ...
    'range',{[1 4],[4 8],[8 13],[13 30]});

bandPower = zeros(1,4);
for k = 1:4
    idx = fw >= bands(k).range(1) & fw < bands(k).range(2);
    if nnz(idx) >= 2
        bandPower(k) = trapz(fw(idx),Pxx(idx));
    else
        bandPower(k) = 0;
    end
end

if sum(bandPower) > 0
    relPower = 100*bandPower/sum(bandPower);
else
    relPower = zeros(1,4);
end
[~,domIdx] = max(bandPower);

% EMG: RMS amplitude and 100-ms moving RMS envelope
RMS_raw = sqrt(mean(emg_raw.^2));
RMS_f = sqrt(mean(emg_f.^2));
win = max(1,round(0.1*fs_emg));
emg_env = sqrt(movmean(emg_f.^2,win));

% Approximate filter output-to-removed-component ratios.
% For synthetic signals, the known-clean-signal SNRs are reported below.
SNR_ecg_est = snr_db(ecg_f,ecg_raw-ecg_f);
SNR_eeg_est = snr_db(eeg_f,eeg_raw-eeg_f);
SNR_emg_est = snr_db(emg_f,emg_raw-emg_f);

%% 6. ECG PLOTS
figure('Name','Experiment 22A: ECG','Color','w', ...
    'Position',[50 50 1100 700]);

subplot(2,2,1);
plot(t_ecg,ecg_raw,'r'); grid on;
title('ECG - Raw waveform');
xlabel('Time (s)'); ylabel('Amplitude (mV)');

subplot(2,2,2);
plot(t_ecg,ecg_f,'b'); hold on;
if ~isempty(loc)
    plot(rPeakTimes,pk,'kv','MarkerFaceColor','g');
end
grid on;
title(sprintf('ECG - Filtered (%.1f-%g Hz) with R-peaks',band_ecg));
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('Filtered','R-peaks');

subplot(2,2,3);
plot(fE,XE_raw,'r'); hold on; plot(fE,XE_filt,'b');
grid on; xlim([0 min(100,fs_ecg/2)]);
title('ECG - Magnitude spectrum');
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
legend('Raw','Filtered');

subplot(2,2,4);
if ~isempty(HR_all)
    stem(rPeakTimes(2:end),HR_all,'filled'); grid on;
    ylim([0 max(HR_all)*1.3]);
else
    text(0.1,0.5,'Not enough R-peaks to calculate heart rate');
    axis off;
end
title(sprintf('Instantaneous heart rate (mean = %.1f bpm)',HR_avg));
xlabel('Time (s)'); ylabel('HR (beats/min)');

%% 7. EEG PLOTS
figure('Name','Experiment 22B: EEG','Color','w', ...
    'Position',[80 80 1100 700]);

subplot(2,2,1);
plot(t_eeg,eeg_raw,'r'); grid on;
title('EEG - Raw waveform');
xlabel('Time (s)'); ylabel('Amplitude (\muV)');

subplot(2,2,2);
plot(t_eeg,eeg_f,'b'); grid on;
title(sprintf('EEG - Filtered (%g-%g Hz)',band_eeg));
xlabel('Time (s)'); ylabel('Amplitude (\muV)');

subplot(2,2,3);
plot(fG,XG_raw,'r'); hold on; plot(fG,XG_filt,'b');
grid on; xlim([0 min(80,fs_eeg/2)]);
title('EEG - Magnitude spectrum');
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
legend('Raw','Filtered');

subplot(2,2,4);
bar(relPower);
set(gca,'XTick',1:4, ...
    'XTickLabel',{'Delta 1-4','Theta 4-8','Alpha 8-13','Beta 13-30'});
grid on; title('EEG - Relative band power');
ylabel('Power (%)');

%% 8. EMG PLOTS
figure('Name','Experiment 22C: EMG','Color','w', ...
    'Position',[110 110 1100 700]);

subplot(2,2,1);
plot(t_emg,emg_raw,'r'); grid on;
title('EMG - Raw waveform');
xlabel('Time (s)'); ylabel('Amplitude (mV)');

subplot(2,2,2);
plot(t_emg,emg_f,'b'); hold on;
plot(t_emg,emg_env,'k','LineWidth',1.5);
grid on; title(sprintf('EMG - Filtered (%g-%g Hz) + RMS envelope',band_emg));
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('Filtered','RMS envelope (100 ms)');

subplot(2,2,3);
plot(fM,XM_raw,'r'); hold on; plot(fM,XM_filt,'b');
grid on; xlim([0 fs_emg/2]);
title('EMG - Magnitude spectrum');
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
legend('Raw','Filtered');

subplot(2,2,4);
bar([RMS_raw RMS_f]);
set(gca,'XTick',1:2,'XTickLabel',{'Raw RMS','Filtered RMS'});
grid on; title('EMG - RMS amplitude'); ylabel('Amplitude');

%% 9. FILTER FREQUENCY RESPONSES
[sE,gE] = bp_sos(band_ecg,fs_ecg,filtOrder);
[sG,gG] = bp_sos(band_eeg,fs_eeg,filtOrder);
[sM,gM] = bp_sos(band_emg,fs_emg,filtOrder);

[hE,wE] = freqz(sE,4096,fs_ecg);
[hG,wG] = freqz(sG,4096,fs_eeg);
[hM,wM] = freqz(sM,4096,fs_emg);

figure('Name','Experiment 22D: Filter Responses','Color','w', ...
    'Position',[140 140 1100 450]);
plot(wE,20*log10(abs(hE)+eps),'r','LineWidth',1.3); hold on;
plot(wG,20*log10(abs(hG)+eps),'g','LineWidth',1.3);
plot(wM,20*log10(abs(hM)+eps),'b','LineWidth',1.3);
grid on; xlim([0 fs_emg/2]); ylim([-80 5]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Band-pass filter magnitude responses (single pass)');
legend('ECG 0.5-40 Hz','EEG 1-40 Hz','EMG 20-450 Hz');

%% 10. OUTPUT PARAMETER TABLE
fprintf('\n============================================================\n');
fprintf(' EXPERIMENT 22: ACQUISITION AND DSP OF ECG, EEG AND EMG\n');
fprintf('============================================================\n');
fprintf('%-28s %-12s %-12s %-12s\n','Parameter','ECG','EEG','EMG');
fprintf('%-28s %-12d %-12d %-12d\n','Sampling frequency (Hz)', ...
    fs_ecg,fs_eeg,fs_emg);
fprintf('%-28s %-12s %-12s %-12s\n','Filter band (Hz)', ...
    sprintf('%g-%g',band_ecg),sprintf('%g-%g',band_eeg), ...
    sprintf('%g-%g',band_emg));
fprintf('%-28s %-12d %-12d %-12d\n','Samples', ...
    numel(ecg_raw),numel(eeg_raw),numel(emg_raw));
fprintf('%-28s %-12.2f %-12.2f %-12.2f\n', ...
    'Filter output/removal ratio (dB)', ...
    SNR_ecg_est,SNR_eeg_est,SNR_emg_est);

fprintf('\n--- ECG features ---\n');
fprintf('Number of R-peaks : %d\n',numel(loc));
if ~isempty(RR)
    fprintf('Mean RR interval : %.3f s\n',mean(RR));
    fprintf('Mean heart rate : %.1f beats/min\n',HR_avg);
    fprintf('HR range (min/max) : %.1f / %.1f beats/min\n', ...
        min(HR_all),max(HR_all));
else
    fprintf('Heart rate : unavailable (insufficient detected R-peaks)\n');
end

fprintf('\n--- EEG band power ---\n');
for k = 1:4
    fprintf('%-6s (%2d-%2d Hz) : %.3f uV^2 (%.1f%%)\n', ...
        bands(k).name,bands(k).range(1),bands(k).range(2), ...
        bandPower(k),relPower(k));
end
fprintf('Dominant rhythm : %s\n',bands(domIdx).name);

fprintf('\n--- EMG features ---\n');
fprintf('RMS (raw) : %.4f mV\n',RMS_raw);
fprintf('RMS (filtered) : %.4f mV\n',RMS_f);
fprintf('Peak RMS envelope : %.4f mV\n',max(emg_env));

if ~useRealData
    fprintf('\n--- SNR using known synthetic clean signals ---\n');
    fprintf('ECG SNR raw : %.2f dB | filtered : %.2f dB\n', ...
        snr_db(ecg_clean,ecg_raw-ecg_clean), ...
        snr_db(ecg_clean,ecg_f-ecg_clean));
    fprintf('EEG SNR raw : %.2f dB | filtered : %.2f dB\n', ...
        snr_db(eeg_clean,eeg_raw-eeg_clean), ...
        snr_db(eeg_clean,eeg_f-eeg_clean));
    fprintf('EMG SNR raw : %.2f dB | filtered : %.2f dB\n', ...
        snr_db(emg_clean,emg_raw-emg_clean), ...
        snr_db(emg_clean,emg_f-emg_clean));
end

fprintf('\nRESULT: ECG, EEG and EMG filtering and feature extraction completed.\n');
fprintf('============================================================\n');

%% LOCAL FUNCTIONS
function [sos,g] = bp_sos(band,fs,order)
    [z,p,k] = butter(order,band/(fs/2),'bandpass');
    [sos,g] = zp2sos(z,p,k);
end

function y = bp_filter(x,band,fs,order)
    [sos,g] = bp_sos(band,fs,order);
    y = filtfilt(sos,g,x);
end

function [f,X] = fft_spectrum(x,fs)
    x = x(:).';
    N = numel(x);
    Xf = fft(x-mean(x));
    X = abs(Xf(1:floor(N/2)+1))/N;
    if numel(X) > 2
        X(2:end-1) = 2*X(2:end-1);
    end
    f = (0:floor(N/2))*fs/N;
end

function value = snr_db(sig,noise)
    value = 10*log10((sum(sig.^2)+eps)/(sum(noise.^2)+eps));
end

function [x,clean] = synth_ecg(fs,T,mains)
    t = 0:1/fs:T-1/fs;
    clean = zeros(size(t));
    beatT = 0.4;

    % Gaussian approximations of P, Q, R, S and T waves.
    % Columns: amplitude, centre offset (s), standard deviation (s).
    waves = [ 0.15 -0.20 0.025;
             -0.15 -0.04 0.010;
              1.20  0.00 0.012;
             -0.25  0.04 0.012;
              0.30  0.25 0.050];

    while beatT < T
        for q = 1:size(waves,1)
            clean = clean + waves(q,1)*exp( ...
                -((t-beatT-waves(q,2)).^2)/(2*waves(q,3)^2));
        end
        beatT = beatT + 0.83 + 0.03*randn;
    end

    baseline = 0.5*sin(2*pi*0.25*t);
    hum = 0.15*sin(2*pi*mains*t);
    noise = 0.03*randn(size(t));
    x = clean + baseline + hum + noise;
end

function [x,clean] = synth_eeg(fs,T,mains)
    t = 0:1/fs:T-1/fs;
    clean = 12*sin(2*pi*2*t+1) + 8*sin(2*pi*6*t+2) + ...
        25*sin(2*pi*10*t) + 6*sin(2*pi*20*t+0.5);
    clean = clean + 3*randn(size(t));

    drift = 40*sin(2*pi*0.1*t);
    hum = 10*sin(2*pi*mains*t);
    x = clean + drift + hum;
end

function [x,clean] = synth_emg(fs,T,mains)
    t = 0:1/fs:T-1/fs;
    noise = randn(size(t));
    [sos,g] = bp_sos([20 450],fs,4);
    noise = filtfilt(sos,g,noise);

    envelope = zeros(size(t));
    for c = [1.5 4.5 7.5]
        envelope = envelope + exp(-((t-c).^2)/(2*0.5^2));
    end

    clean = 0.5*envelope.*noise;
    motion = 0.3*sin(2*pi*3*t);
    hum = 0.1*sin(2*pi*mains*t);
    x = clean + motion + hum + 0.01*randn(size(t));
end
