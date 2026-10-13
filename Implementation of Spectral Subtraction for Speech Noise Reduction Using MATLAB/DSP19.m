
clc; clear; close all;

%% EXPERIMENT 19: SPECTRAL SUBTRACTION FOR SPEECH NOISE REDUCTION

%% 1. SIGNAL PREPARATION AND NOISE CORRUPTION
Fs = 16000;                  % Sampling frequency (Hz)
duration = 1.5;              % Duration (seconds)
t = 0:1/Fs:duration-1/Fs;

f1 = 300;
f2 = 800;
f3 = 2200;

% Generate clean test signal
envelope = sin(pi*t/duration).^2;
x_clean = (sin(2*pi*f1*t) ...
    + 0.6*sin(2*pi*f2*t) ...
    + 0.3*sin(2*pi*f3*t)) .* envelope;

% Add AWGN at the specified input SNR
target_SNR = 10;
speech_pwr = mean(x_clean.^2);
noise_pwr = speech_pwr / (10^(target_SNR/10));

rng(42);
noise = sqrt(noise_pwr) * randn(size(x_clean));
x_noisy = x_clean + noise;

%% 2. SPECTRAL SUBTRACTION PARAMETERS
frame_duration = 0.020;      % 20 ms
hop_duration = 0.010;        % 10 ms

N_frame = round(frame_duration*Fs);  % 320 samples
hop = round(hop_duration*Fs);        % 160 samples
N_fft = 512;

win = hamming(N_frame)';
alpha = 2.0;                 % Over-subtraction factor
beta = 0.02;                 % Spectral floor factor
N_noise_frames = 10;

%% 3. NOISE SPECTRUM ESTIMATION
% Estimate noise spectrum from the initial low-energy frames
noise_spectrum_sum = zeros(1,N_fft);

for k = 1:N_noise_frames
    idx = (k-1)*hop + (1:N_frame);

    frame = x_noisy(idx).*win;
    noise_spectrum_sum = noise_spectrum_sum + abs(fft(frame,N_fft));
end

N_est = noise_spectrum_sum/N_noise_frames;

%% 4. SPECTRAL SUBTRACTION AND OVERLAP-ADD
num_frames = floor((length(x_noisy)-N_frame)/hop)+1;

x_enhanced = zeros(size(x_noisy));
win_gain = zeros(size(x_noisy));

for m = 1:num_frames
    idx = (m-1)*hop + (1:N_frame);

    % Window the noisy signal frame
    frame = x_noisy(idx).*win;

    % Compute FFT
    Y = fft(frame,N_fft);
    Y_mag = abs(Y);
    Y_phase = angle(Y);

    % Spectral magnitude subtraction
    S_hat_mag = max(Y_mag-alpha*N_est, beta*Y_mag);

    % Preserve noisy-signal phase
    S_hat = S_hat_mag.*exp(1i*Y_phase);

    % Inverse FFT
    frame_enhanced = real(ifft(S_hat,N_fft));

    % Overlap-add synthesis
    x_enhanced(idx) = x_enhanced(idx) ...
        + frame_enhanced(1:N_frame).*win;

    win_gain(idx) = win_gain(idx)+win.^2;
end

% Normalize overlap-add output
valid_idx = win_gain > 1e-3;
x_enhanced(valid_idx) = x_enhanced(valid_idx) ...
    ./win_gain(valid_idx);

%% 5. PERFORMANCE ANALYSIS
% Calculate SNR using the clean reference and the actual noise/error
SNR_in = 10*log10(sum(x_clean.^2) ...
    /sum((x_noisy-x_clean).^2));

% Evaluate over samples covered by overlap-add processing
eval_idx = valid_idx & (x_clean ~= 0);

SNR_out = 10*log10(sum(x_clean(eval_idx).^2) ...
    /sum((x_enhanced(eval_idx)-x_clean(eval_idx)).^2));

SNR_improvement = SNR_out-SNR_in;

%% 6. COMMAND WINDOW RESULTS
fprintf('==================================================\n');
fprintf(' EXPERIMENT 19: SPEECH NOISE REDUCTION\n');
fprintf('       USING SPECTRAL SUBTRACTION\n');
fprintf('==================================================\n');

fprintf('Sampling Frequency: %d Hz\n',Fs);
fprintf('Signal Duration: %.2f seconds\n',duration);
fprintf('Original Signal Length: %d samples\n',length(x_clean));
fprintf('Frame Duration: %.0f ms (%d samples)\n', ...
    frame_duration*1000,N_frame);
fprintf('Frame Shift: %.0f ms (%d samples)\n', ...
    hop_duration*1000,hop);
fprintf('FFT Length: %d\n',N_fft);
fprintf('Noise Estimation Frames: %d\n\n',N_noise_frames);

fprintf('PERFORMANCE METRICS:\n');
fprintf('Input SNR: %.2f dB\n',SNR_in);
fprintf('Output SNR: %.2f dB\n',SNR_out);
fprintf('SNR Improvement: %.2f dB\n\n',SNR_improvement);

fprintf('Status: Spectral subtraction processing completed.\n');
fprintf('==================================================\n');

%% 7. TIME-DOMAIN WAVEFORMS
figW = 8.5;
figH = 5.5;

figA = figure('Name','Experiment 19A: Clean Signal', ...
    'NumberTitle','off');
plot(t,x_clean,'b','LineWidth',1.1);
grid on; box on;
title('Clean Input Signal');
xlabel('Time (s)');
ylabel('Amplitude');
exportgraphics(figA,'Exp19_SubfigA_CleanSpeech.png', ...
    'Resolution',600);

figB = figure('Name','Experiment 19B: Noisy Signal', ...
    'NumberTitle','off');
plot(t,x_noisy,'r','LineWidth',1.0);
grid on; box on;
title(sprintf('Noisy Signal (SNR = %.2f dB)',SNR_in));
xlabel('Time (s)');
ylabel('Amplitude');
exportgraphics(figB,'Exp19_SubfigB_NoisySpeech.png', ...
    'Resolution',600);

figC = figure('Name','Experiment 19C: Enhanced Signal', ...
    'NumberTitle','off');
plot(t,x_enhanced,'g','LineWidth',1.1);
grid on; box on;
title(sprintf('Enhanced Signal (SNR = %.2f dB)',SNR_out));
xlabel('Time (s)');
ylabel('Amplitude');
exportgraphics(figC,'Exp19_SubfigC_EnhancedSpeech.png', ...
    'Resolution',600);

%% 8. FREQUENCY-DOMAIN COMPARISON
Nfft_plot = 4096;
f = (0:Nfft_plot/2)*(Fs/Nfft_plot);

X_clean = abs(fft(x_clean,Nfft_plot));
X_noisy = abs(fft(x_noisy,Nfft_plot));
X_enhanced = abs(fft(x_enhanced,Nfft_plot));

figD = figure('Name','Experiment 19D: Spectrum Comparison', ...
    'NumberTitle','off');

plot(f,X_clean(1:Nfft_plot/2+1),'b','LineWidth',1.1);
hold on;
plot(f,X_noisy(1:Nfft_plot/2+1),'r','LineWidth',0.9);
plot(f,X_enhanced(1:Nfft_plot/2+1),'g','LineWidth',1.0);

grid on; box on;
title('Clean, Noisy and Enhanced Signal Spectra');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Clean Signal','Noisy Signal','Enhanced Signal');
xlim([0 Fs/2]);

exportgraphics(figD,'Exp19_SubfigD_SpectrumComparison.png', ...
    'Resolution',600);
