
clc; clear; close all;

%% EXPERIMENT 20: LPC-BASED SPEECH ANALYSIS AND VOCODER IMPLEMENTATION

%% 1. INPUT SIGNAL
Fs = 8000;                 % Sampling frequency
N_frame = 240;             % Frame length (30 ms)
hop = 120;                 % Frame shift (15 ms)
p = 12;                    % LPC prediction order

if isfile('speech_sample.wav')
    [x_speech, Fs_file] = audioread('speech_sample.wav');

    % Convert stereo to mono
    if size(x_speech,2) > 1
        x_speech = mean(x_speech,2);
    end

    % Resample if the input sampling rate is different
    if Fs_file ~= Fs
        x_speech = resample(x_speech, Fs, Fs_file);
    end

    x_speech = x_speech(:)';
else
    % Fallback test signal if the audio file is unavailable
    duration = 0.5;
    t_temp = 0:1/Fs:duration-1/Fs;
    x_speech = sin(2*pi*250*t_temp) ...
             + 0.5*sin(2*pi*700*t_temp);
end

x_speech = x_speech / (max(abs(x_speech)) + eps);
t = (0:length(x_speech)-1)/Fs;

%% 2. FRAME-BY-FRAME LPC ANALYSIS AND VOCODER SYNTHESIS
win = hamming(N_frame)';
num_frames = floor((length(x_speech)-N_frame)/hop)+1;

x_synth = zeros(1,length(x_speech));
win_gain = zeros(1,length(x_speech));

rng(42);

for m = 1:num_frames
    idx = (m-1)*hop + (1:N_frame);

    frame = x_speech(idx).*win;

    % LPC analysis
    [a, E] = lpc(frame,p);

    % Gain from prediction error energy
    G = sqrt(max(E,0));

    % Random-noise excitation for an unvoiced-style vocoder
    excitation = G*randn(1,N_frame);

    % All-pole LPC synthesis filter
    synth_frame = filter(1,a,excitation);

    % Overlap-add reconstruction
    x_synth(idx) = x_synth(idx) + synth_frame.*win;
    win_gain(idx) = win_gain(idx) + win.^2;
end

valid = win_gain > 1e-6;
x_synth(valid) = x_synth(valid)./win_gain(valid);

%% 3. LPC ANALYSIS OF A SINGLE FRAME
frame_start = min(1+floor(length(x_speech)/4), ...
                  length(x_speech)-N_frame+1);

frame_sample = x_speech(frame_start:frame_start+N_frame-1);
frame_sample = frame_sample.*win;

[a_sample,E_sample] = lpc(frame_sample,p);
G_sample = sqrt(max(E_sample,0));

% LPC spectral envelope
[H_lpc,f_lpc] = freqz(G_sample,a_sample,2048,Fs);

% Periodogram estimate of the frame spectrum
[Pxx,f_fft] = pwelch(frame_sample,win,0,2048,Fs);

%% 4. FORMANT FREQUENCY ESTIMATION
r = roots(a_sample);
r = r(imag(r)>0);

formant_freq = angle(r)*Fs/(2*pi);
bandwidth = -Fs/pi*log(max(abs(r),eps));

% Keep plausible formants with limited bandwidth
keep = formant_freq > 90 & ...
       formant_freq < Fs/2 & ...
       bandwidth < 400;

formant_freq = sort(formant_freq(keep));
bandwidth = bandwidth(keep);

%% 5. COMMAND WINDOW RESULTS
fprintf('==================================================\n');
fprintf(' EXPERIMENT 20: LPC SPEECH ANALYSIS AND VOCODER\n');
fprintf('==================================================\n');

fprintf('Input Signal Length: %d samples\n',length(x_speech));
fprintf('Sampling Rate: %d Hz\n',Fs);
fprintf('Frame Length: %d samples (%.1f ms)\n', ...
    N_frame,1000*N_frame/Fs);
fprintf('Frame Shift: %d samples (%.1f ms)\n', ...
    hop,1000*hop/Fs);
fprintf('LPC Prediction Order: %d\n\n',p);

fprintf('LPC Filter Coefficients (a):\n');
fprintf('%.4f ',a_sample);
fprintf('\n\n');

fprintf('Gain (G): %.6f\n',G_sample);
fprintf('Prediction Error Energy (E): %.6f\n\n',E_sample);

fprintf('Estimated Formant Frequencies:\n');
for k = 1:min(3,length(formant_freq))
    fprintf('Formant %d (F%d): %.2f Hz\n', ...
        k,k,formant_freq(k));
end

if isempty(formant_freq)
    fprintf('No reliable formants detected in this frame.\n');
end

fprintf('\nStatus: LPC analysis and vocoder synthesis completed.\n');
fprintf('==================================================\n');

%% 6. FIGURE 20A: INPUT SPEECH WAVEFORM
figA = figure('Name','Experiment 20A: Input Speech', ...
    'NumberTitle','off');

plot(t,x_speech,'b','LineWidth',1.1);
grid on; box on;
title('Input Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');

exportgraphics(figA,'Exp20_SubfigA_InputSpeech.png', ...
    'Resolution',600);

%% 7. FIGURE 20B: LPC SPECTRAL ANALYSIS
figB = figure('Name','Experiment 20B: LPC Spectrum', ...
    'NumberTitle','off');

plot(f_fft,10*log10(Pxx+eps),'b','LineWidth',1.1);
hold on;
plot(f_lpc,20*log10(abs(H_lpc)+eps),'r','LineWidth',1.2);

grid on; box on;
title('Speech Spectrum and LPC Spectral Envelope');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Estimated Speech Spectrum','LPC Envelope');
xlim([0 Fs/2]);

exportgraphics(figB,'Exp20_SubfigB_LPC_Spectrum.png', ...
    'Resolution',600);

%% 8. FIGURE 20C: VOCODER OUTPUT
figC = figure('Name','Experiment 20C: Vocoder Output', ...
    'NumberTitle','off');

plot(t,x_synth,'m','LineWidth',1.1);
grid on; box on;
title('LPC Vocoder Synthesized Signal');
xlabel('Time (s)');
ylabel('Amplitude');

exportgraphics(figC,'Exp20_SubfigC_Vocoder_Output.png', ...
    'Resolution',600);
