
clc; clear; close all;

%% EXPERIMENT 12: MULTIRATE SIGNAL PROCESSING

%% 1. SIGNAL GENERATION
Fs = 8000;                 % Original sampling frequency
f1 = 500;
f2 = 1500;                 % Signal frequencies
M = 2;                     % Decimation factor
L = 2;                     % Interpolation factor
N = 1000;                  % Number of samples

fc_dec = 1800;             % Anti-aliasing cutoff frequency
fc_interp = 3600;          % Anti-imaging cutoff frequency

n = 0:N-1;
t = n/Fs;

x = sin(2*pi*f1*t) + 0.8*sin(2*pi*f2*t);

%% 2. DECIMATION

% Direct downsampling without filtering
x_down_direct = x(1:M:end);
Fs_dec = Fs/M;

% Anti-aliasing FIR filter
b_dec = fir1(60, fc_dec/(Fs/2));
x_dec_filtered = filter(b_dec, 1, x);

% Downsample after filtering
x_dec = x_dec_filtered(1:M:end);

%% 3. INTERPOLATION

% Insert L-1 zeros between input samples
x_up = zeros(1, N*L);
x_up(1:L:end) = x;

Fs_interp = Fs*L;

% Anti-imaging FIR filter designed at interpolated rate
b_interp = L*fir1(60, fc_interp/(Fs_interp/2));

% Filter the upsampled signal
x_interp_filtered = filter(b_interp, 1, x_up);

%% 4. COMMAND WINDOW OUTPUT

fprintf('\n============================================\n');
fprintf(' EXPERIMENT 12: MULTIRATE SIGNAL PROCESSING\n');
fprintf('============================================\n');

fprintf('Original Sampling Frequency (Fs): %d Hz\n', Fs);
fprintf('Decimation Factor (M): %d\n', M);
fprintf('New Sampling Frequency after Decimation: %d Hz\n\n', Fs_dec);

fprintf('Interpolation Factor (L): %d\n', L);
fprintf('New Sampling Frequency after Interpolation: %d Hz\n\n', Fs_interp);

fprintf('Signal Lengths:\n');
fprintf('Original Signal (x): %d samples\n', length(x));
fprintf('Direct Downsampled Signal: %d samples\n', length(x_down_direct));
fprintf('Filtered Decimated Signal: %d samples\n', length(x_dec));
fprintf('Upsampled Signal: %d samples\n', length(x_up));
fprintf('Interpolated Signal: %d samples\n\n', length(x_interp_filtered));

fprintf('Decimation Filter Cutoff: %d Hz\n', fc_dec);
fprintf('Interpolation Filter Cutoff: %d Hz\n\n', fc_interp);

fprintf('Status: Decimation and Interpolation executed successfully.\n');
fprintf('============================================\n');

%% 5. SPECTRUM COMPUTATION

Nfft = 2048;

[X_mag, f_x] = compute_spectrum(x, Fs, Nfft);

[X_down_mag, f_down] = ...
    compute_spectrum(x_down_direct, Fs_dec, Nfft);

[X_dec_mag, f_dec] = ...
    compute_spectrum(x_dec, Fs_dec, Nfft);

[X_up_mag, f_up] = ...
    compute_spectrum(x_up, Fs_interp, Nfft);

[X_interp_mag, f_interp] = ...
    compute_spectrum(x_interp_filtered, Fs_interp, Nfft);

%% 6. ORIGINAL SIGNAL

figure('Name', 'Experiment 12A: Original Signal');

subplot(2,1,1);
plot(t(1:100)*1000, x(1:100), 'b', 'LineWidth', 1.2);
grid on;
title('Original Input Waveform');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_x, X_mag, 'b', 'LineWidth', 1.2);
grid on;
title('Original Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-Fs/2 Fs/2]);

%% 7. DECIMATED SIGNAL

t_dec = (0:length(x_dec)-1)/Fs_dec;

figure('Name', 'Experiment 12B: Decimation');

subplot(2,1,1);
plot(t_dec(1:50)*1000, x_dec(1:50), 'r', 'LineWidth', 1.2);
grid on;
title('Filtered Decimated Waveform (M = 2)');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_down, X_down_mag, 'k--', 'LineWidth', 1);
hold on;
plot(f_dec, X_dec_mag, 'r', 'LineWidth', 1.2);
grid on;
title('Direct and Filtered Decimated Spectra');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Direct Downsampling', 'Filtered Decimation');
xlim([-Fs_dec/2 Fs_dec/2]);

%% 8. INTERPOLATED SIGNAL

t_interp = (0:length(x_interp_filtered)-1)/Fs_interp;

figure('Name', 'Experiment 12C: Interpolation');

subplot(2,1,1);
plot(t_interp(1:200)*1000, ...
    x_interp_filtered(1:200), 'm', 'LineWidth', 1.2);
grid on;
title('Interpolated Waveform (L = 2)');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_up, X_up_mag, 'k--', 'LineWidth', 1);
hold on;
plot(f_interp, X_interp_mag, 'm', 'LineWidth', 1.2);
grid on;
title('Direct and Filtered Interpolated Spectra');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Upsampled Signal', 'Filtered Interpolation');
xlim([-Fs_interp/2 Fs_interp/2]);

%% 9. OUTPUT PARAMETER TABLE

Parameter = {
    'Original signal length'
    'Decimated signal length'
    'Interpolated signal length'
    'Original sampling frequency'
    'Decimated sampling frequency'
    'Interpolated sampling frequency'
    'Decimation factor'
    'Interpolation factor'
};

Value = [
    length(x)
    length(x_dec)
    length(x_interp_filtered)
    Fs
    Fs_dec
    Fs_interp
    M
    L
];

Unit = {
    'Samples'
    'Samples'
    'Samples'
    'Hz'
    'Hz'
    'Hz'
    'Factor'
    'Factor'
};

Output_Table = table(Parameter, Value, Unit);

fprintf('\nOUTPUT PARAMETER TABLE\n');
disp(Output_Table);

%% 10. FINAL RESULT

fprintf('\n============================================\n');
fprintf(' FINAL RESULT\n');
fprintf('============================================\n');
fprintf('1. Multirate signal processing was implemented.\n');
fprintf('2. Decimation was performed with anti-aliasing filtering.\n');
fprintf('3. Interpolation was performed with anti-imaging filtering.\n');
fprintf('4. Time-domain waveforms and spectra were plotted.\n');
fprintf('5. All signal lengths and sampling rates were verified.\n');

%% HELPER FUNCTION

function [X_mag, f] = compute_spectrum(signal, Fs_curr, Nfft)

    X = fft(signal, Nfft);
    X_mag = abs(fftshift(X))/length(signal);

    f = (-Nfft/2:Nfft/2-1)*(Fs_curr/Nfft);

end
