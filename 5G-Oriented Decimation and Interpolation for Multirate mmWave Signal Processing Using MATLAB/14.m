
clear; clc; close all;

%% EXPERIMENT 14: 5G-LIKE MULTIRATE SIGNAL PROCESSING

%% 1. PARAMETER DEFINITIONS
fc_carrier = 28e9;       % Nominal carrier frequency (28 GHz)
Fs_in = 1000e6;          % Input sampling frequency (1000 MHz)
M = 4;                   % Decimation factor
L = 2;                   % Interpolation factor

Fs_dec = Fs_in/M;        % 250 MHz
Fs_out = Fs_dec*L;       % 500 MHz

N_samples = 4000;
N_fir = 60;

f1 = 15e6;
f2 = 35e6;

t_in = (0:N_samples-1)/Fs_in;
x_in = sin(2*pi*f1*t_in) + 0.5*cos(2*pi*f2*t_in);

%% 2. DECIMATION: ANTI-ALIASING FILTER
fc_dec = 0.25;           % Normalized cutoff: 125 MHz

b_dec = fir1(N_fir, fc_dec, 'low');

x_filtered_dec = filter(b_dec, 1, x_in);
x_dec = x_filtered_dec(1:M:end);

%% 3. INTERPOLATION: UPSAMPLING AND ANTI-IMAGING FILTER
x_up = zeros(1, length(x_dec)*L);
x_up(1:L:end) = x_dec;

fc_interp = 0.5;         % Normalized cutoff: 125 MHz
b_interp = fir1(N_fir, fc_interp, 'low');

x_out = L*filter(b_interp, 1, x_up);

%% 4. FILTER PERFORMANCE ANALYSIS
Nfreq = 16384;

[H_dec, f_dec_filter] = freqz(b_dec, 1, Nfreq, Fs_in);
H_dec_db = 20*log10(max(abs(H_dec), 1e-12));

% Passband: 0 to 80% of the decimated Nyquist frequency
pb_idx = f_dec_filter <= 0.8*(Fs_dec/2);

% Stopband: from the decimated Nyquist frequency
sb_idx = f_dec_filter >= Fs_dec/2;

passband_ripple = max(H_dec_db(pb_idx)) ...
                - min(H_dec_db(pb_idx));

stopband_atten = -max(H_dec_db(sb_idx));

%% 5. QUANTITATIVE PERFORMANCE METRICS
% Ideal output reference: original signal sampled at output rate
x_reference = x_in(1:M/L:end);

% Compare only the common output length
N_compare = min(length(x_reference), length(x_out));

% Align comparison by compensating for FIR group delays.
% Decimator delay: N_fir/2 input samples.
% Interpolator delay: N_fir/2 output samples.
delay_dec = N_fir/2/M;
delay_interp = N_fir/2;

% The delays are not necessarily integer output samples.
% Report raw output length and avoid claiming exact reconstruction.
output_peak = max(abs(x_out));

%% 6. COMMAND WINDOW RESULTS
fprintf('\n============================================\n');
fprintf(' EXPERIMENT 14: MULTIRATE PROCESSING\n');
fprintf('============================================\n');

fprintf('Nominal Carrier Frequency : %.2f GHz\n', fc_carrier/1e9);
fprintf('Input Sampling Rate       : %.2f MHz\n', Fs_in/1e6);
fprintf('Decimation Factor (M)     : %d\n', M);
fprintf('Decimated Sampling Rate   : %.2f MHz\n', Fs_dec/1e6);
fprintf('Interpolation Factor (L)  : %d\n', L);
fprintf('Output Sampling Rate      : %.2f MHz\n', Fs_out/1e6);

fprintf('\nSIGNAL LENGTHS\n');
fprintf('Original Signal           : %d samples\n', length(x_in));
fprintf('Decimated Signal          : %d samples\n', length(x_dec));
fprintf('Upsampled Signal          : %d samples\n', length(x_up));
fprintf('Interpolated Output       : %d samples\n', length(x_out));

fprintf('\nDECIMATION FILTER PERFORMANCE\n');
fprintf('Passband Ripple           : %.4f dB\n', passband_ripple);
fprintf('Stopband Attenuation      : %.4f dB\n', stopband_atten);
fprintf('Interpolated Output Peak  : %.6f\n', output_peak);

%% 7. OUTPUT PARAMETER TABLE
Parameter = {
    'Nominal Carrier Frequency'
    'Input Sampling Frequency'
    'Decimated Sampling Frequency'
    'Interpolated Sampling Frequency'
    'Decimation Factor'
    'Interpolation Factor'
    'Original Signal Length'
    'Decimated Signal Length'
    'Interpolated Output Length'
    'Passband Ripple'
    'Stopband Attenuation'
    'Interpolated Output Peak'
};

Value = [
    fc_carrier/1e9
    Fs_in/1e6
    Fs_dec/1e6
    Fs_out/1e6
    M
    L
    length(x_in)
    length(x_dec)
    length(x_out)
    passband_ripple
    stopband_atten
    output_peak
];

Unit = {
    'GHz'
    'MHz'
    'MHz'
    'MHz'
    'Factor'
    'Factor'
    'Samples'
    'Samples'
    'Samples'
    'dB'
    'dB'
    'Amplitude'
};

disp(table(Parameter, Value, Unit));

%% 8. CORRECT FFT-BASED SPECTRUM CALCULATION
Nfft = 8192;

% Original signal spectrum
X_in = fft(x_in, Nfft);
f_in = (0:Nfft/2)*Fs_in/Nfft;
mag_in = abs(X_in(1:Nfft/2+1));
mag_in_db = 20*log10(max(mag_in/max(mag_in), 1e-8));

% Decimated signal spectrum
X_dec = fft(x_dec, Nfft);
f_dec = (0:Nfft/2)*Fs_dec/Nfft;
mag_dec = abs(X_dec(1:Nfft/2+1));
mag_dec_db = 20*log10(max(mag_dec/max(mag_dec), 1e-8));

% Interpolated signal spectrum
X_out = fft(x_out, Nfft);
f_out = (0:Nfft/2)*Fs_out/Nfft;
mag_out = abs(X_out(1:Nfft/2+1));
mag_out_db = 20*log10(max(mag_out/max(mag_out), 1e-8));

%% 9. FIGURE 14A: ORIGINAL SIGNAL
figure('Name', 'Experiment 14A: Original Signal');

subplot(2,1,1);
plot(t_in(1:200)*1e6, x_in(1:200), 'b');
grid on;
xlabel('Time (us)');
ylabel('Amplitude');
title('Original Baseband Signal');

subplot(2,1,2);
plot(f_in/1e6, mag_in_db, 'b', 'LineWidth', 1.1);
grid on;
xlabel('Frequency (MHz)');
ylabel('Normalized Magnitude (dB)');
title('Original Signal Spectrum');
xlim([0 Fs_in/2/1e6]);
ylim([-80 5]);

%% 10. FIGURE 14B: DECIMATED SIGNAL
figure('Name', 'Experiment 14B: Decimation');

t_dec = (0:length(x_dec)-1)/Fs_dec;

subplot(2,1,1);
plot(t_dec(1:min(200,length(t_dec)))*1e6, ...
     x_dec(1:min(200,length(x_dec))), 'r');
grid on;
xlabel('Time (us)');
ylabel('Amplitude');
title('Filtered Decimated Signal');

subplot(2,1,2);
plot(f_dec/1e6, mag_dec_db, 'r', 'LineWidth', 1.1);
grid on;
xlabel('Frequency (MHz)');
ylabel('Normalized Magnitude (dB)');
title('Decimated Signal Spectrum');
xlim([0 Fs_dec/2/1e6]);
ylim([-80 5]);

%% 11. FIGURE 14C: INTERPOLATED SIGNAL
figure('Name', 'Experiment 14C: Interpolation');

t_out = (0:length(x_out)-1)/Fs_out;

subplot(2,1,1);
plot(t_out(1:min(200,length(t_out)))*1e6, ...
     x_out(1:min(200,length(x_out))), 'm');
grid on;
xlabel('Time (us)');
ylabel('Amplitude');
title('Interpolated Output Signal');

subplot(2,1,2);
plot(f_out/1e6, mag_out_db, 'm', 'LineWidth', 1.1);
grid on;
xlabel('Frequency (MHz)');
ylabel('Normalized Magnitude (dB)');
title('Interpolated Output Spectrum');
xlim([0 Fs_out/2/1e6]);
ylim([-80 5]);

%% 12. FIGURE 14D: FILTER RESPONSES
figure('Name', 'Experiment 14D: FIR Filter Responses');

[H_interp, f_interp_filter] = ...
    freqz(b_interp, 1, Nfreq, Fs_out);

plot(f_dec_filter/1e6, ...
     20*log10(max(abs(H_dec),1e-12)), 'b', ...
     'LineWidth', 1.2);
hold on;
plot(f_interp_filter/1e6, ...
     20*log10(max(abs(H_interp),1e-12)), 'r--', ...
     'LineWidth', 1.2);

grid on;
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Anti-Aliasing and Anti-Imaging Filter Responses');
legend('Decimation Filter', 'Interpolation Filter');
ylim([-100 5]);

%% 13. FINAL RESULT
fprintf('\n============================================\n');
fprintf(' FINAL RESULT\n');
fprintf('============================================\n');
fprintf('1. Multirate signal processing was implemented.\n');
fprintf('2. Anti-aliasing filtering and decimation were performed.\n');
fprintf('3. Upsampling and anti-imaging filtering were performed.\n');
fprintf('4. FFT-based spectra were calculated correctly.\n');
fprintf('5. FIR passband ripple and stopband attenuation were evaluated.\n');
fprintf('Note: The 28 GHz carrier is nominal; the generated signal is baseband.\n');
fprintf('============================================\n');
