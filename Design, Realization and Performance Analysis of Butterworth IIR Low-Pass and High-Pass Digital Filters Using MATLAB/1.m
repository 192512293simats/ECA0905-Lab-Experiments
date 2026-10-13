clear; 
clc; 
close all;

% Design Parameters
Fs = 8000;      	% Sampling frequency in Hz
Ap = 1;         	% Passband ripple in dB
As = 40;	        % Stopband attenuation in dB

% 1. Low-Pass Filter Design
fp_lp = 1000;   	% LPF Passband frequency (Hz)
fs_lp = 1500;   	% LPF Stopband frequency (Hz)
wp_lp = fp_lp / (Fs/2);
ws_lp = fs_lp / (Fs/2);

[N_lp, Wn_lp] = buttord(wp_lp, ws_lp, Ap, As);
[b_lp, a_lp] = butter(N_lp, Wn_lp, 'low');

% 2. High-Pass Filter Design
fp_hp = 2500;   	% HPF Passband frequency (Hz)
fs_hp = 1800;   	% HPF Stopband frequency (Hz)
wp_hp = fp_hp / (Fs/2);
ws_hp = fs_hp / (Fs/2);

[N_hp, Wn_hp] = buttord(wp_hp, ws_hp, Ap, As);
[b_hp, a_hp] = butter(N_hp, Wn_hp, 'high');

% ==========================================
% DISPLAY NUMERICAL OUTPUTS IN COMMAND WINDOW
% ==========================================
fprintf('=== EXPERIMENT 1: BUTTERWORTH IIR FILTERS ===\n\n');

fprintf('--- LOW-PASS FILTER ---\n');
fprintf('Calculated Filter Order (N_lp): %d\n', N_lp);
fprintf('3-dB Cutoff Frequency (Wn_lp): %.4f Hz\n', Wn_lp * (Fs/2));
fprintf('Numerator Coefficients (b_lp):\n');
disp(b_lp);
fprintf('Denominator Coefficients (a_lp):\n');
disp(a_lp);

fprintf('--- HIGH-PASS FILTER ---\n');
fprintf('Calculated Filter Order (N_hp): %d\n', N_hp);
fprintf('3-dB Cutoff Frequency (Wn_hp): %.4f Hz\n', Wn_hp * (Fs/2));
fprintf('Numerator Coefficients (b_hp):\n');
disp(b_hp);
fprintf('Denominator Coefficients (a_hp):\n');
disp(a_hp);

% Frequency Response Calculation
[H_lp, f] = freqz(b_lp, a_lp, 1024, Fs);
[H_hp, ~] = freqz(b_hp, a_hp, 1024, Fs);

% Plotting
figure('Name', 'Butterworth IIR LPF and HPF');
subplot(2,2,1);
plot(f, 20*log10(abs(H_lp)), 'b', 'LineWidth', 1.5); grid on;
title(['Butterworth LPF Magnitude (Order N = ', num2str(N_lp), ')']);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); ylim([-80 5]);

subplot(2,2,2);
plot(f, angle(H_lp)*180/pi, 'r', 'LineWidth', 1.5); grid on;
title('Butterworth LPF Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');

subplot(2,2,3);
plot(f, 20*log10(abs(H_hp)), 'b', 'LineWidth', 1.5); grid on;
title(['Butterworth HPF Magnitude (Order N = ', num2str(N_hp), ')']);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); ylim([-80 5]);

subplot(2,2,4);
plot(f, angle(H_hp)*180/pi, 'r', 'LineWidth', 1.5); grid on;
title('Butterworth HPF Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
