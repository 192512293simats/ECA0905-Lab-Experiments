clear; clc; close all;

% Design Parameters
Fs = 10000;         % Sampling frequency in Hz
Ap = 0.5;       	% Passband ripple in dB
As = 30;        	% Stopband attenuation in dB

% 1. Low-Pass Filter Design
fp_lp = 1500;   	% Passband frequency (Hz)
fs_lp = 2000;   	% Stopband frequency (Hz)
wp_lp = fp_lp / (Fs/2);
ws_lp = fs_lp / (Fs/2);

[N_lp, Wn_lp] = cheb1ord(wp_lp, ws_lp, Ap, As);
[b_lp, a_lp] = cheby1(N_lp, Ap, Wn_lp, 'low');

% 2. High-Pass Filter Design
fp_hp = 2500;   	% Passband frequency (Hz)
fs_hp = 1800;   	% Stopband frequency (Hz)
wp_hp = fp_hp / (Fs/2);
ws_hp = fs_hp / (Fs/2);

[N_hp, Wn_hp] = cheb1ord(wp_hp, ws_hp, Ap, As);
[b_hp, a_hp] = cheby1(N_hp, Ap, Wn_hp, 'high');

% ==========================================
% DISPLAY NUMERICAL OUTPUTS IN COMMAND WINDOW
% ==========================================
fprintf('=== EXPERIMENT 2: CHEBYSHEV TYPE-I IIR FILTERS ===\n\n');

fprintf('--- CHEBYSHEV LOW-PASS FILTER ---\n');
fprintf('Calculated Filter Order (N_lp): %d\n', N_lp);
fprintf('Passband Edge Cutoff Frequency (Wn_lp): %.4f Hz\n', Wn_lp * (Fs/2));
fprintf('Numerator Coefficients (b_lp):\n');
disp(b_lp);
fprintf('Denominator Coefficients (a_lp):\n');
disp(a_lp);

fprintf('--- CHEBYSHEV HIGH-PASS FILTER ---\n');
fprintf('Calculated Filter Order (N_hp): %d\n', N_hp);
fprintf('Passband Edge Cutoff Frequency (Wn_hp): %.4f Hz\n', Wn_hp * (Fs/2));
fprintf('Numerator Coefficients (b_hp):\n');
disp(b_hp);
fprintf('Denominator Coefficients (a_hp):\n');
disp(a_hp);

% Frequency Response Calculation
[H_lp, f] = freqz(b_lp, a_lp, 1024, Fs);
[H_hp, ~] = freqz(b_hp, a_hp, 1024, Fs);

% Plotting
figure('Name', 'Chebyshev Type-I LPF and HPF Responses');
subplot(2,2,1);
plot(f, 20*log10(abs(H_lp)), 'b', 'LineWidth', 1.5); grid on;
title(['Chebyshev LPF Magnitude (Order N = ', num2str(N_lp), ')']);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); ylim([-60 5]);

subplot(2,2,2);
plot(f, angle(H_lp)*180/pi, 'r', 'LineWidth', 1.5); grid on;
title('Chebyshev LPF Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');

subplot(2,2,3);
plot(f, 20*log10(abs(H_hp)), 'b', 'LineWidth', 1.5); grid on;
title(['Chebyshev HPF Magnitude (Order N = ', num2str(N_hp), ')']);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); ylim([-60 5]);

subplot(2,2,4);
plot(f, angle(H_hp)*180/pi, 'r', 'LineWidth', 1.5); grid on;
title('Chebyshev HPF Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
