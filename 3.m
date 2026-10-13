clear; clc; close all;

Fs = 5000;
T = 1/Fs;
fc = 1000;
Omega_c = 2 * pi * fc;

% Analog Filter Prototype: H(s) = Omega_c / (s + Omega_c)
b_a = [Omega_c];
a_a = [1, Omega_c];

% 1. Impulse Invariance Method
[b_imp, a_imp] = impinvar(b_a, a_a, Fs);

% 2. Bilinear Transformation Method (with Prewarping)
Omega_p = (2/T) * tan(Omega_c * T / 2); % Prewarped frequency
b_a_pre = [Omega_p];
a_a_pre = [1, Omega_p];
[b_blt, a_blt] = bilinear(b_a_pre, a_a_pre, Fs);

% ==========================================
% DISPLAY NUMERICAL OUTPUTS IN COMMAND WINDOW
% ==========================================
fprintf('=== EXPERIMENT 3: IMPULSE INVARIANCE vs BILINEAR TRANSFORMATION ===\n\n');

fprintf('Analog Cutoff Frequency (Omega_c): %.2f rad/s\n', Omega_c);
fprintf('Prewarped Analog Frequency (Omega_p): %.2f rad/s\n\n', Omega_p);

fprintf('--- IMPULSE INVARIANCE COEFFICIENTS ---\n');
fprintf('Digital Numerator (b_imp):\n');
disp(b_imp);
fprintf('Digital Denominator (a_imp):\n');
disp(a_imp);

fprintf('--- BILINEAR TRANSFORMATION COEFFICIENTS ---\n');
fprintf('Digital Numerator (b_blt):\n');
disp(b_blt);
fprintf('Digital Denominator (a_blt):\n');
disp(a_blt);

% Frequency Responses
[H_imp, f] = freqz(b_imp, a_imp, 1024, Fs);
[H_blt, ~] = freqz(b_blt, a_blt, 1024, Fs);

% Plotting
figure('Name', 'Impulse Invariance vs Bilinear Transformation');
subplot(2,1,1);
plot(f, abs(H_imp), 'b', 'LineWidth', 1.5); hold on;
plot(f, abs(H_blt), 'r--', 'LineWidth', 1.5); grid on;
legend('Impulse Invariance', 'Bilinear Transformation (Prewarped)');
title('Magnitude Response Comparison');
xlabel('Frequency (Hz)'); ylabel('Absolute Magnitude');

subplot(2,1,2);
plot(f, angle(H_imp)*180/pi, 'b', 'LineWidth', 1.5); hold on;
plot(f, angle(H_blt)*180/pi, 'r--', 'LineWidth', 1.5); grid on;
legend('Impulse Invariance', 'Bilinear Transformation');
title('Phase Response Comparison');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
