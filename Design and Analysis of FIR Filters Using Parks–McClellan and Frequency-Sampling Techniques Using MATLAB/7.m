clear;
clc;
close all;

Fs = 10000;
fp = 1500;
fs = 2000;
N = 40;

% Normalized Frequency Bands
f_bands = [0 fp fs Fs/2]/(Fs/2);
a_bands = [1 1 0 0];

%% 1. Parks-McClellan (Equiripple) Design
b_pm = firpm(N,f_bands,a_bands);

%% 2. Frequency-Sampling Design
b_fs = fir2(N,f_bands,a_bands);

%% 3. Display Numerical Outputs
fprintf('=== EXPERIMENT 5: PARKS-MCCLELLAN VS FREQUENCY SAMPLING ===\n\n');

fprintf('Filter Order (N): %d\n',N);
fprintf('Sampling Frequency: %d Hz\n',Fs);
fprintf('Passband Edge: %d Hz\n',fp);
fprintf('Stopband Edge: %d Hz\n\n',fs);

fprintf('--- PARKS-MCCLELLAN COEFFICIENTS (First 10) ---\n');
disp(b_pm(1:10));

fprintf('--- FREQUENCY-SAMPLING COEFFICIENTS (First 10) ---\n');
disp(b_fs(1:10));

%% 4. Frequency Response
[H_pm,f] = freqz(b_pm,1,2048,Fs);
[H_fs,~] = freqz(b_fs,1,2048,Fs);

% Avoid log of zero
Mag_pm = 20*log10(max(abs(H_pm),1e-10));
Mag_fs = 20*log10(max(abs(H_fs),1e-10));

Phase_pm = angle(H_pm)*180/pi;
Phase_fs = angle(H_fs)*180/pi;

%% 5. Magnitude and Phase Response Comparison
figure('Name','Parks-McClellan vs Frequency-Sampling');

subplot(2,1,1);
plot(f,Mag_pm,'b','LineWidth',1.5);
hold on;
plot(f,Mag_fs,'r--','LineWidth',1.5);
grid on;
legend('Parks-McClellan (Equiripple)','Frequency-Sampling', ...
       'Location','best');
title(['FIR Filter Magnitude Response Comparison (Order N = ',num2str(N),')']);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 Fs/2]);
ylim([-80 10]);

subplot(2,1,2);
plot(f,Phase_pm,'b','LineWidth',1.5);
hold on;
plot(f,Phase_fs,'r--','LineWidth',1.5);
grid on;
legend('Parks-McClellan','Frequency-Sampling','Location','best');
title('Phase Response Comparison');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
xlim([0 Fs/2]);
ylim([-200 200]);

%% 6. Display Filter Performance
passband = f <= fp;
stopband = f >= fs;

Rp_pm = max(Mag_pm(passband))-min(Mag_pm(passband));
Rp_fs = max(Mag_fs(passband))-min(Mag_fs(passband));

As_pm = -max(Mag_pm(stopband));
As_fs = -max(Mag_fs(stopband));

fprintf('\n--- PERFORMANCE COMPARISON ---\n');
fprintf('Parks-McClellan Passband Ripple: %.4f dB\n',Rp_pm);
fprintf('Frequency-Sampling Passband Ripple: %.4f dB\n',Rp_fs);
fprintf('Parks-McClellan Stopband Attenuation: %.4f dB\n',As_pm);
fprintf('Frequency-Sampling Stopband Attenuation: %.4f dB\n',As_fs);

fprintf('\nTransition Bandwidth: %d Hz\n',fs-fp);
fprintf('\nComparison completed successfully.\n');