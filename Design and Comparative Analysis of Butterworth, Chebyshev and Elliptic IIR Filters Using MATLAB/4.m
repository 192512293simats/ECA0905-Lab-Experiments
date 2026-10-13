clc;
clear;
close all;

%% Design and Comparative Analysis of IIR Filters

Fs = 1000;
Fp = 100;
Fst = 150;
Rp = 1;
Rs = 40;

Wp = Fp/(Fs/2);
Wst = Fst/(Fs/2);

%% 1. Butterworth Filter

[Nb, Wb] = buttord(Wp, Wst, Rp, Rs);
[bB, aB] = butter(Nb, Wb, 'low');

%% 2. Chebyshev Type-I Filter

[Nc, Wc] = cheb1ord(Wp, Wst, Rp, Rs);
[bC, aC] = cheby1(Nc, Rp, Wc, 'low');

%% 3. Elliptic Filter

[Ne, We] = ellipord(Wp, Wst, Rp, Rs);
[bE, aE] = ellip(Ne, Rp, Rs, We, 'low');

%% Display Filter Orders

fprintf('\nFILTER ORDER COMPARISON\n');
fprintf('Butterworth       = %d\n', Nb);
fprintf('Chebyshev Type-I  = %d\n', Nc);
fprintf('Elliptic          = %d\n', Ne);

%% Display Coefficients

fprintf('\nBUTTERWORTH COEFFICIENTS\n');
disp('Numerator:');
disp(bB);
disp('Denominator:');
disp(aB);

fprintf('\nCHEBYSHEV TYPE-I COEFFICIENTS\n');
disp('Numerator:');
disp(bC);
disp('Denominator:');
disp(aC);

fprintf('\nELLIPTIC COEFFICIENTS\n');
disp('Numerator:');
disp(bE);
disp('Denominator:');
disp(aE);

%% Frequency Response

Nfft = 4096;

[HB, f] = freqz(bB, aB, Nfft, Fs);
[HC, ~] = freqz(bC, aC, Nfft, Fs);
[HE, ~] = freqz(bE, aE, Nfft, Fs);

MB = 20*log10(max(abs(HB), 1e-12));
MC = 20*log10(max(abs(HC), 1e-12));
ME = 20*log10(max(abs(HE), 1e-12));

%% Performance Analysis

P = f <= Fp;
S = f >= Fst;

RippleB = max(MB(P)) - min(MB(P));
RippleC = max(MC(P)) - min(MC(P));
RippleE = max(ME(P)) - min(ME(P));

StopB = -max(MB(S));
StopC = -max(MC(S));
StopE = -max(ME(S));

fprintf('\nPASSBAND RIPPLE (dB)\n');
fprintf('Butterworth       = %.4f\n', RippleB);
fprintf('Chebyshev Type-I  = %.4f\n', RippleC);
fprintf('Elliptic          = %.4f\n', RippleE);

fprintf('\nSTOPBAND ATTENUATION (dB)\n');
fprintf('Butterworth       = %.4f\n', StopB);
fprintf('Chebyshev Type-I  = %.4f\n', StopC);
fprintf('Elliptic          = %.4f\n', StopE);

fprintf('\nTRANSITION BANDWIDTH = %.2f Hz\n', Fst-Fp);

%% Magnitude Response Comparison

figure('Name', 'IIR Magnitude Response Comparison');

plot(f, MB, 'LineWidth', 1.5);
hold on;
plot(f, MC, 'LineWidth', 1.5);
plot(f, ME, 'LineWidth', 1.5);

xline(Fp, '--', 'Fp');
xline(Fst, '--', 'Fst');
yline(-Rp, ':', 'Passband limit');
yline(-Rs, ':', 'Stopband limit');

grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Comparison of IIR Filter Magnitude Responses');
legend('Butterworth', 'Chebyshev Type-I', 'Elliptic', ...
       'Location', 'southwest');
xlim([0 Fs/2]);
ylim([-80 5]);
hold off;

%% Phase Response Comparison

PhaseB = unwrap(angle(HB))*180/pi;
PhaseC = unwrap(angle(HC))*180/pi;
PhaseE = unwrap(angle(HE))*180/pi;

figure('Name', 'IIR Phase Response Comparison');

plot(f, PhaseB, 'LineWidth', 1.3);
hold on;
plot(f, PhaseC, 'LineWidth', 1.3);
plot(f, PhaseE, 'LineWidth', 1.3);

grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Comparison of IIR Filter Phase Responses');
legend('Butterworth', 'Chebyshev Type-I', 'Elliptic');
xlim([0 Fs/2]);
hold off;

%% Filter Order Comparison

figure('Name', 'IIR Filter Order Comparison');

orders = [Nb Nc Ne];
bar(orders);

set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', ...
    {'Butterworth', 'Chebyshev-I', 'Elliptic'});

ylabel('Filter Order');
xlabel('Filter Type');
title('Comparison of IIR Filter Orders');
grid on;

for k = 1:3
    text(k, orders(k), sprintf('%d', orders(k)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end