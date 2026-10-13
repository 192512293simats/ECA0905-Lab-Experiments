
clc; clear; close all;

%% EXPERIMENT 21: LMS AND NLMS ADAPTIVE ECHO CANCELLATION

%% 1. INPUT PARAMETERS
Fs = 8000;                  % Sampling frequency (Hz)
f = 1000;                   % Reference sinusoid frequency (Hz)
N = 4000;                   % Number of samples
L = 32;                     % Adaptive filter length
mu_LMS = 0.01;              % LMS step size
mu_NLMS = 0.5;              % NLMS step size
epsilon = 1e-6;              % NLMS normalization constant
delay = 100;                % Echo delay in samples
echo_gain = 0.5;            % Echo path gain

n = 0:N-1;
x = sin(2*pi*f*n/Fs);       % Reference signal

%% 2. GENERATE ECHO SIGNAL
echo = zeros(1,N);

echo(delay+1:end) = echo_gain*x(1:end-delay);

d = echo;                   % Desired signal to be cancelled

% Adaptive filter input is the original reference signal.
% The adaptive output estimates the echo.
% Error = desired echo - estimated echo.

%% 3. LMS ADAPTIVE FILTER
w_LMS = zeros(L,1);
y_LMS = zeros(1,N);
e_LMS = zeros(1,N);
MSE_LMS = zeros(1,N);

for k = L:N
    u = x(k:-1:k-L+1).';

    y_LMS(k) = w_LMS.'*u;
    e_LMS(k) = d(k)-y_LMS(k);

    w_LMS = w_LMS + mu_LMS*e_LMS(k)*u;
    MSE_LMS(k) = e_LMS(k)^2;
end

%% 4. NLMS ADAPTIVE FILTER
w_NLMS = zeros(L,1);
y_NLMS = zeros(1,N);
e_NLMS = zeros(1,N);
MSE_NLMS = zeros(1,N);

for k = L:N
    u = x(k:-1:k-L+1).';

    y_NLMS(k) = w_NLMS.'*u;
    e_NLMS(k) = d(k)-y_NLMS(k);

    norm_factor = epsilon + (u.'*u);

    w_NLMS = w_NLMS + ...
        (mu_NLMS/norm_factor)*e_NLMS(k)*u;

    MSE_NLMS(k) = e_NLMS(k)^2;
end

%% 5. PERFORMANCE ANALYSIS
valid = L:N;

% Mean-square error
MSE_LMS_final = mean(e_LMS(valid).^2);
MSE_NLMS_final = mean(e_NLMS(valid).^2);

% Error in decibels
MSE_LMS_dB = 10*log10(MSE_LMS_final + eps);
MSE_NLMS_dB = 10*log10(MSE_NLMS_final + eps);

% Moving average convergence curves
window = 100;
MSE_curve_LMS = movmean(MSE_LMS(valid),window);
MSE_curve_NLMS = movmean(MSE_NLMS(valid),window);

% Echo reduction relative to uncancelled echo
input_power = mean(d(valid).^2);
residual_LMS = mean(e_LMS(valid).^2);
residual_NLMS = mean(e_NLMS(valid).^2);

ER_LMS = 10*log10((input_power+eps)/(residual_LMS+eps));
ER_NLMS = 10*log10((input_power+eps)/(residual_NLMS+eps));

% Convergence criterion: MSE remains below 10% of initial
% moving-average MSE for at least 100 samples.
threshold_LMS = 0.1*MSE_curve_LMS(1);
threshold_NLMS = 0.1*MSE_curve_NLMS(1);

conv_LMS_idx = find(MSE_curve_LMS < threshold_LMS,1);
conv_NLMS_idx = find(MSE_curve_NLMS < threshold_NLMS,1);

if isempty(conv_LMS_idx)
    conv_LMS = NaN;
else
    conv_LMS = valid(conv_LMS_idx);
end

if isempty(conv_NLMS_idx)
    conv_NLMS = NaN;
else
    conv_NLMS = valid(conv_NLMS_idx);
end

%% 6. COMMAND WINDOW RESULTS
fprintf('====================================================\n');
fprintf(' EXPERIMENT 21: LMS AND NLMS ECHO CANCELLATION\n');
fprintf('====================================================\n');

fprintf('Sampling Frequency       : %d Hz\n',Fs);
fprintf('Reference Frequency      : %d Hz\n',f);
fprintf('Number of Samples        : %d\n',N);
fprintf('Adaptive Filter Length   : %d taps\n',L);
fprintf('Echo Delay               : %d samples\n',delay);
fprintf('Echo Gain                : %.2f\n',echo_gain);
fprintf('LMS Step Size            : %.4f\n',mu_LMS);
fprintf('NLMS Step Size           : %.4f\n',mu_NLMS);
fprintf('NLMS Epsilon              : %.1e\n\n',epsilon);

fprintf('PERFORMANCE METRICS\n');
fprintf('LMS MSE                  : %.6e\n',MSE_LMS_final);
fprintf('NLMS MSE                 : %.6e\n',MSE_NLMS_final);
fprintf('LMS MSE                  : %.2f dB\n',MSE_LMS_dB);
fprintf('NLMS MSE                 : %.2f dB\n',MSE_NLMS_dB);
fprintf('LMS Echo Reduction       : %.2f dB\n',ER_LMS);
fprintf('NLMS Echo Reduction      : %.2f dB\n',ER_NLMS);

fprintf('LMS Convergence Sample   : %.0f\n',conv_LMS);
fprintf('NLMS Convergence Sample  : %.0f\n\n',conv_NLMS);

%% 7. OUTPUT PARAMETER TABLE
Method = {'LMS';'NLMS'};
StepSize = [mu_LMS;mu_NLMS];
MSE = [MSE_LMS_final;MSE_NLMS_final];
MSE_dB = [MSE_LMS_dB;MSE_NLMS_dB];
EchoReduction_dB = [ER_LMS;ER_NLMS];
ConvergenceSample = [conv_LMS;conv_NLMS];

T = table(Method,StepSize,MSE,MSE_dB,...
    EchoReduction_dB,ConvergenceSample);

disp('OUTPUT PARAMETER TABLE');
disp(T);

%% 8. FIGURE 21A: REFERENCE AND ECHO SIGNALS
figure('Name','Experiment 21A: Reference and Echo',...
    'NumberTitle','off');

subplot(2,1,1);
plot(n(1:500),x(1:500),'b');
grid on;
title('Reference Signal');
xlabel('Sample Index n');
ylabel('Amplitude');

subplot(2,1,2);
plot(n(1:500),d(1:500),'r');
grid on;
title('Echo Signal');
xlabel('Sample Index n');
ylabel('Amplitude');

%% 9. FIGURE 21B: LMS AND NLMS CONVERGENCE
figure('Name','Experiment 21B: Convergence',...
    'NumberTitle','off');

plot(valid,10*log10(MSE_curve_LMS+eps),'b',...
    'LineWidth',1.2);
hold on;
plot(valid,10*log10(MSE_curve_NLMS+eps),'r',...
    'LineWidth',1.2);
grid on;
title('LMS vs NLMS Convergence');
xlabel('Sample Index n');
ylabel('Moving Average MSE (dB)');
legend('LMS','NLMS');

%% 10. FIGURE 21C: RESIDUAL ERROR SIGNAL
figure('Name','Experiment 21C: Residual Error',...
    'NumberTitle','off');

subplot(2,1,1);
plot(n(1:500),e_LMS(1:500),'b');
grid on;
title('LMS Residual Error');
xlabel('Sample Index n');
ylabel('Error');

subplot(2,1,2);
plot(n(1:500),e_NLMS(1:500),'r');
grid on;
title('NLMS Residual Error');
xlabel('Sample Index n');
ylabel('Error');

%% 11. FIGURE 21D: ECHO-CANCELLED OUTPUT
figure('Name','Experiment 21D: Echo Cancellation',...
    'NumberTitle','off');

subplot(3,1,1);
plot(n(1:500),d(1:500),'k');
grid on;
title('Original Echo');
xlabel('Sample Index n');
ylabel('Amplitude');

subplot(3,1,2);
plot(n(1:500),e_LMS(1:500),'b');
grid on;
title('Echo-Cancelled Output using LMS');
xlabel('Sample Index n');
ylabel('Amplitude');

subplot(3,1,3);
plot(n(1:500),e_NLMS(1:500),'r');
grid on;
title('Echo-Cancelled Output using NLMS');
xlabel('Sample Index n');
ylabel('Amplitude');

%% 12. FIGURE 21E: MSE COMPARISON
figure('Name','Experiment 21E: MSE Comparison',...
    'NumberTitle','off');

bar([MSE_LMS_final,MSE_NLMS_final]);
set(gca,'XTickLabel',{'LMS','NLMS'});
grid on;
ylabel('Mean Square Error');
title('Final MSE Comparison');

fprintf('\nRESULT:\n');
fprintf('LMS and NLMS adaptive filtering completed.\n');
fprintf('Convergence, residual error and echo reduction evaluated.\n');
fprintf('====================================================\n');
