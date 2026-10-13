
clc; clear; close all;

%% EXPERIMENT: OVERFLOW, SATURATION, WRAP-AROUND AND SIGNAL SCALING

%% INPUT PARAMETERS
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

A = 0.9;
f = 50;
Fs = 1000;
N = 200;
B = 8;
SF = 0.5;

n = 0:N-1;
x = A*sin(2*pi*f*n/Fs);

%% Q1.7 FIXED-POINT PARAMETERS
S = 2^(B-1);
QMAX = (S-1)/S;
QMIN = -1;

% Signed Q1.7 quantizer with saturation
q = @(z) min(max(round(z*S)/S,QMIN),QMAX);

%% 1. IDEAL AND QUANTIZED RESPONSE
y = filter(b,a,x);
ys = q(y);

Peak_Unscaled = max(abs(y));

%% 2. OVERFLOW DETECTION
ov = (y < QMIN) | (y > QMAX);
No_Overflow = sum(ov);

%% 3. SATURATION
ysat = q(y);
Peak_Saturation = max(abs(ysat));

%% 4. WRAP-AROUND
yint = round(y*S);
ywrap = (mod(yint + S,2*S)-S)/S;

Peak_Wrap = max(abs(ywrap));

%% 5. SIGNAL SCALING
xs = SF*x;
y_scaled_ideal = filter(b,a,xs);
yscale = q(y_scaled_ideal);

Peak_Scaled = max(abs(yscale));

ov_scaled = (y_scaled_ideal < QMIN) | ...
            (y_scaled_ideal > QMAX);
No_Overflow_Scaled = sum(ov_scaled);

%% 6. ERROR CALCULATION
RMSE_Sat = sqrt(mean((y-ysat).^2));
RMSE_Wrap = sqrt(mean((y-ywrap).^2));
RMSE_Scale = sqrt(mean((y_scaled_ideal-yscale).^2));

%% 7. DISPLAY RESULTS
fprintf('\n========================================\n');
fprintf(' FIXED-POINT ARITHMETIC RESULTS\n');
fprintf('========================================\n');

fprintf('Q1.7 Range: [%.7f, %.7f]\n',QMIN,QMAX);
fprintf('Quantization Step Size: %.7f\n',1/S);
fprintf('Unscaled Peak: %.6f\n',Peak_Unscaled);
fprintf('Overflow Count: %d\n',No_Overflow);
fprintf('Saturation Peak: %.6f\n',Peak_Saturation);
fprintf('Wrap-around Peak: %.6f\n',Peak_Wrap);
fprintf('Scaling Factor: %.2f\n',SF);
fprintf('Scaled Peak: %.6f\n',Peak_Scaled);
fprintf('Scaled Overflow Count: %d\n',No_Overflow_Scaled);
fprintf('Saturation RMSE: %.6e\n',RMSE_Sat);
fprintf('Wrap-around RMSE: %.6e\n',RMSE_Wrap);
fprintf('Scaled Output RMSE: %.6e\n',RMSE_Scale);

%% 8. OUTPUT PARAMETER TABLE
Method = {'Unscaled';'Saturation';'Wrap-around';'Scaled'};

Peak = [Peak_Unscaled;
        Peak_Saturation;
        Peak_Wrap;
        Peak_Scaled];

RMSE = [0;
        RMSE_Sat;
        RMSE_Wrap;
        RMSE_Scale];

Overflow = [No_Overflow;
            0;
            0;
            No_Overflow_Scaled];

T = table(Method,Peak,RMSE,Overflow);

fprintf('\nOUTPUT PARAMETER TABLE\n');
disp(T);

%% 9. GRAPH: IDEAL VS SATURATION
figure;
plot(n,y,'k','LineWidth',1.2);
hold on;
plot(n,ysat,'r--','LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Ideal vs Saturation');
legend('Ideal','Saturation','Location','best');

%% 10. GRAPH: IDEAL VS WRAP-AROUND
figure;
plot(n,y,'k','LineWidth',1.2);
hold on;
plot(n,ywrap,'r--','LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Ideal vs Wrap-around');
legend('Ideal','Wrap-around','Location','best');

%% 11. GRAPH: SCALED OUTPUT
figure;
plot(n,y_scaled_ideal,'k','LineWidth',1.2);
hold on;
plot(n,yscale,'r--','LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Scaled Output');
legend('Scaled Ideal','Quantized','Location','best');

%% 12. GRAPH: ERROR COMPARISON
figure;
plot(n,y-ysat,'LineWidth',1.2);
hold on;
plot(n,y-ywrap,'--','LineWidth',1.2);
plot(n,y_scaled_ideal-yscale,':','LineWidth',1.5);
grid on;
xlabel('Sample n');
ylabel('Error');
title('Error Comparison');
legend('Saturation','Wrap-around','Scaling', ...
       'Location','best');

%% 13. GRAPH: PEAK MAGNITUDE COMPARISON
figure;
bar(Peak);
grid on;
set(gca,'XTick',1:4,'XTickLabel',Method);
xlabel('Method');
ylabel('Peak Magnitude');
title('Peak Magnitude Comparison');

%% 14. FINAL RESULT
fprintf('\n========================================\n');
fprintf(' FINAL RESULT\n');
fprintf('========================================\n');
fprintf('Q1.7 fixed-point quantization was implemented.\n');
fprintf('Overflow, saturation and wrap-around were analyzed.\n');
fprintf('Signal scaling and quantization errors were evaluated.\n');
fprintf('All five comparison graphs were generated.\n');
