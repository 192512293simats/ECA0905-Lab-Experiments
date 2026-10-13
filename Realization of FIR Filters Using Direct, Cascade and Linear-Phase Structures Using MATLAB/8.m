clc;
clear;
close all;

%% FIR FILTER PARAMETERS
b = [0.0625 0.2500 0.3750 0.2500 0.0625];
Fs = 1000;
N = 128;

x = [1 zeros(1,N-1)];
n = 0:N-1;

%% 1. DIRECT REALIZATION
yd = filter(b,1,x);

%% 2. CASCADE REALIZATION
% H(z) = (0.25 + 0.5z^-1 + 0.25z^-2)^2
b1 = [0.25 0.5 0.25];
b2 = [0.25 0.5 0.25];

yc = filter(b1,1,x);
yc = filter(b2,1,yc);

%% 3. SYMMETRY AND LINEAR PHASE
symErr = max(abs(b-fliplr(b)));
D = (length(b)-1)/2;

if symErr < 1e-12
    symmetry = 'Symmetric';
else
    symmetry = 'Not Symmetric';
end

%% 4. FREQUENCY, PHASE AND GROUP DELAY
[H,w] = freqz(b,1,1024);
phase = unwrap(angle(H));
gd = grpdelay(b,1,1024);

%% 5. NUMERICAL COMPARISON
err = max(abs(yd-yc));
dev = err/max(abs(yd))*100;

%% 6. DISPLAY RESULTS
fprintf('\n============================================\n');
fprintf(' FIR FILTER REALIZATION RESULTS\n');
fprintf('============================================\n');

fprintf('\nFIR Coefficients:\n');
disp(b);

fprintf('Sampling Frequency = %d Hz\n',Fs);
fprintf('Number of Samples = %d\n',N);
fprintf('Coefficient Symmetry = %s\n',symmetry);
fprintf('Symmetry Error = %.3e\n',symErr);
fprintf('Group Delay = %.2f samples\n',D);
fprintf('Maximum Error = %.3e\n',err);
fprintf('Percentage Deviation = %.3e %%\n',dev);

%% 7. FIRST 10 IMPULSE RESPONSE VALUES
fprintf('\n============================================\n');
fprintf('FIRST 10 IMPULSE RESPONSE VALUES\n');
fprintf('============================================\n');

fprintf('%4s %12s %12s\n','n','Direct','Cascade');

for k = 1:10
    fprintf('%4d %12.6f %12.6f\n', ...
        n(k),yd(k),yc(k));
end

%% 8. CASCADE REALIZATION ZEROS
fprintf('\n============================================\n');
fprintf('CASCADE REALIZATION ZEROS\n');
fprintf('============================================\n');

z = roots(b);
disp(z);

%% 9. OUTPUT PARAMETER TABLE
Parameter = {
    'Direct realization output'
    'Cascade realization output'
    'Impulse response'
    'Phase response'
    'Comparison of realizations'
};

Observation = {
    'Direct FIR filtering'
    'Cascade of first-order sections'
    'Equal to FIR coefficient sequence'
    'Linear phase'
    'Equivalent transfer-function behavior'
};

Value = {
    sprintf('[%.4f %.4f %.4f %.4f %.4f]',yd(1:5))
    sprintf('[%.4f %.4f %.4f %.4f %.4f]',yc(1:5))
    sprintf('[%.4f %.4f %.4f %.4f %.4f]',yd(1:5))
    sprintf('Group delay = %.2f samples',D)
    sprintf('Max error = %.3e; Deviation = %.3e %%',err,dev)
};

Unit_Note = {
    'Amplitude'
    'Amplitude'
    'h[n]'
    'Samples'
    'Numerical precision'
};

Output_Table = table(Parameter,Observation,Value,Unit_Note);

fprintf('\n============================================\n');
fprintf('OUTPUT PARAMETER TABLE\n');
fprintf('============================================\n');

disp(Output_Table);

%% 10. GRAPH: DIRECT VS CASCADE IMPULSE RESPONSE
figure;

stem(n,yd,'filled');
hold on;
stem(n,yc,'--o','MarkerSize',3);

grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Direct and Cascade FIR Impulse Responses');
legend('Direct Realization','Cascade Realization');

%% 11. GRAPH: MAGNITUDE RESPONSE
figure;

plot(w/pi,abs(H),'LineWidth',1.5);
grid on;

xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('|H(e^{j\omega})|');
title('FIR Magnitude Response');
xlim([0 1]);

%% 12. GRAPH: UNWRAPPED PHASE RESPONSE
figure;

plot(w/pi,phase,'LineWidth',1.5);
grid on;

xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');
title('Unwrapped Phase Response');
xlim([0 1]);

%% 13. GRAPH: GROUP DELAY
figure;

plot(w/pi,gd,'LineWidth',1.5);
grid on;

xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Group Delay (samples)');
title('Group Delay');
xlim([0 1]);
ylim([0 3]);

%% 14. GRAPH: NUMERICAL ERROR
figure;

plot(n,abs(yd-yc),'LineWidth',1.5);
grid on;

xlabel('Sample n');
ylabel('Absolute Error');
title('Direct vs Cascade Numerical Error');
xlim([0 N-1]);

%% 15. FINAL RESULT
fprintf('\n============================================\n');
fprintf('FINAL RESULT\n');
fprintf('============================================\n');

fprintf('1. FIR filter realized using Direct structure.\n');
fprintf('2. FIR filter realized using Cascade structure.\n');
fprintf('3. Coefficients are %s.\n',symmetry);
fprintf('4. Linear-phase group delay = %.2f samples.\n',D);
fprintf('5. Maximum realization error = %.3e.\n',err);
fprintf('6. Percentage deviation = %.3e %%.\n',dev);
fprintf('7. Direct and Cascade responses are equivalent.\n');