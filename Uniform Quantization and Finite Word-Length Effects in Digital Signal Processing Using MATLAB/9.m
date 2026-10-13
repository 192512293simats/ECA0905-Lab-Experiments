clc;
clear;
close all;

%% PARAMETERS
A = 1;
f = 50;
Fs = 1000;
N = 1000;

n = 0:N-1;
x = A*sin(2*pi*f*n/Fs);

B = 8;
bits = [4 8 12 16];

%% 1. 8-BIT QUANTIZATION
D = 2*A/2^B;

xr = D*round(x/D);   % Rounding
xt = D*fix(x/D);     % Truncation

er = x - xr;         % Rounding error
et = x - xt;         % Truncation error

SQr = 10*log10(mean(x.^2)/mean(er.^2));
SQt = 10*log10(mean(x.^2)/mean(et.^2));

%% 2. INPUT / OUTPUT QUANTIZATION
xq = xr;
y = x;
yq = D*round(y/D);
eo = y - yq;

SQout = 10*log10(mean(y.^2)/mean(eo.^2));

%% 3. WORD-LENGTH COMPARISON
step = zeros(size(bits));
sqnrR = zeros(size(bits));
sqnrT = zeros(size(bits));
maxER = zeros(size(bits));
maxET = zeros(size(bits));

for k = 1:length(bits)
    b = bits(k);
    d = 2*A/2^b;

    qr = d*round(x/d);
    qt = d*fix(x/d);

    er1 = x - qr;
    et1 = x - qt;

    step(k) = d;
    sqnrR(k) = 10*log10(mean(x.^2)/mean(er1.^2));
    sqnrT(k) = 10*log10(mean(x.^2)/mean(et1.^2));

    maxER(k) = max(abs(er1));
    maxET(k) = max(abs(et1));
end

%% 4. DISPLAY RESULTS
fprintf('\n========== QUANTIZATION RESULTS ==========\n');
fprintf('Amplitude = %.1f | Frequency = %d Hz | Fs = %d Hz\n',A,f,Fs);
fprintf('Samples = %d | Main word length = %d bits\n',N,B);
fprintf('8-bit Step Size = %.7f\n',D);
fprintf('Rounding SQNR = %.4f dB\n',SQr);
fprintf('Truncation SQNR = %.4f dB\n',SQt);
fprintf('Output SQNR = %.4f dB\n',SQout);
fprintf('Max Rounding Error = %.6f\n',max(abs(er)));
fprintf('Max Truncation Error = %.6f\n',max(abs(et)));

%% 5. FIRST 10 SAMPLES TABLE
T1 = table(n(1:10)',x(1:10)',xr(1:10)',er(1:10)', ...
    xt(1:10)',et(1:10)', ...
    'VariableNames',{'n','Original','Rounded','RoundError', ...
    'Truncated','TruncError'});

disp('FIRST 10 SAMPLES');
disp(T1);

%% 6. WORD-LENGTH TABLE
T2 = table(bits',step',sqnrR',sqnrT',maxER',maxET', ...
    'VariableNames',{'Bits','StepSize','SQNR_Round_dB', ...
    'SQNR_Trunc_dB','MaxError_Round','MaxError_Trunc'});

disp('WORD-LENGTH COMPARISON');
disp(T2);

%% 7. OUTPUT PARAMETER TABLE
Parameter = {
    'Quantized signal'
    'Quantization error'
    'SQNR'
    'Word-length comparison'
};

Observation = {
    sprintf('8-bit rounding, Delta = %.7f',D)
    sprintf('Max rounding error = %.6f',max(abs(er)))
    sprintf('Rounding = %.4f dB; Truncation = %.4f dB',SQr,SQt)
    'SQNR increases as word length increases'
};

Unit = {
    'Amplitude'
    'Amplitude'
    'dB'
    'Bits / dB'
};

Output_Table = table(Parameter,Observation,Unit);

disp('OUTPUT PARAMETER TABLE');
disp(Output_Table);

%% 8. ORIGINAL AND QUANTIZED SIGNAL
figure;
plot(n(1:100),x(1:100),'b','LineWidth',1.2);
hold on;
stairs(n(1:100),xr(1:100),'r--','LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Original and 8-bit Quantized Signal');
legend('Original','Quantized');

%% 9. ROUNDING AND TRUNCATION ERROR
figure;
plot(n,er,'b','LineWidth',1.1);
hold on;
plot(n,et,'r--','LineWidth',1.1);
grid on;
xlabel('Sample n');
ylabel('Error');
title('Rounding and Truncation Quantization Error');
legend('Rounding','Truncation');

%% 10. SQNR VS WORD LENGTH
figure;
plot(bits,sqnrR,'o-','LineWidth',1.5);
hold on;
plot(bits,sqnrT,'s--','LineWidth',1.5);
grid on;
xlabel('Word Length (bits)');
ylabel('SQNR (dB)');
title('SQNR vs Word Length');
legend('Rounding','Truncation');

%% 11. QUANTIZATION ERROR VS WORD LENGTH
figure;
semilogy(bits,maxER,'o-','LineWidth',1.5);
hold on;
semilogy(bits,maxET,'s--','LineWidth',1.5);
grid on;
xlabel('Word Length (bits)');
ylabel('Maximum Quantization Error');
title('Quantization Error vs Word Length');
legend('Rounding','Truncation');

%% 12. FINAL RESULT
fprintf('\n========== FINAL RESULT ==========\n');
fprintf('Uniform quantization using rounding and truncation was performed.\n');
fprintf('Quantization error and SQNR were calculated.\n');
fprintf('Increasing word length reduces quantization error and improves SQNR.\n');
fprintf('Finite word-length effects are therefore verified.\n');