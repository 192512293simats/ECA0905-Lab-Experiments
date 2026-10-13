clc; clear; close all;

%% INPUT PARAMETERS
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

Bc = 8;              % Coefficient quantization bits
Bp = 16;             % Product quantization bits
N = 200;

x = [1 zeros(1,N-1)]; 
n = 0:N-1;

%% 1. COEFFICIENT QUANTIZATION
S = 2^(Bc-1)-1;

bq = round(b*S)/S;
aq = round(a*S)/S;
aq = aq/aq(1);

CoeffError = [b-bq, a-aq];

%% 2. IDEAL AND QUANTIZED FILTER RESPONSE
y_ideal = filter(b,a,x);
y_coeff = filter(bq,aq,x);

e_coeff = y_ideal-y_coeff;

MaxCoeffError = max(abs(e_coeff));
SQNR_coeff = 10*log10(sum(y_ideal.^2)/sum(e_coeff.^2));

%% 3. PRODUCT QUANTIZATION
Sp = 2^(Bp-1)-1;

y_prod = zeros(1,N);
xstate = zeros(1,length(b));
ystate = zeros(1,length(a)-1);

for k = 1:N
    xstate = [x(k), xstate(1:end-1)];

    feedforward = round((bq.*xstate)*Sp)/Sp;
    feedback = round((-aq(2:end).*ystate)*Sp)/Sp;

    accumulator = sum(feedforward)+sum(feedback);
    y_prod(k) = round(accumulator*Sp)/Sp;

    ystate = [y_prod(k), ystate(1:end-1)];
end

e_prod = y_coeff-y_prod;

MaxProdError = max(abs(e_prod));
SQNR_prod = 10*log10(sum(y_coeff.^2)/sum(e_prod.^2));

%% 4. ZERO-INPUT LIMIT-CYCLE TEST
y_zero = zeros(1,N);
state = [0.01 0.01 0.01 0.01];

for k = 1:N
    products = round((-aq(2:end).*state)*Sp)/Sp;
    acc = sum(products);

    y_zero(k) = round(acc*Sp)/Sp;
    state = [y_zero(k), state(1:end-1)];
end

LimitAmplitude = max(abs(y_zero));
LimitSamples = find(abs(y_zero)>1e-12);

%% 5. OVERFLOW AND SATURATION TEST
xo = 1.2*ones(1,N);
y_over = filter(bq,aq,xo);

y_sat = min(max(y_over,-1),1);

OverflowCount = sum(abs(y_over)>1);
MaxOverflow = max(abs(y_over));

%% 6. DISPLAY RESULTS
fprintf('\n========================================\n');
fprintf(' FINITE WORD-LENGTH EFFECTS IN IIR FILTERS\n');
fprintf('========================================\n');

fprintf('\n--- QUANTIZATION RESULTS ---\n');
fprintf('Coefficient Word Length = %d bits\n',Bc);
fprintf('Product Word Length = %d bits\n',Bp);
fprintf('Maximum Coefficient Response Error = %.6e\n',MaxCoeffError);
fprintf('Coefficient SQNR = %.2f dB\n',SQNR_coeff);
fprintf('Maximum Product Quantization Error = %.6e\n',MaxProdError);
fprintf('Product SQNR = %.2f dB\n',SQNR_prod);
fprintf('Zero-Input Maximum Amplitude = %.6e\n',LimitAmplitude);
fprintf('Number of Nonzero Limit-Cycle Samples = %d\n',length(LimitSamples));
fprintf('Overflow Samples = %d\n',OverflowCount);
fprintf('Maximum Output Magnitude = %.6f\n',MaxOverflow);

%% 7. OUTPUT PARAMETER TABLE
Test = {
    'Coefficient Quantization'
    'Product Quantization'
    'Zero-Input Test'
    'Overflow Test'
};

MaximumValue = [
    MaxCoeffError
    MaxProdError
    LimitAmplitude
    MaxOverflow
];

SQNR_dB = [
    SQNR_coeff
    SQNR_prod
    NaN
    NaN
];

T = table(Test,MaximumValue,SQNR_dB);
disp(T);

%% 8. DISPLAY COEFFICIENTS
fprintf('\nOriginal Numerator Coefficients:\n');
disp(b);

fprintf('Quantized Numerator Coefficients:\n');
disp(bq);

fprintf('Original Denominator Coefficients:\n');
disp(a);

fprintf('Quantized Denominator Coefficients:\n');
disp(aq);

fprintf('Coefficient Quantization Errors:\n');
disp(CoeffError);

%% 9. IDEAL VS QUANTIZED RESPONSE
figure;
stem(n,y_ideal,'b','filled');
hold on;
stem(n,y_coeff,'r--');
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Ideal vs Coefficient-Quantized Response');
legend('Ideal','Quantized');

%% 10. COEFFICIENT QUANTIZATION ERROR
figure;
plot(n,e_coeff,'LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Error');
title('Coefficient Quantization Noise');

%% 11. PRODUCT QUANTIZATION ERROR
figure;
stem(n,e_prod,'filled');
grid on;
xlabel('Sample n');
ylabel('Error');
title('Product Quantization Noise');

%% 12. ZERO-INPUT LIMIT-CYCLE TEST
figure;
plot(n,y_zero,'LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Zero-Input Limit-Cycle Test');

%% 13. OVERFLOW AND SATURATION
figure;
plot(n,y_over,'b','LineWidth',1.2);
hold on;
plot(n,y_sat,'r--','LineWidth',1.2);
grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Overflow and Saturation Test');
legend('Without Saturation','With Saturation');

%% 14. IDEAL FILTER POLE-ZERO PLOT
figure;
zplane(b,a);
grid on;
title('Ideal Filter Pole-Zero Plot');

%% 15. QUANTIZED FILTER POLE-ZERO PLOT
figure;
zplane(bq,aq);
grid on;
title('Quantized Filter Pole-Zero Plot');

%% 16. FINAL RESULT
fprintf('\n========================================\n');
fprintf(' FINAL RESULT\n');
fprintf('========================================\n');
fprintf('Coefficient quantization was performed.\n');
fprintf('Product quantization was simulated.\n');
fprintf('Zero-input behavior was analyzed.\n');
fprintf('Overflow and saturation were evaluated.\n');
fprintf('Finite word-length effects were demonstrated.\n');