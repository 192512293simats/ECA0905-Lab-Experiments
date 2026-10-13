
clc; clear; close all;

%% EXPERIMENT 13: POLYPHASE DECOMPOSITION OF DECIMATION FILTER

%% 1. INPUT PARAMETERS
N = 24;                  % FIR filter order
M = 2;                   % Decimation factor
Fs = 8000;               % Sampling frequency in Hz
N_samples = 128;         % Number of input samples

% Design FIR prototype filter
b = fir1(N, 0.4);

% Polyphase decomposition
E0 = b(1:M:end);         % Even-indexed polyphase branch
E1 = b(2:M:end);         % Odd-indexed polyphase branch

% Input signal: 500 Hz sinusoid
n = 0:N_samples-1;
t = n/Fs;
x = sin(2*pi*500*t);

%% 2. CONVENTIONAL DECIMATION
% Filter first, then downsample
y_direct_full = filter(b, 1, x);
y_direct = y_direct_full(1:M:end);

%% 3. POLYPHASE DECIMATION
% Downsample input into even and odd branches
x0 = x(1:M:end);

% Delay input by one sample before downsampling
x_delayed = [0, x(1:end-1)];
x1 = x_delayed(1:M:end);

% Filter the two polyphase branches
y0 = filter(E0, 1, x0);
y1 = filter(E1, 1, x1);

% Combine branch outputs
y_poly = y0 + y1;

%% 4. VERIFICATION
rec_error = max(abs(y_direct - y_poly));

if rec_error < 1e-12
    verification = 'PASS';
else
    verification = 'FAIL';
end

%% 5. COMPUTATIONAL COMPLEXITY
N_taps = length(b);
N_out = length(y_direct);

% Approximate MAC counts
MAC_direct = N_samples * N_taps;
MAC_poly = N_out * (length(E0) + length(E1));

comp_saving = (1 - MAC_poly/MAC_direct)*100;

%% 6. COMMAND WINDOW OUTPUT
fprintf('\n============================================\n');
fprintf(' EXPERIMENT 13: POLYPHASE DECIMATION\n');
fprintf('============================================\n');

fprintf('Sampling Frequency       : %d Hz\n', Fs);
fprintf('Decimation Factor (M)    : %d\n', M);
fprintf('FIR Filter Order         : %d\n', N);
fprintf('Number of FIR Taps       : %d\n', N_taps);
fprintf('Input Samples            : %d\n', N_samples);
fprintf('Output Samples           : %d\n', N_out);
fprintf('Output Sampling Frequency: %d Hz\n', Fs/M);

fprintf('\nPOLYPHASE COEFFICIENTS\n');
fprintf('E0 has %d coefficients:\n', length(E0));
disp(E0);

fprintf('E1 has %d coefficients:\n', length(E1));
disp(E1);

fprintf('\nCOMPUTATIONAL COMPLEXITY\n');
fprintf('Direct MAC Operations  : %d\n', MAC_direct);
fprintf('Polyphase MAC Operations: %d\n', MAC_poly);
fprintf('Computational Savings  : %.2f%%\n', comp_saving);

fprintf('\nOUTPUT VERIFICATION\n');
fprintf('Maximum Absolute Difference: %.6e\n', rec_error);
fprintf('Verification Status: %s\n', verification);

%% 7. OUTPUT PARAMETER TABLE
Parameter = {
    'FIR Filter Order'
    'Number of FIR Taps'
    'Decimation Factor'
    'Input Signal Length'
    'Output Signal Length'
    'Output Sampling Frequency'
    'Direct MAC Operations'
    'Polyphase MAC Operations'
    'Computational Savings'
    'Maximum Absolute Difference'
};

Value = [
    N
    N_taps
    M
    length(x)
    N_out
    Fs/M
    MAC_direct
    MAC_poly
    comp_saving
    rec_error
];

Unit = {
    'Order'
    'Taps'
    'Factor'
    'Samples'
    'Samples'
    'Hz'
    'MACs'
    'MACs'
    'Percent'
    'Amplitude'
};

Output_Table = table(Parameter, Value, Unit);
disp(Output_Table);

%% 8. FIGURE 13A: PROTOTYPE IMPULSE RESPONSE
figure('Name', 'Experiment 13A: Prototype FIR Filter');

stem(0:N, b, 'filled', 'LineWidth', 1.2);
grid on;
box on;
xlabel('Sample Index (n)');
ylabel('Coefficient Value');
title('Prototype FIR Filter Impulse Response');

%% 9. FIGURE 13B: POLYPHASE COEFFICIENTS
figure('Name', 'Experiment 13B: Polyphase Coefficients');

subplot(2,1,1);
stem(0:length(E0)-1, E0, 'b', 'filled');
grid on;
xlabel('Branch Sample Index');
ylabel('Coefficient');
title('Polyphase Branch E0');

subplot(2,1,2);
stem(0:length(E1)-1, E1, 'r', 'filled');
grid on;
xlabel('Branch Sample Index');
ylabel('Coefficient');
title('Polyphase Branch E1');

%% 10. FIGURE 13C: OUTPUT COMPARISON
figure('Name', 'Experiment 13C: Output Comparison');

n_out = 0:N_out-1;

plot(n_out, y_direct, 'b-', 'LineWidth', 1.5);
hold on;
plot(n_out, y_poly, 'r--', 'LineWidth', 1.5);

grid on;
box on;
xlabel('Output Sample Index');
ylabel('Amplitude');
title('Conventional vs Polyphase Decimation');
legend('Direct Decimation', 'Polyphase Decimation', ...
       'Location', 'best');

%% 11. FIGURE 13D: NUMERICAL ERROR
figure('Name', 'Experiment 13D: Numerical Error');

stem(n_out, abs(y_direct-y_poly), 'filled');
grid on;
xlabel('Output Sample Index');
ylabel('Absolute Error');
title('Direct vs Polyphase Numerical Error');

%% 12. FINAL RESULT
fprintf('\n============================================\n');
fprintf(' FINAL RESULT\n');
fprintf('============================================\n');
fprintf('1. FIR prototype filter designed successfully.\n');
fprintf('2. Polyphase decomposition completed.\n');
fprintf('3. Conventional and polyphase outputs compared.\n');
fprintf('4. Computational savings calculated.\n');

if rec_error < 1e-12
    fprintf('5. Equivalent output verified within tolerance.\n');
else
    fprintf('5. Output mismatch detected; check branch alignment.\n');
end

fprintf('============================================\n');
