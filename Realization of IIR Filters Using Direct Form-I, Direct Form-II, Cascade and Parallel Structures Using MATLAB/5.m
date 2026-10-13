clc;
clear;
close all;
format short;

%% IIR FILTER REALIZATION
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

N = 50;
x = [1 zeros(1,N-1)];
n = 0:N-1;

%% 1. DIRECT FORM-I
y1 = zeros(1,N);

for m = 1:N
    for i = 1:min(m,length(b))
        y1(m) = y1(m) + b(i)*x(m-i+1);
    end

    for i = 2:min(m,length(a))
        y1(m) = y1(m) - a(i)*y1(m-i+1);
    end
end

%% 2. DIRECT FORM-II TRANSPOSED
y2 = filter(b,a,x);

%% 3. CASCADE REALIZATION
[sos,g] = tf2sos(b,a);

yc = g*x;

for i = 1:size(sos,1)
    yc = filter(sos(i,1:3),sos(i,4:6),yc);
end

%% 4. PARALLEL REALIZATION
[r,p,k] = residuez(b,a);

yp = filter(k,1,x);

for i = 1:length(r)
    yp = yp + filter(r(i),[1 -p(i)],x);
end

yp = real(yp);

%% 5. DISPLAY COEFFICIENTS
fprintf('\n========================================\n');
fprintf('       IIR FILTER REALIZATION\n');
fprintf('========================================\n');

fprintf('\nNumerator coefficients b:\n');
disp(b);

fprintf('Denominator coefficients a:\n');
disp(a);

fprintf('\nTransfer Function:\n');
disp(tf(b,a,1));

fprintf('\n========================================\n');
fprintf('       DIRECT FORM-I\n');
fprintf('========================================\n');
disp('DF-I numerator coefficients:');
disp(b);
disp('DF-I denominator coefficients:');
disp(a);

fprintf('\n========================================\n');
fprintf('       DIRECT FORM-II\n');
fprintf('========================================\n');
disp('DF-II numerator coefficients:');
disp(b);
disp('DF-II denominator coefficients:');
disp(a);

fprintf('\n========================================\n');
fprintf('       CASCADE REALIZATION\n');
fprintf('========================================\n');
disp('Second-Order Sections (SOS):');
disp(sos);
fprintf('Overall gain:\n');
disp(g);

fprintf('\n========================================\n');
fprintf('       PARALLEL REALIZATION\n');
fprintf('========================================\n');
disp('Residues:');
disp(r);
disp('Poles:');
disp(p);
disp('Direct term:');
disp(k);

%% 6. REFERENCE IMPULSE RESPONSE
yref = filter(b,a,x);

fprintf('\n========================================\n');
fprintf('       REFERENCE IMPULSE RESPONSE\n');
fprintf('========================================\n');
disp(yref(1:min(10,N)));

%% 7. NUMERICAL ERROR AND DEVIATION
E = [max(abs(yref-y1)); ...
     max(abs(yref-y2)); ...
     max(abs(yref-yc)); ...
     max(abs(yref-yp))];

peakValue = max(abs(yref));

if peakValue > 0
    D = (E/peakValue)*100;
else
    D = zeros(size(E));
end

Structure = {'Direct Form-I'; ...
             'Direct Form-II'; ...
             'Cascade'; ...
             'Parallel'};

Comparison = table(Structure,E,D, ...
    'VariableNames', ...
    {'Structure','Maximum_Error','Deviation_percent'});

fprintf('\n========================================\n');
fprintf('       NUMERICAL COMPARISON\n');
fprintf('========================================\n');
disp(Comparison);

fprintf('\n========================================\n');
fprintf('       PERCENTAGE DEVIATION\n');
fprintf('========================================\n');

fprintf('DF-I     = %.6e %%\n',D(1));
fprintf('DF-II    = %.6e %%\n',D(2));
fprintf('Cascade  = %.6e %%\n',D(3));
fprintf('Parallel = %.6e %%\n',D(4));

%% 8. IMPULSE RESPONSE COMPARISON
figure('Name','Impulse Response of IIR Filter Realizations');

stem(n,y1,'o-','LineWidth',1.2);
hold on;
stem(n,y2,'--o','LineWidth',1.0);
stem(n,yc,':o','LineWidth',1.3);
stem(n,yp,'-.o','LineWidth',1.0);

grid on;
xlabel('n');
ylabel('Amplitude');
title('Impulse Response of IIR Filter Realizations');
legend('Direct Form-I','Direct Form-II', ...
       'Cascade','Parallel','Location','best');
xlim([0 N-1]);
hold off;

%% 9. REFERENCE IMPULSE RESPONSE
figure('Name','Reference Impulse Response');

stem(n,yref,'o','LineWidth',1.3);
grid on;
xlabel('Sample n');
ylabel('h[n]');
title('Impulse Response of IIR Filter');
xlim([0 N-1]);

%% 10. NUMERICAL ERROR PLOT
figure('Name','Numerical Error Between Realizations');

plot(n,abs(yref-y1),'LineWidth',1.2);
hold on;
plot(n,abs(yref-y2),'--','LineWidth',1.2);
plot(n,abs(yref-yc),':','LineWidth',1.5);
plot(n,abs(yref-yp),'-.','LineWidth',1.2);

grid on;
xlabel('Sample n');
ylabel('Absolute Error');
title('Numerical Error Between Realizations');
legend('DF-I Error','DF-II Error', ...
       'Cascade Error','Parallel Error','Location','best');
xlim([0 N-1]);
hold off;

%% 11. FREQUENCY RESPONSE
figure('Name','Frequency Response of IIR Filter');

freqz(b,a,1024);
sgtitle('Frequency Response of IIR Filter');

%% 12. POLE-ZERO PLOT
figure('Name','Pole-Zero Plot');

zplane(b,a);
title('Pole-Zero Plot');

%% 13. COMPLETION MESSAGE
fprintf('\n========================================\n');
fprintf('IIR FILTER REALIZATION COMPLETED\n');
fprintf('All five output figures generated.\n');
fprintf('========================================\n');