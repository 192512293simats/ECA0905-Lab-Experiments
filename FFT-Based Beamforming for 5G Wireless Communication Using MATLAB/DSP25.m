%% FFT-Based Beamforming for 5G (28 GHz, 8-element ULA)
% Base MATLAB implementation. No Phased Array System Toolbox required.
% Carrier frequency: 28 GHz | Elements: 8 | Spacing: lambda/2
% Desired direction: +20 deg | Interferers: -30 deg and +45 deg
% Save as DSP25.m and run as a script.
%
% Notes:
% - Uses spatial FFT and a Capon/MVDR reference spectrum.
% - The local functions are placed at the end of the script.
% - xline is replaced by plot-based reference lines for compatibility.

clc; clear; close all;
rng(1);

%% 1. INPUT PARAMETERS
fc = 28e9;
c = 3e8;
lambda = c/fc;
N = 8;
d = lambda/2;
thD = 20;
thI = [-30 45];
steerCases = [-30 0 20 45];
K = 2000;
SNRdB = 10;
INRdB = 10;
noiseVar = 1;
Nfft = 2048;

n = (0:N-1).';
a = @(t) exp(1j*n*(2*pi*d/lambda)*sind(t));
th = -90:0.05:90;
A = a(th);

%% 2. ARRAY DATA MODEL
Pd = noiseVar*10^(SNRdB/10);
Pi = noiseVar*10^(INRdB/10);
cg = @(M,L) (randn(M,L) + 1j*randn(M,L))/sqrt(2);

Xd = a(thD)*(sqrt(Pd)*cg(1,K));
Xi = a(thI)*(sqrt(Pi)*cg(numel(thI),K));
Xn = sqrt(noiseVar)*cg(N,K);
X = Xd + Xi + Xn;

%% 3. FFT SPATIAL SPECTRUM AND CAPON REFERENCE
% FFT of the array dimension estimates spatial frequency.
S = fftshift(mean(abs(fft(X,Nfft,1)).^2,2));
f = (-Nfft/2:Nfft/2-1).'/Nfft;
u = f*lambda/d;
u = max(-1,min(1,u));             % protect asind against round-off
thS = asind(u);
SdB = 10*log10(S/max(S) + eps);

R = (X*X')/K;
Rloaded = R + 1e-3*trace(R)/N*eye(N);
Rinv = Rloaded \ eye(N);
denom = real(sum(conj(A).*(Rinv*A),1));
Pcap = 1./max(denom,eps);
CdB = 10*log10(Pcap/max(Pcap) + eps);

pk = top_peaks(SdB,3);
doaFFT = sort(thS(pk)).';

%% 4. BEAMFORMER WEIGHTS
aD = a(thD);
w_conv = conj(aD)/N;

taper = 0.54 - 0.46*cos(2*pi*n/(N-1));
w_hamm = taper.*conj(aD)/sum(taper);

% MVDR weights, conjugated because output/pattern below use non-conjugate transpose.
w_std = (Rloaded\aD)/(aD'*(Rloaded\aD));
w_mvdr = conj(w_std);

AF_conv = 20*log10(abs(w_conv.'*A) + 1e-6);
AF_hamm = 20*log10(abs(w_hamm.'*A) + 1e-6);
AF_mvdr = 20*log10(abs(w_mvdr.'*A) + 1e-6);

% FFT-based evaluation of conventional array factor.
AF_fft = fftshift(ifft(w_conv,Nfft)*Nfft);
AFfftdB = 20*log10(abs(AF_fft) + 1e-6);

%% 5. STEERING CASES
AFsteer = zeros(numel(steerCases),numel(th));
peakAng = zeros(1,numel(steerCases));
for i = 1:numel(steerCases)
    w = conj(a(steerCases(i)))/N;
    AFsteer(i,:) = 20*log10(abs(w.'*A) + 1e-6);
    [~,ix] = max(AFsteer(i,:));
    peakAng(i) = th(ix);
end

%% 6. PERFORMANCE METRICS
mC = beam_metrics(th,AF_conv);
mH = beam_metrics(th,AF_hamm);
[~,ipk] = max(AFfftdB);
steerFFT = thS(ipk);

gainAt = @(w,t) 20*log10(abs(w.'*a(t)) + eps);
W = {w_conv,w_hamm,w_mvdr};
names = {'Conventional (FFT)','Hamming taper','MVDR (adaptive)'};

SINR_in = 10*log10(mean(abs(Xd(:)).^2) / ...
    (mean(abs(Xi(:)).^2) + mean(abs(Xn(:)).^2)));
SINR_out = zeros(1,3);
for i = 1:3
    SINR_out(i) = sinr_db(W{i},Xd,Xi,Xn);
end

%% 7. FIGURE 1: SPATIAL SPECTRUM
figure('Name','Experiment 25A: Spatial Spectrum','Color','w', ...
    'Position',[40 40 950 500]);
plot(thS,SdB,'b','LineWidth',1.5); hold on;
plot(th,CdB,'r--','LineWidth',1.1);
yl = [-40 3];
plot([thD thD],yl,'g-','LineWidth',1.2);
plot([thI(1) thI(1)],yl,'k:','LineWidth',1.2);
plot([thI(2) thI(2)],yl,'k:','LineWidth',1.2);
plot(thS(pk),SdB(pk),'bv','MarkerFaceColor','b','MarkerSize',7);
grid on; xlim([-90 90]); ylim(yl);
xlabel('Arrival angle \theta (deg)');
ylabel('Normalised power (dB)');
title('FFT-based spatial spectrum of the 8-element array data');
legend('FFT spatial spectrum','Capon reference','Desired +20 deg', ...
    'Interferer -30 deg','Interferer +45 deg','Detected FFT peaks', ...
    'Location','southoutside');

%% 8. FIGURE 2: BEAM PATTERN
figure('Name','Experiment 25B: Beam Pattern','Color','w', ...
    'Position',[80 80 1150 480]);
subplot(1,2,1);
plot(th,AF_conv,'b','LineWidth',1.5); hold on;
plot(thS(1:25:end),AFfftdB(1:25:end),'ro','MarkerSize',3);
yl = [-50 5];
plot([thD thD],yl,'g-','LineWidth',1.2);
plot([thI(1) thI(1)],yl,'k:','LineWidth',1.2);
plot([thI(2) thI(2)],yl,'k:','LineWidth',1.2);
grid on; xlim([-90 90]); ylim(yl);
xlabel('\theta (deg)'); ylabel('|AF| (dB)');
title(sprintf('Beam steered to %d deg (FFT peak %.1f deg)',thD,steerFFT));
legend('Direct array-factor formula','FFT-based array factor', ...
    'Desired direction','Interference directions','Location','southoutside');

subplot(1,2,2);
polarplot(deg2rad(th),max(AF_conv,-50),'b','LineWidth',1.4); hold on;
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
thetalim([-90 90]); rlim([-50 0]);
polarplot(deg2rad([thD thD]),[-50 0],'g','LineWidth',1.2);
polarplot(deg2rad([thI(1) thI(1)]),[-50 0],'k:','LineWidth',1.2);
polarplot(deg2rad([thI(2) thI(2)]),[-50 0],'k:','LineWidth',1.2);
title('Polar beam pattern (dB)');

%% 9. FIGURE 3: BEAM STEERING COMPARISON
figure('Name','Experiment 25C: Beam Steering','Color','w', ...
    'Position',[120 120 950 480]);
plot(th,AFsteer.','LineWidth',1.5);
grid on; xlim([-90 90]); ylim([-50 5]);
xlabel('\theta (deg)'); ylabel('|AF| (dB)');
title('Beam steering using phase-shift weights');
legend(arrayfun(@(s)sprintf('Steer %d deg',s),steerCases, ...
    'UniformOutput',false),'Location','southoutside');

%% 10. FIGURE 4: INTERFERENCE SUPPRESSION AND SINR
figure('Name','Experiment 25D: Interference Suppression','Color','w', ...
    'Position',[160 160 1150 480]);
subplot(1,2,1);
plot(th,AF_conv,'b','LineWidth',1.4); hold on;
plot(th,AF_hamm,'m','LineWidth',1.4);
plot(th,AF_mvdr,'r','LineWidth',1.4);
yl = [-70 5];
plot([thD thD],yl,'g-','LineWidth',1.2);
plot([thI(1) thI(1)],yl,'k:','LineWidth',1.2);
plot([thI(2) thI(2)],yl,'k:','LineWidth',1.2);
grid on; xlim([-90 90]); ylim(yl);
xlabel('\theta (deg)'); ylabel('|AF| (dB)');
title('Beam patterns of three beamformers');
legend(names,'Location','southoutside');

subplot(1,2,2);
bar([SINR_in SINR_out]); grid on;
set(gca,'XTick',1:4,'XTickLabel', ...
    {'Input','Conventional','Hamming','MVDR'});
ylabel('SINR (dB)'); title('Input and output SINR');
v = [SINR_in SINR_out];
for i = 1:4
    text(i,v(i)+0.6,sprintf('%.1f',v(i)), ...
        'HorizontalAlignment','center');
end

%% 11. OUTPUT PARAMETER TABLE
fprintf('\n================ SETUP ================\n');
fprintf('Carrier frequency       : %.0f GHz\n',fc/1e9);
fprintf('Wavelength              : %.3f mm\n',lambda*1e3);
fprintf('Element spacing         : %.3f mm\n',d*1e3);
fprintf('Number of elements      : %d\n',N);
fprintf('Desired direction       : %d deg\n',thD);
fprintf('Interferer directions   : %d deg and %d deg\n',thI(1),thI(2));
fprintf('Desired SNR             : %d dB\n',SNRdB);
fprintf('Interferer INR          : %d dB each\n',INRdB);

fprintf('\n=========== OUTPUT PARAMETER TABLE ===========\n');
fprintf('FFT detected peak angles (deg): %s\n',mat2str(round(doaFFT,1)));
fprintf('True source directions (deg)  : %s\n',mat2str(sort([thD thI])));
fprintf('Conventional beam peak        : %.2f deg (target %d deg)\n',mC.peak,thD);
fprintf('Peak from FFT beam pattern    : %.2f deg\n',steerFFT);

fprintf('\n--- Beam pattern characteristics ---\n');
fprintf('%-22s %-12s %-24s %-12s\n', ...
    'Beamformer','HPBW (deg)','Approx. nulls (deg)','PSLL (dB)');
fprintf('%-22s %-12.2f [%6.2f, %6.2f] %-12.2f\n', ...
    names{1},mC.hpbw,mC.nulls(1),mC.nulls(2),mC.psl);
fprintf('%-22s %-12.2f [%6.2f, %6.2f] %-12.2f\n', ...
    names{2},mH.hpbw,mH.nulls(1),mH.nulls(2),mH.psl);

fprintf('\n--- Steering cases ---\n');
for i = 1:numel(steerCases)
    fprintf('Command %4d deg -> peak at %6.2f deg\n', ...
        steerCases(i),peakAng(i));
end

fprintf('\n--- Array response (dB) at source directions ---\n');
fprintf('%-22s %-14s %-14s %-14s\n', ...
    'Beamformer','Desired +20','Interferer -30','Interferer +45');
for i = 1:3
    fprintf('%-22s %-14.2f %-14.2f %-14.2f\n',names{i}, ...
        gainAt(W{i},thD),gainAt(W{i},thI(1)),gainAt(W{i},thI(2)));
end

fprintf('\n--- SINR results ---\n');
fprintf('Input SINR (per element): %.2f dB\n',SINR_in);
for i = 1:3
    fprintf('%-22s: %.2f dB (change %.2f dB)\n', ...
        names{i},SINR_out(i),SINR_out(i)-SINR_in);
end
fprintf('\nRESULT: FFT-based spatial spectrum, beam steering, beamformer comparison, and SINR analysis completed.\n');

%% LOCAL FUNCTIONS
function idx = top_peaks(y,num)
    y = y(:);
    if numel(y) < 3
        idx = (1:numel(y)).';
        return;
    end
    loc = find(y(2:end-1) > y(1:end-2) & ...
               y(2:end-1) >= y(3:end)) + 1;
    if isempty(loc)
        [~,order] = sort(y,'descend');
        idx = order(1:min(num,numel(order)));
        return;
    end
    [~,order] = sort(y(loc),'descend');
    idx = loc(order(1:min(num,numel(loc))));
end

function m = beam_metrics(th,dB)
    dB = dB - max(dB);
    [~,ip] = max(dB);

    % Half-power beamwidth, using the -3 dB crossings.
    left = ip;
    while left > 1 && dB(left) >= -3
        left = left - 1;
    end
    right = ip;
    while right < numel(th) && dB(right) >= -3
        right = right + 1;
    end
    m.hpbw = th(right) - th(left);

    % Approximate first local minima on either side of the main peak.
    ln = ip;
    while ln > 2 && dB(ln-1) <= dB(ln)
        ln = ln - 1;
    end
    rn = ip;
    while rn < numel(th)-1 && dB(rn+1) <= dB(rn)
        rn = rn + 1;
    end
    m.nulls = [th(ln) th(rn)];

    side = dB;
    side(max(1,ln):min(numel(th),rn)) = -Inf;
    m.psl = max(side);
    m.peak = th(ip);
end

function s = sinr_db(w,Xd,Xi,Xn)
    Ps = mean(abs(w.'*Xd).^2);
    Pint = mean(abs(w.'*Xi).^2);
    Pn = mean(abs(w.'*Xn).^2);
    s = 10*log10(Ps/(Pint + Pn + eps));
end
