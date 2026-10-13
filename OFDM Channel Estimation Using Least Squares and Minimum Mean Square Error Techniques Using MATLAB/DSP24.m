%% OFDM Channel Estimation using LS and MMSE Techniques
% 64 subcarriers, QPSK, 4-tap Rayleigh fading channel
% Corrected standalone MATLAB script (no Communications Toolbox required).
% Save this file as DSP24.m and press Run.
clc; clear; close all;
rng(1);

%% 1. INPUT PARAMETERS
P.N = 64;
P.pilots = [8 22 44 58];             % MATLAB 1-based indices
P.Lcp = 16;
P.Xp = [1+1j; 1-1j; -1+1j; -1-1j]/sqrt(2);
P.dataIdx = setdiff(1:P.N, P.pilots);

SNR_dB = 20;
NsymMain = 1000;                     % Increase to 3000 for longer simulation
NsymSweep = 1500;                    % Increase to 4000 for smoother BER curves
snrSweep = 0:5:30;

% Exponential 4-tap power-delay profile, normalized to unit total power
pdp = exp(-(0:3)/1.5);
pdp = pdp/sum(pdp);
P.pdp = pdp(:);

% Frequency-domain channel correlation matrix
dk = (0:P.N-1).' - (0:P.N-1);
P.R = zeros(P.N);
for ell = 0:numel(pdp)-1
    P.R = P.R + pdp(ell+1)*exp(-1j*2*pi*dk*ell/P.N);
end

%% 2. MAIN RUN AT 20 dB
out = run_ofdm(SNR_dB, NsymMain, P);

%% 3. BER AND MSE SWEEP
berLS = zeros(size(snrSweep));
berMMSE = zeros(size(snrSweep));
berPerf = zeros(size(snrSweep));
mseLS = zeros(size(snrSweep));
mseMMSE = zeros(size(snrSweep));

for n = 1:numel(snrSweep)
    o = run_ofdm(snrSweep(n), NsymSweep, P);
    berLS(n) = o.ber(1);
    berMMSE(n) = o.ber(2);
    berPerf(n) = o.ber(3);
    mseLS(n) = o.mse(1);
    mseMMSE(n) = o.mse(2);
end

% Approximate QPSK BER in a Rayleigh channel with ideal channel knowledge.
% SNR_dB here represents Es/N0, so Eb/N0 = Es/N0 / 2 for QPSK.
gam = 10.^(snrSweep/10)/2;
berTh = 0.5*(1 - sqrt(gam./(1+gam)));

%% 4. FIGURE 1: CHANNEL MAGNITUDE AND PHASE
k = (1:P.N).';
figure('Name','Experiment 24A: Channel Estimates','Color','w');
subplot(1,2,1);
plot(k,abs(out.H1(:,1)),'k','LineWidth',1.5); hold on;
plot(k,abs(out.H1(:,2)),'r--','LineWidth',1.2);
plot(k,abs(out.H1(:,3)),'b-.','LineWidth',1.2);
plot(P.pilots,abs(out.Hpil1),'go','MarkerFaceColor','g','MarkerSize',6);
grid on; xlabel('Subcarrier index'); ylabel('|H(k)|');
title(sprintf('Channel magnitude (SNR = %d dB)',SNR_dB));
legend('Actual','LS (interpolated)','MMSE','LS at pilots','Location','best');

subplot(1,2,2);
plot(k,angle(out.H1(:,1)),'k','LineWidth',1.5); hold on;
plot(k,angle(out.H1(:,2)),'r--','LineWidth',1.2);
plot(k,angle(out.H1(:,3)),'b-.','LineWidth',1.2);
grid on; xlabel('Subcarrier index'); ylabel('Phase (rad)');
title('Channel phase');
legend('Actual','LS','MMSE','Location','best');

%% 5. FIGURE 2: ESTIMATION ERROR
figure('Name','Experiment 24B: Estimation Error','Color','w');
subplot(1,2,1);
semilogy(k,max(out.mse_k(:,1),eps),'r-o','MarkerSize',3); hold on;
semilogy(k,max(out.mse_k(:,2),eps),'b-s','MarkerSize',3);
for p = P.pilots
    xline(p,':k');
end
grid on; xlabel('Subcarrier index'); ylabel('MSE = E|H-H_{est}|^2');
title(sprintf('Per-subcarrier error (SNR = %d dB)',SNR_dB));
legend('LS','MMSE','Location','best');

subplot(1,2,2);
semilogy(snrSweep,max(mseLS,eps),'r-o','LineWidth',1.4); hold on;
semilogy(snrSweep,max(mseMMSE,eps),'b-s','LineWidth',1.4);
grid on; xlabel('SNR (dB)'); ylabel('Channel estimation MSE');
title('Estimation MSE versus SNR');
legend('LS','MMSE','Location','best');

%% 6. FIGURE 3: EQUALIZED CONSTELLATIONS
figure('Name','Experiment 24C: Equalized Constellations','Color','w');
titles = {'Perfect CSI','LS estimate','MMSE estimate'};
names = {'perfect','ls','mmse'};
berIndices = [3 1 2];  % out.ber is [LS MMSE Perfect]
ideal = [1+1j 1-1j -1+1j -1-1j]/sqrt(2);

for ii = 1:3
    subplot(1,3,ii);
    z = out.const.(names{ii});
    plot(real(z),imag(z),'.','MarkerSize',4); hold on;
    plot(real(ideal),imag(ideal),'r+','MarkerSize',10,'LineWidth',1.8);
    grid on; axis equal; axis([-2 2 -2 2]);
    xlabel('In-phase'); ylabel('Quadrature');
    title(sprintf('%s (BER = %.2e)',titles{ii},out.ber(berIndices(ii))));
end

%% 7. FIGURE 4: BER VERSUS SNR
figure('Name','Experiment 24D: BER Comparison','Color','w');
semilogy(snrSweep,max(berLS,1e-5),'r-o','LineWidth',1.4); hold on;
semilogy(snrSweep,max(berMMSE,1e-5),'b-s','LineWidth',1.4);
semilogy(snrSweep,max(berPerf,1e-5),'k-^','LineWidth',1.4);
semilogy(snrSweep,max(berTh,1e-5),'g--','LineWidth',1.2);
grid on; xlabel('SNR (dB)'); ylabel('Bit error rate (BER)');
title('BER: LS versus MMSE channel estimation');
legend('LS','MMSE','Perfect CSI','Theory (Rayleigh, ideal CSI)','Location','southwest');
ylim([1e-5 1]);

%% 8. PRINT OUTPUT TABLES
fprintf('\n========== OUTPUT PARAMETER TABLE (SNR = %d dB) ==========\n',SNR_dB);
fprintf('Channel at pilot subcarriers (first OFDM symbol)\n');
fprintf('%-8s %-22s %-22s %-22s\n','Pilot','Actual H','LS estimate','MMSE estimate');
for ii = 1:numel(P.pilots)
    pp = P.pilots(ii);
    fprintf('%-8d %-22s %-22s %-22s\n',pp, ...
        complex_string(out.H1(pp,1)), ...
        complex_string(out.H1(pp,2)), ...
        complex_string(out.H1(pp,3)));
end

fprintf('\n%-30s %-14s %-14s\n','Metric','LS','MMSE');
fprintf('%-30s %-14.4e %-14.4e\n','Channel estimation MSE',out.mse(1),out.mse(2));
fprintf('%-30s %-14.4e %-14.4e\n','Bit error rate (BER)',out.ber(1),out.ber(2));
fprintf('Perfect-CSI BER: %.4e\n',out.ber(3));
fprintf('Bits simulated at %d dB: %d\n',SNR_dB,2*numel(P.dataIdx)*NsymMain);

fprintf('\n%-10s %-12s %-12s %-12s %-12s %-12s\n', ...
    'SNR (dB)','BER LS','BER MMSE','BER perfect','MSE LS','MSE MMSE');
for n = 1:numel(snrSweep)
    fprintf('%-10d %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e\n', ...
        snrSweep(n),berLS(n),berMMSE(n),berPerf(n),mseLS(n),mseMMSE(n));
end

fprintf('\nRESULT: OFDM LS and MMSE channel estimation, equalization, BER and MSE analysis completed.\n');

%% LOCAL FUNCTIONS (keep these at the end of this script)
function s = complex_string(z)
    s = sprintf('%.3f%+.3fj',real(z),imag(z));
end

function out = run_ofdm(snrdB,Nsym,P)
    N = P.N;
    pil = P.pilots(:);
    Np = numel(pil);
    dIdx = P.dataIdx(:);
    Nd = numel(dIdx);
    L = numel(P.pdp);
    sig2 = 10^(-snrdB/10);   % complex noise variance for unit average symbol power

    % Transmitter: random QPSK data and known pilot symbols
    b1 = randi([0 1],Nd,Nsym);
    b2 = randi([0 1],Nd,Nsym);
    dataSym = ((1-2*b1) + 1j*(1-2*b2))/sqrt(2);
    X = zeros(N,Nsym);
    X(dIdx,:) = dataSym;
    X(pil,:) = repmat(P.Xp(:),1,Nsym);

    % OFDM modulation and cyclic prefix
    x = sqrt(N)*ifft(X,N,1);
    xcp = [x(N-P.Lcp+1:N,:); x];

    % Independent 4-tap Rayleigh channel for each OFDM symbol
    h = sqrt(P.pdp/2).*(randn(L,Nsym)+1j*randn(L,Nsym));
    y = zeros(N+P.Lcp,Nsym);
    for ss = 1:Nsym
        y(:,ss) = filter(h(:,ss),1,xcp(:,ss));
    end
    y = y + sqrt(sig2/2)*(randn(size(y))+1j*randn(size(y)));

    % Remove cyclic prefix and transform to frequency domain
    Y = fft(y(P.Lcp+1:end,:),N,1)/sqrt(N);
    Hact = fft(h,N,1);

    % LS channel estimate at pilot subcarriers, then linear interpolation
    Hls_p = Y(pil,:)./repmat(P.Xp(:),1,Nsym);
    Hls = interp1(pil,Hls_p,(1:N).','linear','extrap');

    % LMMSE estimate using the known channel correlation matrix
    Rpp = P.R(pil,pil);
    Rfp = P.R(:,pil);
    W = Rfp/(Rpp + sig2*eye(Np));
    Hmmse = W*Hls_p;

    % One-tap zero-forcing equalization on data subcarriers
    Yd = Y(dIdx,:);
    Xe.perfect = Yd./Hact(dIdx,:);
    Xe.ls = Yd./Hls(dIdx,:);
    Xe.mmse = Yd./Hmmse(dIdx,:);

    % BER calculation for QPSK: negative real/imaginary part means bit 1
    nbits = 2*Nd*Nsym;
    ber = zeros(1,3);
    fieldNames = {'perfect','ls','mmse'};
    for jj = 1:numel(fieldNames)
        z = Xe.(fieldNames{jj});
        bitErrors = sum(sum((real(z)<0) ~= (b1==1))) + ...
                    sum(sum((imag(z)<0) ~= (b2==1)));
        ber(jj) = bitErrors/nbits;
    end
    out.ber = [ber(2),ber(3),ber(1)];  % [LS MMSE Perfect]

    % Channel estimation MSE
    eLS = abs(Hact-Hls).^2;
    eMMSE = abs(Hact-Hmmse).^2;
    out.mse = [mean(eLS(:)),mean(eMMSE(:))];
    out.mse_k = [mean(eLS,2),mean(eMMSE,2)];

    % Save first-symbol channel curves and pilot estimates for plots
    out.H1 = [Hact(:,1),Hls(:,1),Hmmse(:,1)];
    out.Hpil1 = Hls_p(:,1);

    % Store a subset of equalized symbols for constellation plots
    nPlot = min(300,Nsym);
    out.const.perfect = reshape(Xe.perfect(:,1:nPlot),[],1);
    out.const.ls = reshape(Xe.ls(:,1:nPlot),[],1);
    out.const.mmse = reshape(Xe.mmse(:,1:nPlot),[],1);
end
