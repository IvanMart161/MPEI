clear; close all; clc;

%% Параметры варианта Мартыненко И.А. (ЭРэ-19-22)
Pfa = 1e-5;
nTrial = 1e6;
SNRdB_set = 1:1:30;
N = 16;
nGuardCell = 1;
NM = 3;

%% CFAR Environment Generation
SNR_set = db2pow(SNRdB_set); % Signal to noise ratio
nSNR = length(SNRdB_set); % Number of samples in SNR set
noisePwr = 1; % Noise Power
nWindow = N + 2*nGuardCell + 1;
iCUT = N/2 + nGuardCell + 1;
iMasking = [iCUT - nGuardCell - 1, ...
    iCUT + nGuardCell + 1, iCUT + nGuardCell + 2];
nMasking = NM;
iTarget = [iMasking iCUT];
nTarget = length(iTarget);
pwrSignal = exprnd(noisePwr,nTrial,nWindow);

%% Simulation
threshold = zeros(nTrial,1);
Pd = zeros(nSNR,1);
h = waitbar(0,'Calculating...');

for iSNR = 1:nSNR
    
    waitbar(iSNR/nSNR);
    SNR = SNR_set(iSNR);
    targetPwr = noisePwr*SNR;
    pwrSignal(:,iTarget) = exprnd(targetPwr + noisePwr,nTrial,nTarget);
    
    threshold(:,1) = SOCA_CFAR1(Pfa,N,nGuardCell,pwrSignal);
    Pd(iSNR) = sum(pwrSignal(:,iCUT) > threshold(:,1))/nTrial;
end
close(h); % close progress bar

%% Plotting
SNRdB_set = SNRdB_set(:);

%% PROBABILITY OF DETECTION CURVE
figure
    set(gca,'FontName','Times New Roman','FontSize',24);
    plot(SNRdB_set,Pd_inf(Pfa,SNRdB_set),'k-','linewidth',1.5);
    hold on
    plot(SNRdB_set,Pd,'linewidth',1.5)
    legend('Нейман—Пирсон','СУЛТ-МСЗ','Location','northwest')
    ylabel('Вероятность правильного обнаружения')
    xlabel('ОСШ, дБ')
    title(['Обнаружение близкорасположенных целей: NM = ' num2str(NM)])
    grid on
    box on
    set(gca,'FontName','Times New Roman','FontSize',18)
    saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'HMTE_pd.png'));
% CFAR LOSS CURVE

SNRdB_set = SNRdB_set(:);
CFARloss = zeros(nSNR,1);

figure
    set(gca,'FontName','Times New Roman','FontSize',24);
    hold on
    
    SNRoptimal = OptSNR(Pd,Pfa);
    CFARloss = SNRdB_set - SNRoptimal;
    plot(SNRdB_set,CFARloss,'linewidth',1.5)
    legend('СУЛТ-МСЗ','Location','northwest')
    ylabel('Потери обнаружения, дБ')
    xlabel('ОСШ, дБ')
    xMin = min(SNRdB_set);
    xMax = max(SNRdB_set);
    yMin = 0;
    yMax = ceil(2*max(max(CFARloss)))/2;
    axis([xMin xMax yMin yMax]);
    title('Потери при обнаружении близкорасположенных целей')
    grid on;
    box on;
    set(gca,'FontName','Times New Roman','FontSize',18)
    saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'HMTE_pd_loss.png'));