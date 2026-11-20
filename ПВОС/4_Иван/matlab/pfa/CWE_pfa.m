clear; close all; clc

%% Параметры варианта Мартыненко И.А. (ЭРэ-19-22)
INR_dB = 15;
Pfa = 1e-5;
nTrial = 1e5;
N = 16;
nGuardCell = 1;

%% Environment Generation
INR = db2pow(INR_dB);
noisePwr = 1; % Noise Power
clutterPwr = noisePwr*INR;
nWindow = N + 2*nGuardCell + 1;
pwrNoise = exprnd(noisePwr,nTrial,nWindow);
pwrClutter = exprnd(clutterPwr,nTrial,nWindow);
iCUT = N/2 + nGuardCell + 1;
iLead = 1:N/2;
iTest = N/2 + nGuardCell + 1;
iLag = N/2 + 2*nGuardCell + 2:nWindow;
iClutterWall = [0 iLead iTest iLag];
nClutterWall = length(iClutterWall);

%% Simulation
threshold = zeros(nTrial,1);
expPfa = zeros(nClutterWall,1);
h = waitbar(0,'Calculating...');

for i = 1:nClutterWall
    waitbar(i/nClutterWall);
    clutterWall = iClutterWall(i);
    if (clutterWall == 0)
        pwrSignal = pwrNoise;
        pwrCUT = noisePwr;
    elseif (clutterWall == nWindow)
        pwrSignal = pwrClutter;
        pwrCUT = clutterPwr;
    else
        iClutter = 1:clutterWall;
        iNoise = clutterWall + 1:nWindow;
        pwrSignal(:,iNoise) = pwrNoise(:,iNoise);
        pwrSignal(:,iClutter) = pwrClutter(:,iClutter);
        if (clutterWall >= iCUT)
            pwrCUT = clutterPwr;
        else
            pwrCUT = noisePwr;
        end
    end

    threshold(:,1) = SOCA_CFAR1(Pfa,N,nGuardCell,pwrSignal);

    expPfa(i,:) = monteCarloFalseAlarm(threshold,pwrCUT);
end
close(h) % close progress bar

%% Plotting %
cellLead = 1:N/2 + 1;
iLead = 0:N/2;
cellLag = N/2 + 2:N + 2;
iLag = N/2:N;
xMin = 0;
xMax = N;
yMin = 10^-8;
yMax = 10^-0;

figure
semilogy(iLead,expPfa(cellLead,1),'linewidth',1.5)
hold on
semilogy(iLag,expPfa(cellLag,1),'linewidth',1.5)
line([xMin xMax],[Pfa Pfa],'color','m','linestyle',':','linewidth',1.5)
line([N/2 N/2],[yMin yMax],'color','b','linestyle',':','linewidth',1.5)
ylabel('Вероятность ложной тревоги')
xlabel('Число опорных ячеек в помехе')
title('Вероятность ложной тревоги МСЗ на границе помех')
axis([xMin xMax yMin yMax])
grid on
legend('СУЛТ-МСЗ','Заданная ВЛТ','Граница окна')
hold off
set(gca,'FontName','Times New Roman','FontSize',18)
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'CWE_pfa.png'));