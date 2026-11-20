clear; close all; clc

%% Параметры варианта Мартыненко И.А. (ЭРэ-19-22)
SNR_dB = 20;
INR_dB = 15;
Pfa = 1e-5;
nTrial = 1e5;
N = 16;
nGuardCell = 1;

%% Формирование окна обнаружителя
SNR = db2pow(SNR_dB);
INR = db2pow(INR_dB);
noisePwr = 1;
clutterPwr = noisePwr*INR;
nWindow = N + 2*nGuardCell + 1;
targetPwr = noisePwr*SNR;
pwrNoise = exprnd(noisePwr,nTrial,nWindow);
pwrClutter = exprnd(clutterPwr,nTrial,nWindow);
iCUT = N/2 + nGuardCell + 1;
iLead = 1:N/2;
iTest = iCUT;
iLag = N/2 + 2*nGuardCell + 2:nWindow;
iClutterWall = [0 iLead iTest iLag];
nClutterWall = length(iClutterWall);

%% Моделирование
threshold = zeros(nTrial,1);
expPd = zeros(nClutterWall,1);
h = waitbar(0,'Вычисление...');

for i = 1:nClutterWall
    clutterWall = iClutterWall(i);
    waitbar(i/nClutterWall,h);
    if clutterWall == 0
        pwrSignal = pwrNoise;
        pwrCUT = noisePwr;
    elseif clutterWall == nWindow
        pwrSignal = pwrClutter;
        pwrCUT = clutterPwr;
    else
        iClutter = 1:clutterWall;
        iNoise = clutterWall + 1:nWindow;
        pwrSignal(:,iNoise) = pwrNoise(:,iNoise);
        pwrSignal(:,iClutter) = pwrClutter(:,iClutter);
        if clutterWall >= iCUT
            pwrCUT = clutterPwr;
        else
            pwrCUT = noisePwr;
        end
    end

    threshold(:,1) = SOCA_CFAR1(Pfa,N,nGuardCell,pwrSignal);
    if i <= nClutterWall/2
        pwrTarget = exprnd(noisePwr + targetPwr,nTrial,1);
    else
        pwrTarget = exprnd(clutterPwr + targetPwr,nTrial,1);
    end
    expPd(i) = sum(pwrTarget > threshold(:,1))/nTrial;
end
close(h)

%% График вероятности правильного обнаружения
cellLead = 1:N/2 + 1;
iLead = 0:N/2;
cellLag = N/2 + 2:N + 2;
iLag = N/2:N;
xMin = 0;
xMax = N;
yMin = 0;
yMax = 1;

figure
hold on
plot(iLead,expPd(cellLead),'linewidth',1.5)
line([N/2 N/2],[yMin yMax],'color','m','linestyle',':','linewidth',1.5)
plot(iLag,expPd(cellLag),'linewidth',1.5)
legend('СУЛТ-МСЗ','Граница окна','Location','best')
title('Вероятность правильного обнаружения на границе помех')
ylabel('Вероятность правильного обнаружения')
xlabel('Число опорных ячеек в помехе')
axis([xMin xMax yMin yMax])
grid on
box on
set(gca,'FontName','Times New Roman','FontSize',18)
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'CWE_pd.png'));