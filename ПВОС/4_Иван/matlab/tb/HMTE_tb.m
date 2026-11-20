clear; close all; clc;

%% Параметры варианта Мартыненко И.А. (ЭРэ-19-22)
nRangeGate = 1000;
N = 16;
nGuardCell = 1;
NM = 3;
INR_dB = 0;
SINR_dB = 20;
Pfa = 1e-5;

%% Три маскирующие цели и проверяемая цель
iCUT = nRangeGate/2;
iMasking = iCUT + [-10 10 20];
iSignalPosition = [iMasking iCUT];
data = simRangeData(INR_dB,nRangeGate,SINR_dB*ones(1,NM+1),iSignalPosition);
detector = phased.CFARDetector('NumTrainingCells',N, ...
    'NumGuardCells',2*nGuardCell, ...
    'ProbabilityFalseAlarm',Pfa, ...
    'ThresholdOutputPort',true,'Method','SOCA');
[detection,threshold] = detector(data,1:nRangeGate);

figure(1)
plot(pow2db(data),'linewidth',1.5)
hold on
plot(pow2db(threshold),'linewidth',1.5)
stairs(detection*5 - 40,'linewidth',1.5,'color','k')
grid on
axis tight
xlabel('Дальность, отсчёты')
ylabel('Мощность, дБ')
legend('Данные после квадратичного детектора','Порог МСЗ','Обнаружение')
hold off
ylim([-40 SINR_dB + 5])
set(gca,'FontName','Times New Roman','FontSize',18);
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'HMTE_tb.png'));