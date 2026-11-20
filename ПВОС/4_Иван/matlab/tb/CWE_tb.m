clear; close all; clc;

%% Параметры варианта Мартыненко И.А. (ЭРэ-19-22)
nRangeGate = 1000;
iWallPosition = [nRangeGate/2 nRangeGate];
iSignalPosition = nRangeGate/2 - 50;
INR_dB = [15 35];
SINR_dB = 20;
Pfa = 1e-5;
N = 16;
nGuardCell = 1;

%% Моделирование обнаружения на границе помех
data = simRangeData(INR_dB,iWallPosition,SINR_dB,iSignalPosition);
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
ylim([-40 max(INR_dB) + 20])
set(gca,'FontName','Times New Roman','FontSize',18);
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'CWE_tb.png'));
    