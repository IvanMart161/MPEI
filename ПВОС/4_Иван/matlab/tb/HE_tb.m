clear; close all; clc;

%% Параметры варианта Мартыненко И.А. (ЭРэ-19-22)
SINR_dB = 20;
Pfa = 1e-5;
nRangeGate = 1000;
iSignalPosition = nRangeGate/2;
N = 16;
nGuardCell = 1;

%% Моделирование однородного фона
data = simRangeData(0,nRangeGate,SINR_dB,iSignalPosition);
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
saveas(gcf,fullfile(fileparts(mfilename('fullpath')),'HE_tb.png'));