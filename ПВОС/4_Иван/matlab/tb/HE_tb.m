clear all; close all; clc;

%% Исходные параметры (Вариант: Мартыненко И.А. — СУЛТ-МСЗ)
INR_dB = 0;               % Отношение помеха/шум (гомогенные условия — только собственный шум)
SINR_dB = 20;             % Отношение сигнал/(шум+помеха), дБ
Pfa = 1e-5;               % Заданная вероятность ложной тревоги
nRangeGate = 1000;        % Число элементов разрешения по дальности
iSignalPosition = nRangeGate/2; % Положение цели (500-й отсчет)
N = 16;                   % Общее число обучающих ячеек скользящего окна
nGuardCell = 1;           % Число защитных ячеек с каждой стороны

%% Моделирование
% Формирование данных на выходе квадратичного детектора (КД)
data = simRangeData(INR_dB, nRangeGate, SINR_dB, iSignalPosition);

% Создание обнаружителя с СУЛТ-МСЗ (метод 'SOCA')
detector = phased.CFARDetector('NumTrainingCells', N, ...
    'NumGuardCells', 2*nGuardCell, 'ProbabilityFalseAlarm', Pfa, ...
    'ThresholdOutputPort', true, 'Method', 'SOCA');

% Обнаружение и расчет порога
[detection, threshold] = detector(data, 1:nRangeGate);

%% Построение и сохранение графика
fig = figure(1);
plot(pow2db(data), 'LineWidth', 1.2)
hold on
plot(pow2db(threshold), 'LineWidth', 1.5)
stairs(detection*5 - 40, 'LineWidth', 1.5, 'Color', 'k')
grid on
axis tight
title('Обнаружение в гомогенных условиях в ОБН с СУЛТ')
xlabel('Отсчеты дальности')
ylabel('Амплитуда, дБ')
legend('сигнал после КД', 'порог ОБН с СУЛТ-МСЗ', 'решение об обнаружении')
ylim([-40 SINR_dB + 10])
hold off

% Автосохранение в общую папку export_figs в корне отчета
out_dir = '../../export_figs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
set(fig, 'Position', [100, 100, 900, 450]);
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(out_dir, 'fig_he_tb.png'), '-dpng', '-r300');