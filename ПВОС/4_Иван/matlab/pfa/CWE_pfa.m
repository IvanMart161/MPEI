clear all; close all; clc

%% Исходные параметры (Вариант: Мартыненко И.А. — СУЛТ-МСЗ)
INR_dB_set = [15, 35];          % Отношение помеха/шум по варианту, дБ
Pfa = 1e-5;                     % Заданная вероятность ложной тревоги
nTrial = 10^6;                  % Число испытаний Монте-Карло (увеличено для Pfa=10^-5)
N = 16;                         % Общее число обучающих ячеек скользящего окна
nGuardCell = 1;                 % Число защитных ячеек с каждой стороны

noisePwr = 1;                   % Дисперсия собственного шума
nWindow = N + 2*nGuardCell + 1;
iCUT = N/2 + nGuardCell + 1;
iLead = 1:N/2;
iTest = N/2 + nGuardCell + 1;
iLag = N/2 + 2*nGuardCell + 2:nWindow;
iClutterWall = [0 iLead iTest iLag];
nClutterWall = length(iClutterWall);

expPfa_all = zeros(nClutterWall, length(INR_dB_set));
h = waitbar(0, 'Моделирование ВЛТ на границе помех...');

for idx_inr = 1:length(INR_dB_set)
    INR = db2pow(INR_dB_set(idx_inr));
    clutterPwr = noisePwr * INR;
    pwrNoise = exprnd(noisePwr, nTrial, nWindow);
    pwrClutter = exprnd(clutterPwr, nTrial, nWindow);
    
    for i = 1:nClutterWall
        waitbar(((idx_inr - 1)*nClutterWall + i) / (length(INR_dB_set)*nClutterWall));
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
            pwrSignal(:, iNoise) = pwrNoise(:, iNoise);
            pwrSignal(:, iClutter) = pwrClutter(:, iClutter);
            if (clutterWall >= iCUT)
                pwrCUT = clutterPwr;
            else
                pwrCUT = noisePwr;
            end
        end
        
        threshold = SOCA_CFAR1(Pfa, N, nGuardCell, pwrSignal);
        expPfa_all(i, idx_inr) = monteCarloFalseAlarm(threshold, pwrCUT);
    end
end
close(h);

%% Построение и сохранение графика
cellLead = 1:N/2 + 1;
xLead = 0:N/2;
cellLag = N/2 + 2:N + 2;
xLag = N/2:N;
xMin = 0; xMax = N;
yMin = 10^-8; yMax = 10^0;

fig = figure(1);
clr1 = [0.494 0.184 0.556]; % фиолетовый для ОПШ = 15 дБ
clr2 = [0.850 0.325 0.098]; % оранжевый для ОПШ = 35 дБ

% Кривая для ОПШ = 15 дБ
h1 = semilogy(xLead, expPfa_all(cellLead, 1), '-', 'Color', clr1, 'LineWidth', 1.8);
hold on
semilogy(xLag, expPfa_all(cellLag, 1), '-', 'Color', clr1, 'LineWidth', 1.8);

% Кривая для ОПШ = 35 дБ
h2 = semilogy(xLead, expPfa_all(cellLead, 2), '--', 'Color', clr2, 'LineWidth', 1.8);
semilogy(xLag, expPfa_all(cellLag, 2), '--', 'Color', clr2, 'LineWidth', 1.8);

% Опорные линии заданной ВЛТ и середины окна
h3 = line([xMin xMax], [Pfa Pfa], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 1.5);
h4 = line([N/2 N/2], [yMin yMax], 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5);

ylabel('Вероятность ложной тревоги P_{ЛТ}');
xlabel('Расположение границы помех в скользящем окне ОБН с СУЛТ');
title('Зависимость P_{ЛТ} ОБН с СУЛТ-МСЗ на границе помех')
axis([xMin xMax yMin yMax]);
grid on
legend([h1, h2, h3, h4], ...
    'ОБН с СУЛТ-МСЗ (ОПШ = 15 дБ)', ...
    'ОБН с СУЛТ-МСЗ (ОПШ = 35 дБ)', ...
    'Заданный уровень P_{ЛТ} = 10^{-5}', ...
    'Граница помех (N/2 = 8)', ...
    'Location', 'southeast')
hold off

% Автосохранение в общую папку export_figs
out_dir = '../../export_figs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
set(fig, 'Position', [100, 100, 900, 450]);
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(out_dir, 'fig_cwe_pfa.png'), '-dpng', '-r300');