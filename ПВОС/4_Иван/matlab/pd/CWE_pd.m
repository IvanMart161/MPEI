clear all; close all; clc

%% Исходные параметры (Вариант: Мартыненко И.А. — СУЛТ-МСЗ)
SNR_dB = 20;                    % Отношение сигнал/шум, дБ
INR_dB_set = [15, 35];          % Отношение помеха/шум по варианту, дБ
Pfa = 1e-5;                     % Заданная вероятность ложной тревоги
nTrial = 2*10^5;                % Число испытаний Монте-Карло
N = 16;                         % Общее число обучающих ячеек скользящего окна
nGuardCell = 1;                 % Число защитных ячеек с каждой стороны

SNR = db2pow(SNR_dB);
noisePwr = 1;
targetPwr = noisePwr * SNR;
nWindow = N + 2*nGuardCell + 1;
iCUT = N/2 + nGuardCell + 1;
iLead = 1:N/2;
iTest = N/2 + nGuardCell + 1;
iLag = N/2 + 2*nGuardCell + 2:nWindow;
iClutterWall = [0 iLead iTest iLag];
nClutterWall = length(iClutterWall);

expPd_all = zeros(nClutterWall, length(INR_dB_set));
h = waitbar(0, 'Моделирование ВПО на границе помех...');

for idx_inr = 1:length(INR_dB_set)
    INR = db2pow(INR_dB_set(idx_inr));
    clutterPwr = noisePwr * INR;
    pwrNoise = exprnd(noisePwr, nTrial, nWindow);
    pwrClutter = exprnd(clutterPwr, nTrial, nWindow);
    
    for i = 1:nClutterWall
        clutterWall = iClutterWall(i);
        waitbar(((idx_inr - 1)*nClutterWall + i) / (length(INR_dB_set)*nClutterWall));
        if (clutterWall == 0)
            pwrSignal = pwrNoise;
        elseif (clutterWall == nWindow)
            pwrSignal = pwrClutter;
        else
            iClutter = 1:clutterWall;
            iNoise = clutterWall + 1:nWindow;
            pwrSignal(:, iNoise) = pwrNoise(:, iNoise);
            pwrSignal(:, iClutter) = pwrClutter(:, iClutter);
        end
        
        threshold = SOCA_CFAR1(Pfa, N, pwrSignal);
        
        if (i <= (nClutterWall/2))
            pwrTarget = exprnd(noisePwr + targetPwr, nTrial, 1);
        else
            pwrTarget = exprnd(clutterPwr + targetPwr, nTrial, 1);
        end
        
        iDetection = pwrTarget > threshold;
        expPd_all(i, idx_inr) = sum(iDetection) / nTrial;
    end
end
close(h);

%% Построение и сохранение графика
cellLead = 1:N/2 + 1;
xLead = 0:N/2;
cellLag = N/2 + 2:N + 2;
xLag = N/2:N;
xMin = 0; xMax = N;
yMin = 0; yMax = 1;

fig = figure(1);
clr1 = [0.494 0.184 0.556]; % фиолетовый для ОПШ = 15 дБ
clr2 = [0.850 0.325 0.098]; % оранжевый для ОПШ = 35 дБ

hold on
h1 = plot(xLead, expPd_all(cellLead, 1), '-', 'Color', clr1, 'LineWidth', 1.8);
plot(xLag, expPd_all(cellLag, 1), '-', 'Color', clr1, 'LineWidth', 1.8);

h2 = plot(xLead, expPd_all(cellLead, 2), '--', 'Color', clr2, 'LineWidth', 1.8);
plot(xLag, expPd_all(cellLag, 2), '--', 'Color', clr2, 'LineWidth', 1.8);

h3 = line([N/2 N/2], [yMin yMax], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 1.5);

legend([h1, h2, h3], ...
    'ОБН с СУЛТ-МСЗ (ОПШ = 15 дБ)', ...
    'ОБН с СУЛТ-МСЗ (ОПШ = 35 дБ)', ...
    'Граница помех (N/2 = 8)', ...
    'Location', 'northeast')
title('Зависимость P_{ПО} на границе помех ОБН с СУЛТ-МСЗ')
ylabel('Вероятность правильного обнаружения P_{ПО}');
xlabel('Расположение границы помех в скользящем окне ОБН с СУЛТ');
axis([xMin xMax yMin yMax]);
grid on
box on
hold off

% Автосохранение в общую папку export_figs
out_dir = '../../export_figs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
set(fig, 'Position', [100, 100, 900, 450]);
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(out_dir, 'fig_cwe_pd.png'), '-dpng', '-r300');