clear all; close all; clc;

%% Исходные параметры (Вариант: Мартыненко И.А. — СУЛТ-МСЗ)
iMasking = [1 2 3];             % Положения NM = 3 мешающих целей в отстающем полуокне
Pfa = 1e-5;                     % Заданная вероятность ложной тревоги
nTrial = 10^6;                  % Число испытаний Монте-Карло
SNRdB_set = 1:1:30;             % Диапазон ОСШ, дБ

N = 16;                         % Общее число обучающих ячеек скользящего окна
nGuardCell = 1;                 % Число защитных ячеек с каждой стороны

%% Генерация помехо-сигнальной обстановки
SNR_set = db2pow(SNRdB_set);
nSNR = length(SNRdB_set);
noisePwr = 1;
nWindow = N + 2*nGuardCell + 1;
iCUT = N/2 + nGuardCell + 1;
nMasking = length(iMasking);
iTarget = [iMasking iCUT];
nTarget = length(iTarget);
pwrSignal = exprnd(noisePwr, nTrial, nWindow);

%% Моделирование методом Монте-Карло
Pd = zeros(nSNR, 1);
h = waitbar(0, 'Моделирование ВПО и потерь на фоне мешающих целей...');

for iSNR = 1:nSNR
    waitbar(iSNR/nSNR);
    SNR = SNR_set(iSNR);
    targetPwr = noisePwr * SNR;
    pwrSignal(:, iTarget) = exprnd(targetPwr + noisePwr, nTrial, nTarget);
    
    threshold = SOCA_CFAR1(Pfa, N, nGuardCell, pwrSignal);
    
    nDetection = sum(pwrSignal(:, iCUT) > threshold);
    Pd(iSNR) = nDetection / nTrial;
end
close(h);

SNRdB_set = SNRdB_set(:);
out_dir = '../../export_figs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% График 1: Зависимость вероятности правильного обнаружения от ОСШ
fig1 = figure(1);
plot(SNRdB_set, Pd_inf(Pfa, SNRdB_set), 'k-', 'LineWidth', 1.8); % Критерий НП
hold on
plot(SNRdB_set, Pd, '-', 'Color', [0.494 0.184 0.556], 'LineWidth', 1.8);
legend('ОБН с НП', 'ОБН с СУЛТ-МСЗ', 'Location', 'northwest');
ylabel('Вероятность правильного обнаружения P_{ПО}');
xlabel('ОСШ, дБ');
title(['Зависимость P_{ПО}(ОСШ) на фоне ' num2str(nMasking) ' близкорасположенных целей']);
grid on;
box on;
hold off

set(fig1, 'Position', [100, 100, 900, 450]);
set(fig1, 'PaperPositionMode', 'auto');
print(fig1, fullfile(out_dir, 'fig_hmte_pd.png'), '-dpng', '-r300');

%% График 2: Зависимость потерь на обнаружение от ОСШ
SNRoptimal = OptSNR(Pd, Pfa);
CFARloss = SNRdB_set - SNRoptimal;

fig2 = figure(2);
plot(SNRdB_set, CFARloss, '-', 'Color', [0.494 0.184 0.556], 'LineWidth', 1.8);
legend('ОБН с СУЛТ-МСЗ', 'Location', 'northwest');
ylabel('Потери на обнаружение L, дБ');
xlabel('ОСШ, дБ');
title(['Зависимость L(ОСШ) в ОБН с СУЛТ на фоне ' num2str(nMasking) ' близкорасположенных целей']);
xMin = min(SNRdB_set);
xMax = max(SNRdB_set);
yMin = 0;
yMax = max(4, ceil(2*max(CFARloss))/2);
axis([xMin xMax yMin yMax]);
grid on;
box on;

set(fig2, 'Position', [150, 150, 900, 450]);
set(fig2, 'PaperPositionMode', 'auto');
print(fig2, fullfile(out_dir, 'fig_hmte_loss.png'), '-dpng', '-r300');

% Вывод в командное окно значений при заданном по варианту ОСШ = 20 дБ
idx_20 = find(SNRdB_set == 20);
fprintf('При ОСШ = 20 дБ (NM = %d): P_ПО = %.4f, потери L = %.2f дБ\n', ...
    nMasking, Pd(idx_20), CFARloss(idx_20));