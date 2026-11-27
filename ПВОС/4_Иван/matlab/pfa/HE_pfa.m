clear all; close all; clc;

%% Исходные параметры (Вариант: Мартыненко И.А. — СУЛТ-МСЗ)
Pfa_set = 10.^-(1:0.25:7);      % Диапазон заданных вероятностей ложной тревоги
nPfa = length(Pfa_set);
nRangeGate = 10^6;              % Число испытаний Монте-Карло
N = 16;                         % Общее число обучающих ячеек скользящего окна
nGuardCell = 1;                 % Число защитных ячеек с каждой стороны
noisePwr = 1;                   % Дисперсия собственного шума

%% Теоретический расчет порогового множителя для СУЛТ-МСЗ (SOCA)
pwrSignal = exprnd(noisePwr, nRangeGate, 1);
alpha_vec = zeros(size(Pfa_set));

a_set = 0.1:0.1:250;
P = zeros(length(a_set), 1);
for i = 1:length(a_set)
    a_val = a_set(i);
    sum_term = zeros(N/2, 1);
    for k = 0:N/2 - 1
        sum_term(k + 1, 1) = nchoosek(N/2 - 1 + k, k) * (2 + a_val/(N/2))^(-k);
    end
    P(i, 1) = 2 * (2 + a_val/(N/2))^(-N/2) * sum(sum_term);
end

for i = 1:length(Pfa_set)
    alpha_vec(i) = interp1(P, a_set, Pfa_set(i));
end

%% Моделирование методом Монте-Карло
threshold = zeros(nRangeGate, nPfa);
h = waitbar(0, 'Вычисление порога...');

for i = 1:nPfa
    waitbar(i/nPfa);
    [threshold(:, i), iCFAR] = SOCA_CFAR1(Pfa_set(i), N, nGuardCell, pwrSignal);
end

close(h);
expPfa = monteCarloFalseAlarm(threshold(iCFAR, :), noisePwr);

%% Построение и сохранение графика
fig = figure(1);
semilogy(alpha_vec, Pfa_set, '-', 'Color', [0.494 0.184 0.556], 'LineWidth', 1.8);
hold on
semilogy(alpha_vec, expPfa, 'o', 'Color', [0.494 0.184 0.556], 'LineWidth', 1.5, 'MarkerSize', 6);
grid on
title('Зависимость P_{ЛТ}(\alpha) ОБН с СУЛТ-МСЗ в гомогенных условиях')
xlabel('Пороговый множитель \alpha')
ylabel('Вероятность ложной тревоги P_{ЛТ}')
legend('Теоретическая зависимость', 'Эксперимент (Монте-Карло)')
hold off

% Вывод в командное окно значения alpha для заданного Pfa = 1e-5
idx_target = find(abs(Pfa_set - 1e-5) < 1e-12);
fprintf('Для P_ЛТ = 10^-5 (N = %d): alpha_МСЗ = %.2f\n', N, alpha_vec(idx_target));

% Автосохранение в общую папку export_figs
out_dir = '../../export_figs';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
set(fig, 'Position', [100, 100, 900, 450]);
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(out_dir, 'fig_he_pfa.png'), '-dpng', '-r300');