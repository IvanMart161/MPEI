%% ЛАБОРАТОРНАЯ РАБОТА № 3: ДОМАШНЯЯ ПОДГОТОВКА
% Кафедра РЭКС НИУ "МЭИ"
% Тема: Измерение основных параметров ВЧ-усилителей
% Скрипт численного моделирования и решения задач домашней подготовки

clear; clc; close all;

fprintf('=\n');
fprintf('  РАСЧЁТ ЗАДАЧ ДОМАШНЕЙ ПОДГОТОВКИ К ЛР №3 (РЭКС)     \n');
fprintf('=\n\n');

%% Задача 2: Мощность шума на выходе при входном шуме 0 дБм
P_in_noise_dBm = 0;           % Входной шум, дБм
Kp_task2_dB = 10;             % Коэффициент усиления, дБ
Ksh_task2_dB = 3;             % Коэффициент шума, дБ
Pi_task2_Hz = 100e6;          % Полоса, 100 МГц

% Тепловой шум согласованной нагрузки
kT0_dBm_Hz = -174;            % дБм/Гц
P_therm_dBm = kT0_dBm_Hz + 10*log10(Pi_task2_Hz); % -94 дБм

% Мощность на выходе:
% Так как входной шум (0 дБм) >> теплового шума (-94 дБм),
% вклад собственных шумов усилителя пренебрежимо мал.
P_out_noise_dBm = P_in_noise_dBm + Kp_task2_dB;

fprintf('--- Задача 2 ---\n');
fprintf('Входной шум: %.1f дБм в полосе 100 МГц\n', P_in_noise_dBm);
fprintf('Тепловой порог: %.1f дБм\n', P_therm_dBm);
fprintf('Выходная мощность шума: %.2f дБм (%.2f мВт)\n\n', ...
P_out_noise_dBm, 10^(P_out_noise_dBm/10));

%% Задача 3: Выходное отношение сигнал/шум
Pin_sig_dBm = -80;            % Мощность входного сигнала, дБм
Pi_task3_Hz = 10e6;           % Полоса, 10 МГц
NF_task3_dB = 6;              % Коэффициент шума, дБ

P_noise_in_dBm = kT0_dBm_Hz + 10*log10(Pi_task3_Hz); % -104 дБм
SNR_in_dB = Pin_sig_dBm - P_noise_in_dBm;            % 24 дБ
SNR_out_dB = SNR_in_dB - NF_task3_dB;                 % 18 дБ

fprintf('--- Задача 3 ---\n');
fprintf('Мощность входного теплового шума: %.1f дБм\n', P_noise_in_dBm);
fprintf('Входное отношение сигнал/шум (SNR_вх): %.2f дБ\n', SNR_in_dB);
fprintf('Выходное отношение сигнал/шум (SNR_вых): %.2f дБ\n\n', SNR_out_dB);

%% Задача 11: Входной аттенюатор 6 дБ и параметры IP3
Att_dB = 6;
fprintf('--- Задача 11 ---\n');
fprintf('Включение аттенюатора %.1f дБ на вход усилителя:\n', Att_dB);
fprintf('  IIP3 увеличивается на: +%.1f дБ (IIP3_new = IIP3_0 + %.1f дБ)\n', Att_dB, Att_dB);
fprintf('  OIP3 остаётся неизменным: 0.0 дБ (OIP3_new = OIP3_0)\n\n');

%% Задача 12: Полоса интермодуляции 3-го порядка
Bandwidth_sig_MHz = 10;
Bandwidth_IM3_MHz = 3 * Bandwidth_sig_MHz;

fprintf('--- Задача 12 ---\n');
fprintf('Полоса входного сигнала: %.1f МГц\n', Bandwidth_sig_MHz);
fprintf('Ширина спектра продуктов IM3: %.1f МГц\n\n', Bandwidth_IM3_MHz);

%% Задача 13: Сравнение каскадов аттенюатор-усилитель
Kp_amp_dB = 20;               Kp_amp = 10^(Kp_amp_dB/10);       % 100
NF_amp_dB = 3;                Ksh_amp = 10^(NF_amp_dB/10);      % 1.9953
Att_task13_dB = 6;

Kp_att = 10^(-Att_task13_dB/10);                                 % 0.25
Ksh_att = 10^(Att_task13_dB/10);                                 % 3.9811

% Схема А: Аттенюатор на входе
Ksh_total_A = Ksh_att + (Ksh_amp - 1) / Kp_att;
NF_total_A_dB = 10*log10(Ksh_total_A);

% Схема Б: Аттенюатор на выходе
Ksh_total_B = Ksh_amp + (Ksh_att - 1) / Kp_amp;
NF_total_B_dB = 10*log10(Ksh_total_B);

fprintf('--- Задача 13 ---\n');
fprintf('Схема А (Атт на входе):  К_ш = %.3f (NF = %.2f дБ)\n', Ksh_total_A, NF_total_A_dB);
fprintf('Схема Б (Атт на выходе): К_ш = %.3f (NF = %.2f дБ)\n', Ksh_total_B, NF_total_B_dB);
fprintf('Выигрыш по коэффициенту шума при переносе аттенюатора: %.2f дБ\n\n', ...
NF_total_A_dB - NF_total_B_dB);

%% МОДЕЛИРОВАНИЕ: АМПЛИТУДНАЯ ХАРАКТЕРИСТИКА, P1dB И IP3
% Зададим модель усилителя: y(t) = K1x + K2x^2 + K3*x^3
% Линейное усиление Kp = +15 дБ, входное согласованное сопротивление R = 50 Ом
Kp_sim_dB = 15;
K1 = 10^(Kp_sim_dB/20); % Коэффициент передачи по напряжению
K3 = -0.08 * K1;        % Нелинейный коэффициент сжатия (компрессии)

% Диапазон входных мощностей
Pin_dBm_vec = -30:0.5:15;
Pin_W_vec = 1e-3 * 10.^(Pin_dBm_vec / 10);
R_load = 50;
Vin_amp_vec = sqrt(2 * R_load * Pin_W_vec);

% Выходная амплитуда основного тона при однотоновом сигнале:
% y(t) = K1Acos(wt) + K3*(Acos(wt))^3 -> первая гармоника: (K1A + 3/4K3A^3)
Vout_amp_vec = abs(K1 * Vin_amp_vec + 0.75 * K3 * (Vin_amp_vec.^3));
Pout_W_vec = (Vout_amp_vec.^2) / (2 * R_load);
Pout_dBm_vec = 10 * log10(Pout_W_vec / 1e-3);

% Линейная экстраполяция
Pout_lin_dBm_vec = Pin_dBm_vec + Kp_sim_dB;

% Нахождение точки 1-дБ компрессии (P1dB)
delta_comp = Pout_lin_dBm_vec - Pout_dBm_vec;
idx_1dB = find(delta_comp >= 1.0, 1, 'first');
IP1dB_val = Pin_dBm_vec(idx_1dB);
OP1dB_val = Pout_dBm_vec(idx_1dB);

% Расчёт IP3
% При двухтоновом воздействии C: основная гармоника V1 = K1C, интермодуляция V_IM3 = 3/4K3C^3
% OIP3 = K1^2 * C_IP3^2 / (2R), где K1C = 3/4|K3|C^3 => C_IP3 = sqrt(4K1 / (3*|K3|))
C_IP3 = sqrt(4 * K1 / (3 * abs(K3)));
IIP3_mW = (C_IP3^2 / (2 * R_load)) / 1e-3;
IIP3_val = 10*log10(IIP3_mW);
OIP3_val = IIP3_val + Kp_sim_dB;

%% ПОСТРОЕНИЕ ГРАФИКОВ
figure('Name', 'Лабораторная работа №3: Нелинейные характеристики', ...
'Color', 'w', 'Position', [100, 100, 950, 700]);

% График 1: Точка компрессии P1dB и IIP3
subplot(2, 1, 1);
plot(Pin_dBm_vec, Pout_dBm_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'Реальная амплитудная хар-ка');
hold on; grid on;
plot(Pin_dBm_vec, Pout_lin_dBm_vec, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Линейная асимптота');

% Отметка P1dB
plot(IP1dB_val, OP1dB_val, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
'DisplayName', sprintf('P_{1dB} (IP=%.1f дБм, OP=%.1f дБм)', IP1dB_val, OP1dB_val));
xline(IP1dB_val, 'r:');
yline(OP1dB_val, 'r:');

xlabel('Входная мощность P_{вх}, дБм', 'FontSize', 11);
ylabel('Выходная мощность P_{вых}, дБм', 'FontSize', 11);
title('Определение точки децибельной компрессии P_{1dB}', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
ylim([min(Pout_dBm_vec) max(Pout_dBm_vec)+5]);

% График 2: Спектр двухтонового сигнала на выходе нелинейного усилителя
subplot(2, 1, 2);
f1 = 1000; f2 = 1001; % МГц
f_IM3_1 = 2*f1 - f2;  % 999 МГц
f_IM3_2 = 2*f2 - f1;  % 1002 МГц

% Уровни сигналов для типовой точки ниже P1dB на 10 дБ
P_fund_out = OP1dB_val - 10; % Мощность основного тона
IMD3_val = 2 * (OIP3_val - P_fund_out);
P_IM3_out = P_fund_out - IMD3_val;

freqs = [f_IM3_1, f1, f2, f_IM3_2];
powers = [P_IM3_out, P_fund_out, P_fund_out, P_IM3_out];

stem(freqs, powers, 'filled', 'LineWidth', 2, 'MarkerSize', 7, 'Color', [0.85 0.33 0.1]);
hold on; grid on;
xlabel('Частота, МГц', 'FontSize', 11);
ylabel('Мощность, дБм', 'FontSize', 11);
title(sprintf('Выходной двухтоновый спектр (\Delta f = 1 МГц, IMD_3 = %.1f дБ)', IMD3_val), ...
'FontSize', 12, 'FontWeight', 'bold');
xlim([997 1004]);
ylim([P_IM3_out-15 P_fund_out+10]);

% Текстовые аннотации спектра
text(f1, P_fund_out + 3, sprintf('f_1 (%.1f дБм)', P_fund_out), 'HorizontalAlignment', 'center');
text(f2, P_fund_out + 3, sprintf('f_2 (%.1f дБм)', P_fund_out), 'HorizontalAlignment', 'center');
text(f_IM3_1, P_IM3_out + 3, sprintf('2f_1-f_2\n(%.1f дБм)', P_IM3_out), 'HorizontalAlignment', 'center');
text(f_IM3_2, P_IM3_out + 3, sprintf('2f_2-f_1\n(%.1f дБм)', P_IM3_out), 'HorizontalAlignment', 'center');

fprintf('=\n');
fprintf('  Моделирование выполнено успешно!\n');
fprintf('  Точка компрессии:  IP1dB = %.2f дБм, OP1dB = %.2f дБм\n', IP1dB_val, OP1dB_val);
fprintf('  Точка пересечения: IIP3  = %.2f дБм, OIP3  = %.2f дБм\n', IIP3_val, OIP3_val);
fprintf('=\n');