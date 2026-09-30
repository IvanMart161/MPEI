clc
% Скрипт для расчета номиналов и коэффициентов отражения 
% для П, Т-образных аттенюаторов и делителя мощности с округлением до ряда E24.

Z0 = 50; % Волновое сопротивление, Ом
A = 6;   % Ослабление, дБ
K = 10^(A / 20); % Линейный коэффициент по напряжению

% Вспомогательные функции
par = @(R1, R2) (R1 .* R2) ./ (R1 + R2);
calc_G = @(Zin) (Zin - Z0) / (Zin + Z0);
calc_G_dB = @(G) 20 * log10(max(abs(G), 1e-10));

fprintf('=== РАСЧЕТ ХАРАКТЕРИСТИК (Ослабление: %d дБ, Z0 = %d Ом) ===\n\n', A, Z0);

%% 1. П-образный аттенюатор
R_par_pi = Z0 * (K + 1) / (K - 1);
R_ser_pi = Z0 * (K^2 - 1) / (2 * K);

R_par_pi_e24 = round_to_e24(R_par_pi);
R_ser_pi_e24 = round_to_e24(R_ser_pi);

Zin_pi_50_e24 = par(R_par_pi_e24, R_ser_pi_e24 + par(R_par_pi_e24, 50));
Zin_pi_25_e24 = par(R_par_pi_e24, R_ser_pi_e24 + par(R_par_pi_e24, 25));

fprintf('1. П-ОБРАЗНЫЙ АТТЕНЮАТОР\n');
fprintf('   Идеальные: R_паралл = %.2f Ом, R_послед = %.2f Ом\n', R_par_pi, R_ser_pi);
fprintf('   Ряд Е24:   R_паралл = %.2f Ом, R_послед = %.2f Ом\n', R_par_pi_e24, R_ser_pi_e24);
fprintf('   --- Коэффициенты отражения (по Е24) ---\n');
fprintf('   Нагрузка 50 Ом: Zin = %5.2f Ом | Г = %7.4f (%6.2f дБ)\n', Zin_pi_50_e24, calc_G(Zin_pi_50_e24), calc_G_dB(calc_G(Zin_pi_50_e24)));
fprintf('   Нагрузка 25 Ом: Zin = %5.2f Ом | Г = %7.4f (%6.2f дБ)\n\n', Zin_pi_25_e24, calc_G(Zin_pi_25_e24), calc_G_dB(calc_G(Zin_pi_25_e24)));

%% 2. Т-образный аттенюатор
R_ser_t = Z0 * (K - 1) / (K + 1);
R_par_t = Z0 * (2 * K) / (K^2 - 1);

R_ser_t_e24 = round_to_e24(R_ser_t);
R_par_t_e24 = round_to_e24(R_par_t);

Zin_t_50_e24 = R_ser_t_e24 + par(R_par_t_e24, R_ser_t_e24 + 50);
Zin_t_25_e24 = R_ser_t_e24 + par(R_par_t_e24, R_ser_t_e24 + 25);

fprintf('2. Т-ОБРАЗНЫЙ АТТЕНЮАТОР\n');
fprintf('   Идеальные: R_послед = %.2f Ом, R_паралл = %.2f Ом\n', R_ser_t, R_par_t);
fprintf('   Ряд Е24:   R_послед = %.2f Ом, R_паралл = %.2f Ом\n', R_ser_t_e24, R_par_t_e24);
fprintf('   --- Коэффициенты отражения (по Е24) ---\n');
fprintf('   Нагрузка 50 Ом: Zin = %5.2f Ом | Г = %7.4f (%6.2f дБ)\n', Zin_t_50_e24, calc_G(Zin_t_50_e24), calc_G_dB(calc_G(Zin_t_50_e24)));
fprintf('   Нагрузка 25 Ом: Zin = %5.2f Ом | Г = %7.4f (%6.2f дБ)\n\n', Zin_t_25_e24, calc_G(Zin_t_25_e24), calc_G_dB(calc_G(Zin_t_25_e24)));

%% 3. Резистивный делитель мощности 1:2
R_split = Z0 / 3;
R_split_e24 = round_to_e24(R_split);

Zin_split_50_e24 = R_split_e24 + par(R_split_e24 + 50, R_split_e24 + 50);
Zin_split_25_e24 = R_split_e24 + par(R_split_e24 + 25, R_split_e24 + 50);

fprintf('3. РЕЗИСТИВНЫЙ ДЕЛИТЕЛЬ МОЩНОСТИ 1:2\n');
fprintf('   Идеальные: R1 = R2 = R3 = %.2f Ом\n', R_split);
fprintf('   Ряд Е24:   R1 = R2 = R3 = %.2f Ом\n', R_split_e24);
fprintf('   --- Коэффициенты отражения (по Е24) ---\n');
fprintf('   Нагрузки 50 Ом и 50 Ом: Zin = %5.2f Ом | Г = %7.4f (%6.2f дБ)\n', Zin_split_50_e24, calc_G(Zin_split_50_e24), calc_G_dB(calc_G(Zin_split_50_e24)));
fprintf('   Нагрузки 25 Ом и 50 Ом: Zin = %5.2f Ом | Г = %7.4f (%6.2f дБ)\n', Zin_split_25_e24, calc_G(Zin_split_25_e24), calc_G_dB(calc_G(Zin_split_25_e24)));

%% Локальная функция округления до ряда E24
function R_out = round_to_e24(R_in)
    if R_in == 0
        R_out = 0;
        return;
    end
    
    e24 = [1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, ...
           3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1];
           
    exponent = floor(log10(R_in));
    mantissa = R_in / (10^exponent);
    
    % Находим ближайшее значение в массиве
    [~, idx] = min(abs(e24 - mantissa));
    R_out = e24(idx) * (10^exponent);
end