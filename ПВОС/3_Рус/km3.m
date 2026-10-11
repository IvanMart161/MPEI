% Главный скрипт
clear; close all; clc;

%% 1. Исходные параметры (Вариант 4)
N = 4;

fc = (130 + N * 10) * 1e6;       % ПЧ входного сигнала: 170 МГц
fs = (100 + N * 10) * 1e6;       % Частота дискретизации АЦП: 140 МГц
fG = (130 + N * 10) * 1e6;       % Центральная частота ЦГ: 170 МГц
n_a = 5 + N;                     % Разрядность АЦП: 9 бит
n_g = 5 + N;                     % Разрядность ЦГ: 9 бит

pulse_base = 512;                % База сигнала
pulse_bw = 20e6;                 % Полоса сигнала (проверить по Таблице 1)
pulse_width = pulse_base / pulse_bw; % Длительность импульса
prf = 1 / pulse_width;           % Частота повторения

fsv = 40e6;                      % Частота дискретизации видеосигнала
decimation_factor = floor(fs/fsv); % Коэффициент децимации (3)

%% 2. Генерация тестового воздействия (ЛЧМ)
% Генерируем комплексную огибающую (как в pvos_practice_3)
waveform_baseband = generate_lfm_baseband(fsv, pulse_width, prf, pulse_bw);

% Перенос на ПЧ и имитация АЦП (здесь можно задействовать твои вспомогательные функции для АЦП)
t_rf = (0:length(waveform_baseband)*decimation_factor-1)' / fs;
signal_rf = real(rectpuls(t_rf - pulse_width/2, pulse_width) .* ...
            exp(1i * 2 * pi * fc * t_rf + 1i * pi * (pulse_bw/pulse_width) * (t_rf - pulse_width/2).^2));
signal_rf = signal_rf + 0.1 * randn(size(signal_rf)); % Добавление шума
signal_rf = quantize_signal(signal_rf, n_a); % Квантование АЦП

%% 3. ЦДМ: Цифровой гетеродин
nco_signal = exp(-1i * 2 * pi * fG * t_rf);
nco_signal = quantize_signal(nco_signal, n_g); % Квантование ЦГ (все вычисления в double) 
demodulated_signal = signal_rf .* nco_signal;

%% 4. ЦФ: Цифровая фильтрация и децимация
% Проектирование ФНЧ для подавления гармоник ЦДМ
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', pulse_bw/2, ...
         'StopbandFrequency', pulse_bw, 'SampleRate', fs);
filtered_signal = filter(lpFilt, demodulated_signal);

% Децимация
baseband_rx = filtered_signal(1:decimation_factor:end);

%% 5. Согласованная фильтрация (fftfilt)
[mf_out, mf_out_w] = matched_filter_fftfilt(baseband_rx, waveform_baseband);

%% 6. Визуализация 
% График 1: Сигнал после АЦП (Фрагмент, чтобы было видно заполнение)
figure;
plot(t_rf(1:2000)*1e6, signal_rf(1:2000), 'LineWidth', 1.5);
grid on;
title('Сигнал на выходе АЦП (фрагмент)');
xlabel('Время, мкс'); ylabel('Амплитуда');

% График 2: Сигнал после цифрового фильтра и децимации (на видеочастоте)
figure;
plot(real(baseband_rx), 'LineWidth', 1.5); hold on;
plot(imag(baseband_rx), 'LineWidth', 1.5);
grid on; axis tight;
title('Комплексная огибающая после ЦДМ, ЦФ и децимации');
xlabel('Отсчеты'); ylabel('Амплитуда');
legend('Синфазная (I)', 'Квадратурная (Q)');

% График 3: Выход согласованного фильтра (как в image_9d435c.png)
figure;
plot(mag2db(abs(mf_out)), 'LineWidth', 2); hold on;
plot(mag2db(abs(mf_out_w)), 'LineWidth', 2);
grid on; axis tight;
title('Выход согласованного фильтра');
legend('Без окна', 'С окном Хэмминга');
xlabel('Отсчеты'); ylabel('Амплитуда, дБ');
%% --- Вложенные функции модулей ---
function wave = generate_lfm_baseband(fsv, pulse_width, prf, pulse_bw)
    % Использование встроенного объекта phased.LinearFMWaveform 
    hwav = phased.LinearFMWaveform('SampleRate', fsv, 'PulseWidth', pulse_width, ...
        'PRF', prf, 'SweepBandwidth', pulse_bw, 'SweepInterval', 'Symmetric', ...
        'Envelope', 'Rectangular'); 
    wave = step(hwav); 
end

function sig_q = quantize_signal(sig, bits)
    % Имитация ограничения разрядности (без изменения типа данных double) 
    levels = 2^bits - 1;
    sig_q = round(sig * (levels/2)) / (levels/2);
end

function [mf_out, mf_out_windowed] = matched_filter_fftfilt(in_sig, ref_wave)
    % Создание импульсной характеристики (перевернутая комплексно-сопряженная копия) 
    replica = flip(conj(ref_wave)); 
    
    % Метод сегментированной свертки fftfilt без окна 
    mf_out = fftfilt(replica, in_sig); 
    
    % Метод с окном Хэмминга 
    win = hamming(length(replica)); 
    replica_w = replica .* win; 
    mf_out_windowed = fftfilt(replica_w, in_sig); 
end