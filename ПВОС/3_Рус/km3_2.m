% Главный скрипт: Пространственно-временная обработка сигналов
clear all; close all; clc;
rng(0); % Для воспроизводимости шумовой дорожки

%% ===================== ПАРАМЕТРЫ ВАРИАНТА 4 =====================
N_var = 4;

fc      = (130 + N_var*10)*1e6;   % ПЧ входного сигнала: 170 МГц
n_a     = 5 + N_var;              % Разрядность АЦП: 9 бит
fs      = (100 + N_var*10)*1e6;   % Частота дискретизации АЦП: 140 МГц
n_g     = 5 + N_var;              % Разрядность ЦГ: 9 бит
f_G     = (130 + N_var*10)*1e6;   % Центральная частота ЦГ: 170 МГц

B_base  = 512;                    % База сигнала
Delta_f = 20e6;                   % Полоса ЛЧМ (20 МГц)
pulse_width = B_base/Delta_f;     % Длительность импульса: 25.6 мкс

length_test   = 2^14;             % Длина тестового кадра: 16384 отсчета
PRI           = length_test/fs;   % Период повторения
prf           = 1/PRI;            % Частота повторения
delay_samples = length_test/2;    % Задержка импульса: 8192 отсчета

noise_power  = 0.5;               % Мощность шума

fprintf('Вариант %d:\n', N_var);
fprintf('fc=%.1f MHz, fs=%.1f MHz, f_G=%.1f MHz\n', fc/1e6, fs/1e6, f_G/1e6);
fprintf('Разрядность: АЦП=%d бит, ЦГ=%d бит\n', n_a, n_g);

%% ===================== ГЕНЕРАЦИЯ ЛЧМ (через phased) =====================
% Согласно заданию, используем встроенную функцию phased.LinearFMWaveform
hwav = phased.LinearFMWaveform('SampleRate', fs, 'PulseWidth', pulse_width, ...
    'PRF', prf, 'SweepBandwidth', Delta_f, 'SweepInterval', 'Symmetric', ...
    'Envelope', 'Rectangular');

x_lfm_base = step(hwav);
N_imp = round(pulse_width * fs); % Количество отсчетов в самом импульсе (3584)
x_lfm = x_lfm_base(1:N_imp);     % Отрезаем нули до конца периода

%% ===================== ФОРМИРОВАНИЕ ВХОДНОГО СИГНАЛА ====================
n = (0:length_test-1).';
t = n/fs;

% Формируем комплексную огибающую с задержкой
signal_bb = zeros(length_test, 1);
idx = delay_samples + (1:length(x_lfm));
signal_bb(idx) = x_lfm;

% Перенос на промежуточную частоту (ПЧ) и добавление шума
signal_rf = real(signal_bb .* exp(1i*2*pi*fc*t));
noise_rf  = sqrt(noise_power) * randn(length_test, 1);
in_rf     = signal_rf + noise_rf;

%% ===================== АЦП (Субдискретизация) =====================
in_adc = quantize_real(in_rf, n_a);

%% ===================== ЦДМ + ЦФ =====================
% Демодуляция (перенос алиасной частоты на 0)
[x_bb, ~] = digital_demod(in_adc, fs, f_G, n_g);
% Фильтрация побочных продуктов (ФНЧ)
x_bb = digital_lpf(x_bb, fs, Delta_f, 128); 

% Примечание: Децимация (D>1) здесь не применяется для сохранения 
% максимального временного разрешения при сравнении графиков СФ.

%% ===================== РЕПЛИКА СФ =====================
replica = flipud(conj(x_lfm)); % Зеркальная комплексно-сопряженная копия

%% ===================== СФ: ПРЯМАЯ СВЕРТКА (FIR) =====================
mf_fir_nw = matched_filter_fir(x_bb, replica, false);
mf_fir_w  = matched_filter_fir(x_bb, replica, true);

%% ===================== СФ: СПЕКТРАЛЬНАЯ (FFT/IFFT) ==================
mf_fft_nw = matched_filter_fft(x_bb, replica, [], false);
mf_fft_w  = matched_filter_fft(x_bb, replica, [], true);

%% ===================== СФ: СЕГМЕНТИРОВАННАЯ (fftfilt) ===============
mf_fftfilt_nw = matched_filter_fftfilt(x_bb, replica, false);
mf_fftfilt_w  = matched_filter_fftfilt(x_bb, replica, true);


%% ===================== ВИЗУАЛИЗАЦИЯ (с ЗУМОМ) =======================
% Вспомогательная функция для нормировки графиков к 0 дБ
norm_db = @(x) mag2db(abs(x) / max(abs(x)));

% Определение границ для масштабирования графика вокруг пика
% Пик находится на отсчете: Задержка + Длина импульса = 8192 + 3584 = 11776
peak_idx = delay_samples + N_imp;
zoom_lims = [peak_idx - 150, peak_idx + 150];

% 1. Огибающая видеосигнала после ЦДМ и ЦФ
figure('Name','Комплексная огибающая', 'Position', [100, 100, 800, 400]);
plot(abs(x_bb), 'LineWidth', 1.5);
grid on; axis tight;
title('Огибающая комплексного видеосигнала на входе СФ');
xlabel('Отсчет'); ylabel('|x|, Амплитуда');

% 2. СФ: fftfilt (Требуемый по заданию метод, крупный план)
figure('Name','СФ fftfilt (Зум)', 'Position', [150, 150, 800, 400]);
plot(norm_db(mf_fftfilt_nw), 'LineWidth', 1.5); hold on;
plot(norm_db(mf_fftfilt_w),  'LineWidth', 1.5);
grid on;
xlim(zoom_lims);
ylim([-80 5]);
title('Отклик СФ (fftfilt) - Сравнение оконной обработки');
xlabel('Отсчет'); ylabel('|y|, дБ (нормиров.)');
legend('Без окна (Прямоугольное)', 'С окном Хэмминга');

% 3. Сравнение всех трех методов СФ на одном графике
figure('Name','Сравнение методов СФ', 'Position', [200, 200, 800, 400]);
plot(norm_db(mf_fir_nw),     'b-',  'LineWidth', 4); hold on;
plot(norm_db(mf_fft_nw),     'r--', 'LineWidth', 2);
plot(norm_db(mf_fftfilt_nw), 'g-.', 'LineWidth', 1);
grid on;
xlim(zoom_lims);
ylim([-80 5]);
title('Тождественность алгоритмов согласованной фильтрации');
xlabel('Отсчет'); ylabel('|y|, дБ (нормиров.)');
legend('FIR (прямая свертка)', 'FFT/IFFT (спектральная)', 'fftfilt (Overlap-Add)');


%% ===================== ЛОКАЛЬНЫЕ ВЛОЖЕННЫЕ ФУНКЦИИ =====================

function y = quantize_real(x, bits)
    full_scale = max(abs(x(:)))*1.01;
    if full_scale == 0, full_scale = 1; end
    LSB = 2*full_scale/(2^bits - 1);
    y = round(x/LSB)*LSB;
    y = max(min(y, full_scale), -full_scale);
end

function y = quantize_complex(x, bits)
    full_scale = max(abs(x(:)))*1.01;
    y = quantize_real(real(x), bits) + 1i*quantize_real(imag(x), bits);
end

function [y, lo] = digital_demod(x, fs, f_G, n_g)
    n = (0:numel(x)-1).';
    lo = exp(-1i*2*pi*f_G*n/fs);
    lo = quantize_complex(lo, n_g);
    y = x .* lo;
end

function y = digital_lpf(x, fs, f_cut, order)
    if nargin < 4, order = 128; end
    h = fir1(order, f_cut/(fs/2));
    y = filter(h, 1, x);
end

function y = matched_filter_fir(in, replica, use_window)
    if nargin < 3, use_window = false; end
    if use_window
        replica = replica .* hamming(length(replica));
    end
    y = filter(replica, 1, in);
end

function y = matched_filter_fft(in, replica, Nfft, use_window)
    if nargin < 3 || isempty(Nfft)
        Nfft = 2^nextpow2(length(in) + length(replica) - 1);
    end
    if nargin < 4, use_window = false; end
    if use_window
        replica = replica .* hamming(length(replica));
    end
    y_full = ifft(fft(in, Nfft) .* fft(replica, Nfft), Nfft);
    y = y_full(1:length(in));
end

function y = matched_filter_fftfilt(in, replica, use_window)
    if nargin < 3, use_window = false; end
    if use_window
        replica = replica .* hamming(length(replica));
    end
    y = fftfilt(replica, in);
end