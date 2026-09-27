clear all; close all; clc;

%% 1. Настройка директорий экспорта
dest_png  = 'export_figs';
dest_tiff = fullfile(dest_png, 'pic_tiff');
dest_fig  = fullfile(dest_png, 'pic_fig');

if ~exist(dest_png, 'dir'),  mkdir(dest_png);  end
if ~exist(dest_tiff, 'dir'), mkdir(dest_tiff); end
if ~exist(dest_fig, 'dir'),  mkdir(dest_fig);  end

%% 2. Параметры Варианта 8
N_var = 8;
fs = (100 + N_var*10) * 1e6;         
fc = (130 + N_var*10) * 1e6 + 0.1e6; 
fg = (130 + N_var*10) * 1e6;         
noba = 5 + N_var;                    
nobg = 5 + N_var;                    

R = 10;
fsv = fs / R;                        
bandwidth = fsv / 2 / 1.2;           

length_test = 24e3;
t = (0:length_test - 1).' / fs;

snr = 70;
noise_power = 1;
signal_power = noise_power * db2pow(snr);
noise_vec = sqrt(noise_power) * randn(length_test, 1);

%% 3. Математическая инициализация фильтров (БЕЗ зависающих dsp-объектов)
disp('Настройка фильтров и гетеродина...');

% 3.1. Синтез импульсной характеристики CIC-фильтра (каскад из N_cic окон)
N_cic = 8;
hCIC = ones(1, R);
for k = 1:(N_cic - 1)
    hCIC = conv(hCIC, ones(1, R));
end
cic_gain = (R^N_cic) / 2;

% 3.2. Синтез FIR-компенсатора
Fpass = bandwidth / 2;
Fstop = Fpass + 500e3;
FIR_design = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass, ...
    'StopbandFrequency', Fstop, ...
    'PassbandRipple', 0.1, ...
    'StopbandAttenuation', 60, ...
    'SampleRate', fsv, ...
    'DesignMethod', 'equiripple');
hFIR = FIR_design.Coefficients;
hFIR = hFIR ./ max(hFIR);

% 3.3. Гетеродин
nco_real = cos(2*pi*fg*t);
nco_imag = sin(2*pi*fg*t);
nco_sig = single(nco_real - 1i * nco_imag);

%% =========================================================================
%% ПУНКТ 5: Сигнал с линейно возрастающей амплитудой (LA)
%% =========================================================================
disp('--- Выполняется Пункт 5 (Линейно возрастающая амплитуда) ---');
disp(' 1. Генерация сигнала...');
max_adc_val = 2^(noba - 1) - 1;
envelope_la = linspace(0, 1, length_test).';
signal_la = (max_adc_val * envelope_la) .* sin(2*pi*fc*t);
in_la = signal_la + noise_vec;

disp(' 2. Оцифровка и демодуляция...');
adc_la = round(in_la);
adc_la(adc_la > max_adc_val) = max_adc_val;
adc_la(adc_la < -max_adc_val) = -max_adc_val;
dem_la = double(adc_la) .* double(nco_sig);

disp(' 3. Мгновенная математическая фильтрация...');
% CIC фильтрация (простая свертка) и прореживание (1:R:end)
dem_filtered_5 = filter(hCIC, 1, dem_la);
cic_la = dem_filtered_5(1:R:end) / cic_gain;
% FIR фильтрация
fir_la = filter(hFIR, 1, cic_la);

disp(' 4. Отрисовка и экспорт (в фоне)...');
h5 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 800]);
subplot(2, 2, 1); plot(adc_la, 'LineWidth', 1.2); grid on; axis tight;
title('Signal LA on IF'); xlabel('Samples'); ylabel('LSB');
subplot(2, 2, 2); plot(real(fir_la), 'LineWidth', 1.2); hold on; plot(imag(fir_la), 'LineWidth', 1.2);
grid on; axis tight; title('Signal LA on IF = 0'); xlabel('Samples'); ylabel('LSB');

[f_la_in, s_la_in] = get_spectrum(adc_la, fs, 1);
subplot(2, 2, 3); plot(f_la_in/1e6, mag2db(abs(s_la_in)), 'LineWidth', 1.2); grid on; axis tight;
title('Spectrum LA on IF'); xlabel('f, MHz'); ylabel('|A|, dB');

[f_la_out, s_la_out] = get_spectrum(fir_la, fsv, 1);
subplot(2, 2, 4); plot(f_la_out/1e6, mag2db(abs(s_la_out)), 'LineWidth', 1.2); grid on; axis tight;
title('Spectrum LA on IF = 0'); xlabel('f, MHz'); ylabel('|A|, dB');

apply_white_theme(h5);
save_all_formats(h5, 'p5_linear_amplitude', dest_png, dest_tiff, dest_fig);
disp('Пункт 5 завершен успешно!');

%% =========================================================================
%% ПУНКТ 6: ЛЧМ-сигнал
%% =========================================================================
disp('--- Выполняется Пункт 6 (ЛЧМ-сигнал) ---');
disp(' 1. Генерация видео-ЛЧМ сигнала...');

B_target = 512;
sweep_bw = bandwidth;
pulse_width = B_target / sweep_bw;
num_lfm_samples = round(pulse_width * fs);

t_pulse = (0:num_lfm_samples-1).' / fs;
t_sym = t_pulse - pulse_width/2;
lfm_pulse = exp(1i * pi * (sweep_bw / pulse_width) * t_sym.^2); 

lfm_video = zeros(length_test, 1);
start_idx = round((length_test - num_lfm_samples) / 2) + 1;
lfm_video(start_idx : start_idx + num_lfm_samples - 1) = lfm_pulse;

disp(' 2. Перенос на ПЧ и оцифровка...');
lfm_if = max_adc_val * real(lfm_video .* exp(1i * 2*pi*fc*t));
in_lfm = lfm_if + noise_vec;

adc_lfm = round(in_lfm);
adc_lfm(adc_lfm > max_adc_val) = max_adc_val;
adc_lfm(adc_lfm < -max_adc_val) = -max_adc_val;

dem_lfm = double(adc_lfm) .* double(nco_sig);

disp(' 3. Мгновенная математическая фильтрация...');
dem_filtered_6 = filter(hCIC, 1, dem_lfm);
cic_lfm = dem_filtered_6(1:R:end) / cic_gain;
fir_lfm = filter(hFIR, 1, cic_lfm);

disp(' 4. Отрисовка входных графиков (П.6 часть 1)...');
h6_1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 800]);
subplot(2, 2, 1); plot(real(lfm_video), 'LineWidth', 1.2); hold on; plot(imag(lfm_video), 'LineWidth', 1.2);
grid on; axis tight; title('Signal LFM video'); xlabel('Samples'); ylabel('LSB');
subplot(2, 2, 2); plot(lfm_if, 'LineWidth', 1.2); grid on; axis tight;
title('Signal LFM on IF'); xlabel('Samples'); ylabel('LSB');

[f_lfm_v, s_lfm_v] = get_spectrum(lfm_video, fs, 0);
subplot(2, 2, 3); plot(f_lfm_v/1e6, mag2db(abs(s_lfm_v)), 'LineWidth', 1.2); grid on; axis tight;
title('Spectrum LFM video'); xlabel('f, MHz'); ylabel('|A|, dB');

[f_lfm_if, s_lfm_if] = get_spectrum(lfm_if, fs, 0);
subplot(2, 2, 4); plot(f_lfm_if/1e6, mag2db(abs(s_lfm_if)), 'LineWidth', 1.2); grid on; axis tight;
title('Spectrum LFM on IF'); xlabel('f, MHz'); ylabel('|A|, dB');

apply_white_theme(h6_1);
save_all_formats(h6_1, 'p6_lfm_input', dest_png, dest_tiff, dest_fig);

disp(' 5. Отрисовка выходных графиков (П.6 часть 2)...');
h6_2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 800]);
subplot(2, 2, 1); plot(adc_lfm, 'LineWidth', 1.2); grid on; axis tight;
title('Signal LFM on IF'); xlabel('Samples'); ylabel('LSB');
subplot(2, 2, 2); plot(real(fir_lfm), 'LineWidth', 1.2); hold on; plot(imag(fir_lfm), 'LineWidth', 1.2);
grid on; axis tight; title('Signal LFM on IF = 0'); xlabel('Samples'); ylabel('LSB');

[f_lfm_adc, s_lfm_adc] = get_spectrum(adc_lfm, fs, 1);
subplot(2, 2, 3); plot(f_lfm_adc/1e6, mag2db(abs(s_lfm_adc)), 'LineWidth', 1.2); grid on; axis tight;
title('Spectrum LFM on IF'); xlabel('f, MHz'); ylabel('|A|, dB');

[f_lfm_out, s_lfm_out] = get_spectrum(fir_lfm, fsv, 1);
subplot(2, 2, 4); plot(f_lfm_out/1e6, mag2db(abs(s_lfm_out)), 'LineWidth', 1.2); grid on; axis tight;
title('Spectrum LFM on IF = 0'); xlabel('f, MHz'); ylabel('|A|, dB');

apply_white_theme(h6_2);
save_all_formats(h6_2, 'p6_lfm_output', dest_png, dest_tiff, dest_fig);
disp('Пункт 6 завершен успешно!');

%% =========================================================================
%% ПУНКТ 7: Одиночный КИХ-дециматор
%% =========================================================================
disp('--- Выполняется Пункт 7 (Эквивалентный одиночный КИХ) ---');

single_fir_out = decimate(double(dem_la), R, 'fir');

h7 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 500]);
[f_s_fir, s_s_fir] = get_spectrum(single_fir_out, fsv, 1);

plot(f_s_fir/1e6, mag2db(abs(s_s_fir)), 'b', 'LineWidth', 2); hold on;
plot(f_la_out/1e6, mag2db(abs(s_la_out)), '--r', 'LineWidth', 1.5);
grid on; axis tight;
legend('Одиночный КИХ (Single FIR)', 'Каскад CIC + FIR', 'Location', 'South');
title('Сравнение спектров: Одиночный КИХ vs Каскад CIC+FIR');
xlabel('f, MHz'); ylabel('|A|, dB');

apply_white_theme(h7);
save_all_formats(h7, 'p7_single_fir_comparison', dest_png, dest_tiff, dest_fig);

disp('Пункт 7 завершен!');
disp('=== ВСЕ ПУНКТЫ УСПЕШНО ВЫПОЛНЕНЫ! ГРАФИКИ В ПАПКЕ export_figs ===');

%% =========================================================================
%% Вспомогательные функции
%% =========================================================================
function [f, s] = get_spectrum(in, fs, weighting)
    f = linspace(-fs/2, fs/2, length(in));
    if weighting
        s = fftshift(fft(double(in) .* hann(length(in))));
    else
        s = fftshift(fft(double(in)));
    end
end

function save_all_formats(h, filename, p_png, p_tiff, p_fig)
    print(h, '-dpng',  '-r300', fullfile(p_png,  [filename '.png']));
    print(h, '-dtiff', '-r300', fullfile(p_tiff, [filename '.tiff']));
    savefig(h, fullfile(p_fig, [filename '.fig']));
end

function apply_white_theme(h_fig)
    axs = findobj(h_fig, 'Type', 'axes');
    for idx = 1:length(axs)
        ax = axs(idx);
        set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
            'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.5, ...
            'FontSize', 12, 'FontName', 'Times New Roman');
        t = get(ax, 'Title'); set(t, 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 14);
        x = get(ax, 'XLabel'); set(x, 'Color', 'k', 'FontSize', 12);
        y = get(ax, 'YLabel'); set(y, 'Color', 'k', 'FontSize', 12);
    end
    lg = findobj(h_fig, 'Type', 'Legend');
    if ~isempty(lg)
        set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    end
end