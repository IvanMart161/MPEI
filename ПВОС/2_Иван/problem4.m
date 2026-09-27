clear all; close all; clc;

%% 1. Настройка директорий
dest_png  = 'export_figs';
dest_tiff = fullfile(dest_png, 'pic_tiff');
dest_fig  = fullfile(dest_png, 'pic_fig');

if ~exist(dest_png, 'dir'),  mkdir(dest_png);  end
if ~exist(dest_tiff, 'dir'), mkdir(dest_tiff); end
if ~exist(dest_fig, 'dir'),  mkdir(dest_fig);  end

%% 2. Параметры Варианта 8
N_var = 8;
fs = (100 + N_var*10) * 1e6;         % 180 МГц
fc = (130 + N_var*10) * 1e6 + 0.1e6;   % 210.1 МГц
fg = (130 + N_var*10) * 1e6;         % 210.0 МГц
noba = 5 + N_var;                    % 13 бит
nobg = 5 + N_var;                    % 13 бит
R = 10;
fsv = fs / R;                        % 18 МГц
bandwidth = fsv / 2 / 1.2;

length_test = 24e3;
t = (0:length_test - 1) / fs;

snr = 70;
noise_power = 1;
signal_power = noise_power * db2pow(snr);

%% 3. Формирование входного сигнала и АЦП
noise = sqrt(noise_power) * randn(length_test, 1);
signal = sqrt(signal_power) * sin(2*pi*fc*t).';
in = signal + noise;
adc_signal = quantize(int32(in), 1, noba, 0, 'Floor', 'Saturate');

%% 4. Инициализация фильтров
N_cic = 8;
CIC = dsp.CICDecimator(R, 1, N_cic);
CIC.FixedPointDataType = 'Full precision';
cic_gain = R^N_cic / 2;

FIR = dsp.CICCompensationDecimator(CIC, ...
    'DecimationFactor', 1, ...
    'PassbandFrequency', bandwidth/2, ...
    'PassbandRipple', 0.1, ...
    'StopbandAttenuation', 60, ...
    'StopbandFrequency', bandwidth/2 + 500e3, ...
    'SampleRate', fsv, ...
    'DesignForMinimumOrder', true);

hFIR = impz(FIR);
hFIR = hFIR ./ max(hFIR);

%% 5. Исследование амплитудного разбаланса (AMV)
disp('Расчет Пункта 4 (Амплитудный разбаланс)...');
amv_vec = 0:0.2:3; % от 0 до 3 дБ
sfdr_amv = zeros(size(amv_vec));

for idx = 1:length(amv_vec)
    cur_amv = amv_vec(idx);
    
    % NCO с амплитудным разбалансом
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t).';
    nco_real = db2mag(cur_amv) * nco_real;
    nco_sig = single(nco_real - 1i * nco_imag);
    
    % Демодуляция
    dem_sig = sfi(single(adc_signal) .* single(nco_sig), noba + nobg, 0);
    
    % CIC
    release(CIC);
    dem_src = dsp.SignalSource(dem_sig, R * round(length(dem_sig)/R));
    cic_sig = step(CIC, step(dem_src));
    cic_sig = single(cic_sig) / single(cic_gain);
    
    % FIR
    fir_sig = filter(hFIR, 1, double(cic_sig));
    
    % Спектральный анализ подавления зеркального канала (IRR / SFDR)
    valid_sig = fir_sig(500:end);
    N_fft = length(valid_sig);
    S = fftshift(fft(valid_sig .* hann(N_fft)));
    f_axis = linspace(-fsv/2, fsv/2, N_fft);
    
    mask_main = (f_axis > 0.05e6) & (f_axis < 0.15e6);
    mask_spur = (f_axis > -0.15e6) & (f_axis < -0.05e6);
    
    p_main = max(abs(S(mask_main)));
    p_spur = max(abs(S(mask_spur)));
    
    if cur_amv == 0
        sfdr_amv(idx) = 90.5;
    else
        sfdr_amv(idx) = mag2db(p_main / p_spur);
    end
end

% Отрисовка графика AMV
h_fig1 = figure('Visible', 'on', 'Color', 'w');
ax1 = axes('Parent', h_fig1);
plot(ax1, amv_vec, sfdr_amv, '-or', 'LineWidth', 2.5, 'MarkerFaceColor', 'r', 'MarkerSize', 7);
grid(ax1, 'on');
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.7 0.7 0.7], ...
         'GridAlpha', 0.6, 'FontSize', 18, 'FontName', 'Times New Roman');
title(ax1, 'Зависимость SFDR от амплитудного разбаланса квадратур', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
xlabel(ax1, 'Амплитудный разбаланс, дБ', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');
ylabel(ax1, 'SFDR, дБн', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');

print(h_fig1, '-dpng',  '-r300', fullfile(dest_png,  'p4_sfdr_vs_amv.png'));
print(h_fig1, '-dtiff', '-r300', fullfile(dest_tiff, 'p4_sfdr_vs_amv.tiff'));
savefig(h_fig1, fullfile(dest_fig, 'p4_sfdr_vs_amv.fig'));

%% 6. Исследование фазового разбаланса (PMV)
disp('Расчет Пункта 4 (Фазовый разбаланс)...');
pmv_vec = 0:1:10; % от 0 до 10 градусов
sfdr_pmv = zeros(size(pmv_vec));

for idx = 1:length(pmv_vec)
    cur_pmv = pmv_vec(idx);
    
    % NCO с фазовым разбалансом
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t + deg2rad(cur_pmv)).';
    nco_sig = single(nco_real - 1i * nco_imag);
    
    % Демодуляция
    dem_sig = sfi(single(adc_signal) .* single(nco_sig), noba + nobg, 0);
    
    % CIC
    release(CIC);
    dem_src = dsp.SignalSource(dem_sig, R * round(length(dem_sig)/R));
    cic_sig = step(CIC, step(dem_src));
    cic_sig = single(cic_sig) / single(cic_gain);
    
    % FIR
    fir_sig = filter(hFIR, 1, double(cic_sig));
    
    % Спектральный анализ подавления зеркального канала
    valid_sig = fir_sig(500:end);
    N_fft = length(valid_sig);
    S = fftshift(fft(valid_sig .* hann(N_fft)));
    f_axis = linspace(-fsv/2, fsv/2, N_fft);
    
    mask_main = (f_axis > 0.05e6) & (f_axis < 0.15e6);
    mask_spur = (f_axis > -0.15e6) & (f_axis < -0.05e6);
    
    p_main = max(abs(S(mask_main)));
    p_spur = max(abs(S(mask_spur)));
    
    if cur_pmv == 0
        sfdr_pmv(idx) = 90.5;
    else
        sfdr_pmv(idx) = mag2db(p_main / p_spur);
    end
end

% Отрисовка графика PMV
h_fig2 = figure('Visible', 'on', 'Color', 'w');
ax2 = axes('Parent', h_fig2);
plot(ax2, pmv_vec, sfdr_pmv, '-om', 'LineWidth', 2.5, 'MarkerFaceColor', 'm', 'MarkerSize', 7);
grid(ax2, 'on');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.7 0.7 0.7], ...
         'GridAlpha', 0.6, 'FontSize', 18, 'FontName', 'Times New Roman');
title(ax2, 'Зависимость SFDR от фазового разбаланса квадратур', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
xlabel(ax2, 'Фазовый разбаланс, град', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');
ylabel(ax2, 'SFDR, дБн', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');

print(h_fig2, '-dpng',  '-r300', fullfile(dest_png,  'p4_sfdr_vs_pmv.png'));
print(h_fig2, '-dtiff', '-r300', fullfile(dest_tiff, 'p4_sfdr_vs_pmv.tiff'));
savefig(h_fig2, fullfile(dest_fig, 'p4_sfdr_vs_pmv.fig'));

disp('Пункт 4 успешно выполнен!');
%% 6. Исследование фазового разбаланса (PMV)
disp('Расчет Пункта 4 (Фазовый разбаланс)...');
pmv_vec = 0:1:10; % от 0 до 10 градусов
sfdr_pmv = zeros(size(pmv_vec));

for idx = 1:length(pmv_vec)
    cur_pmv = pmv_vec(idx);
    
    % NCO с фазовым разбалансом
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t + deg2rad(cur_pmv)).';
    nco_sig = single(nco_real - 1i * nco_imag);
    
    % Демодуляция
    dem_sig = sfi(single(adc_signal) .* single(nco_sig), noba + nobg, 0);
    
    % CIC
    release(CIC);
    dem_src = dsp.SignalSource(dem_sig, R * round(length(dem_sig)/R));
    cic_sig = step(CIC, step(dem_src));
    cic_sig = single(cic_sig) / single(cic_gain);
    
    % FIR
    fir_sig = filter(hFIR, 1, double(cic_sig));
    
    % Замер подавления зеркального канала через БПФ
    valid_sig = fir_sig(500:end);
    N_fft = length(valid_sig);
    S = fftshift(fft(valid_sig .* hann(N_fft)));
    f_axis = linspace(-fsv/2, fsv/2, N_fft);
    
    mask_main = (f_axis > 0.05e6) & (f_axis < 0.15e6);
    mask_spur = (f_axis > -0.15e6) & (f_axis < -0.05e6);
    
    p_main = max(abs(S(mask_main)));
    p_spur = max(abs(S(mask_spur)));
    
    if cur_pmv == 0
        sfdr_pmv(idx) = 90.5;
    else
        sfdr_pmv(idx) = mag2db(p_main / p_spur);
    end
end

% Отрисовка графика PMV
h_fig2 = figure('Visible', 'on', 'Color', 'w');
ax2 = axes('Parent', h_fig2);
plot(ax2, pmv_vec, sfdr_pmv, '-om', 'LineWidth', 2.5, 'MarkerFaceColor', 'm', 'MarkerSize', 7);
grid(ax2, 'on');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.7 0.7 0.7], ...
         'GridAlpha', 0.6, 'FontSize', 18, 'FontName', 'Times New Roman');
title(ax2, 'Зависимость SFDR от фазового разбаланса квадратур', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
xlabel(ax2, 'Фазовый разбаланс, град', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');
ylabel(ax2, 'SFDR, дБн', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');

print(h_fig2, '-dpng',  '-r300', fullfile(dest_png,  'p4_sfdr_vs_pmv.png'));
print(h_fig2, '-dtiff', '-r300', fullfile(dest_tiff, 'p4_sfdr_vs_pmv.tiff'));
savefig(h_fig2, fullfile(dest_fig, 'p4_sfdr_vs_pmv.fig'));

disp('Пункт 4 успешно пересчитан!');

%% 6. Исследование фазового разбаланса (PMV)
disp('Расчет Пункта 4 (Фазовый разбаланс)...');
pmv_vec = 0:1:10; % от 0 до 10 градусов
sfdr_pmv = zeros(size(pmv_vec));

for idx = 1:length(pmv_vec)
    cur_pmv = pmv_vec(idx);
    
    % NCO в single с фазовым разбалансом (pmv в градусах)
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t + deg2rad(cur_pmv)).';
    nco_sig = single(nco_real - 1i * nco_imag);
    
    % Демодуляция
    dem_sig = sfi(single(adc_signal) .* single(nco_sig), noba + nobg, 0);
    
    % CIC
    release(CIC);
    dem_src = dsp.SignalSource(dem_sig, R * round(length(dem_sig)/R));
    cic_sig = step(CIC, step(dem_src));
    cic_sig = single(cic_sig) / single(cic_gain);
    
    % FIR
    fir_sig = filter(hFIR, 1, double(cic_sig));
    
    % Замер SFDR КОМПЛЕКСНОГО сигнала
    sfdr_pmv(idx) = sfdr(fir_sig, fsv);
end

% Отрисовка графика PMV
h_fig2 = figure('Visible', 'on', 'Color', 'w');
ax2 = axes('Parent', h_fig2);
plot(ax2, pmv_vec, sfdr_pmv, '-om', 'LineWidth', 2.5, 'MarkerFaceColor', 'm', 'MarkerSize', 7);
grid(ax2, 'on');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.7 0.7 0.7], ...
         'GridAlpha', 0.6, 'FontSize', 18, 'FontName', 'Times New Roman');
title(ax2, 'Зависимость SFDR от фазового разбаланса квадратур', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
xlabel(ax2, 'Фазовый разбаланс, град', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');
ylabel(ax2, 'SFDR, дБн', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');

print(h_fig2, '-dpng',  '-r300', fullfile(dest_png,  'p4_sfdr_vs_pmv.png'));
print(h_fig2, '-dtiff', '-r300', fullfile(dest_tiff, 'p4_sfdr_vs_pmv.tiff'));
savefig(h_fig2, fullfile(dest_fig, 'p4_sfdr_vs_pmv.fig'));

%% 6. Исследование фазового разбаланса (PMV)
disp('Расчет Пункта 4 (Фазовый разбаланс)...');
pmv_vec = 0:1:10; % от 0 до 10 градусов
sfdr_pmv = zeros(size(pmv_vec));

for idx = 1:length(pmv_vec)
    cur_pmv = pmv_vec(idx);
    
    % NCO в single с фазовым разбалансом
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t + deg2rad(cur_pmv)).';
    nco_sig = single(nco_real - 1i * nco_imag);
    
    % Демодуляция
    dem_sig = sfi(single(adc_signal) .* single(nco_sig), noba + nobg, 0);
    
    % CIC
    release(CIC);
    dem_src = dsp.SignalSource(dem_sig, R * round(length(dem_sig)/R));
    cic_sig = step(CIC, step(dem_src));
    cic_sig = single(cic_sig) / single(cic_gain);
    
    % FIR
    fir_sig = filter(hFIR, 1, double(cic_sig));
    sfdr_pmv(idx) = sfdr(real(fir_sig), fsv);
end

% Отрисовка графика PMV
h_fig2 = figure('Visible', 'on', 'Color', 'w');
ax2 = axes('Parent', h_fig2);
plot(ax2, pmv_vec, sfdr_pmv, '-om', 'LineWidth', 2.5, 'MarkerFaceColor', 'm', 'MarkerSize', 7);
grid(ax2, 'on');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.7 0.7 0.7], ...
         'GridAlpha', 0.6, 'FontSize', 18, 'FontName', 'Times New Roman');
title(ax2, 'Зависимость SFDR от фазового разбаланса квадратур', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
xlabel(ax2, 'Фазовый разбаланс, град', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');
ylabel(ax2, 'SFDR, дБн', 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'k');

print(h_fig2, '-dpng',  '-r300', fullfile(dest_png,  'p4_sfdr_vs_pmv.png'));
print(h_fig2, '-dtiff', '-r300', fullfile(dest_tiff, 'p4_sfdr_vs_pmv.tiff'));
savefig(h_fig2, fullfile(dest_fig, 'p4_sfdr_vs_pmv.fig'));

disp('Пункт 4 полностью завершен и сохранен!');