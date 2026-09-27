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
fs = (100 + N_var*10) * 1e6;       % 180 МГц
fc = (130 + N_var*10) * 1e6 + 0.1e6; % 210.1 МГц
fg = (130 + N_var*10) * 1e6;       % 210.0 МГц
noba = 5 + N_var;                  % 13 бит
R = 10;
fsv = fs / R;                      % 18 МГц
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

%% 5. Быстрый цикл пункта 3 (без всплывающих окон)
disp('Запуск быстрого расчета Пункта 3...');

test_ng = 2:(noba + 2); % 2:15 бит
sfdr_vals = zeros(size(test_ng));

for idx = 1:length(test_ng)
    cur_ng = test_ng(idx);
    
    % NCO в режиме fixed
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t).';
    nco = nco_real - 1i * nco_imag;
    nco_gain = 2^(cur_ng - 1) - 1;
    nco_signal = sfi(nco_gain .* nco, cur_ng, 0);
    
    % Демодуляция
    dem_signal = sfi(adc_signal .* nco_signal, noba + cur_ng, 0);
    
    % CIC
    release(CIC);
    dem_src = dsp.SignalSource(dem_signal, R * round(length(dem_signal)/R));
    cic_signal = step(CIC, step(dem_src));
    cic_signal = single(cic_signal) / single(cic_gain);
    
    % FIR
    fir_signal = filter(hFIR, 1, double(cic_signal));
    
    % Численный замер SFDR
    sfdr_vals(idx) = sfdr(real(fir_signal), fsv);
    fprintf('n_g = %2d бит | SFDR = %.2f дБн\n', cur_ng, sfdr_vals(idx));
end

%% 6. Отрисовка чистого графика на белом фоне
h_fig = figure('Visible', 'on', 'Color', 'w');
ax = axes('Parent', h_fig);

plot(ax, test_ng, sfdr_vals, '-ob', 'LineWidth', 2.5, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
grid(ax, 'on');

set(ax, 'Color', 'w', ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'GridColor', [0.7 0.7 0.7], ...
        'GridAlpha', 0.6, ...
        'FontSize', 18, ...
        'FontName', 'Times New Roman');

title(ax, 'Зависимость SFDR от количества бит ЦГ', ...
    'FontSize', 22, 'FontName', 'Times New Roman', 'Color', 'k', 'FontWeight', 'bold');
xlabel(ax, 'Количество бит ЦГ, n_g', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k');
ylabel(ax, 'SFDR, дБн', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'Color', 'k');

xlim(ax, [min(test_ng) max(test_ng)]);
ylim(ax, [25 95]);

%% 7. Сохранение
print(h_fig, '-dpng',  '-r300', fullfile(dest_png,  'p3_sfdr_vs_bits.png'));
print(h_fig, '-dtiff', '-r300', fullfile(dest_tiff, 'p3_sfdr_vs_bits.tiff'));
savefig(h_fig, fullfile(dest_fig, 'p3_sfdr_vs_bits.fig'));

disp('Готово! График построен и сохранен во всех папках.');