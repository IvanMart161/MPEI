clear all; close all; clc;
length_test = 24e3;
destdirectory_FIG = 'export_figs'; 
if ~exist(destdirectory_FIG, 'dir')
    mkdir(destdirectory_FIG);
end

%% Параметры устройства (Вариант 8)
fs = (100 + 8*10) * 1e6;         % 180 МГц
fc = (130 + 8*10) * 1e6 + 0.1e6; % 210.1 МГц
fg = (130 + 8*10) * 1e6;         % 210.0 МГц
noba = 5 + 8;                    % 13 бит
nobg = 5 + 8;                    % 13 бит
nobo = 13;

R = 10;                          % Коэффициент децимации
fsv = fs / R;                    % 18 МГц
bandwidth = fsv/2/1.2;           % Рабочая полоса

snr = 70; 
noise_power = 1; 
signal_power = noise_power*db2pow(snr);

%% Формирование входного сигнала
t = (0:length_test - 1)/fs;
noise = sqrt(noise_power) * randn(length_test, 1);
signal = sqrt(signal_power)*sin(2*pi*fc*t).';
in = signal + noise;

%% Входной сигнал и спектр
zoom_pts = 150; 
figure
plot(t*1e3, in, 'LineWidth', 1.5)
grid on; axis tight;
xlim([0, t(zoom_pts)*1e3]);
title('Signal before ADC'); xlabel('t, ms'); ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

[fin, sin_spec] = get_spectrum(in, fs, 1);
figure
plot(fin/1e6, mag2db(abs(sin_spec)), 'LineWidth', 1.2)
grid on; axis tight; ylim([-20 150]);
title('Spectrum of signal before ADC'); xlabel('f, MHz'); ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

%% АЦП
adc_signal = get_adc(in, noba);
figure
stairs(adc_signal, 'LineWidth', 1.5)
grid on; axis tight;
xlim([0, zoom_pts]);
title('Signal after ADC'); xlabel('t, samples'); ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
    
[fadc, sadc] = get_spectrum(adc_signal, fs, 1);
figure
plot(fadc/1e6, mag2db(abs(sadc)), 'LineWidth', 1.2)
grid on; axis tight; ylim([-20 150]);
title('Spectrum of signal after ADC'); xlabel('f, MHz'); ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

%% Гетеродин (NCO)
nco_type = 'single';
amv = 0;
pmv = 0;
[nco_signal, nco_gain] = get_nco(nobg, fg, t, nco_type, amv, pmv);

figure
plot(real(nco_signal), 'LineWidth', 1.5)
hold on; plot(imag(nco_signal), 'LineWidth', 1.5)
grid on; axis tight;
xlim([0, zoom_pts]);
title('Signal after NCO'); xlabel('t, samples'); ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

%% Спектр гетеродина (NCO)
[fnco, snco] = get_spectrum(nco_signal, fs, 1);
snco_db = mag2db(abs(snco));
figure
plot(fnco/1e6, snco_db, 'LineWidth', 1.2) 
grid on; axis tight; ylim([-200 120]);
title('Spectrum of signal after NCO'); xlabel('f, MHz'); ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
hold on;
plot(-30.0, max(snco_db), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(-30.0, max(snco_db) + 20, 'f_{G,alias} = -30.0 MHz', ...
    'FontName', 'Times New Roman', 'FontSize', 28, ...
    'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');
%% Демодулятор (Смеситель)
dem_signal = get_dem(adc_signal, nco_signal, noba, nobg, 'single');

figure
plot(real(dem_signal), 'LineWidth', 1.5)
hold on; plot(imag(dem_signal), 'LineWidth', 1.5)
grid on; axis tight;
xlim([0, zoom_pts]);
title('Signal after DEM'); xlabel('t, samples'); ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

[fdem, sdem] = get_spectrum(dem_signal, fs, 1);
sdem_db = mag2db(abs(sdem));
figure
plot(fdem/1e6, sdem_db, 'LineWidth', 1.2) 
grid on; axis tight; 
xlim([-90 90]);
ylim([-20 200]); % Подняли потолок, чтобы 2 строки 28 кегля встали идеально
title('Spectrum of signal after DEM'); xlabel('f, MHz'); ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

[~, idx_if]  = min(abs(fdem/1e6 - 0.1));
[~, idx_sum] = min(abs(fdem/1e6 - (-60.1)));
hold on;

% Метка +0.1 МГц в 2 строки
plot(fdem(idx_if)/1e6, sdem_db(idx_if), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(fdem(idx_if)/1e6, sdem_db(idx_if) + 22, {'f_{IF} =', '+0.1 MHz'}, ...
    'FontName', 'Times New Roman', 'FontSize', 28, ...
    'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');

% Метка -60.1 МГц в 2 строки
plot(fdem(idx_sum)/1e6, sdem_db(idx_sum), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(fdem(idx_sum)/1e6, sdem_db(idx_sum) + 22, {'f_{sum} =', '-60.1 MHz'}, ...
    'FontName', 'Times New Roman', 'FontSize', 28, ...
    'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');
%% CIC фильтр
N = 8;
CIC = dsp.CICDecimator(R, 1, N);
CIC.FixedPointDataType = 'Full precision';
cic_gain = R^N/2;

dem_src = dsp.SignalSource(dem_signal, R*round(length(dem_signal)/R));    
cic_signal = step(CIC, step(dem_src));    
cic_signal = single(cic_signal) / single(cic_gain); 

figure
plot(real(cic_signal), 'LineWidth', 2)
hold on; plot(imag(cic_signal), 'LineWidth', 2)
grid on; axis tight;
xlim([0 450]); 
title('Signal after CIC'); xlabel('t, samples'); ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

[fcic, scic] = get_spectrum(cic_signal, fsv, 1);
figure
plot(fcic/1e6, mag2db(abs(scic)), 'LineWidth', 1.5) 
grid on; axis tight;
title('Spectrum of signal after CIC'); xlabel('f, MHz'); ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

%% КИХ-компенсатор (FIR)
fPass = bandwidth/2;
passbandRipple = 0.1;
fStop = bandwidth/2 + 500e3;
stopbandAttenuation = 60;

FIR = dsp.CICCompensationDecimator(CIC, ...
    'DecimationFactor', 1, ...
    'PassbandFrequency', fPass, ...
    'PassbandRipple', passbandRipple, ...
    'StopbandAttenuation', stopbandAttenuation, ...
    'StopbandFrequency', fStop, ...
    'SampleRate', fsv, ...
    'DesignForMinimumOrder', true);

hFIR = impz(FIR);
hFIR = hFIR ./ max(hFIR);
fir_signal = filter(hFIR, 1, double(cic_signal));

figure
plot(real(fir_signal), 'LineWidth', 2)
hold on; plot(imag(fir_signal), 'LineWidth', 2)
grid on; axis tight;
xlim([100 550]); 
title('Signal after FIR'); xlabel('t, samples'); ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

[ffir, sfir] = get_spectrum(fir_signal, fsv, 1);
figure
plot(ffir/1e6, mag2db(abs(sfir)), 'LineWidth', 1.5) 
grid on; axis tight;
title('Spectrum of signal after FIR'); xlabel('f, MHz'); ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')

%% Автоматический экспорт всех графиков (Пункт 2)
disp('Перекраска и сохранение графиков...');
figs = findobj('Type', 'figure'); 
for i = 1:length(figs)
    fig = figs(i);
    set(fig, 'Color', 'w');
    
    axs = findobj(fig, 'Type', 'axes');
    for j = 1:length(axs)
        ax = axs(j);
        set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
            'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5);
        
        t_label = get(ax, 'Title');  set(t_label, 'Color', 'k');
        x_label = get(ax, 'XLabel'); set(x_label, 'Color', 'k');
        y_label = get(ax, 'YLabel'); set(y_label, 'Color', 'k');
    end
    
    filename = fullfile(destdirectory_FIG, sprintf('p2_fig_%02d.png', fig.Number));
    print(fig, '-dpng', '-r300', filename);
end
disp('Все 12 графиков успешно сохранены в export_figs!');

%% Функции (строго в конце файла)
function out = get_adc(in, bits)
    out = quantize(int32(in), 1, bits, 0, 'Floor', 'Saturate');
end

function [nco_signal, nco_gain] = get_nco(nobg, fg, t, nco_type, amv, pmv)
    nco_real = cos(2*pi*fg*t).';
    nco_imag = sin(2*pi*fg*t + deg2rad(pmv)).';
    nco_real = max(real(nco_real))*db2mag(amv)*nco_real;
    nco = nco_real - 1i*nco_imag;
    nco_gain = 2^(nobg - 1) - 1;
    switch nco_type
        case 'fixed'
            nco = nco_gain.*nco;
            nco_signal = sfi(nco, nobg, 0);
        case 'single'
            nco_signal = single(nco);
        case 'double'
            nco_signal = double(nco);
    end
end

function dem_signal = get_dem(adc_signal, nco_signal, noba, nobg, dem_type)
    switch dem_type
        case 'fixed'
            dem_signal = adc_signal.*nco_signal;
            dem_signal = sfi(dem_signal, noba + nobg, 0);
        case 'single'
            dem_signal = single(adc_signal).*single(nco_signal);
            dem_signal = sfi(dem_signal, noba + nobg, 0);
        case 'double'
            dem_signal = double(adc_signal).*double(nco_signal);
            dem_signal = sfi(dem_signal, noba + nobg, 0);
    end
end

function [f, s] = get_spectrum(in, fs, weighting)
    f = linspace(-fs/2, fs/2, length(in));
    if weighting
        s = fftshift(fft(double(in).*hann(length(in))));
    else
        s = fftshift(fft(double(in)));
    end
end