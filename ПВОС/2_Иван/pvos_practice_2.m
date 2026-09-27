clear all; close all; clc;

length_test = 24e3;

destdirectory_FIG = 'export_figs'; 
if ~exist(destdirectory_FIG, 'dir')
    mkdir(destdirectory_FIG);
end

%% Параметры устройства (Вариант 8)
% N = 8
fs = (100 + 8*10) * 1e6;       % 180 МГц[cite: 12]
fc = (130 + 8*10) * 1e6 + 0.1e6; % 210.1 МГц (с отстройкой 0.1 МГц, чтобы сигнал не лег в 0)[cite: 12]
fg = (130 + 8*10) * 1e6;       % 210.0 МГц[cite: 12]
noba = 5 + 8;                  % 13 бит[cite: 12]
nobg = 5 + 8;                  % 13 бит[cite: 12]
nobo = 13;

R = 10;                        % Коэффициент децимации
fsv = fs / R;                  % 18 МГц
bandwidth = fsv/2/1.2;         % Рабочая полоса

snr = 70; 
noise_power = 1; 
signal_power = noise_power*db2pow(snr);

%% Формирование входного сигнала
t = (0:length_test - 1)/fs;
noise = sqrt(noise_power) * randn(length_test, 1);
signal = sqrt(signal_power)*sin(2*pi*fc*t).';
in = signal + noise;

%% Входной сигнал и спектр
figure
    plot(t*1e3, in, 'LineWidth', 2)
    grid on; axis tight;
    title('Signal before ADC'); xlabel('t, ms'); ylabel('LSB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

[fin, sin_spec] = get_spectrum(in, fs, 1);
figure
    plot(fin/1e6, mag2db(abs(sin_spec)), 'LineWidth', 2)
    grid on; axis tight;
    title('Spectrum of signal before ADC'); xlabel('f, MHz'); ylabel('dB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

%% АЦП
adc_signal = get_adc(in, noba);

figure
    stairs(adc_signal, 'LineWidth', 2)
    grid on; axis tight;
    title('Signal after ADC'); xlabel('t, samples'); ylabel('LSB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')
    
[fadc, sadc] = get_spectrum(adc_signal, fs, 1);
figure
    plot(fadc/1e6, mag2db(abs(sadc)), 'LineWidth', 2)
    grid on; axis tight;
    title('Spectrum of signal after ADC'); xlabel('f, MHz'); ylabel('dB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

%% Гетеродин (NCO) и Демодулятор
nco_type = 'single';
amv = 0;
pmv = 0;
[nco_signal, nco_gain] = get_nco(nobg, fg, t, nco_type, amv, pmv);

figure
    plot(real(nco_signal), 'LineWidth', 2)
    hold on; plot(imag(nco_signal), 'LineWidth', 2)
    grid on; axis tight;
    title('Signal after NCO'); xlabel('t, samples'); ylabel('LSB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

[fnco, snco] = get_spectrum(nco_signal, fs, 1);
figure
    plot(fnco/1e6, mag2db(abs(snco)), 'LineWidth', 2) 
    grid on; axis tight;
    title('Spectrum of signal after NCO'); xlabel('f, MHz'); ylabel('dB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

dem_signal = get_dem(adc_signal, nco_signal, noba, nobg, 'single');

figure
    plot(real(dem_signal), 'LineWidth', 2)
    hold on; plot(imag(dem_signal), 'LineWidth', 2)
    grid on; axis tight;
    title('Signal after DEM'); xlabel('t, samples'); ylabel('LSB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

[fdem, sdem] = get_spectrum(dem_signal, fs, 1);
figure
    plot(fdem/1e6, mag2db(abs(sdem)), 'LineWidth', 2) 
    grid on; axis tight;
    title('Spectrum of signal after DEM'); xlabel('f, MHz'); ylabel('dB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

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
    title('Signal after CIC'); xlabel('t, samples'); ylabel('LSB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

[fcic, scic] = get_spectrum(cic_signal, fsv, 1);
figure
    plot(fcic/1e6, mag2db(abs(scic)), 'LineWidth', 2) 
    grid on; axis tight;
    title('Spectrum of signal after CIC'); xlabel('f, MHz'); ylabel('dB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

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
    title('Signal after FIR'); xlabel('t, samples'); ylabel('LSB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')

[ffir, sfir] = get_spectrum(fir_signal, fsv, 1);
figure
    plot(ffir/1e6, mag2db(abs(sfir)), 'LineWidth', 2) 
    grid on; axis tight;
    title('Spectrum of signal after FIR'); xlabel('f, MHz'); ylabel('dB')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')    

% АЧХ каскада
FC = dsp.FilterCascade(CIC, FIR);                                                  
RESP = fvtool(CIC, FIR, FC, 'Fs', [fs fsv fs]);
RESP.NormalizeMagnitudeto1 = 'on';
щ                                                                                                                 й
%% Пункт 3: Зависимость SFDR от количества бит ЦГ (nco_type = 'fixed')
% Построение итогового графика зависимости на чистом белом фоне
h_fig = figure('Color', 'w');
ax = axes('Parent', h_fig);

plot(ax, test_ng, sfdr_vals, '-ob', 'LineWidth', 2.5, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
grid(ax, 'on');

% Принудительное назначение цветов осей, сетки и фона
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

% Сохранение графиков
print(h_fig, '-dpng', '-r300', fullfile(destdirectory_FIG, 'p3_sfdr_vs_bits.png'));
print(h_fig, '-dtiff', '-r300', fullfile(fullfile(destdirectory_FIG, 'pic_tiff'), 'p3_sfdr_vs_bits.tiff'));
savefig(h_fig, fullfile(fullfile(destdirectory_FIG, 'pic_fig'), 'p3_sfdr_vs_bits.fig'));

%% Функции
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