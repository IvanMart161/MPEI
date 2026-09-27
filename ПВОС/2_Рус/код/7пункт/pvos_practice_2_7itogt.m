clear all; close all; clc;

length_test = 24e3;

%% Some properties
snr = 70; % signal to noise ratio, [dB]
noise_power = 1; % power of noise
signal_power = noise_power*db2pow(snr); % power of signal

fs = 140e6; % sampling frequency before decimation, [Hz] (Частота дискретизации АЦП)
fsv = 14e6; % sampling frequency after decimation, [Hz] (Сохраняем коэффициент децимации R=10)

fc = 170.1e6; % if mono = 0 fc = 170e6; (Промежуточная частота входного сигнала)
fg = 170e6; % center demodulation geterodin freq, [Hz] (Центральная частота ЦГ)

noba = 9; % number of bits after analog to digital conversion (Количество бит АЦП)
nobg = 9; % number of bits for digital geterodin (Количество бит ЦГ)
nobo = 14; % number of bits for output

bandwidth = fsv/2/1.2; % maximum bandwith of signal, [Hz]
%% Some code here
t = (0:length_test - 1)/fs; % time vector, s
noise = sqrt(noise_power)*(normrnd(0,1,[length_test 1])); % noise generate

signal = sqrt(signal_power)*sin(2*pi*fc*t).';

in = signal + noise; % input generate

%% INPUT
% Plot results
figure
    plot(t*1e3,in,'LineWidth',2)
    grid on
    axis tight
    title('Signal before ADC')
    xlabel('t, ms')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

% Spectrum
[fin,sin] = get_spectrum(in,fs,1);

% Plot results
figure
    plot(fin/1e6,mag2db(abs(sin)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal before ADC')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

%% ADC
adc_signal = get_adc(in,noba);

% Plot results
figure
    stairs(adc_signal,'LineWidth',2)
    grid on
    axis tight
    title('Signal after ADC')
    xlabel('t, samples')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')
    
% Spectrum ADC
[fadc,sadc] = get_spectrum(adc_signal,fs,1);

% Plot results
figure
    plot(fadc/1e6,mag2db(abs(sadc)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal after ADC')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

%% NCO
nco_type = 'single';

amv = 0;
pmv = 0;

[nco_signal,nco_gain] = get_nco(nobg,fg,t,nco_type,amv,pmv);

% Plot results
figure
    plot(real(nco_signal),'LineWidth',2)
    hold on
    plot(imag(nco_signal),'LineWidth',2)
    grid on
    axis tight
    title('Signal after NCO')
    xlabel('t, samples')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

% Spectrum NCO
[fnco,snco] = get_spectrum(nco_signal,fs,1);

% Plot results
figure
    plot(fnco/1e6,mag2db(abs(snco)),'LineWidth',2) 
    grid on
    axis tight
    title('Spectrum of signal after NCO')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

dem_signal = get_dem(adc_signal,nco_signal,noba,nobg,'single');

% Plot results
figure
    plot(real(dem_signal),'LineWidth',2)
    hold on
    plot(imag(dem_signal),'LineWidth',2)
    grid on
    axis tight
    title('Signal after NCO')
    xlabel('t, samples')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

% Spectrum NCO
[fdem,sdem] = get_spectrum(dem_signal,fs,1);

% Plot results
figure
    plot(fdem/1e6,mag2db(abs(sdem)),'LineWidth',2) 
    grid on
    axis tight
    title('Spectrum of signal after NCO')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')
    
%% Practice 2 start here

R = fs/fsv; % decimation factor
N = 8; % number of section

%% GET CIC FILTER

% Create object
% 
% CIC = dsp.CICDecimator(R,1,N); % cic filter object
% 
% CIC.FixedPointDataType = 'Full precision';
% 
% % Gain calculation
% cic_gain = R^N/2;
% 
% % Input generate
% dem_signal = dsp.SignalSource(dem_signal, R*round(length(dem_signal)/R));    
% 
% % CIC filter output generate
% cic_signal = step(CIC,step(dem_signal));    
% 
% % Normscale by gain
% cic_signal = single(cic_signal)/single(cic_gain); 
% 
% % Plot results
% figure
%     plot(real(cic_signal),'LineWidth',2)
%     hold on
%     plot(imag(cic_signal),'LineWidth',2)
%     grid on
%     axis tight
%     title('Signal after CIC')
%     xlabel('t, samples')
%     ylabel('LSB')
%     set(gca,'Fontsize',28,'Fontname','Times New Roman')
% 
% % Spectrum CIC
% [fcic,scic] = get_spectrum(cic_signal,fsv,1);
% 
% % Plot results
% figure
%     plot(fcic/1e6,mag2db(abs(scic)),'LineWidth',2) 
%     grid on
%     axis tight
%     title('Spectrum of signal after CIC')
%     xlabel('f, MHz')
%     ylabel('dB')
%     set(gca,'Fontsize',28,'Fontname','Times New Roman')

% % FIR
% % FIR properties
% fPass = bandwidth/2;
% passbandRipple = 0.1;
% fStop = bandwidth/2 + 500e3;
% stopbandAttenuation = 60;
% 
% % Create object
% FIR = dsp.CICCompensationDecimator(CIC, ...
%     'DecimationFactor',1,...
%     'PassbandFrequency',fPass,...
%     'PassbandRipple',passbandRipple,...
%     'StopbandAttenuation',stopbandAttenuation,...
%     'StopbandFrequency',fStop,...
%     'SampleRate',fsv,...
%     'DesignForMinimumOrder',true);
% 
% % Impulse responce generate
% hFIR = impz(FIR);
% hFIR = hFIR./max(hFIR);
% 
% % Generate signal
% fir_signal = filter(hFIR,1,double(cic_signal)); % ./nco_gain;
% 
% % fir_signal = sfi(fir_signal,nobo,0);
% 
% % Plot results
% figure
%     plot(real(fir_signal),'LineWidth',2)
%     hold on
%     plot(imag(fir_signal),'LineWidth',2)
%     grid on
%     axis tight
%     title('Signal after FIR')
%     xlabel('t, samples')
%     ylabel('LSB')
%     set(gca,'Fontsize',28,'Fontname','Times New Roman')
% 
% % Spectrum CIC
% [ffir,sfir] = get_spectrum(fir_signal,fsv,1);
% 
% % Plot results
% figure
%     plot(ffir/1e6,mag2db(abs(sfir)),'LineWidth',2) 
%     grid on
%     axis tight
%     title('Spectrum of signal after FIR')
%     xlabel('f, MHz')
%     ylabel('dB')
%     set(gca,'Fontsize',28,'Fontname','Times New Roman')    
% 
% % Calculate cascade filter impulse responce
% FC = dsp.FilterCascade(CIC,FIR);

%% PRACTICE 2: FIR DECIMATOR (Пункт 7)

R = fs/fsv; % Коэффициент децимации (R = 10 для варианта 4)

% Параметры КИХ-фильтра из задания
fPass = bandwidth/2;
passbandRipple = 0.1;
fStop = bandwidth/2 + 500e3;
stopbandAttenuation = 60;

% 1. Расчет коэффициентов ФНЧ на исходной частоте дискретизации (fs)
lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', fPass, ...
    'StopbandFrequency', fStop, ...
    'PassbandRipple', passbandRipple, ...
    'StopbandAttenuation', stopbandAttenuation, ...
    'SampleRate', fs);

% 2. Создание объекта КИХ-дециматора
% Используем рассчитанные коэффициенты (Numerator) и задаем децимацию (R)
FIRDecim = dsp.FIRDecimator('DecimationFactor', R, 'Numerator', lpFilt.Coefficients);

% 3. Подготовка входного сигнала
% Обрезаем сигнал до длины, кратной коэффициенту децимации
len_to_use = R * floor(length(dem_signal)/R);
dem_signal_adj = dem_signal(1:len_to_use);

% 4. Фильтрация и децимация отсчетов в заданное число раз
fir_signal = step(FIRDecim, double(dem_signal_adj));

% --- Вывод результатов ---

% Временная диаграмма после КИХ-дециматора
figure
    plot(real(fir_signal), 'LineWidth', 2)
    hold on
    plot(imag(fir_signal), 'LineWidth', 2)
    grid on
    axis tight
    title('Signal after FIR Decimator')
    xlabel('t, samples')
    ylabel('LSB')
    set(gca, 'Fontsize', 28, 'Fontname', 'Times New Roman')

% Спектр после КИХ-дециматора
[ffir, sfir] = get_spectrum(fir_signal, fsv, 1);

figure
    plot(ffir/1e6, mag2db(abs(sfir)), 'LineWidth', 2) 
    grid on
    axis tight
    title('Spectrum of signal after FIR Decimator')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca, 'Fontsize', 28, 'Fontname', 'Times New Roman')    

% Просмотр АЧХ и ФЧХ спроектированного фильтра-дециматора
fvtool(FIRDecim, 'Fs', fs);

%% ПУНКТ 3: Исследование влияния разрядности ЦГ на SFDR

% Исходные данные для пункта 3
nco_type_p3 = 'fixed'; % Фиксированная точка для гетеродина
amv_p3 = 0; % Без амплитудных искажений
pmv_p3 = 0; % Без фазовых искажений

% Диапазон изменения бит: от 2 до na + 2 (9 + 2 = 11)
nobg_array = 2:(noba + 2); 
sfdr_values = zeros(size(nobg_array));

% Подготовка длины сигнала (чтобы избежать ошибок размерности при децимации)
len_to_use = R * floor(length(adc_signal)/R);
adc_signal_adj = adc_signal(1:len_to_use);
t_adj = t(1:len_to_use);

for k = 1:length(nobg_array)
    current_nobg = nobg_array(k);
    
    % 1. Генерация сигнала ЦГ с текущей разрядностью
    [nco_sig_loop, nco_gain_loop] = get_nco(current_nobg, fg, t_adj, nco_type_p3, amv_p3, pmv_p3);
    
    % 2. Демодуляция (умножение сигнала АЦП на ЦГ)
    dem_sig_loop = get_dem(adc_signal_adj, nco_sig_loop, noba, current_nobg, nco_type_p3);
    
    % 3. Фильтрация и децимация через единый КИХ-фильтр
    % Обязательно сбрасываем состояния фильтра перед каждой новой итерацией
    reset(FIRDecim); 
    fir_sig_loop = step(FIRDecim, double(dem_sig_loop));
    
    % 4. Измерение SFDR
    % Анализируем динамический диапазон свободный от помех (SFDR)
    % Измерение проводится для синфазной (вещественной) составляющей
    sfdr_values(k) = sfdr(real(fir_sig_loop), fsv);
end

% Построение графика зависимости SFDR от разрядности
figure
    plot(nobg_array, sfdr_values, '-o', 'LineWidth', 2, 'MarkerSize', 8)
    grid on
    title('Зависимость SFDR от разрядности ЦГ')
    xlabel('Количество бит ЦГ, n_g')
    ylabel('SFDR, дБн')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')
    
% Автоматическое сохранение результатов моделирования (п. 2 задания)
saveas(gcf, 'SFDR_vs_NCO_bits.fig');
saveas(gcf, 'SFDR_vs_NCO_bits.tiff');

%% ПУНКТ 4: Исследование амплитудных и фазовых искажений

% Исходные данные для пункта 4
nco_type_p4 = 'single'; % Одинарная точность для анализа искажений

% 1. Исследование фазовых искажений (pmv)
pmv_array = 0:1:15; % Изменение фазы от 0 до 15 градусов
sfdr_pmv = zeros(size(pmv_array));
amv_fixed = 0; % Амплитудные искажения отключены

for k = 1:length(pmv_array)
    current_pmv = pmv_array(k);
    
    % Генерация ЦГ, демодуляция, фильтрация
    [nco_sig_loop, ~] = get_nco(nobg, fg, t_adj, nco_type_p4, amv_fixed, current_pmv);
    dem_sig_loop = get_dem(adc_signal_adj, nco_sig_loop, noba, nobg, nco_type_p4);
    
    reset(FIRDecim); % Сброс состояний КИХ-фильтра
    fir_sig_loop = step(FIRDecim, double(dem_sig_loop));
    
    sfdr_pmv(k) = sfdr(real(fir_sig_loop), fsv);
end

% Построение графика для фазовых искажений
figure
    plot(pmv_array, sfdr_pmv, '-o', 'LineWidth', 2, 'MarkerSize', 8)
    grid on
    title('Влияние фазовых искажений на SFDR')
    xlabel('Фазовое искажение (pmv), градусы')
    ylabel('SFDR, дБн')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')
    
saveas(gcf, 'SFDR_vs_Phase_Mismatch.fig');
saveas(gcf, 'SFDR_vs_Phase_Mismatch.tiff');

% 2. Исследование амплитудных искажений (amv)
amv_array = 0:0.2:3; % Изменение амплитуды от 0 до 3 дБ
sfdr_amv = zeros(size(amv_array));
pmv_fixed = 0; % Фазовые искажения отключены

for k = 1:length(amv_array)
    current_amv = amv_array(k);
    
    % Генерация ЦГ, демодуляция, фильтрация
    [nco_sig_loop, ~] = get_nco(nobg, fg, t_adj, nco_type_p4, current_amv, pmv_fixed);
    dem_sig_loop = get_dem(adc_signal_adj, nco_sig_loop, noba, nobg, nco_type_p4);
    
    reset(FIRDecim); % Сброс состояний КИХ-фильтра
    fir_sig_loop = step(FIRDecim, double(dem_sig_loop));
    
    sfdr_amv(k) = sfdr(real(fir_sig_loop), fsv);
end

% Построение графика для амплитудных искажений
figure
    plot(amv_array, sfdr_amv, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', '#D95319')
    grid on
    title('Влияние амплитудных искажений на SFDR')
    xlabel('Амплитудное искажение (amv), дБ')
    ylabel('SFDR, дБн')
    set(gca, 'Fontsize', 20, 'Fontname', 'Times New Roman')
    
saveas(gcf, 'SFDR_vs_Amplitude_Mismatch.fig');
saveas(gcf, 'SFDR_vs_Amplitude_Mismatch.tiff');

%% ПУНКТ 5: Монохроматический сигнал с линейно возрастающей амплитудой (LA)

% 1. Настройки частоты и динамического диапазона
fc_la = fg + 1e6; % Смещение на 1 МГц для попадания в рабочую полосу
max_amp = 2^(noba - 1) - 1; % Макс. амплитуда для АЦП (при 9 битах = 255)

% 2. Формирование линейно возрастающей огибающей
A_la = linspace(0, max_amp, length_test).'; 

% 3. Генерация входного тестового сигнала (чистый сигнал без шума)
%in_la = A_la .* builtin('sin', 2*pi*fc_la*t).';
noise_la = sqrt(noise_power) * normrnd(0, 1, [length_test 1]); % Генерируем вектор шума
in_la = A_la .* builtin('sin', 2*pi*fc_la*t).' + noise_la;

% 4. Прохождение тракта: АЦП -> ЦГ -> Демодулятор
adc_signal_la = get_adc(in_la, noba);
[nco_sig_la, ~] = get_nco(nobg, fg, t, 'single', 0, 0);
dem_sig_la = get_dem(adc_signal_la, nco_sig_la, noba, nobg, 'single');

% 5. Фильтрация и децимация (используем ранее созданный КИХ-фильтр)
len_to_use_la = R * floor(length(dem_sig_la)/R);
dem_sig_la_adj = dem_sig_la(1:len_to_use_la);

reset(FIRDecim); % Очистка памяти фильтра
fir_sig_la = step(FIRDecim, double(dem_sig_la_adj));

% --- ВИЗУАЛИЗАЦИЯ И СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ---

% График 1: Временная диаграмма входного сигнала (Signal LA on IF)
fig1 = figure;
plot(in_la, 'LineWidth', 1)
grid on; axis tight;
title('Signal LA on IF')
xlabel('Samples')
ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig1, 'P5_Signal_LA_on_IF.fig');
saveas(fig1, 'P5_Signal_LA_on_IF.tiff');

% График 2: Спектр входного сигнала (Spectrum LA on IF)
[fin_la, sin_la] = get_spectrum(in_la, fs, 1);
fig2 = figure;
plot(fin_la/1e6, mag2db(abs(sin_la)), 'LineWidth', 1)
grid on; axis tight;
title('Spectrum LA on IF')
xlabel('f, MHz')
ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig2, 'P5_Spectrum_LA_on_IF.fig');
saveas(fig2, 'P5_Spectrum_LA_on_IF.tiff');

% График 3: Временная диаграмма на выходе КИХ (Signal LA on IF = 0)
fig3 = figure;
plot(real(fir_sig_la), 'LineWidth', 1.5)
hold on
plot(imag(fir_sig_la), 'LineWidth', 1.5)
grid on; axis tight;
title('Signal LA on IF = 0')
xlabel('Samples')
ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig3, 'P5_Signal_LA_out.fig');
saveas(fig3, 'P5_Signal_LA_out.tiff');

% График 4: Спектр на выходе КИХ (Spectrum LA on IF = 0)
[ffir_la, sfir_la] = get_spectrum(fir_sig_la, fsv, 1);
fig4 = figure;
plot(ffir_la/1e6, mag2db(abs(sfir_la)), 'LineWidth', 1)
grid on; axis tight;
title('Spectrum LA on IF = 0')
xlabel('f, MHz')
ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig4, 'P5_Spectrum_LA_out.fig');
saveas(fig4, 'P5_Spectrum_LA_out.tiff');

%% ПУНКТ 6: Симметричный ЛЧМ сигнал (Финальная отлаженная версия)

clear sin; % Очистка переменной спектра от прошлых конфликтов

% 1. Расчет параметров ЛЧМ
B_lfm = 512; % База сигнала
tau_lfm = B_lfm / bandwidth; % Длительность импульса
prf_lfm = fs / length_test; % Период повторения (один импульс на все окно)

% 2. Генерация базового ЛЧМ-сигнала
lfm_waveform = phased.LinearFMWaveform( ...
    'SampleRate', fs, ...
    'PulseWidth', tau_lfm, ...
    'PRF', prf_lfm, ...
    'SweepBandwidth', bandwidth, ...
    'SweepDirection', 'Up', ...
    'SweepInterval', 'Symmetric', ... % Центрируем спектр вокруг нуля
    'OutputFormat', 'Samples', ...
    'NumSamples', length_test);

lfm_baseband = step(lfm_waveform); 

% Идеальное центрирование импульса во временном окне (сдвиг на половину пустого пространства)
delay_samples = round((length_test - (tau_lfm * fs)) / 2); 
lfm_baseband = circshift(lfm_baseband, delay_samples);

% 3. Перенос на несущую (fc) и добавление шума
max_amp_lfm = 2^(noba - 1) - 1; 
t_row = (0:length_test-1)/fs; % Гарантированный вектор-строка

% Формируем радиосигнал (при fc = fg = 170 МГц получим условие IF = 0)
in_lfm_clean = max_amp_lfm * real(lfm_baseband .* exp(1i * 2 * pi * fc * t_row.'));

noise_lfm = sqrt(noise_power) * normrnd(0, 1, [length_test 1]);
in_lfm = in_lfm_clean + noise_lfm;

% 4. Обработка в тракте
adc_signal_lfm = get_adc(in_lfm, noba);
[nco_sig_lfm, ~] = get_nco(nobg, fg, t_row, 'single', 0, 0);
dem_sig_lfm = get_dem(adc_signal_lfm, nco_sig_lfm, noba, nobg, 'single');

% 5. КИХ-фильтрация и децимация
len_to_use_lfm = R * floor(length(dem_sig_lfm)/R);
dem_sig_lfm_adj = dem_sig_lfm(1:len_to_use_lfm);
dem_sig_lfm_adj = dem_sig_lfm_adj(:); % Строгий вектор-столбец для КИХ

reset(FIRDecim); 
fir_sig_lfm = step(FIRDecim, double(dem_sig_lfm_adj));


% --- ВИЗУАЛИЗАЦИЯ (П.6) ---

% График 1: Signal LFM on IF
fig5 = figure;
plot(in_lfm, 'LineWidth', 1)
grid on; axis tight;
title('Signal LFM on IF')
xlabel('Samples')
ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig5, 'P6_Signal_LFM_on_IF.fig');
saveas(fig5, 'P6_Signal_LFM_on_IF.tiff');

% График 2: Spectrum LFM on IF
% Отключаем окно Ханна (0), чтобы видеть прямоугольный спектр
[fin_lfm, sin_lfm] = get_spectrum(in_lfm, fs, 0);
fig6 = figure;
plot(fin_lfm/1e6, mag2db(abs(sin_lfm)), 'LineWidth', 1.5)
grid on; axis tight;
title('Spectrum LFM on IF')
xlabel('f, MHz')
ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig6, 'P6_Spectrum_LFM_on_IF.fig');
saveas(fig6, 'P6_Spectrum_LFM_on_IF.tiff');

% График 3: Signal LFM on IF = 0
fig7 = figure;
plot(real(fir_sig_lfm), 'LineWidth', 1)
hold on
plot(imag(fir_sig_lfm), 'LineWidth', 1)
grid on; axis tight;
title('Signal LFM on IF = 0')
xlabel('Samples')
ylabel('LSB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig7, 'P6_Signal_LFM_video.fig');
saveas(fig7, 'P6_Signal_LFM_video.tiff');

% График 4: Spectrum LFM on IF = 0
% Отключаем окно Ханна (0), чтобы видеть прямоугольный спектр
[ffir_lfm, sfir_lfm] = get_spectrum(fir_sig_lfm, fsv, 0);
fig8 = figure;
plot(ffir_lfm/1e6, mag2db(abs(sfir_lfm)), 'LineWidth', 1.5)
grid on; axis tight;
title('Spectrum LFM on IF = 0')
xlabel('f, MHz')
ylabel('dB')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(fig8, 'P6_Spectrum_LFM_video.fig');
saveas(fig8, 'P6_Spectrum_LFM_video.tiff');

%% ДОПОЛНЕНИЕ К ПУНКТУ 6 (Видеосигнал с КИХ-фильтрацией и нормировкой)

% 1. Берем чистый ЛЧМ-сигнал (без шума), чтобы четко видеть зоны нулевой задержки
in_lfm_vid = in_lfm_clean;

% 2. Нормируем входной радиосигнал для графика (диапазон [-1, 1])
in_lfm_vid_norm = in_lfm_vid / max(abs(in_lfm_vid));

% 3. Пропускаем чистый сигнал через АЦП и перемножаем с гетеродином (демодуляция)
adc_vid = get_adc(in_lfm_vid, noba);
dem_vid = get_dem(adc_vid, nco_sig_lfm, noba, nobg, 'single');

% 4. Добавляем фильтрацию и децимацию (КИХ-фильтр)
len_to_use_vid = R * floor(length(dem_vid)/R);
dem_vid_adj = dem_vid(1:len_to_use_vid);
dem_vid_adj = dem_vid_adj(:); % Гарантируем вектор-столбец для КИХ

reset(FIRDecim); 
fir_vid = step(FIRDecim, double(dem_vid_adj)); % Пропускаем через фильтр

% 5. Нормируем I и Q каналы отфильтрованного видеосигнала
fir_vid_real = real(fir_vid) / max(abs(real(fir_vid)));
fir_vid_imag = imag(fir_vid) / max(abs(imag(fir_vid)));


% --- График 1: Временная диаграмма ЛЧМ на ПЧ (без шума) ---
figure;
plot(in_lfm_vid_norm, 'LineWidth', 1, 'Color', '#D95319')
grid on; axis tight;
title('Signal LFM on IF')
xlabel('Samples')
ylabel('LSB')
ylim([-1.2 1.2])
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(gcf, 'P6_Add_Signal_LFM_on_IF.tiff');

% --- График 2: Спектр ЛЧМ на ПЧ (без окна Ханна) ---
[fin_vid, sin_vid] = get_spectrum(in_lfm_vid_norm, fs, 0);
figure;
plot(fin_vid/1e6, mag2db(abs(sin_vid)), 'LineWidth', 1.5, 'Color', '#0072BD')
grid on; axis tight;
title('Spectrum LFM on IF')
xlabel('f, MHz')
ylabel('|A|')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(gcf, 'P6_Add_Spectrum_LFM_on_IF.tiff');

% --- График 3: Временная диаграмма видеосигнала (ПОСЛЕ КИХ-фильтра) ---
figure;
plot(fir_vid_real, 'LineWidth', 1, 'Color', '#0072BD') % Синий (I-канал)
hold on;
plot(fir_vid_imag, 'LineWidth', 1, 'Color', '#D95319') % Оранжевый (Q-канал)
grid on; axis tight;
title('Signal LFM video')
xlabel('Samples') % Количество отсчетов здесь будет в 10 раз меньше (децимация)
ylabel('LSB')
ylim([-1.2 1.2])
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(gcf, 'P6_Add_Signal_LFM_video.tiff');

% --- График 4: Спектр видеосигнала (ПОСЛЕ КИХ и децимации) ---
% ВАЖНО: используем fsv (14 МГц), так как сигнал уже децимирован
[fvid_out, svid_out] = get_spectrum(fir_vid, fsv, 0);
figure;
plot(fvid_out/1e6, mag2db(abs(svid_out)), 'LineWidth', 1.5, 'Color', '#0072BD')
grid on; axis tight;
title('Spectrum LFM video')
xlabel('f, MHz')
ylabel('|A|')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
saveas(gcf, 'P6_Add_Spectrum_LFM_video.tiff');