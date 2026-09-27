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