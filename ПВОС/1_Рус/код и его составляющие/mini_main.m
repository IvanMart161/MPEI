clear all; close all; clc;

length_test = 100e3;

%% Some properties
snr = 75; % signal to noise ratio, [dB]
noise_power = 1; % power of noise
signal_power = noise_power*db2pow(snr); % power of signal

N = 4; % номер в журнале

% Индивидуальные параметры из Таблицы 1
fs = (100 + N * 10) * 1e6;       % Частота дискретизации АЦП (f_S) (sampling frequency before decimation, [Hz])
fc = (130 + N * 10) * 1e6 + 0.1e6; % Промежуточная частота (f_C) + отстройка 0,1 МГц (if mono = 0 fc = 312.6e6)
fg = (130 + N * 10) * 1e6;       % Центральная частота ЦГ (f_G) (center demodulation geterodin freq, [Hz])

na = 14; % number of bits after analog to digital conversion
ng = 8; % number of bits for digital geterodin

weighting = 1; % on or off weigth for spectrum

%% Some code here
t = (0:length_test - 1)/fs; % time vector, s
noise = sqrt(noise_power).*(normrnd(0,1,[length_test 1])); % noise generate
signal = sqrt(signal_power).*sin(2*pi*fc*t).';
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
[fin,sin] = get_spectrum(in,fs,weighting);

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
adc_signal = get_adc(in,na);

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
[fadc,sadc] = get_spectrum(adc_signal,fs,weighting);

% Plot results
figure
    plot(fadc/1e6,mag2db(abs(sadc)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal after ADC')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')
  

%% 3. Исследование зависимости от количества бит ЦГ
nco_type = 'fixed'; % Требование задания для пункта 3
amv = 0; % Амплитудные искажения отключены
pmv = 0; % Фазовые искажения отключены
dem_type = 'single'; % Тип демодулятора из исходного кода

% Диапазон изменения бит гетеродина от 2 до na + 2
ng_range = 2:(na + 2); 

% Создаем окно графиков до начала цикла
fig_nco = figure('Name', 'Зависимость спектра DEM от ng');

for i = 1:length(ng_range)
    ng_current = ng_range(i);
    
    % 1. Формируем сигнал ЦГ с текущим количеством бит
    [nco_signal, nco_gain] = get_nco(ng_current, fg, t, nco_type, amv, pmv);
    
    % 2. Проводим демодуляцию
    dem_signal = get_dem(adc_signal, nco_signal, na, ng_current, dem_type);
    
    % 3. Рассчитываем спектр
    [fdem, sdem] = get_spectrum(dem_signal, fs, weighting);
    
    % 4. Отрисовка графика
    figure(fig_nco); % Переключаем фокус на наше окно
    plot(fdem/1e6, mag2db(abs(sdem)), 'LineWidth', 2);
    grid on;
    axis tight;
    
    % Динамический заголовок для удобства отслеживания
    title(sprintf('Спектр после DEM: ng = %d бит', ng_current));
    xlabel('f, МГц');
    ylabel('дБ');
    set(gca, 'Fontsize', 24, 'Fontname', 'Times New Roman');
    
    % Остановка цикла для сохранения графика
    fprintf('Текущее значение ng = %d. Сохраните график (если нужно) и нажмите любую клавишу...\n', ng_current);
    pause; 
end
fprintf('Цикл исследования для пункта 3 завершен.\n');

%% 4. Имитационное моделирование: контурный график искажений (как в лекции)
nco_type = 'single'; % Требование задания
dem_type = 'single'; % Тип демодулятора

% Задаем двумерную сетку значений (20x20 точек для оптимальной скорости расчета)
pmv_vec = linspace(0, 4, 20);  % Фаза от 0 до 4 градусов
amv_vec = linspace(0, 0.6, 20); % Амплитуда от 0 до 0.6 дБ

[PMV, AMV] = meshgrid(pmv_vec, amv_vec);
P_rel = zeros(size(PMV)); % Матрица для сохранения результатов (дБ)

h = waitbar(0, 'Моделирование сетки искажений. Пожалуйста, подождите...');
total_iters = numel(PMV);

for i = 1:total_iters
    amv_current = AMV(i);
    pmv_current = PMV(i);
    
    % 1. Формируем сигнал ЦГ
    [nco_signal, ~] = get_nco(ng, fg, t, nco_type, amv_current, pmv_current);
    
    % 2. Проводим демодуляцию
    dem_signal = get_dem(adc_signal, nco_signal, na, ng, dem_type);
    
    % 3. Рассчитываем спектр
    [fdem, sdem] = get_spectrum(dem_signal, fs, weighting);
    sdem_db = mag2db(abs(sdem));
    
    % 4. Умный поиск пиков в спектре
    % Исключаем нулевую частоту (DC) из поиска
    valid_idx = find(abs(fdem) > 0.05e6);
    
    % Находим глобальный максимум - это наш полезный сигнал
    [pwr_main, local_max_idx] = max(sdem_db(valid_idx));
    f_main = fdem(valid_idx(local_max_idx));
    
    % Зеркальный канал возникает на противоположной частоте (-f_main)
    f_image_target = -f_main;
    search_range = 0.05e6; % Окно поиска 50 кГц
    
    % Ищем мощность зеркального канала
    idx_image = find(fdem > (f_image_target - search_range) & fdem < (f_image_target + search_range));
    pwr_image = max(sdem_db(idx_image));
    
    % 5. Сохраняем разницу мощностей
    P_rel(i) = pwr_image - pwr_main;
    
    % Обновление прогресс-бара
    if mod(i, 20) == 0
        waitbar(i/total_iters, h);
    end
end
close(h);

% 6. Отрисовка контурного графика
figure('Name', 'Моделирование: Относительная мощность зеркального канала');
levels = [-50, -45, -40, -35, -30]; % Требуемые уровни изолиний
[C, h_contour] = contour(PMV, AMV, P_rel, levels, 'k', 'LineWidth', 1.5);
clabel(C, h_contour, 'FontSize', 12, 'LabelSpacing', 200);

grid on;
title('Относительная мощность зеркального канала (Моделирование)');
xlabel('Phase mismatch (degrees)');
ylabel('Amplitude mismatch (dB)');
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman');
%% NCO
nco_type = 'single';
amv = 0; % amplitude mismatch IQ, dB
pmv = 0; % phase mismatch IQ, deg
[nco_signal,nco_gain] = get_nco(ng,fg,t,nco_type,amv,pmv);

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
[fnco,snco] = get_spectrum(nco_signal,fs,weighting);

% Plot results
figure
    plot(fnco/1e6,mag2db(abs(snco)),'LineWidth',2) 
    grid on
    axis tight
    title('Spectrum of signal after NCO')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

%% DEM
dem_type = 'single';
dem_signal = get_dem(adc_signal,nco_signal,na,ng,dem_type);

% Plot results
figure
    plot(real(dem_signal),'LineWidth',2)
    hold on
    plot(imag(dem_signal),'LineWidth',2)
    grid on
    axis tight
    title('Signal after DEM')
    xlabel('t, samples')
    ylabel('LSB')
    legend('I','Q')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')
    
% Spectrum DEM
[fdem,sdem] = get_spectrum(dem_signal,fs,weighting);

figure
    plot(fdem/1e6,mag2db(abs(sdem)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal after DEM')
    xlabel('f, MHz')
    ylabel('dB')    
    set(gca,'Fontsize',28,'Fontname','Times New Roman')  