clear all; close all; clc;

length_test = 24e3;

%% Some properties
snr = 70; % signal to noise ratio, [dB]
noise_power = 1; % power of noise
signal_power = noise_power*db2pow(snr); % power of signal

fs = 140e6; % sampling frequency before decimation, [Hz] (Частота дискретизации АЦП)
fsv = 14e6; % sampling frequency after decimation, [Hz] (Сохраняем коэффициент децимации R=10)

fc = 170e6; % if mono = 0 fc = 170e6; (Промежуточная частота входного сигнала)
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
CIC = dsp.CICDecimator(R,1,N); % cic filter object
CIC.FixedPointDataType = 'Full precision';

% Gain calculation
cic_gain = R^N/2;

% Input generate (используем новую переменную src)
src = dsp.SignalSource(dem_signal, R*round(length(dem_signal)/R));    

% CIC filter output generate
cic_signal = step(CIC,step(src));    

% Normscale by gain
cic_signal = single(cic_signal)/single(cic_gain); 

% Plot results
figure
    plot(real(cic_signal),'LineWidth',2)
    hold on
    plot(imag(cic_signal),'LineWidth',2)
    grid on
    axis tight
    title('Signal after CIC')
    xlabel('t, samples')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

% Spectrum CIC
[fcic,scic] = get_spectrum(cic_signal,fsv,1);

% Plot results
figure
    plot(fcic/1e6,mag2db(abs(scic)),'LineWidth',2) 
    grid on
    axis tight
    title('Spectrum of signal after CIC')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

%% FIR Compensator
% FIR properties
fPass = bandwidth/2;
passbandRipple = 0.1;
fStop = bandwidth/2 + 500e3;
stopbandAttenuation = 60;

% Create object
FIR = dsp.CICCompensationDecimator(CIC, ...
    'DecimationFactor',1,...
    'PassbandFrequency',fPass,...
    'PassbandRipple',passbandRipple,...
    'StopbandAttenuation',stopbandAttenuation,...
    'StopbandFrequency',fStop,...
    'SampleRate',fsv,...
    'DesignForMinimumOrder',true);

% Impulse responce generate
hFIR = impz(FIR);
hFIR = hFIR./max(hFIR);

% Generate signal
fir_signal = filter(hFIR,1,double(cic_signal)); 

% Plot results
figure
    plot(real(fir_signal),'LineWidth',2)
    hold on
    plot(imag(fir_signal),'LineWidth',2)
    grid on
    axis tight
    title('Signal after FIR')
    xlabel('t, samples')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

% Spectrum FIR
[ffir,sfir] = get_spectrum(fir_signal,fsv,1);

% Plot results
figure
    plot(ffir/1e6,mag2db(abs(sfir)),'LineWidth',2) 
    grid on
    axis tight
    title('Spectrum of signal after FIR')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')    

% Calculate cascade filter impulse responce
FC = dsp.FilterCascade(CIC,FIR);

% Create and plot object
RESP = fvtool(CIC, FIR, FC,'Fs',[fs fsv fs]);
RESP.NormalizeMagnitudeto1 = 'on';

%% Пункт 3: Исследование зависимости побочных составляющих от разрядности ЦГ
% Для варианта N=4: na = 9, следовательно ng варьируется от 2 до 11
nobg_array = 2:(noba + 2); 
sfdr_values = zeros(size(nobg_array));

nco_type_p3 = 'fixed'; % Согласно заданию используем fixed

for i = 1:length(nobg_array)
    current_nobg = nobg_array(i);
    
    % 1. Генерируем сигнал гетеродина с текущей разрядностью
    [nco_signal_p3, nco_gain_p3] = get_nco(current_nobg, fg, t, nco_type_p3, 0, 0);
    
    % 2. Демодуляция (используем уже имеющийся adc_signal)
    dem_signal_p3 = get_dem(adc_signal, nco_signal_p3, noba, current_nobg, nco_type_p3);
    
    % 3. Фильтрация (пропускаем через CIC и КИХ-компенсатор)
    src_p3 = dsp.SignalSource(dem_signal_p3, R*round(length(dem_signal_p3)/R));
    
    release(CIC); % Разблокируем объект перед сменой типа данных
    cic_signal_p3 = step(CIC, step(src_p3));
    cic_signal_p3 = single(cic_signal_p3)/single(cic_gain); 
    
    % Сигнал на выходе устройства (после КИХ)
    fir_signal_p3 = filter(hFIR, 1, double(cic_signal_p3));
    
    % 4. Измерение SFDR
    signal_for_sfdr = real(fir_signal_p3);
    sfdr_values(i) = sfdr(signal_for_sfdr, fsv); 
    
    % Выводим и безопасно сохраняем показательные спектры для 2, 6 и 11 бит
    if current_nobg == 2 || current_nobg == 6 || current_nobg == 11
        fig_spectrum = figure;
        sfdr(signal_for_sfdr, fsv); 
        title(sprintf('Спектр на выходе (ЦГ: %d бит)', current_nobg));
        set(gca,'Fontsize',18,'Fontname','Times New Roman');
        
        filename_base = sprintf('p3_spectrum_%d_bit', current_nobg);
        filename_fig = [filename_base, '.fig'];
        filename_tif = [filename_base, '.tif'];
        
        if ~isfile(filename_fig)
            savefig(fig_spectrum, filename_fig);
        else
            disp(['Файл ', filename_fig, ' уже существует. Пропуск...']);
        end
        
        if ~isfile(filename_tif)
            print(fig_spectrum, filename_base, '-dtiff', '-r300');
        else
            disp(['Файл ', filename_tif, ' уже существует. Пропуск...']);
        end
    end
end

% 5. Построение итогового графика зависимости SFDR от разрядности ЦГ
fig_sfdr = figure; 
plot(nobg_array, sfdr_values, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
title('Зависимость SFDR от количества бит ЦГ');
xlabel('Количество бит ЦГ, n_g');
ylabel('SFDR, дБн');
set(gca,'Fontsize',24,'Fontname','Times New Roman');

% Безопасное сохранение итогового графика
if ~isfile('p3_sfdr_vs_bits.fig')
    savefig(fig_sfdr, 'p3_sfdr_vs_bits.fig');
else
    disp('Файл p3_sfdr_vs_bits.fig уже существует. Пропуск...');
end

if ~isfile('p3_sfdr_vs_bits.tif')
    print(fig_sfdr, 'p3_sfdr_vs_bits', '-dtiff', '-r300');
else
    disp('Файл p3_sfdr_vs_bits.tif уже существует. Пропуск...');
end

%% Пункт 4: Исследование квадратурных искажений (с проверкой наличия файлов)
% Возвращаем разрядность ЦГ к исходному значению для N=4
current_nobg = 9; 
nco_type_p4 = 'single'; 

% Задаем диапазоны изменения искажений
amv_array = 0:0.02:0.2; % Амплитудный дисбаланс (от 0 до 0.2)
pmv_array = 0:2:20;     % Фазовый дисбаланс (от 0 до 20)

sfdr_amv = zeros(size(amv_array));
sfdr_pmv = zeros(size(pmv_array));

%% 4.1 Исследование амплитудных искажений (amv) при pmv = 0
for i = 1:length(amv_array)
    current_amv = amv_array(i);
    
    % Генерируем гетеродин с амплитудным перекосом квадратур
    [nco_signal_p4, ~] = get_nco(current_nobg, fg, t, nco_type_p4, current_amv, 0);
    dem_signal_p4 = get_dem(adc_signal, nco_signal_p4, noba, current_nobg, nco_type_p4);
    
    % Фильтрация
    src_p4 = dsp.SignalSource(dem_signal_p4, R*round(length(dem_signal_p4)/R));
    release(CIC);
    cic_signal_p4 = step(CIC, step(src_p4));
    cic_signal_p4 = single(cic_signal_p4)/single(cic_gain); 
    fir_signal_p4 = filter(hFIR, 1, double(cic_signal_p4));
    
    % Измерение SFDR
    signal_for_sfdr = real(fir_signal_p4);
    sfdr_amv(i) = sfdr(signal_for_sfdr, fsv);
    
    % Построение и безопасное сохранение показательного спектра
    if abs(current_amv - 0.1) < 1e-5
        fig_amv_spec = figure;
        sfdr(signal_for_sfdr, fsv);
        title(sprintf('Спектр при амплитудном искажении amv = %.2f', current_amv));
        set(gca,'Fontsize',18,'Fontname','Times New Roman');
        
        if ~isfile('p4_spectrum_amv_0.1.fig')
            savefig(fig_amv_spec, 'p4_spectrum_amv_0.1.fig');
        else
            disp('Файл p4_spectrum_amv_0.1.fig уже существует. Пропуск...');
        end
        
        if ~isfile('p4_spectrum_amv_0.1.tif')
            print(fig_amv_spec, 'p4_spectrum_amv_0.1', '-dtiff', '-r300');
        else
            disp('Файл p4_spectrum_amv_0.1.tif уже существует. Пропуск...');
        end
    end
end

%% 4.2 Исследование фазовых искажений (pmv) при amv = 0
for i = 1:length(pmv_array)
    current_pmv = pmv_array(i);
    
    % Генерируем гетеродин с фазовым перекосом квадратур
    [nco_signal_p4, ~] = get_nco(current_nobg, fg, t, nco_type_p4, 0, current_pmv);
    dem_signal_p4 = get_dem(adc_signal, nco_signal_p4, noba, current_nobg, nco_type_p4);
    
    % Фильтрация
    src_p4 = dsp.SignalSource(dem_signal_p4, R*round(length(dem_signal_p4)/R));
    release(CIC);
    cic_signal_p4 = step(CIC, step(src_p4));
    cic_signal_p4 = single(cic_signal_p4)/single(cic_gain); 
    fir_signal_p4 = filter(hFIR, 1, double(cic_signal_p4));
    
    % Измерение SFDR
    signal_for_sfdr = real(fir_signal_p4);
    sfdr_pmv(i) = sfdr(signal_for_sfdr, fsv);
    
    % Построение и безопасное сохранение показательного спектра
    if abs(current_pmv - 10) < 1e-5
        fig_pmv_spec = figure;
        sfdr(signal_for_sfdr, fsv);
        title(sprintf('Спектр при фазовом искажении pmv = %d', current_pmv));
        set(gca,'Fontsize',18,'Fontname','Times New Roman');
        
        if ~isfile('p4_spectrum_pmv_10.fig')
            savefig(fig_pmv_spec, 'p4_spectrum_pmv_10.fig');
        else
            disp('Файл p4_spectrum_pmv_10.fig уже существует. Пропуск...');
        end
        
        if ~isfile('p4_spectrum_pmv_10.tif')
            print(fig_pmv_spec, 'p4_spectrum_pmv_10', '-dtiff', '-r300');
        else
            disp('Файл p4_spectrum_pmv_10.tif уже существует. Пропуск...');
        end
    end
end

%% 4.3 Построение итоговых графиков с безопасным сохранением
fig_amv = figure;
plot(amv_array, sfdr_amv, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
title('Зависимость SFDR от амплитудных искажений (amv)');
xlabel('Амплитудный дисбаланс, amv');
ylabel('SFDR, дБн');
set(gca,'Fontsize',24,'Fontname','Times New Roman');

if ~isfile('p4_sfdr_vs_amv.fig')
    savefig(fig_amv, 'p4_sfdr_vs_amv.fig');
end
if ~isfile('p4_sfdr_vs_amv.tif')
    print(fig_amv, 'p4_sfdr_vs_amv', '-dtiff', '-r300');
end

fig_pmv = figure;
plot(pmv_array, sfdr_pmv, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
grid on;
title('Зависимость SFDR от фазовых искажений (pmv)');
xlabel('Фазовый дисбаланс, pmv');
ylabel('SFDR, дБн');
set(gca,'Fontsize',24,'Fontname','Times New Roman');

if ~isfile('p4_sfdr_vs_pmv.fig')
    savefig(fig_pmv, 'p4_sfdr_vs_pmv.fig');
end
if ~isfile('p4_sfdr_vs_pmv.tif')
    print(fig_pmv, 'p4_sfdr_vs_pmv', '-dtiff', '-r300');
end

%% Пункт 5: Тестовое воздействие с линейно нарастающей амплитудой (LA)
% 1. Задаем параметры сигнала
f_offset = 1e6; % Смещение на 1 МГц, чтобы сигнал на выходе был отличен от нуля
fc_la = fc + f_offset; % Новая частота входного сигнала (171 МГц)

% Формируем линейно нарастающую огибающую в пределах динамического диапазона АЦП
max_amp = 2^(noba-1) - 1; % Для 9 бит это 255
envelope = linspace(0, max_amp, length_test).'; 

clear sin; % <--- ДОБАВЬ ЭТУ СТРОКУ, чтобы MATLAB снова понял, что sin - это функция

% Генерируем тестовый сигнал и добавляем шум
signal_la = envelope .* sin(2*pi*fc_la*t).';
in_la = signal_la + noise;

% 2. Пропускаем через АЦП
adc_signal_la = get_adc(in_la, noba);

% Получаем спектр после АЦП для первого графика
[f_adc_la, s_adc_la] = get_spectrum(adc_signal_la, fs, 1);

% 3. Цифровой гетеродин и демодуляция (используем стандартные параметры без искажений)
[nco_signal_la, nco_gain_la] = get_nco(nobg, fg, t, 'fixed', 0, 0);
dem_signal_la = get_dem(adc_signal_la, nco_signal_la, noba, nobg, 'fixed');

% 4. Фильтрация
src_la = dsp.SignalSource(dem_signal_la, R*round(length(dem_signal_la)/R));
release(CIC); % Не забываем разблокировать объект
cic_signal_la = step(CIC, step(src_la));
cic_signal_la = single(cic_signal_la)/single(cic_gain); 
fir_signal_la = filter(hFIR, 1, double(cic_signal_la));

% Получаем спектр на выходе фильтра
[f_fir_la, s_fir_la] = get_spectrum(fir_signal_la, fsv, 1);

% 5. Построение графиков (2x2 как в методичке)
fig_p5 = figure('Position', [100, 100, 1200, 800]);

% График 1: Сигнал до фильтрации (после АЦП)
subplot(2, 2, 1);
plot(adc_signal_la, 'LineWidth', 1); 
axis tight;
grid on;
title('Signal LA on IF');
xlabel('Samples');
ylabel('LSB');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% График 2: Сигнал на выходе устройства
subplot(2, 2, 2);
plot(real(fir_signal_la), 'LineWidth', 1.5);
hold on;
plot(imag(fir_signal_la), 'LineWidth', 1.5);
axis tight;
grid on;
title('Signal LA on IF = 0');
xlabel('Samples');
ylabel('LSB');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% График 3: Спектр до фильтрации
subplot(2, 2, 3);
plot(f_adc_la/1e6, mag2db(abs(s_adc_la)), 'LineWidth', 1.5);
axis tight;
grid on;
title('Spectrum LA on IF');
xlabel('f, MHz');
ylabel('|A|');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% График 4: Спектр на выходе устройства
subplot(2, 2, 4);
plot(f_fir_la/1e6, mag2db(abs(s_fir_la)), 'LineWidth', 1.5);
axis tight;
grid on;
title('Spectrum LA on IF = 0');
xlabel('f, MHz');
ylabel('|A|');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% 6. Безопасное сохранение
if ~isfile('p5_linear_amplitude.fig')
    savefig(fig_p5, 'p5_linear_amplitude.fig');
else
    disp('Файл p5_linear_amplitude.fig уже существует. Пропуск...');
end

if ~isfile('p5_linear_amplitude.tif')
    print(fig_p5, 'p5_linear_amplitude', '-dtiff', '-r300');
else
    disp('Файл p5_linear_amplitude.tif уже существует. Пропуск...');
end

%% Пункт 6: Симметричный ЛЧМ сигнал с прямоугольной огибающей
% 1. Расчет параметров ЛЧМ согласно Таблице 2
B_lfm = 512;
delta_f_lfm = bandwidth; % Согласовано со значением из скрипта (около 5.83 МГц)
T_pulse_lfm = B_lfm / delta_f_lfm;
T_rep_lfm = length_test / fs; % Период повторения исходя из количества отсчетов
prf_lfm = 1 / T_rep_lfm;      % Частота повторения импульсов

% 2. Генерация ЛЧМ видеосигнала
hwav = phased.LinearFMWaveform('SampleRate', fs, ...
    'PulseWidth', T_pulse_lfm, ...
    'PRF', prf_lfm, ...
    'SweepBandwidth', delta_f_lfm, ...
    'SweepInterval', 'Symmetric', ...
    'OutputFormat', 'Samples', ...
    'NumSamples', length_test);

lfm_base = hwav();

% Центрируем импульс для наглядности (как на эталонном графике)
pulse_samples = round(T_pulse_lfm * fs);
shift_amount = round((length_test - pulse_samples) / 2);
lfm_base = circshift(lfm_base, shift_amount);

% 3. Перенос на промежуточную частоту fc и добавление шума
% max_amp уже рассчитан в 5 пункте (для 9 бит = 255)
signal_lfm = max_amp * real(lfm_base .* exp(1j*2*pi*fc*t(:)));
in_lfm = signal_lfm + noise;

% 4. Пропускаем через АЦП
adc_signal_lfm = get_adc(in_lfm, noba);
[f_adc_lfm, s_adc_lfm] = get_spectrum(adc_signal_lfm, fs, 1);

% 5. Цифровой гетеродин и демодуляция
[nco_signal_lfm, nco_gain_lfm] = get_nco(nobg, fg, t, 'fixed', 0, 0);
dem_signal_lfm = get_dem(adc_signal_lfm, nco_signal_lfm, noba, nobg, 'fixed');

% 6. Фильтрация
src_lfm = dsp.SignalSource(dem_signal_lfm, R*round(length(dem_signal_lfm)/R));
release(CIC); % Сброс состояния фильтра
cic_signal_lfm = step(CIC, step(src_lfm));
cic_signal_lfm = single(cic_signal_lfm)/single(cic_gain); 
fir_signal_lfm = filter(hFIR, 1, double(cic_signal_lfm));

[f_fir_lfm, s_fir_lfm] = get_spectrum(fir_signal_lfm, fsv, 1);

% 7. Построение графиков (2x2 как в методичке)
fig_p6 = figure('Position', [150, 150, 1200, 800]);

% График 1: Сигнал до фильтрации
subplot(2, 2, 1);
plot(adc_signal_lfm, 'LineWidth', 1); 
axis tight;
grid on;
title('Signal LFM on IF');
xlabel('Samples');
ylabel('LSB');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% График 2: Сигнал на выходе устройства
subplot(2, 2, 2);
plot(real(fir_signal_lfm), 'LineWidth', 1.5);
hold on;
plot(imag(fir_signal_lfm), 'LineWidth', 1.5);
axis tight;
grid on;
title('Signal LFM on IF = 0');
xlabel('Samples');
ylabel('LSB');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% График 3: Спектр до фильтрации
subplot(2, 2, 3);
plot(f_adc_lfm/1e6, mag2db(abs(s_adc_lfm)), 'LineWidth', 1.5);
axis tight;
grid on;
title('Spectrum LFM on IF');
xlabel('f, MHz');
ylabel('|A|');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% График 4: Спектр на выходе устройства (в полосе пропускания)
subplot(2, 2, 4);
plot(f_fir_lfm/1e6, mag2db(abs(s_fir_lfm)), 'LineWidth', 1.5);
axis tight;
grid on;
title('Spectrum LFM on IF = 0');
xlabel('f, MHz');
ylabel('|A|');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% 8. Безопасное сохранение
if ~isfile('p6_lfm_signal.fig')
    savefig(fig_p6, 'p6_lfm_signal.fig');
else
    disp('Файл p6_lfm_signal.fig уже существует. Пропуск...');
end

if ~isfile('p6_lfm_signal.tif')
    print(fig_p6, 'p6_lfm_signal', '-dtiff', '-r300');
else
    disp('Файл p6_lfm_signal.tif уже существует. Пропуск...');
end

%% Дополнение к пункту 6: Огибающая видеосигнала ЛЧМ
% В качестве видеосигнала выступает базовая комплексная огибающая lfm_base
% Сформируем нормированный ЛЧМ-сигнал на ПЧ (с амплитудой от -1 до 1, как на скриншоте)
signal_lfm_norm = real(lfm_base .* exp(1j*2*pi*fc*t(:)));

% Вычисляем спектры (передаем fs, т.к. дискретизация исходная)
[f_video, s_video] = get_spectrum(lfm_base, fs, 1);
[f_if_norm, s_if_norm] = get_spectrum(signal_lfm_norm, fs, 1);

% Создаем окно графиков
fig_p6_video = figure('Position', [200, 200, 1200, 800]);

% 1. Временная область: Signal LFM video
subplot(2, 2, 1);
plot(real(lfm_base), 'LineWidth', 1);
hold on;
plot(imag(lfm_base), 'LineWidth', 1);
axis tight;
grid on;
title('Signal LFM video');
xlabel('Samples');
ylabel('LSB'); 
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% 2. Временная область: Signal LFM on IF
subplot(2, 2, 2);
% Рисуем оранжевым цветом, чтобы соответствовать визуалу твоего примера
plot(signal_lfm_norm, 'LineWidth', 1, 'Color', '#D95319'); 
axis tight;
grid on;
title('Signal LFM on IF');
xlabel('Samples');
ylabel('LSB');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% 3. Спектр: Spectrum LFM video
subplot(2, 2, 3);
plot(f_video/1e6, mag2db(abs(s_video)), 'LineWidth', 1.5);
axis tight;
grid on;
title('Spectrum LFM video');
xlabel('f, MHz');
ylabel('|A|');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% 4. Спектр: Spectrum LFM on IF
subplot(2, 2, 4);
plot(f_if_norm/1e6, mag2db(abs(s_if_norm)), 'LineWidth', 1.5);
axis tight;
grid on;
title('Spectrum LFM on IF');
xlabel('f, MHz');
ylabel('|A|');
set(gca, 'Fontsize', 14, 'Fontname', 'Times New Roman');

% Безопасное сохранение
if ~isfile('p6_lfm_video_supplement.fig')
    savefig(fig_p6_video, 'p6_lfm_video_supplement.fig');
else
    disp('Файл p6_lfm_video_supplement.fig уже существует. Пропуск...');
end

if ~isfile('p6_lfm_video_supplement.tif')
    print(fig_p6_video, 'p6_lfm_video_supplement', '-dtiff', '-r300');
else
    disp('Файл p6_lfm_video_supplement.tif уже существует. Пропуск...');
end