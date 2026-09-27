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

% Create and plot object
RESP = fvtool(CIC, FIR, FC,'Fs',[fs fsv fs]);

RESP.NormalizeMagnitudeto1 = 'on';   
    