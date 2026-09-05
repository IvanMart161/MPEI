clear all; close all; clc;

length_test = 100e3;

%% Some properties
snr = 75; % signal to noise ratio, [dB]
noise_power = 1; % power of noise
signal_power = noise_power*db2pow(snr); % power of signal

fs = 210e6; % sampling frequency before decimation, [Hz]
fc = 180.6e6; % if mono = 0 fc = 312.6e6; 
fg = 210.6e6; % center demodulation geterodin freq, [Hz]

na = 14; % number of bits after analog to digital conversion
ng = 14; % number of bits for digital geterodin

weighting = 1; % on or off weigth for spectrum

%% Some code here
t = (0:length_test - 1)/fs; % time vector, s
noise = sqrt(noise_power).*(normrnd(0,1,[length_test 1])); % noise generate
signal = sqrt(signal_power).*sin(2*pi*fc*t).';
in = signal + noise; % input generate

%% INPUT
% Plot results
figure(1)
    plot(t*1e3,in,'LineWidth',2)
    grid on
    axis tight
    title('Signal before ADC')
    xlabel('t, ms')
    ylabel('LSB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')


% Spectrum INPUT
[fin, sin_spec] = get_spectrum(in, fs, weighting);

% Plot results
figure(2)
    plot(fin/1e6,mag2db(abs(sin_spec)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal before ADC')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

%% ADC
adc_signal = get_adc(in,na);

% Plot results
figure(3)
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
figure(4)
    plot(fadc/1e6,mag2db(abs(sadc)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal after ADC')
    xlabel('f, MHz')
    ylabel('dB')
    set(gca,'Fontsize',28,'Fontname','Times New Roman')

%% NCO
nco_type = 'single';
amv = 0; % amplitude mismatch IQ, dB
pmv = 2; % phase mismatch IQ, deg
[nco_signal,nco_gain] = get_nco(ng,fg,t,nco_type,amv,pmv);

% Plot results
figure(5)
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
figure(6)
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
figure(7)
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

figure(8)
    plot(fdem/1e6,mag2db(abs(sdem)),'LineWidth',2)
    grid on
    axis tight
    title('Spectrum of signal after DEM')
    xlabel('f, MHz')
    ylabel('dB')    
    set(gca,'Fontsize',28,'Fontname','Times New Roman')  

function out = get_adc(in,bits)
    out = quantize(int32(in),1,bits,0,'Floor','Saturate');
end

function [nco_signal,nco_gain] = get_nco(nobg,fg,t,nco_type,amv,pmv)

    % Get signal nco
    nco_real = cos(2*pi*fg*t).';
    % Get Phase mismatch IQ
    nco_imag = sin(2*pi*fg*t + deg2rad(pmv)).';
    % Get am[lide mismatch IQ
    nco_real = max(real(nco_real))*db2mag(amv)*nco_real;
    % Calulate nco
    nco = nco_real - 1i*nco_imag;
    
    % Calcilate gain nco
    nco_gain = 2^(nobg - 1) - 1;

    % Fromat conversion
    switch nco_type
        case 'fixed'
            nco = nco_gain.*nco;
            nco_signal = sfi(nco,nobg,0);
        case 'single'
            nco_signal = single(nco);
        case 'double'
            nco_signal = double(nco);
    end
end

function dem_signal = get_dem(adc_signal,nco_signal,noba,nobg,dem_type)
    switch dem_type
        case 'fixed'
            dem_signal = adc_signal.*nco_signal;
            dem_signal = sfi(dem_signal,noba + nobg,0);
        case 'single'
            dem_signal = single(adc_signal).*single(nco_signal);
            dem_signal = sfi(dem_signal,noba + nobg,0);
        case 'double'
            dem_signal = double(adc_signal).*double(nco_signal);
            dem_signal = sfi(dem_signal,noba + nobg,0);  
    end
end

function [f,s] = get_spectrum(in,fs,weighting)
    f = linspace(-fs/2,fs/2,length(in));
    if (weighting)
        s = fftshift(fft(double(in).*hann(length(in))));
    else
        s = fftshift(fft(double(in)));
    end
end

