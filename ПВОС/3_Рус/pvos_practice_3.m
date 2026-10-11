clear all; close all; clc;

length_test = 2^14; % lenght of test 2^N for matched filtering

%% Some properties
fsv = 40e6; % sampling frequency after decimation, [Hz]

pulse_base = 512; % pulse base
pulse_bw = 20e6; % pulse bandwith, [Hz]
pulse_width = pulse_base/pulse_bw; % pulse width, [s]
prf = 1/pulse_width; % pulse repetition frequency, [Hz]

delay_samples = length_test/2; % samples for position of echo 

noise_power = 1; % power of noise
SNR = 50; % signal to noise ratio, [dB]
signal_power = noise_power*10.^(SNR/10); % power of signal

t = (0:length_test - 1)/fsv;

% Generate waveform 
waveform = phased.LinearFMWaveform(...
    'SampleRate',fsv,...
    'PulseWidth',pulse_width,...
    'PRF',prf,...
    'SweepBandwidth',pulse_bw,...
    'SweepInterval','Symmetric',...
    'Envelope','Rectangular');

waveform = step(waveform); % waveform generate at video

% Generate complex signal
signal = sqrt(signal_power)*[zeros(delay_samples,1);...
    waveform;...
    zeros(length_test - length(waveform) - delay_samples,1)];

% Generate complex noise
noise = sqrt(noise_power).*(normrnd(0,1,[length_test 1]))...
    + 1i*sqrt(noise_power).*(normrnd(0,1,[length_test 1])); % noise generate

% Generate input data for matched filtering
in = signal + noise;

% Plot input data
figure
    plot(real(in),'LineWidth',2)
    hold on
    plot(imag(in),'LineWidth',2)
    grid on
    axis tight
    
%% First method of matched filtering
% Use FIR
% Without windowing
replica_time = flip(conj(waveform));
mf_signal_one = filter(replica_time,1,in);

% With windowing
window_time = hamming(length(replica_time));
replica_time_w = replica_time.*window_time;
mf_signal_one_w = filter(replica_time_w,1,in);

figure
    plot(mag2db(abs(mf_signal_one)),'LineWidth',2)
    hold on
    plot(mag2db(abs(mf_signal_one_w)),'LineWidth',2)
    grid on
    axis tight

%% Second method of matched filtering
% Use FFT and IFFT
% Without windowing
in_fft = fft(in,length_test); % get FFT

replica_fft = fft(replica_time,length_test); % get FFT
replica_fft_w = fft(replica_time_w,length_test); % get FFT

mult_fft = in_fft.*replica_fft;
mult_fft_w = in_fft.*replica_fft_w;

mf_signal_two = ifft(mult_fft,length_test); 
mf_signal_two_w = ifft(mult_fft_w,length_test);

figure
    plot(mag2db(abs(mf_signal_two)),'LineWidth',2)
    hold on
    plot(mag2db(abs(mf_signal_two_w)),'LineWidth',2)
       
%% Third method of matched filtering
% Use overlap_save or overlap-add method
% Without windowing
mf_signal_three = fftfilt(replica_time,in);

% With windowing
mf_signal_three_w = fftfilt(replica_time_w,in);

% Plot results
figure
    plot(mag2db(abs(mf_signal_three)),'LineWidth',2)
    hold on
    plot(mag2db(abs(mf_signal_three_w)),'LineWidth',2)
    axis tight
    grid on