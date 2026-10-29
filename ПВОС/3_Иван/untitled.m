clear all; close all; clc;
rng(0);

%% ===================== ПАРАМЕТРЫ ВАРИАНТА 8 =====================
N_PZh = 8;

fc      = (130 + N_PZh*10)*1e6;   % 210 МГц
n_a     = 5 + N_PZh;              % 13 бит АЦП
fs      = (100 + N_PZh*10)*1e6;   % 180 МГц
n_g     = 5 + N_PZh;              % 13 бит ЦГ
f_G     = (130 + N_PZh*10)*1e6;   % 210 МГц

B_base  = 512;                    % база сигнала
Delta_f = 20e6;                   % полоса ЛЧМ (20 МГц)
pulse_width = B_base/Delta_f;     % 25.6 мкс

length_test   = 2^14;             % 16384 отсчета
PRI           = length_test/fs;   % период повторения
prf           = 1/PRI;
delay_samples = length_test/2;    % задержка (8192)

noise_power  = 1;
SNR_dB       = 50;
signal_power = noise_power*10^(SNR_dB/10);

fprintf('Вариант 8:\n');
fprintf('fc=%.1f MHz, fs=%.1f MHz, f_G=%.1f MHz\n', fc/1e6, fs/1e6, f_G/1e6);
fprintf('n_a=%d, n_g=%d, pulse_width=%.3f us, N_imp=%d отсчетов\n', ...
        n_a, n_g, pulse_width*1e6, round(pulse_width*fs));

N_lfm = round(pulse_width * fs); % 4608 отсчетов

t_lfm = linspace(-pulse_width/2, pulse_width/2, N_lfm).';

chirp_slope = Delta_f / pulse_width;

x_lfm = exp(1i * pi * chirp_slope * t_lfm.^2);

fprintf('Длина одного ЛЧМ-импульса: %d отсчетов\n', length(x_lfm));

n = (0:length_test-1).';
t = n/fs;

signal_bb = zeros(length_test,1);
idx = delay_samples + (1:length(x_lfm));
signal_bb(idx) = sqrt(signal_power)*x_lfm;

signal_if = real(signal_bb .* exp(1i*2*pi*fc*t));
noise_if  = sqrt(noise_power)*randn(length_test,1);
in_if     = signal_if + noise_if;

in_adc = quantize_real(in_if, n_a);

[x_bb, lo] = digital_demod(in_adc, fs, f_G, n_g);
x_bb = digital_lpf(x_bb, fs, 20e6, 128);

replica = flipud(conj(x_lfm));

mf_fir_nw = matched_filter_fir(x_bb, replica, false);
mf_fir_w  = matched_filter_fir(x_bb, replica, true);

mf_fft_nw = matched_filter_fft(x_bb, replica, [], false);
mf_fft_w  = matched_filter_fft(x_bb, replica, [], true);

mf_fftfilt_nw = matched_filter_fftfilt(x_bb, replica, false);
mf_fftfilt_w  = matched_filter_fftfilt(x_bb, replica, true);

fig1 = figure('Name','Вход после ЦДМ и ЦФ');
plot(abs(x_bb), 'LineWidth', 1.5);
grid on; axis tight;
xlabel('Отсчет'); ylabel('|x|, Амплитуда (огибающая)');
save_fig(fig1, 'fig_input_bb');

norm_db = @(x) mag2db(abs(x) / max(abs(x)));

fig2 = figure('Name','СФ FIR');
plot(norm_db(mf_fir_nw), 'LineWidth', 1.5); hold on;
plot(norm_db(mf_fir_w),  'LineWidth', 1.5);
grid on;
xlim([11500 13000]); % Зум вокруг пика (пик на 12288)
ylim([-100 5]);      % Ограничение по Y, чтобы видеть шум
xlabel('Отсчет'); ylabel('|y|, дБ (норм.)');
legend('без окна','с окном');
save_fig(fig2, 'fig_mf_fir');

fig3 = figure('Name','СФ FFT/IFFT');
plot(norm_db(mf_fft_nw), 'LineWidth', 1.5); hold on;
plot(norm_db(mf_fft_w),  'LineWidth', 1.5);
grid on;
xlim([11500 13000]);
ylim([-100 5]);
xlabel('Отсчет'); ylabel('|y|, дБ (норм.)');
legend('без окна','с окном');
save_fig(fig3, 'fig_mf_fft');

fig4 = figure('Name','СФ fftfilt');
plot(norm_db(mf_fftfilt_nw), 'LineWidth', 1.5); hold on;
plot(norm_db(mf_fftfilt_w),  'LineWidth', 1.5);
grid on;
xlim([11500 13000]);
ylim([-100 5]);
xlabel('Отсчет'); ylabel('|y|, дБ (норм.)');
legend('без окна','с окном');
save_fig(fig4, 'fig_mf_fftfilt');

fig5 = figure('Name','Сравнение методов СФ');
plot(norm_db(mf_fir_nw),     'b-', 'LineWidth', 4); hold on;
plot(norm_db(mf_fft_nw),     'r--', 'LineWidth', 2);
plot(norm_db(mf_fftfilt_nw), 'g-.', 'LineWidth', 1);

grid on; axis tight;
ylim([-100 5]);
xlabel('Отсчет'); ylabel('|y|, дБ (норм.)');
legend('FIR (толстая)', 'FFT/IFFT (штриховая)', 'fftfilt (штрихпунктирная)');
save_fig(fig5, 'fig_mf_compare');

function y = quantize_real(x, bits, full_scale)
    if nargin < 3 || isempty(full_scale)
        full_scale = max(abs(x(:)))*1.01;
        if full_scale == 0, full_scale = 1; end
    end
    LSB = 2*full_scale/(2^bits - 1);
    y = round(x/LSB)*LSB;
    y = max(min(y, full_scale), -full_scale);
end

function y = quantize_complex(x, bits, full_scale)
    if nargin < 3 || isempty(full_scale)
        full_scale = max(abs(x(:)))*1.01;
    end
    y = quantize_real(real(x), bits, full_scale) + ...
        1i*quantize_real(imag(x), bits, full_scale);
end

function [y, lo] = digital_demod(x, fs, f_G, n_g)
    n = (0:numel(x)-1).';
    lo = exp(-1i*2*pi*f_G*n/fs);
    lo = quantize_complex(lo, n_g, 1.0);
    y = x .* lo;
end

function y = digital_lpf(x, fs, f_cut, order)
    if nargin < 4, order = 128; end
    h = fir1(order, f_cut/(fs/2));
    y = filter(h, 1, x);
end

function y = matched_filter_fir(in, replica, use_window)
    if nargin < 3, use_window = false; end
    if use_window
        replica = replica .* hamming(length(replica));
    end
    y = filter(replica, 1, in);
end

function y = matched_filter_fft(in, replica, Nfft, use_window)
    if nargin < 3 || isempty(Nfft)
        Nfft = 2^nextpow2(length(in) + length(replica) - 1);
    end
    if nargin < 4, use_window = false; end
    if use_window
        replica = replica .* hamming(length(replica));
    end
    y_full = ifft(fft(in, Nfft) .* fft(replica, Nfft), Nfft);
    y = y_full(1:length(in));
end

function y = matched_filter_fftfilt(in, replica, use_window)
    if nargin < 3, use_window = false; end
    if use_window
        replica = replica .* hamming(length(replica));
    end
    y = fftfilt(replica, in);
end

function save_fig(figHandle, fileName)
    output_dir = 'export_figs';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    set(figHandle, 'Position', [100, 100, 900, 450]);
    
    set(figHandle, 'PaperPositionMode', 'auto');
    
    full_path = fullfile(output_dir, [fileName '.png']);
    print(figHandle, full_path, '-dpng', '-r300');
    fprintf('График сохранен: %s\n', full_path);
end
