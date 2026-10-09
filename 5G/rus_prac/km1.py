import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import psd # Пользовательский модуль psd.py должен быть в той же папке

# Загрузка данных
mat = sio.loadmat('test1_25p0.mat')
x = mat['pdin']
y = mat['pdout']

# Извлечение частот и выравнивание массивов (перенесено из шаблона вверх)
tx_freq = float(np.reshape(mat['tx_freq'], 1))
rx_freq = float(np.reshape(mat['rx_freq'], 1))
x = np.reshape(x, x.size)
y = np.reshape(y, y.size)
Fs = 245.76E6

# STEP 1. Remove ADC DC
# Вычитание постоянной составляющей приемника
y = y - np.mean(y)

# STEP 2. Frequency shift compensation
# Компенсация частотного смещения между rx_freq и tx_freq
delta_f = rx_freq - tx_freq
t_y = np.arange(y.size) / Fs
y = y * np.exp(1j * 2 * np.pi * delta_f * t_y)

# STEP 3. Time shift compensation
# Взаимная корреляция для поиска задержки (используем быструю свертку через БПФ)
corr = signal.correlate(y, x, mode='valid', method='fft')
delay = np.argmax(np.abs(corr))
print('delay = %d samples' % delay)

# Выравниваем выходной сигнал по длине входного с учетом найденной задержки
y_aligned = y[delay : delay + x.size]

# STEP 4. Gain control calculation
# Расчет комплексного коэффициента выравнивания по амплитуде и фазе
# Формула: g = (y^H * x) / (y^H * y)
numerator = np.sum(np.conj(y_aligned) * x)
denominator = np.sum(np.conj(y_aligned) * y_aligned)
g = numerator / denominator

mag_g_dB = 20 * np.log10(np.abs(g))
phase_g_rad = np.angle(g)
print('Complex Gain magnitude = %.4f dB' % mag_g_dB)
print('Complex Gain phase = %.4f rad' % phase_g_rad)

# Применяем коэффициент к выровненному сигналу
y_aligned_g = y_aligned * g

# STEP 5. Error calculation and PSD plotting
# Расчет вектора искажений в трактах
e = x - y_aligned_g

# Расчет спектральной плотности мощности (СПМ) методом Уэлча
(X, frq) = psd.psd_welch(x, 2048, signal.windows.blackmanharris(2048), 1024, Fs)
(Y, frq) = psd.psd_welch(y_aligned_g, 2048, signal.windows.blackmanharris(2048), 1024, Fs)
(E, frq) = psd.psd_welch(e, 2048, signal.windows.blackmanharris(2048), 1024, Fs)

# Отрисовка графиков
plt.figure(figsize=(10, 6))
plt.plot(frq / 1E6, X, label='Входной сигнал (x)')
plt.plot(frq / 1E6, Y, label='Выходной сигнал с компенсацией (y)')
plt.plot(frq / 1E6, E, label='Сигнал ошибки (e)')
plt.xlabel('freq, MHz')
plt.ylabel('PSD, dB/Hz')
plt.grid(True)
plt.legend()
plt.title('Спектральная плотность мощности (PSD)')
plt.show()