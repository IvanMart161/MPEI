import numpy as np
from scipy.signal.windows import blackmanharris
import scipy.io as sio
import scipy.signal as signal
import matplotlib.pyplot as plt
import psd

filename = 'test1_25p0.mat'
mat = sio.loadmat(filename)

x = mat['pdin'].flatten()
y = mat['pdout'].flatten()
tx_freq = float(mat['tx_freq'].squeeze())
rx_freq = float(mat['rx_freq'].squeeze())

Fs = 245.76e6  

win = blackmanharris(2048)



X_init, frq = psd.psd_welch(x, 2048, win, 1024, Fs)
Y_init, _   = psd.psd_welch(y, 2048, win, 1024, Fs)

plt.figure(1, figsize=(9, 6))
plt.plot(frq / 1e6, X_init, label='pdin (x)')
plt.plot(frq / 1e6, Y_init, label='pdout (y)')
plt.title('Figure 1: Исходные спектры сигналов')
plt.xlabel('freq, MHz')
plt.ylabel('PSD, dB/Hz')
plt.grid(True)
plt.legend()

y_no_dc = y - np.mean(y)



Y_no_dc, _ = psd.psd_welch(y_no_dc, 2048, win, 1024, Fs)

plt.figure(2, figsize=(9, 6))
plt.plot(frq / 1e6, X_init, label='pdin (x)')
plt.plot(frq / 1e6, Y_no_dc, label='pdout без DC (y)')
plt.title('Figure 2 (Step 1): Удаление постоянной составляющей')
plt.xlabel('freq, MHz')
plt.ylabel('PSD, dB/Hz')
plt.grid(True)
plt.legend()



delta_f = rx_freq - tx_freq
t = np.arange(len(y_no_dc)) / Fs
y_freq_comp = y_no_dc * np.exp(1j * 2 * np.pi * delta_f * t)



Y_freq_comp, _ = psd.psd_welch(y_freq_comp, 2048, win, 1024, Fs)

plt.figure(3, figsize=(9, 6))
plt.plot(frq / 1e6, X_init, label='pdin (x)')
plt.plot(frq / 1e6, Y_freq_comp, label='pdout скомпенсированный по частоте')
plt.title('Figure 2 (Step 2): Компенсация частотного смещения rx_freq и tx_freq')
plt.xlabel('freq, MHz')
plt.ylabel('PSD, dB/Hz')
plt.grid(True)
plt.legend()



corr = signal.correlate(y_freq_comp, x, mode='valid')
delay = int(np.argmax(np.abs(corr)))
print(f'delay = {delay} samples')

y_sync = y_freq_comp[delay:delay + len(x)]



fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), num=4)

ax1.plot(np.abs(x), label='|x| (pdin, длина 43840)', alpha=0.8)
ax1.plot(np.abs(y_freq_comp), label='|y| (pdout, длина 65536)', alpha=0.7)
ax1.set_title(f'Figure 3: Выравнивание во времени (delay = {delay} отсчетов)')
ax1.set_xlabel('Отсчеты')
ax1.set_ylabel('Амплитуда')
ax1.grid(True)
ax1.legend()

zoom_idx = np.arange(2600, 2900)
ax2.plot(zoom_idx, np.abs(x[zoom_idx]), label='|x| (вход)', color='C0')
ax2.plot(zoom_idx, np.abs(y_sync[zoom_idx]), label='|y_sync| (выход после сдвига)', color='C1')
ax2.set_title('Зум синхронизированного участка (отсчеты 2600-2900)')
ax2.set_xlabel('Отсчеты')
ax2.set_ylabel('Амплитуда')
ax2.grid(True)
ax2.legend()
plt.tight_layout()



g = np.vdot(y_sync, x) / np.vdot(y_sync, y_sync)
gain_mag_db = 20 * np.log10(np.abs(g))
gain_phase_rad = np.angle(g)

print(f'Complex Gain magnitude = {gain_mag_db:.2f} dB')
print(f'Complex Gain phase = {gain_phase_rad:.4f} rad')

y_aligned = y_sync * g

e = x - y_aligned

Y_aligned_psd, _ = psd.psd_welch(y_aligned, 2048, win, 1024, Fs)
E_psd, _         = psd.psd_welch(e, 2048, win, 1024, Fs)

plt.figure(5, figsize=(9, 6))
plt.plot(frq / 1e6, X_init, label='x (входной сигнал)')
plt.plot(frq / 1e6, Y_aligned_psd, label='y * g (выход скомпенсированный)')
plt.plot(frq / 1e6, E_psd, label='e = x - y*g (сигнал ошибки)', color='green')
plt.title('Figure 4: Спектральная плотность мощности и искажения в тракте')
plt.xlabel('freq, MHz')
plt.ylabel('PSD, dB/Hz')
plt.grid(True)
plt.legend()

plt.show()
