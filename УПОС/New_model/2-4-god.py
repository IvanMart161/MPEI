import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------------------
# Параметры моделирования для режима Uc = 0.5 В
# -------------------------------------------------------------
f0 = 50.0e3           # Центральная частота 50 кГц
df_sh = 20.0e3        # Полоса шума 20 кГц
U_sh_in = 0.1         # СКО шума 0.1 В
Uc = 0.5              # Амплитуда сигнала 0.5 В
fs = 400.0e3          # Частота дискретизации

N = int(fs * 0.005)   # 5 мс
t = np.arange(N) / fs

# Прямоугольный спектр шума
np.random.seed(31)
w = np.random.normal(0, 1.0, N)
W = np.fft.rfft(w)
freqs = np.fft.rfftfreq(N, 1.0 / fs)
mask = (freqs >= f0 - df_sh/2.0) & (freqs <= f0 + df_sh/2.0)
W[~mask] = 0.0
noise = np.fft.irfft(W, n=N)
noise = noise * (U_sh_in / np.std(noise))

# Квадратуры шума (аналитический сигнал, сдвинутый на f0)
an_noise = hilbert(noise)
env_base = an_noise * np.exp(-1j * 2.0 * np.pi * f0 * t)
uc_noise = np.real(env_base)
us_noise = np.imag(env_base)

# Сигнал и смесь
signal = Uc * np.cos(2.0 * np.pi * f0 * t)
mix = signal + noise
mix_env = np.abs(hilbert(mix))

# Интервал наблюдения [0.2; 0.3] мс
mask_int = (t >= 0.2e-3) & (t <= 0.3e-3)
t_ms = t[mask_int] * 1e3
mix_seg = mix[mask_int]
env_seg = mix_env[mask_int]

# Квадратуры смеси: косинусная смещается на Uc
Uc_mix = uc_noise[mask_int] + Uc
Us_mix = us_noise[mask_int]

# -------------------------------------------------------------
# Отрисовка триптиха стенда (3 графика)
# -------------------------------------------------------------
fig = plt.figure(figsize=(13, 7.5))

# 1. Реализация смеси сигнала и шума (слева вверху)
ax1 = fig.add_axes([0.08, 0.56, 0.40, 0.36])
ax1.plot(t_ms, mix_seg, color='#0011cc', lw=1.2, label='Смесь')
ax1.plot(t_ms, env_seg, color='red', lw=1.2, label='Огибающая')
ax1.plot(t_ms[-1], env_seg[-1], 'ro', markersize=4)
ax1.set_title("Реализация смеси сигнала и шума", fontsize=10)
ax1.set_xlabel(r"$t$ (мс)", fontsize=9)
ax1.set_xlim(0.2, 0.3)
ax1.set_ylim(-1.0, 1.0)
ax1.grid(True, ls='--', color='gray', alpha=0.6)

# 2. Квадратурные составляющие шума (справа вверху)
ax2 = fig.add_axes([0.55, 0.56, 0.40, 0.36])
ax2.plot(t_ms, Uc_mix, color='red', lw=1.3, label=r"$U^c$")
ax2.plot(t_ms, Us_mix, color='#0011cc', lw=1.3, label=r"$U^s$")
ax2.plot(t_ms[-1], Uc_mix[-1], 'ro', markersize=4)
ax2.plot(t_ms[-1], Us_mix[-1], 'bo', markersize=4)
ax2.axhline(0, color='black', lw=0.6)
ax2.set_title("Квадратурные составляющие шума", fontsize=10)
ax2.set_xlabel(r"$t$ (мс)", fontsize=9)
ax2.set_xlim(0.2, 0.3)
ax2.set_ylim(-0.8, 0.8)
ax2.grid(True, ls='--', color='gray', alpha=0.6)
ax2.legend(loc='upper right', framealpha=0.9)

# 3. Годограф комплексной огибающей (снизу по центру)
ax3 = fig.add_axes([0.31, 0.08, 0.42, 0.40])
ax3.plot(Uc_mix, Us_mix, color='#0011cc', lw=1.4)
ax3.axhline(0, color='black', lw=0.7)
ax3.axvline(0, color='black', lw=0.7)

# Зеленый вектор сигнала Uc
ax3.annotate('', xy=(Uc, 0), xytext=(0, 0),
             arrowprops=dict(edgecolor='#00bb00', facecolor='#00bb00', width=2.5, headwidth=6))

# Точка в конце реализации
pt_end = -1
# Вектор шума (желтый/оранжевый от конца Uc к траектории)
ax3.plot([Uc, Uc_mix[pt_end]], [0, Us_mix[pt_end]], color='#ffaa00', lw=2)
# Маркер текущей точки
ax3.plot(Uc_mix[pt_end], Us_mix[pt_end], 'ko', markersize=4)

ax3.set_title("Годограф комплексной огибающей", fontsize=10)
ax3.set_xlabel(r"$U^c$, В", fontsize=9)
ax3.set_ylabel(r"$U^s$, В", fontsize=9)
ax3.set_xlim(-0.5, 1.1)
ax3.set_ylim(-0.8, 0.8)
ax3.grid(True, ls='--', color='gray', alpha=0.6)

save_path = os.path.join(output_dir, "pic_godograf_stand_uc05.png")
plt.savefig(save_path, dpi=200)
plt.close()
print(f"График сохранен: {save_path}")
