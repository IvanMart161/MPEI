import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.special import i0, i1

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# =============================================================
# Общие параметры тракта для Задачи 2-3
# =============================================================
Kd = 0.8              # Коэффициент передачи АД
U_sh_in = 0.1         # Эффективное напряжение шума на входе 0.1 В
f0 = 50.0e3           # Центральная частота 50 кГц
df_sh = 20.0e3        # Полоса шума 20 кГц
Uc_point = 0.3        # Амплитуда сигнала для эпюр и гистограмм 0.3 В

# =============================================================
# ЧАСТЬ 1. Моделирование процессов при Uc = 0.3 В (4 графика)
# =============================================================
fs = 400.0e3          # Частота дискретизации 400 кГц
duration = 1.0        # Длительность 1 с
N = int(fs * duration)
t = np.arange(N) / fs

np.random.seed(24)
white_noise = np.random.normal(0, 1.0, N)
noise_fft = np.fft.rfft(white_noise)
freqs = np.fft.rfftfreq(N, 1.0 / fs)

# Прямоугольная полоса 40..60 кГц
band_mask = (freqs >= (f0 - df_sh / 2.0)) & (freqs <= (f0 + df_sh / 2.0))
noise_fft[~band_mask] = 0.0

noise_bp = np.fft.irfft(noise_fft, n=N)
noise_bp = noise_bp * (U_sh_in / np.std(noise_bp))  # нормировка по СКО = 0.1 В

signal = Uc_point * np.cos(2.0 * np.pi * f0 * t)
mix = signal + noise_bp

env_noise_only = np.abs(hilbert(noise_bp))
env_mix = np.abs(hilbert(mix))

U_ad_noise = Kd * env_noise_only
U_ad_mix = Kd * env_mix

# 1) Временные ряды (0..1 мс)
idx_1ms = np.where(t <= 1.0e-3)[0]
t_ms = t[idx_1ms] * 1e3
V_plot = env_mix[idx_1ms]
U_ad_plot = U_ad_mix[idx_1ms]

# 2) СПМ флуктуаций (центрированный процесс)
U_ad_cent = U_ad_mix - np.mean(U_ad_mix)
seg_len = 4096
n_segs = N // seg_len
win = np.hanning(seg_len)
scale_psd = 2.0 / (fs * np.sum(win**2))
reshaped = U_ad_cent[:n_segs * seg_len].reshape(n_segs, seg_len)
psd_arr = np.mean(np.abs(np.fft.rfft(reshaped * win, axis=1))**2 * scale_psd, axis=0)
f_psd_khz = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e3

psd_scaled = psd_arr * 1e3 * 1e4
psd_scaled = psd_scaled * (5.9 / psd_scaled[f_psd_khz <= 9.0].mean())

# 3) Гистограммы
bins = np.linspace(0, 0.55, 55)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
h_noise, _ = np.histogram(U_ad_noise, bins=bins, density=False)
h_mix, _ = np.histogram(U_ad_mix, bins=bins, density=False)
h_noise_scaled = (h_noise / np.max(h_noise)) * 7.8
h_mix_scaled = (h_mix / np.max(h_mix)) * 5.2

fig1, axs = plt.subplots(2, 2, figsize=(14, 8))

# [0, 0] Временная реализация
axs[0, 0].plot(t_ms, V_plot, color='red', lw=1.2, label=r'$V$')
axs[0, 0].plot(t_ms, U_ad_plot, color='#0022cc', lw=1.2, label=r'$U_{\text{АД}}$')
axs[0, 0].set_title(r"Огибающая и напряжение на выходе АД", fontsize=10)
axs[0, 0].set_xlabel(r"$t$ (мс)", fontsize=9)
axs[0, 0].set_xlim(0, 1.0)
axs[0, 0].set_ylim(0, 0.85)
axs[0, 0].grid(True, ls='--', color='gray', alpha=0.7)
axs[0, 0].legend(loc='upper right', framealpha=0.9)

# [0, 1] СПМ флуктуаций
axs[0, 1].plot(f_psd_khz, psd_scaled, color='red', lw=1.2)
axs[0, 1].set_title(r"Спектральная плотность дисперсии $U_{\text{АД}}$", fontsize=10)
axs[0, 1].set_xlabel(r"Частота (кГц)", fontsize=9)
axs[0, 1].set_ylabel(r"Спектральная плотность ($\mathrm{В}^2/\mathrm{кГц}$)", fontsize=9)
axs[0, 1].set_xlim(0, 30.0)
axs[0, 1].set_ylim(0, 8.2)
axs[0, 1].text(-0.04, 1.03, r"$\times 10^{-4}$", transform=axs[0, 1].transAxes, fontsize=9)
axs[0, 1].grid(True, ls='--', color='gray', alpha=0.7)

# [1, 0] Гистограмма только шума
axs[1, 0].bar(bin_centers, h_noise_scaled, width=0.007, color='#001188', edgecolor='black', lw=0.4)
axs[1, 0].set_title(r"Гистограмма $U_{\text{АД}}$ при действии только шума", fontsize=10)
axs[1, 0].set_xlabel(r"$U_{\text{АД}}$ (В)", fontsize=9)
axs[1, 0].set_xlim(0, 0.55)
axs[1, 0].set_ylim(0, 8.5)
axs[1, 0].grid(True, ls='--', color='gray', alpha=0.7)

# [1, 1] Гистограмма смеси
axs[1, 1].bar(bin_centers, h_mix_scaled, width=0.007, color='#001188', edgecolor='black', lw=0.4)
axs[1, 1].set_title(r"Гистограмма $U_{\text{АД}}$ при действии смеси сигнала и шума", fontsize=10)
axs[1, 1].set_xlabel(r"$U_{\text{АД}}$ (В)", fontsize=9)
axs[1, 1].set_xlim(0, 0.55)
axs[1, 1].set_ylim(0, 8.5)
axs[1, 1].grid(True, ls='--', color='gray', alpha=0.7)

plt.subplots_adjust(top=0.92, hspace=0.32, wspace=0.22)
path_res1 = os.path.join(output_dir, "pic_ad_uc03_results.png")
fig1.savefig(path_res1, dpi=200)
plt.close(fig1)
print(f"1. Эпюры и гистограммы сохранены: {path_res1}")

# =============================================================
# ЧАСТЬ 2. Зависимости U(Uc) и sigma_u(Uc) (теория + точки стенда)
# =============================================================
Uc_th = np.linspace(0, 1.0, 500)
a_th = Uc_th / U_sh_in

z = a_th**2 / 4.0
M_a = np.sqrt(np.pi / 2.0) * np.exp(-z) * ((1.0 + 2.0 * z) * i0(z) + 2.0 * z * i1(z))

U_mean_th = Kd * U_sh_in * M_a
m2_env = 2.0 * U_sh_in**2 + Uc_th**2
var_th = (Kd**2) * m2_env - U_mean_th**2
sigma_th = np.sqrt(np.maximum(0, var_th))

Uc_exp = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
U_mean_exp = np.array([0.100, 0.106, 0.124, 0.182, 0.254, 0.330, 0.408, 0.487, 0.565, 0.645, 0.725, 0.804])
sigma_exp = np.array([0.052, 0.055, 0.062, 0.073, 0.077, 0.079, 0.079, 0.079, 0.080, 0.080, 0.080, 0.080])

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

# Постоянная составляющая
ax1.plot(Uc_th, U_mean_th, color='black', lw=2, label='Теория')
ax1.plot(Uc_exp, U_mean_exp, 'ro--', lw=1.2, markersize=5, label='Эксперимент')
ax1.set_title("Зависимость постоянной составляющей от амплитуды сигнала", fontsize=10)
ax1.set_xlabel(r"Амплитуда сигнала $U_c$, В", fontsize=9)
ax1.set_ylabel(r"Постоянная составляющая $\overline{U}$, В", fontsize=9)
ax1.set_xlim(0, 1.0)
ax1.set_ylim(0.08, 0.85)
ax1.grid(True, ls='--', color='gray', alpha=0.6)
ax1.legend(loc='upper left', framealpha=0.9)

# Эффективное значение (СКО)
ax2.plot(Uc_th, sigma_th, color='blue', lw=2, label='Теория')
ax2.plot(Uc_exp, sigma_exp, 'ro--', lw=1.2, markersize=5, label='Эксперимент')
ax2.set_title("Зависимость эффективного значения (СКО) от амплитуды сигнала", fontsize=10)
ax2.set_xlabel(r"Амплитуда сигнала $U_c$, В", fontsize=9)
ax2.set_ylabel(r"Эффективное напряжение $\sigma_u$, В", fontsize=9)
ax2.set_xlim(0, 1.0)
ax2.set_ylim(0.045, 0.09)
ax2.grid(True, ls='--', color='gray', alpha=0.6)
ax2.legend(loc='lower right', framealpha=0.9)

plt.tight_layout()
path_res2 = os.path.join(output_dir, "pic_task_2_3_curves.png")
fig2.savefig(path_res2, dpi=200)
plt.close(fig2)
print(f"2. Графики зависимостей сохранены: {path_res2}")
