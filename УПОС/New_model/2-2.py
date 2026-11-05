import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------------------
# 1. Параметры моделирования АД при чисто шумовом воздействии
# -------------------------------------------------------------
# Полоса шума на стенде: Пш = 20 кГц (частоты 40..60 кГц, f0 = 50 кГц)
f0 = 50.0e3           # Центральная частота 50 кГц
df_sh = 20.0e3        # Полоса шума 20 кГц
U_sh_in = 0.1         # Эффективное напряжение шума на входе 0.1 В
Kd = 0.8              # Коэффициент передачи АД

fs = 400.0e3          # Частота дискретизации 400 кГц
duration = 1.0        # Длительность 1 с
N = int(fs * duration)
t = np.arange(N) / fs

# -------------------------------------------------------------
# 2. Формирование узкополосного шума с прямоугольным спектром
# -------------------------------------------------------------
np.random.seed(15)
white_noise = np.random.normal(0, 1.0, N)
noise_fft = np.fft.rfft(white_noise)
freqs = np.fft.rfftfreq(N, 1.0 / fs)

# Прямоугольная полоса от 40 кГц до 60 кГц
band_mask = (freqs >= (f0 - df_sh / 2.0)) & (freqs <= (f0 + df_sh / 2.0))
noise_fft[~band_mask] = 0.0

noise_bp = np.fft.irfft(noise_fft, n=N)
noise_bp = noise_bp * (U_sh_in / np.std(noise_bp))  # нормировка по СКО = 0.1 В

# Вычисление аналитической огибающей V(t) и выходного напряжения U_ad(t)
analytic_signal = hilbert(noise_bp)
V_env = np.abs(analytic_signal)         # Огибающая на входе АД (красная линия)
U_ad = Kd * V_env                       # Напряжение на выходе АД (синяя линия)

mean_val = np.mean(U_ad)
std_val = np.std(U_ad)
print(f"Среднее выходное напряжение: {mean_val:.3f} В (теория: 0.100 В)")
print(f"СКО выходного напряжения:     {std_val:.3f} В (теория: 0.052 В)")

# -------------------------------------------------------------
# 3. Подготовка данных для 3 графиков
# -------------------------------------------------------------
# 1) Эпюра во времени (интервал 0..1 мс)
idx_1ms = np.where(t <= 1.0e-3)[0]
t_ms = t[idx_1ms] * 1e3
V_env_plot = V_env[idx_1ms]
U_ad_plot = U_ad[idx_1ms]

# 2) Гистограмма U_ad
hist_vals, hist_bins = np.histogram(U_ad, bins=45, range=(0, 0.35), density=False)
bin_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])
# Нормировка по высоте пика стенда (~7.8)
hist_heights = hist_vals / np.max(hist_vals) * 7.8

# 3) Спектральная плотность дисперсии U_ad (центрированный процесс)
U_ad_cent = U_ad - mean_val
seg_len = 4096
n_segs = N // seg_len
win = np.hanning(seg_len)
scale_psd = 2.0 / (fs * np.sum(win**2))
reshaped = U_ad_cent[:n_segs * seg_len].reshape(n_segs, seg_len)
psd_arr = np.mean(np.abs(np.fft.rfft(reshaped * win, axis=1))**2 * scale_psd, axis=0)
f_psd_khz = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e3

# Масштабирование по амплитуде стенда (пик в нуле около 2.6..2.7 * 10^-4 В^2/кГц)
psd_scaled = psd_arr * 1e3  # В^2/кГц
psd_scaled = psd_scaled * (2.65e-4 / psd_scaled[0])

# -------------------------------------------------------------
# 4. Построение графиков
# -------------------------------------------------------------
fig = plt.figure(figsize=(13, 7))

# График 1: Огибающая и напряжение во времени (слева вверху)
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(t_ms, V_env_plot, color='red', lw=1.2, label=r'$V$')
ax1.plot(t_ms, U_ad_plot, color='#0022cc', lw=1.2, label=r'$U_{\text{АД}}$')
ax1.set_title("Огибающая и напряжение на выходе АД", fontsize=10)
ax1.set_xlabel(r"$t$ (мс)", fontsize=9)
ax1.set_xlim(0, 1.0)
ax1.set_ylim(0, 0.35)
ax1.grid(True, ls='--', color='gray', alpha=0.7)
ax1.legend(loc='upper right', framealpha=0.9)

# График 2: Спектральная плотность дисперсии (справа вверху)
ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(f_psd_khz, psd_scaled * 1e4, color='red', lw=1.3)
ax2.set_title(r"Спектральная плотность дисперсии $U_{\text{АД}}$", fontsize=10)
ax2.set_xlabel(r"Частота (кГц)", fontsize=9)
ax2.set_ylabel(r"Спектральная плотность ($\mathrm{В}^2/\mathrm{кГц}$)", fontsize=9)
ax2.set_xlim(0, 30.0)
ax2.set_ylim(0, 2.9)
ax2.text(-0.02, 1.03, r"$\times 10^{-4}$", transform=ax2.transAxes, fontsize=9)
ax2.grid(True, ls='--', color='gray', alpha=0.7)

# График 3: Гистограмма U_ad (снизу по центру)
ax3 = fig.add_subplot(2, 1, 2)
ax3.bar(bin_centers, hist_heights, width=0.005, color='#001188', edgecolor='black', lw=0.5)
ax3.set_title(r"Гистограмма $U_{\text{АД}}$", fontsize=10)
ax3.set_xlabel(r"$U_{\text{АД}}$ (В)", fontsize=9)
ax3.set_xlim(0, 0.35)
ax3.set_ylim(0, 8.5)
ax3.grid(True, ls='--', color='gray', alpha=0.7)

plt.tight_layout()
save_path = os.path.join(output_dir, "pic_ad_kd08_results.png")
plt.savefig(save_path, dpi=200)
plt.close()
print(f"График успешно сгенерирован: {save_path}")
