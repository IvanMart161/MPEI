import numpy as np
import matplotlib.pyplot as plt
import os

# Создаем папку для экспорта графиков, если её нет
os.makedirs("exports_fig", exist_ok=True)

# ==========================================
# 1. Параметры тракта БВЧ (Рис. 2 отчета)
# ==========================================
fs = 200.0e6  # 200 МГц
N = 1_920_000
t = np.arange(N) / fs

np.random.seed(42)

# 1) Входной шум
u_in_std = 0.45e-5
w_in = np.random.normal(0, u_in_std, N)

# 2) Шум на выходе преселектора (f0 = 50 МГц, П = 3 МГц)
freqs = np.fft.rfftfreq(N, 1.0 / fs)
x_p = 2.0 * (freqs - 50.0e6) / 3.0e6
H_pres = 1.0 / np.sqrt(1.0 + x_p**4)
w_pres_raw = np.fft.irfft(np.fft.rfft(np.random.normal(0, 1.0, N)) * H_pres, n=N)
w_pres = w_pres_raw * (0.28e-6 / np.std(w_pres_raw))

# 3) Шум на выходе УПЧ (f_p = 10 МГц, П = 1 МГц, СКО = 6.36 мВ)
u_upch_rms = 6.36e-3
mask_upch = (freqs >= 9.5e6) & (freqs <= 10.5e6)
H_upch = np.zeros_like(freqs)
H_upch[mask_upch] = 1.0
w_upch_raw = np.fft.irfft(np.fft.rfft(np.random.normal(0, 1.0, N)) * H_upch, n=N)
w_upch = w_upch_raw * (u_upch_rms / np.std(w_upch_raw))

# ==========================================
# 2. Вычисление спектров
# ==========================================
seg_len = 8192
n_segs = N // seg_len
win = np.hanning(seg_len)
scale_psd = 2.0 / (fs * np.sum(win**2))

def get_psd(sig):
    reshaped = sig[:n_segs * seg_len].reshape(n_segs, seg_len)
    fft_res = np.fft.rfft(reshaped * win, axis=1)
    return np.mean(np.abs(fft_res)**2 * scale_psd, axis=0)

f_axis = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e6  # МГц

psd_in_curve = (get_psd(w_in) / 50.0) * 1e6
psd_in_curve = psd_in_curve * (4.26e-14 / np.mean(psd_in_curve[(f_axis >= 10) & (f_axis <= 90)]))

psd_pres_raw = (get_psd(w_pres) / 50.0) * 1e6
psd_pres_curve = psd_pres_raw * (5.5e-14 / np.max(psd_pres_raw))

psd_upch_raw = np.sqrt(get_psd(w_upch) * 1e6)
psd_upch_curve = psd_upch_raw * (6.6e-3 / np.max(psd_upch_raw))

# ==========================================
# 3. Отрисовка и сохранение в exports_fig/
# ==========================================
fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))

idx_10us = np.where((t >= 0) & (t <= 10.0e-6))[0]
t_us = t[idx_10us] * 1e6

# [0, 0] Входной шум
axs[0, 0].plot(t_us, w_in[idx_10us] * 1e5, color='#0011dd', lw=0.35)
axs[0, 0].set_title("Шум, приведённый ко входу\nБВЧ (в полосе 500 МГц)", fontsize=9)
axs[0, 0].set_xlabel("t, мкс", fontsize=8)
axs[0, 0].set_xlim(0, 10)
axs[0, 0].set_ylim(-2.2, 2.2)
axs[0, 0].text(-0.02, 1.03, r"$\times 10^{-5}$", transform=axs[0, 0].transAxes, fontsize=8)
axs[0, 0].grid(True, linestyle='--', color='lightgray')

# [0, 1] Выход преселектора
axs[0, 1].plot(t_us, w_pres[idx_10us] * 1e6, color='#0011dd', lw=0.45)
axs[0, 1].set_title("Шум на выходе преселектора", fontsize=9)
axs[0, 1].set_xlabel("t, мкс", fontsize=8)
axs[0, 1].set_xlim(0, 10)
axs[0, 1].set_ylim(-1.2, 1.2)
axs[0, 1].text(-0.02, 1.03, r"$\times 10^{-6}$", transform=axs[0, 1].transAxes, fontsize=8)
axs[0, 1].grid(True, linestyle='--', color='lightgray')

# [0, 2] Выход УПЧ
axs[0, 2].plot(t_us, w_upch[idx_10us], color='#0011dd', lw=0.5)
axs[0, 2].set_title("Шум на выходе УПЧ", fontsize=9)
axs[0, 2].set_xlabel("t, мкс", fontsize=8)
axs[0, 2].set_xlim(0, 10)
axs[0, 2].set_ylim(-0.015, 0.015)
axs[0, 2].grid(True, linestyle='--', color='lightgray')

# [1, 0] СПМ приведенного шума
axs[1, 0].plot(f_axis, psd_in_curve * 1e14, color='#cc0000', lw=1.1)
axs[1, 0].set_title("Спектральная плотность мощности\nприведенного шума", fontsize=9)
axs[1, 0].set_xlabel("f, МГц", fontsize=8)
axs[1, 0].set_ylabel("Спектральная плотность (Вт/МГц)", fontsize=8)
axs[1, 0].set_xlim(0, 100)
axs[1, 0].set_ylim(0, 6.0)
axs[1, 0].text(-0.02, 1.03, r"$\times 10^{-14}$", transform=axs[1, 0].transAxes, fontsize=8)
axs[1, 0].grid(True, linestyle='--', color='lightgray')

# [1, 1] СПМ преселектора
axs[1, 1].plot(f_axis, psd_pres_curve * 1e14, color='#cc0000', lw=1.1)
axs[1, 1].set_title("Спектральная плотность мощности шума\nна выходе преселектора", fontsize=9)
axs[1, 1].set_xlabel("f, МГц", fontsize=8)
axs[1, 1].set_xlim(40, 60)
axs[1, 1].set_ylim(0, 8.0)
axs[1, 1].text(-0.02, 1.03, r"$\times 10^{-14}$", transform=axs[1, 1].transAxes, fontsize=8)
axs[1, 1].grid(True, linestyle='--', color='lightgray')

# [1, 2] СП напряжения УПЧ
axs[1, 2].plot(f_axis, psd_upch_curve * 1e3, color='#cc0000', lw=1.1)
axs[1, 2].set_title("Спектральная плотность напряжения\nшума на выходе УПЧ", fontsize=9)
axs[1, 2].set_xlabel("f, МГц", fontsize=8)
axs[1, 2].set_ylabel(r"Спектральная плотность (В/$\sqrt{\mathrm{МГц}}$)", fontsize=8)
axs[1, 2].set_xlim(8.0, 12.0)
axs[1, 2].set_ylim(0, 10.0)
axs[1, 2].text(-0.02, 1.03, r"$\times 10^{-3}$", transform=axs[1, 2].transAxes, fontsize=8)
axs[1, 2].grid(True, linestyle='--', color='lightgray')
# Путь с тильдой, который развернется в полный путь к домашней директории
output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

output_path = os.path.join(output_dir, "pic1.png")

# Сохранение графика
plt.savefig(output_path, dpi=200)
print(f"График успешно сохранен: {output_path}")
