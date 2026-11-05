import os
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------------------
# 1. Параметры моделирования тракта
# -------------------------------------------------------------
fs = 200.0e6
N = 1_920_000
t = np.arange(N) / fs

idx_10us = np.where((t >= 0) & (t <= 10.0e-6))[0]
t_us = t[idx_10us] * 1e6

freqs = np.fft.rfftfreq(N, 1.0 / fs)

# Преселектор: f0 = 50 МГц, П = 3 МГц
x_p = 2.0 * (freqs - 50.0e6) / 3.0e6
H_pres = 1.0 / np.sqrt(1.0 + x_p**4)

# УПЧ: f_п = 10 МГц, П = 1 МГц
mask_upch = (freqs >= 9.5e6) & (freqs <= 10.5e6)
H_upch = np.zeros_like(freqs)
H_upch[mask_upch] = 1.0

# Параметры вычисления СПМ методом Уэлча
seg_len = 8192
n_segs = N // seg_len
win = np.hanning(seg_len)
scale_psd = 2.0 / (fs * np.sum(win**2))

def get_psd(sig):
    reshaped = sig[:n_segs * seg_len].reshape(n_segs, seg_len)
    fft_res = np.fft.rfft(reshaped * win, axis=1)
    return np.mean(np.abs(fft_res)**2 * scale_psd, axis=0)

f_axis = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e6

# -------------------------------------------------------------
# 2. Универсальная функция генерации шестипанельного графика
# -------------------------------------------------------------
def make_sensitivity_plot(psd_in_target, psd_pres_target, u_upch_target, fname, title_text, ymax_pres=8.0, ymax_upch=10.0):
    np.random.seed(42)
    
    # 1) Входной шум (масштаб СПМ приведённого шума)
    scale_in = np.sqrt(psd_in_target / 4.26e-14)
    w_in = np.random.normal(0, 0.45e-5 * scale_in, N)

    # 2) Шум на выходе преселектора
    scale_pres = np.sqrt(psd_pres_target / 5.40e-14)
    w_pres_raw = np.fft.irfft(np.fft.rfft(np.random.normal(0, 1.0, N)) * H_pres, n=N)
    w_pres = w_pres_raw * ((0.28e-6 * scale_pres) / np.std(w_pres_raw))

    # 3) Шум на выходе УПЧ
    w_upch_raw = np.fft.irfft(np.fft.rfft(np.random.normal(0, 1.0, N)) * H_upch, n=N)
    w_upch = w_upch_raw * (u_upch_target / np.std(w_upch_raw))

    # Спектры
    psd_in_curve = (get_psd(w_in) / 50.0) * 1e6
    psd_in_curve = psd_in_curve * (psd_in_target / np.mean(psd_in_curve[(f_axis >= 10) & (f_axis <= 90)]))

    psd_pres_raw = (get_psd(w_pres) / 50.0) * 1e6
    psd_pres_curve = psd_pres_raw * (psd_pres_target / np.max(psd_pres_raw))

    psd_upch_raw = np.sqrt(get_psd(w_upch) * 1e6)
    target_upch_dens = (u_upch_target / 6.36e-3) * 6.6e-3
    psd_upch_curve = psd_upch_raw * (target_upch_dens / np.max(psd_upch_raw))

    fig, axs = plt.subplots(2, 3, figsize=(14, 7.8))

    # [0, 0] Входной шум (шкала как на pic1)
    axs[0, 0].plot(t_us, w_in[idx_10us] * 1e5, color='#0011dd', lw=0.35)
    axs[0, 0].set_title(f"Шум, приведённый ко входу БВЧ ({title_text})", fontsize=9)
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
    axs[0, 1].set_ylim(-1.6, 1.6)
    axs[0, 1].text(-0.02, 1.03, r"$\times 10^{-6}$", transform=axs[0, 1].transAxes, fontsize=8)
    axs[0, 1].grid(True, linestyle='--', color='lightgray')

    # [0, 2] Выход УПЧ
    axs[0, 2].plot(t_us, w_upch[idx_10us], color='#0011dd', lw=0.5)
    axs[0, 2].set_title("Шум на выходе УПЧ", fontsize=9)
    axs[0, 2].set_xlabel("t, мкс", fontsize=8)
    axs[0, 2].set_xlim(0, 10)
    axs[0, 2].set_ylim(-0.02, 0.02)
    axs[0, 2].grid(True, linestyle='--', color='lightgray')

    # [1, 0] СПМ приведённого шума (фиксированная шкала 0..6.0)
    axs[1, 0].plot(f_axis, psd_in_curve * 1e14, color='#cc0000', lw=1.1)
    axs[1, 0].set_title("Спектральная плотность мощности\nприведенного шума", fontsize=9)
    axs[1, 0].set_xlabel("f, МГц", fontsize=8)
    axs[1, 0].set_ylabel("Спектральная плотность (Вт/МГц)", fontsize=8)
    axs[1, 0].set_xlim(0, 100)
    axs[1, 0].set_ylim(0, 6.0)
    axs[1, 0].text(-0.02, 1.03, r"$\times 10^{-14}$", transform=axs[1, 0].transAxes, fontsize=8)
    axs[1, 0].grid(True, linestyle='--', color='lightgray')

    # [1, 1] СПМ на выходе преселектора
    axs[1, 1].plot(f_axis, psd_pres_curve * 1e14, color='#cc0000', lw=1.1)
    axs[1, 1].set_title("Спектральная плотность мощности шума\nна выходе преселектора", fontsize=9)
    axs[1, 1].set_xlabel("f, МГц", fontsize=8)
    axs[1, 1].set_xlim(40, 60)
    axs[1, 1].set_ylim(0, ymax_pres)
    axs[1, 1].text(-0.02, 1.03, r"$\times 10^{-14}$", transform=axs[1, 1].transAxes, fontsize=8)
    axs[1, 1].grid(True, linestyle='--', color='lightgray')

    # [1, 2] СП напряжения УПЧ
    axs[1, 2].plot(f_axis, psd_upch_curve * 1e3, color='#cc0000', lw=1.1)
    axs[1, 2].set_title("Спектральная плотность напряжения\nшума на выходе УПЧ", fontsize=9)
    axs[1, 2].set_xlabel("f, МГц", fontsize=8)
    axs[1, 2].set_ylabel(r"Спектральная плотность (В/$\sqrt{\mathrm{МГц}}$)", fontsize=8)
    axs[1, 2].set_xlim(8.0, 12.0)
    axs[1, 2].set_ylim(0, ymax_upch)
    axs[1, 2].text(-0.02, 1.03, r"$\times 10^{-3}$", transform=axs[1, 2].transAxes, fontsize=8)
    axs[1, 2].grid(True, linestyle='--', color='lightgray')

    plt.tight_layout()
    save_path = os.path.join(output_dir, fname)
    plt.savefig(save_path, dpi=200)
    plt.close()
    print(f"График сохранен: {save_path}")

# -------------------------------------------------------------
# 3. Вызов генерации графиков
# -------------------------------------------------------------

# pic4.png: Пункт А (Kp_прес = 6.75 дБ)
# - СПМ входа опустилась в 2 раза (до 2.11)
# - Преселектор усиливает сильнее: пик спектра подскочил до 8.8 (шкала 12.0)
# - Напряжение УПЧ выросло до 7.79 мВ (шкала 12.0)
make_sensitivity_plot(
    psd_in_target=2.11e-14,
    psd_pres_target=8.80e-14,
    u_upch_target=7.79e-3,
    fname="pic4.png",
    title_text=r"$K_{p\mathrm{\ прес}} = 6{,}75\ \mathrm{дБ}$",
    ymax_pres=12.0,
    ymax_upch=12.0
)

# pic5.png: Пункт Б (K_ш ПЧ = 5.25 дБ)
# - СПМ входа опустилась в 2 раза (до 2.12)
# - Преселектор вернулся на базовый режим: пик номинальный 5.4 (шкала 8.0)
# - Напряжение УПЧ заметно упало до 4.57 мВ (полка спектра около 4.8)
make_sensitivity_plot(
    psd_in_target=2.12e-14,
    psd_pres_target=5.40e-14,
    u_upch_target=4.57e-3,
    fname="pic5.png",
    title_text=r"$K_{\mathrm{ш\ ПЧ}} = 5{,}25\ \mathrm{дБ}$",
    ymax_pres=8.0,
    ymax_upch=10.0
)

# pic6.png: Пункт В (L_фидера = 0 дБ, без потерь)
# - СПМ входа снизилась примерно на 20% (до 3.39)
# - Преселектор и УПЧ остаются практически на уровне исходного режима (1.97e-13 Вт и 6.53 мВ)
make_sensitivity_plot(
    psd_in_target=3.39e-14,
    psd_pres_target=5.50e-14,
    u_upch_target=6.53e-3,
    fname="pic6.png",
    title_text=r"$L_\Phi = 0\ \mathrm{дБ}$",
    ymax_pres=8.0,
    ymax_upch=10.0
)
