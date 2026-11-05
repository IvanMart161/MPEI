import os
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

fs = 200.0e6
N = 1_920_000
t = np.arange(N) / fs

idx_10us = np.where((t >= 0) & (t <= 10.0e-6))[0]
t_us = t[idx_10us] * 1e6

freqs = np.fft.rfftfreq(N, 1.0 / fs)
x_p = 2.0 * (freqs - 50.0e6) / 3.0e6
H_pres = 1.0 / np.sqrt(1.0 + x_p**4)

mask_upch = (freqs >= 9.5e6) & (freqs <= 10.5e6)
H_upch = np.zeros_like(freqs)
H_upch[mask_upch] = 1.0

seg_len = 8192
n_segs = N // seg_len
win = np.hanning(seg_len)
scale_psd = 2.0 / (fs * np.sum(win**2))

def get_psd(sig):
    reshaped = sig[:n_segs * seg_len].reshape(n_segs, seg_len)
    fft_res = np.fft.rfft(reshaped * win, axis=1)
    return np.mean(np.abs(fft_res)**2 * scale_psd, axis=0)

f_axis = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e6

def generate_bvch_plot(t_ant, psd_in_target, psd_pres_target, u_upch_target, fname, title_ant_desc):
    np.random.seed(42)

    # 1) Входной шум (масштабируется под целевую СПМ)
    scale_factor = np.sqrt(psd_in_target / 4.26e-14)
    w_in = np.random.normal(0, 0.45e-5 * scale_factor, N)

    # 2) Шум на выходе преселектора
    w_pres_raw = np.fft.irfft(np.fft.rfft(np.random.normal(0, 1.0, N)) * H_pres, n=N)
    target_pres_rms = 0.28e-6 * scale_factor
    w_pres = w_pres_raw * (target_pres_rms / np.std(w_pres_raw))

    # 3) Шум на выходе УПЧ
    w_upch_raw = np.fft.irfft(np.fft.rfft(np.random.normal(0, 1.0, N)) * H_upch, n=N)
    w_upch = w_upch_raw * (u_upch_target / np.std(w_upch_raw))

    # Спектры
    psd_in_curve = (get_psd(w_in) / 50.0) * 1e6
    psd_in_curve = psd_in_curve * (psd_in_target / np.mean(psd_in_curve[(f_axis >= 10) & (f_axis <= 90)]))

    psd_pres_raw = (get_psd(w_pres) / 50.0) * 1e6
    psd_pres_curve = psd_pres_raw * (psd_pres_target / np.max(psd_pres_raw))

    psd_upch_raw = np.sqrt(get_psd(w_upch) * 1e6)
    target_upch_dens = 6.6e-3 * scale_factor
    psd_upch_curve = psd_upch_raw * (target_upch_dens / np.max(psd_upch_raw))

    fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))

    # [0, 0] Входной шум
    axs[0, 0].plot(t_us, w_in[idx_10us] * 1e5, color='#0011dd', lw=0.35)
    axs[0, 0].set_title(f"Шум на входе БВЧ ({title_ant_desc})", fontsize=9)
    axs[0, 0].set_xlabel("t, мкс", fontsize=8)
    axs[0, 0].set_xlim(0, 10)
    ylim_in = 2.2 * scale_factor
    axs[0, 0].set_ylim(-ylim_in, ylim_in)
    axs[0, 0].text(-0.02, 1.03, r"$\times 10^{-5}$", transform=axs[0, 0].transAxes, fontsize=8)
    axs[0, 0].grid(True, linestyle='--', color='lightgray')

    # [0, 1] Выход преселектора
    axs[0, 1].plot(t_us, w_pres[idx_10us] * 1e6, color='#0011dd', lw=0.45)
    axs[0, 1].set_title("Шум на выходе преселектора", fontsize=9)
    axs[0, 1].set_xlabel("t, мкс", fontsize=8)
    axs[0, 1].set_xlim(0, 10)
    ylim_pres = 1.2 * scale_factor
    axs[0, 1].set_ylim(-ylim_pres, ylim_pres)
    axs[0, 1].text(-0.02, 1.03, r"$\times 10^{-6}$", transform=axs[0, 1].transAxes, fontsize=8)
    axs[0, 1].grid(True, linestyle='--', color='lightgray')

    # [0, 2] Выход УПЧ
    axs[0, 2].plot(t_us, w_upch[idx_10us], color='#0011dd', lw=0.5)
    axs[0, 2].set_title("Шум на выходе УПЧ", fontsize=9)
    axs[0, 2].set_xlabel("t, мкс", fontsize=8)
    axs[0, 2].set_xlim(0, 10)
    ylim_upch = 0.015 * scale_factor
    axs[0, 2].set_ylim(-ylim_upch, ylim_upch)
    axs[0, 2].grid(True, linestyle='--', color='lightgray')

    # [1, 0] СПМ приведенного шума
    axs[1, 0].plot(f_axis, psd_in_curve * 1e14, color='#cc0000', lw=1.1)
    axs[1, 0].set_title("СПМ приведённого шума", fontsize=9)
    axs[1, 0].set_xlabel("f, МГц", fontsize=8)
    axs[1, 0].set_ylabel("Спектральная плотность (Вт/МГц)", fontsize=8)
    axs[1, 0].set_xlim(0, 100)
    ylim_psd_in = 12.0 if scale_factor > 1.2 else 6.0
    axs[1, 0].set_ylim(0, ylim_psd_in)
    axs[1, 0].text(-0.02, 1.03, r"$\times 10^{-14}$", transform=axs[1, 0].transAxes, fontsize=8)
    axs[1, 0].grid(True, linestyle='--', color='lightgray')

    # [1, 1] СПМ преселектора
    axs[1, 1].plot(f_axis, psd_pres_curve * 1e14, color='#cc0000', lw=1.1)
    axs[1, 1].set_title("СПМ шума на выходе преселектора", fontsize=9)
    axs[1, 1].set_xlabel("f, МГц", fontsize=8)
    axs[1, 1].set_xlim(40, 60)
    ylim_psd_pres = 15.0 if scale_factor > 1.2 else 8.0
    axs[1, 1].set_ylim(0, ylim_psd_pres)
    axs[1, 1].text(-0.02, 1.03, r"$\times 10^{-14}$", transform=axs[1, 1].transAxes, fontsize=8)
    axs[1, 1].grid(True, linestyle='--', color='lightgray')

    # [1, 2] СП напряжения УПЧ
    axs[1, 2].plot(f_axis, psd_upch_curve * 1e3, color='#cc0000', lw=1.1)
    axs[1, 2].set_title("СП напряжения шума на выходе УПЧ", fontsize=9)
    axs[1, 2].set_xlabel("f, МГц", fontsize=8)
    axs[1, 2].set_ylabel(r"Спектральная плотность (В/$\sqrt{\mathrm{МГц}}$)", fontsize=8)
    axs[1, 2].set_xlim(8.0, 12.0)
    ylim_psd_upch = 15.0 if scale_factor > 1.2 else 10.0
    axs[1, 2].set_ylim(0, ylim_psd_upch)
    axs[1, 2].text(-0.02, 1.03, r"$\times 10^{-3}$", transform=axs[1, 2].transAxes, fontsize=8)
    axs[1, 2].grid(True, linestyle='--', color='lightgray')

    plt.tight_layout()
    save_path = os.path.join(output_dir, fname)
    plt.savefig(save_path, dpi=200)
    plt.close()
    print(f"График сохранён: {save_path}")

# Генерация pic2.png (T_A = 293 К, исходная СПМ 4.24e-14 Вт/МГц)
generate_bvch_plot(
    t_ant=293.0,
    psd_in_target=4.24e-14,
    psd_pres_target=5.47e-14,
    u_upch_target=6.35e-3,
    fname="pic2.png",
    title_ant_desc=r"$T_A = 293\ \mathrm{K}$"
)

# Генерация pic3.png (T_A = 3371 К, удвоенная СПМ 8.48e-14 Вт/МГц)
generate_bvch_plot(
    t_ant=3371.0,
    psd_in_target=8.48e-14,
    psd_pres_target=10.94e-14,
    u_upch_target=8.98e-3,
    fname="pic3.png",
    title_ant_desc=r"$T_A = 3371\ \mathrm{K}$"
)
