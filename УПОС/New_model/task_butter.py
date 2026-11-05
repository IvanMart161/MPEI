import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, butter, filtfilt

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------------------
# Параметры тракта (Бригада №8: фильтр Баттерворта)
# -------------------------------------------------------------
f0 = 50.0e3           # Центральная частота f0 = 50 кГц
df_sh = 20.0e3        # Шумовая полоса 20 кГц
U_sh_in = 0.1         # Эффективное напряжение шума Uш = 0.1 В
fs = 2.0e6            # Высокая частота дискретизации для гладких синусоид
duration = 0.003      # 3 мс
N = int(fs * duration)
t = np.arange(N) / fs

# Полосовой фильтр Баттерворта 2-го порядка (эквивалент ФСС стенда)
low = (f0 - df_sh / 2.0) / (fs / 2.0)
high = (f0 + df_sh / 2.0) / (fs / 2.0)
b_bw, a_bw = butter(2, [low, high], btype='bandpass')

# Подбор зерна для отчетливого провала огибающей к нулю около 0.25 мс
np.random.seed(108)
raw_noise = np.random.normal(0, 1.0, N)
noise_bp = filtfilt(b_bw, a_bw, raw_noise)
noise_bp = noise_bp * (U_sh_in / np.std(noise_bp))  # нормировка по СКО = 0.1 В

# Комплексная огибающая и квадратурные составляющие
analytic_signal = hilbert(noise_bp)
env_base = analytic_signal * np.exp(-1j * 2.0 * np.pi * f0 * t)
Uc = np.real(env_base)
Us = np.imag(env_base)
V_env = np.abs(analytic_signal)

def draw_stand_window(t_start, t_end, filename, title_suffix):
    mask = (t >= t_start) & (t <= t_end)
    t_plot = (t[mask] - t_start + (0.2e-3 if t_start > 0 else 0)) * 1e3
    
    # Исходная временная шкала для осей
    t_axis_ms = t[mask] * 1e3
    mix_seg = noise_bp[mask]
    env_seg = V_env[mask]
    uc_seg = Uc[mask]
    us_seg = Us[mask]

    fig = plt.figure(figsize=(13, 7.5))

    # 1. Реализация шума и огибающей (слева вверху)
    ax1 = fig.add_axes([0.08, 0.56, 0.40, 0.36])
    ax1.plot(t_axis_ms, mix_seg, color='#0011cc', lw=0.9, label='Шум')
    ax1.plot(t_axis_ms, env_seg, color='red', lw=1.2, label='Огибающая')
    ax1.plot(t_axis_ms[-1], env_seg[-1], 'ko', markersize=4)
    ax1.set_title(f"Реализация шума ({title_suffix})", fontsize=10)
    ax1.set_xlabel(r"$t$ (мс)", fontsize=9)
    ax1.set_xlim(t_axis_ms[0], t_axis_ms[-1])
    ax1.set_ylim(-0.4, 0.4)
    ax1.grid(True, ls='--', color='gray', alpha=0.6)

    # 2. Квадратурные составляющие (справа вверху)
    ax2 = fig.add_axes([0.55, 0.56, 0.40, 0.36])
    ax2.plot(t_axis_ms, uc_seg, color='red', lw=1.3, label=r"$U_{\text{ш}}^c$")
    ax2.plot(t_axis_ms, us_seg, color='#0011cc', lw=1.3, label=r"$U_{\text{ш}}^s$")
    ax2.plot(t_axis_ms[-1], uc_seg[-1], 'ro', markersize=4)
    ax2.plot(t_axis_ms[-1], us_seg[-1], 'bo', markersize=4)
    ax2.axhline(0, color='black', lw=0.6)
    ax2.set_title("Квадратурные составляющие шума", fontsize=10)
    ax2.set_xlabel(r"$t$ (мс)", fontsize=9)
    ax2.set_xlim(t_axis_ms[0], t_axis_ms[-1])
    ax2.set_ylim(-0.4, 0.4)
    ax2.grid(True, ls='--', color='gray', alpha=0.6)
    ax2.legend(loc='upper right', framealpha=0.9)

    # 3. Годограф комплексной огибающей (снизу по центру)
    ax3 = fig.add_axes([0.31, 0.08, 0.42, 0.40])
    ax3.plot(uc_seg, us_seg, color='#0011cc', lw=1.1)
    ax3.axhline(0, color='black', lw=0.7)
    ax3.axvline(0, color='black', lw=0.7)

    # Вектор в текущую точку конца реализации
    ax3.plot([0, uc_seg[-1]], [0, us_seg[-1]], color='#ffaa00', lw=2)
    ax3.plot(uc_seg[-1], us_seg[-1], 'ko', markersize=4)

    # Если интервал с провалом — подсвечиваем точку прохождения через нуль
    if t_start > 0:
        min_idx = np.argmin(env_seg)
        ax3.plot(uc_seg[min_idx], us_seg[min_idx], 'ro', markersize=5, label='Точка провала')
        ax3.legend(loc='upper right', fontsize=8)

    ax3.set_title("Годограф комплексной огибающей", fontsize=10)
    ax3.set_xlabel(r"$U_{\text{ш}}^c$, В", fontsize=9)
    ax3.set_ylabel(r"$U_{\text{ш}}^s$, В", fontsize=9)
    ax3.set_xlim(-0.4, 0.4)
    ax3.set_ylim(-0.4, 0.4)
    ax3.grid(True, ls='--', color='gray', alpha=0.6)

    save_path = os.path.join(output_dir, filename)
    plt.savefig(save_path, dpi=200)
    plt.close()
    print(f"Сохранён график: {save_path}")

# Построение для интервала [0; 1] мс
draw_stand_window(0.0, 1.0e-3, "pic_stand_butterworth_full.png", "[0; 1] мс")

# Построение для детального интервала [0.2; 0.3] мс (провал огибающей)
draw_stand_window(0.2e-3, 0.3e-3, "pic_stand_butterworth_dip.png", "[0.2; 0.3] мс")
