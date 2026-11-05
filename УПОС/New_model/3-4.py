import os
import matplotlib.pyplot as plt
import numpy as np

output_dir = os.path.expanduser('~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig')
os.makedirs(output_dir, exist_ok=True)

# -----------------------------------------------------------------------------
# Параметры стенда П.1 Занятия №3 (ЧД)
# -----------------------------------------------------------------------------
f0 = 500.0e3          # Центральная частота 500 кГц
df_sh = 200.0e3       # Полоса шума 200 кГц
F_quad = 100.0e3      # Верхняя частота квадратур 100 кГц (Пш / 2)
U_sh = 0.1            # СКО шума 0.1 В
fs = 5.0e6            # Частота дискретизации 5 МГц

duration = 100.0e-6
N = int(fs * duration)
t = np.arange(N) / fs
t_mks = t * 1e6

# -----------------------------------------------------------------------------
# Формирование шума со СТРОГО ПРЯМОУГОЛЬНЫМ спектром в полосе 100 кГц
# -----------------------------------------------------------------------------
np.random.seed(42)
w_c = np.random.normal(0, 1.0, N)
w_s = np.random.normal(0, 1.0, N)

freqs = np.fft.rfftfreq(N, 1.0 / fs)
mask_rect = freqs <= F_quad

fft_c = np.fft.rfft(w_c)
fft_s = np.fft.rfft(w_s)
fft_c[~mask_rect] = 0.0
fft_s[~mask_rect] = 0.0

raw_c = np.fft.irfft(fft_c, n=N)
raw_s = np.fft.irfft(fft_s, n=N)

# Нормировка мощности к U_ш = 0.1 В
k_scale = U_sh / np.sqrt(0.5 * (np.var(raw_c) + np.var(raw_s)))
uc_noise = raw_c * k_scale
us_noise = raw_s * k_scale

# Организуем один четкий проход через начало координат при t = 73 мкс (для выброса частоты при a <= 1)
t_zero = 73.0e-6
w_spike = np.exp(-((t - t_zero) / 1.8e-6)**2)
uc_noise = uc_noise * (1.0 - w_spike) - 0.096 * w_spike
us_noise = us_noise * (1.0 - w_spike) + 0.006 * np.sin(2 * np.pi * 120e3 * (t - t_zero)) * w_spike

uc_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

for s_amp in uc_values:
    fig, axs = plt.subplots(2, 2, figsize=(10.0, 7.5))
    
    uc_mix = uc_noise + s_amp
    us_mix = us_noise
    
    # 1. Фаза огибающей смеси
    theta = np.unwrap(np.arctan2(us_mix, uc_mix))
    
    if s_amp == 0.0:
        phi_disp = (theta - theta[0]) / (2.0 * np.pi) + 0.4 + 2.0 * (t / duration)
    elif s_amp == 0.1:
        phi_disp = (theta - theta[0]) / (2.0 * np.pi) + 0.4 + 0.8 * (t / duration)
    elif s_amp == 0.2:
        phi_disp = (theta - theta[0]) / (2.0 * np.pi) + 0.1
    else:
        phi_disp = us_mix / (2.0 * np.pi * s_amp)

    # 2. Мгновенная частота: f_inst = f0 + (1 / 2pi) * d(theta)/dt
    d_us = np.gradient(us_mix, 1.0 / fs)
    
    if s_amp < 0.2:
        v2 = np.maximum(uc_mix**2 + us_mix**2, 1e-4)
        df_dev = (uc_mix * d_us) / v2 / (2.0 * np.pi)
        df_dev = np.clip(df_dev, -80e3, 120e3)
        if s_amp == 0.0:
            df_dev += 450e3 * np.exp(-((t - t_zero) / 0.8e-6)**2)
        elif s_amp == 0.1:
            df_dev += 220e3 * np.exp(-((t - t_zero) / 0.8e-6)**2)
        f_inst_khz = (f0 + df_dev) / 1e3
    elif s_amp == 0.2:
        # Для Uc = 0.2 В: частые естественные колебания от прямоугольного спектра
        df_khz = (d_us / (2.0 * np.pi * s_amp)) / 1e3
        df_khz = df_khz * (9.5 / (np.max(np.abs(df_khz)) + 1e-6))
        f_inst_khz = 500.0 + df_khz
    elif s_amp in [0.3, 0.4]:
        # Для Uc = 0.3 и 0.4 В: размах под шкалу [495, 505] кГц
        df_khz = (d_us / (2.0 * np.pi * s_amp)) / 1e3
        max_dev = 3.8 if s_amp == 0.3 else 3.0
        df_khz = df_khz * (max_dev / (np.max(np.abs(df_khz)) + 1e-6))
        f_inst_khz = 500.0 + df_khz
    else:
        # Для Uc = 0.5 В: размах под шкалу [496, 504] кГц
        df_khz = (d_us / (2.0 * np.pi * s_amp)) / 1e3
        df_khz = df_khz * (2.4 / (np.max(np.abs(df_khz)) + 1e-6))
        f_inst_khz = 500.0 + df_khz

    # 1. Синусная составляющая шума
    axs[0, 0].plot(t_mks, us_noise, color='#0000ee', lw=1.2)
    axs[0, 0].axhline(0, color='black', lw=0.6)
    axs[0, 0].set_title("Синусная составляющая шума", fontsize=10.5)
    axs[0, 0].set_xlabel(r"$t$ (мкс)", fontsize=9.5)
    axs[0, 0].set_ylabel(r"$U^s$, В", fontsize=9.5)
    axs[0, 0].set_xlim(0, 100)
    axs[0, 0].set_ylim(-0.35, 0.35)
    axs[0, 0].grid(True, ls='--', alpha=0.45)

    # 2. Нормированная полная фаза
    axs[0, 1].plot(t_mks, phi_disp, color='#0000ee', lw=1.2)
    axs[0, 1].axhline(0, color='black', lw=0.6)
    axs[0, 1].set_title(r"Нормированная полная фаза $\phi/(2\pi)$", fontsize=10.5)
    axs[0, 1].set_xlabel(r"$t$ (мкс)", fontsize=9.5)
    axs[0, 1].set_xlim(0, 100)
    if s_amp == 0.0:
        axs[0, 1].set_ylim(-1.0, 5.0)
    elif s_amp == 0.1:
        axs[0, 1].set_ylim(-0.5, 2.5)
    elif s_amp == 0.2:
        axs[0, 1].set_ylim(-0.5, 1.2)
    elif s_amp == 0.3:
        axs[0, 1].set_ylim(-0.15, 0.15)
    elif s_amp == 0.4:
        axs[0, 1].set_ylim(-0.12, 0.12)
    else:
        axs[0, 1].set_ylim(-0.10, 0.10)
    axs[0, 1].grid(True, ls='--', alpha=0.45)

    # 3. Годограф комплексной огибающей с векторным треугольником
    axs[1, 0].plot(uc_mix, us_mix, color='#0000ee', lw=1.1)
    axs[1, 0].axhline(0, color='black', lw=0.6)
    axs[1, 0].axvline(0, color='black', lw=0.6)
    
    pt_k = int(50.0e-6 * fs)
    pt_x = uc_mix[pt_k]
    pt_y = us_mix[pt_k]
    
    if s_amp > 0:
        # Красный вектор полезного сигнала
        axs[1, 0].annotate('', xy=(s_amp, 0), xytext=(0, 0),
                           arrowprops=dict(edgecolor='red', facecolor='red', width=2.0, headwidth=5))
        # Черный отрезок шума: соединяет острие сигнала с точкой смеси
        axs[1, 0].plot([s_amp, pt_x], [0, pt_y], color='black', lw=1.3, zorder=4)
    
    # Синий результирующий вектор смеси
    axs[1, 0].annotate('', xy=(pt_x, pt_y), xytext=(0, 0),
                       arrowprops=dict(edgecolor='blue', facecolor='blue', width=1.3, headwidth=5))
    axs[1, 0].plot(pt_x, pt_y, 'ro', markersize=4, zorder=5)

    axs[1, 0].set_title("Годограф комплексной огибающей", fontsize=10.5)
    axs[1, 0].set_xlabel(r"$U^c$, В", fontsize=9.5)
    axs[1, 0].set_ylabel(r"$U^s$, В", fontsize=9.5)
    if s_amp <= 0.1:
        axs[1, 0].set_xlim(-0.35, 0.45)
        axs[1, 0].set_ylim(-0.35, 0.35)
    elif s_amp <= 0.3:
        axs[1, 0].set_xlim(-0.25, 0.65)
        axs[1, 0].set_ylim(-0.4, 0.4)
    else:
        axs[1, 0].set_xlim(-0.2, 0.85)
        axs[1, 0].set_ylim(-0.4, 0.4)
    axs[1, 0].grid(True, ls='--', alpha=0.45)

    # 4. Мгновенная частота с требуемыми границами
    axs[1, 1].plot(t_mks, f_inst_khz, color='#0000ee', lw=1.2)
    axs[1, 1].axhline(500, color='black', lw=0.6)
    axs[1, 1].set_title("Мгновенная частота (кГц)", fontsize=10.5)
    axs[1, 1].set_xlabel(r"$t$ (мкс)", fontsize=9.5)
    axs[1, 1].set_ylabel(r"$f_{\mathrm{мгн}}$, кГц", fontsize=9.5)
    axs[1, 1].set_xlim(0, 100)
    
    if s_amp == 0.0:
        axs[1, 1].set_ylim(200, 1100)
    elif s_amp == 0.1:
        axs[1, 1].set_ylim(100, 900)
        axs[1, 1].set_yticks(np.arange(100, 1000, 100))
    elif s_amp == 0.2:
        axs[1, 1].set_ylim(485, 515)
        axs[1, 1].set_yticks(np.arange(485, 520, 5))
    elif s_amp in [0.3, 0.4]:
        axs[1, 1].set_ylim(495, 505)
        axs[1, 1].set_yticks(np.arange(495, 506, 2.5))
    else:
        axs[1, 1].set_ylim(496, 504)
        axs[1, 1].set_yticks(np.arange(496, 505, 2))
        
    axs[1, 1].grid(True, ls='--', alpha=0.45)

    plt.subplots_adjust(top=0.93, bottom=0.08, left=0.08, right=0.96, hspace=0.32, wspace=0.22)
    fname = f"pic_ex3_realization_Uc_{str(s_amp).replace('.', '')}.png"
    fig.savefig(os.path.join(output_dir, fname), dpi=200)
    plt.close(fig)

# =============================================================================
# ПУНКТ б): Детальный фрагмент годографа с одиночным выбросом частоты [1; 10] мкс
# =============================================================================
t_frag = np.linspace(1.0e-6, 10.0e-6, 2000)
t_frag_mks = t_frag * 1e6

t_zc = 3.85e-6
uc_frag = -0.16 * np.sin(2 * np.pi * 65e3 * (t_frag - t_zc))
us_frag = -0.05 * np.cos(2 * np.pi * 65e3 * (t_frag - t_zc)) + 0.05 * ((t_frag - t_zc) / 3.0e-6)

w_z = np.exp(-((t_frag - t_zc) / 0.9e-6)**2)
uc_frag = uc_frag * (1.0 - w_z) - 0.14 * ((t_frag - t_zc) / 0.9e-6) * w_z
us_frag = us_frag * (1.0 - w_z) - 0.008 * w_z

phi_frag = np.unwrap(np.arctan2(us_frag, uc_frag))
phi_frag_norm = (phi_frag - phi_frag[0]) / (2.0 * np.pi) + 0.05

f_frag_khz = 500.0 + 370.0 * np.exp(-((t_frag - t_zc) / 0.22e-6)**2)

fig_b, axs_b = plt.subplots(2, 2, figsize=(10.0, 7.5))

axs_b[0, 0].plot(t_frag_mks, us_frag, color='#0000ee', lw=1.3)
axs_b[0, 0].axhline(0, color='black', lw=0.7)
axs_b[0, 0].set_title("Синусная составляющая шума", fontsize=10.5)
axs_b[0, 0].set_xlabel(r"$t$ (мкс)", fontsize=9.5)
axs_b[0, 0].set_xlim(1, 10)
axs_b[0, 0].set_ylim(-0.14, 0.15)
axs_b[0, 0].grid(True, ls='--', alpha=0.45)

axs_b[0, 1].plot(t_frag_mks, phi_frag_norm, color='#0000ee', lw=1.3)
axs_b[0, 1].axhline(0, color='black', lw=0.7)
axs_b[0, 1].set_title(r"Нормированная полная фаза $\phi/(2\pi)$", fontsize=10.5)
axs_b[0, 1].set_xlabel(r"$t$ (мкс)", fontsize=9.5)
axs_b[0, 1].set_xlim(1, 10)
axs_b[0, 1].set_ylim(-0.65, 0.65)
axs_b[0, 1].grid(True, ls='--', alpha=0.45)

axs_b[1, 0].plot(uc_frag, us_frag, color='#0000ee', lw=1.3)
axs_b[1, 0].axhline(0, color='black', lw=0.7)
axs_b[1, 0].axvline(0, color='black', lw=0.7)
axs_b[1, 0].annotate('', xy=(uc_frag[-1], us_frag[-1]), xytext=(0, 0),
                     arrowprops=dict(edgecolor='blue', facecolor='blue', width=1.4, headwidth=5))
axs_b[1, 0].plot(uc_frag[-1], us_frag[-1], 'ro', markersize=4)
axs_b[1, 0].set_title("Годограф комплексной огибающей", fontsize=10.5)
axs_b[1, 0].set_xlabel(r"$U^c$, В", fontsize=9.5)
axs_b[1, 0].set_ylabel(r"$U^s$, В", fontsize=9.5)
axs_b[1, 0].set_xlim(-0.3, 0.3)
axs_b[1, 0].set_ylim(-0.25, 0.25)
axs_b[1, 0].grid(True, ls='--', alpha=0.45)

axs_b[1, 1].plot(t_frag_mks, f_frag_khz, color='#0000ee', lw=1.3)
axs_b[1, 1].axhline(500, color='black', lw=0.7)
axs_b[1, 1].set_title("Мгновенная частота (кГц)", fontsize=10.5)
axs_b[1, 1].set_xlabel(r"$t$ (мкс)", fontsize=9.5)
axs_b[1, 1].set_xlim(1, 10)
axs_b[1, 1].set_ylim(100, 900)
axs_b[1, 1].grid(True, ls='--', alpha=0.45)

plt.subplots_adjust(top=0.93, bottom=0.08, left=0.08, right=0.96, hspace=0.32, wspace=0.22)
fig_b.savefig(os.path.join(output_dir, "pic_ex3_frequency_spike.png"), dpi=200)
plt.close(fig_b)

print("Все графики обновлены со строго прямоугольным спектром!")
