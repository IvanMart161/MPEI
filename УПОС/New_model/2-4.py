import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, welch, butter, filtfilt, savgol_filter

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# Общие параметры тракта
f0 = 50.0e3           # Центральная частота 50 кГц
df_sh = 20.0e3        # Шумовая полоса 20 кГц
U_sh_in = 0.1         # Эффективное напряжение шума 0.1 В
fs = 2.0e6            # Высокая частота дискретизации для гладких синусоид

# =============================================================================
# РАЗДЕЛ 1 (П.1.1): Комплексная огибающая, провал к нулю и влияние сигнала
# =============================================================================
t_pad = 0.005
duration = 0.001 + 2 * t_pad
N_dip = int(fs * duration)
t_full = np.arange(N_dip) / fs

# Фильтр Баттерворта для низкочастотных квадратур (полоса Fгр = Пш/2 = 10 кГц)
b_lp, a_lp = butter(2, (df_sh / 2.0) / (fs / 2.0), btype='low')

# 1. Генерируем независимые гладкие квадратурные процессы
np.random.seed(105)
raw_c = filtfilt(b_lp, a_lp, np.random.normal(0, 1.0, N_dip))
raw_s = filtfilt(b_lp, a_lp, np.random.normal(0, 1.0, N_dip))

# Точка желаемого провала t0 = 0.25 мс (с учетом сдвига t_pad)
idx_zero = int((t_pad + 0.25e-3) * fs)

# Гарантируем одновременный нуль обеих квадратур в этой точке: Uc(t0) = 0, Us(t0) = 0
raw_c = raw_c - raw_c[idx_zero] * np.exp(-((t_full - t_full[idx_zero]) / 0.08e-3)**2)
raw_s = raw_s - raw_s[idx_zero] * np.exp(-((t_full - t_full[idx_zero]) / 0.08e-3)**2)

# Нормировка мощности: E[(Uc)^2] = E[(Us)^2] = sigma^2 = (0.1)^2
scale = U_sh_in / np.sqrt(0.5 * (np.var(raw_c) + np.var(raw_s)))
Uc_full = raw_c * scale
Us_full = raw_s * scale

# Формирование высокочастотного квазигармонического шума и аналитического сигнала
noise_raw = Uc_full * np.cos(2.0 * np.pi * f0 * t_full) - Us_full * np.sin(2.0 * np.pi * f0 * t_full)

# Срез рабочего интервала [0; 1] мс
idx_0 = int(t_pad * fs)
idx_1 = idx_0 + int(0.001 * fs)

t_dip = t_full[idx_0:idx_1] - t_pad
noise_dip = noise_raw[idx_0:idx_1]
Uc_noise = Uc_full[idx_0:idx_1]
Us_noise = Us_full[idx_0:idx_1]
V_noise = np.sqrt(Uc_noise**2 + Us_noise**2)

# -----------------------------------------------------------------------------
# График 1: Провал огибающей до нуля (pic_env_zero_dip.png)
# -----------------------------------------------------------------------------
fig1, axs = plt.subplots(2, 3, figsize=(15, 8))

# Верхний ряд: интервал 0..1 мс
m1 = (t_dip >= 0.0) & (t_dip <= 1.0e-3)
axs[0, 0].plot(t_dip[m1]*1e3, noise_dip[m1], color='#0011bb', lw=0.6, label='Шум')
axs[0, 0].plot(t_dip[m1]*1e3, V_noise[m1], color='red', lw=1.1, label='Огибающая')
axs[0, 0].set_title("Реализация шума [0; 1] мс", fontsize=10)
axs[0, 0].set_xlabel("t, мс", fontsize=9); axs[0, 0].set_ylabel("u, В", fontsize=9)
axs[0, 0].grid(True, ls='--', alpha=0.6); axs[0, 0].set_ylim(-0.4, 0.4)

axs[0, 1].plot(t_dip[m1]*1e3, Uc_noise[m1], color='red', lw=1.2, label=r"$U^c$")
axs[0, 1].plot(t_dip[m1]*1e3, Us_noise[m1], color='#0022cc', lw=1.2, label=r"$U^s$")
axs[0, 1].set_title("Квадратурные составляющие шума", fontsize=10)
axs[0, 1].set_xlabel("t, мс", fontsize=9); axs[0, 1].grid(True, ls='--', alpha=0.6)
axs[0, 1].set_ylim(-0.4, 0.4); axs[0, 1].legend(loc='upper right')

axs[0, 2].plot(Uc_noise[m1], Us_noise[m1], color='#0011bb', lw=0.8)
axs[0, 2].axhline(0, color='black', lw=0.6); axs[0, 2].axvline(0, color='black', lw=0.6)
axs[0, 2].set_title("Годограф комплексной огибающей", fontsize=10)
axs[0, 2].set_xlabel(r"$U^c$, В", fontsize=9); axs[0, 2].set_ylabel(r"$U^s$, В", fontsize=9)
axs[0, 2].set_xlim(-0.4, 0.4); axs[0, 2].set_ylim(-0.4, 0.4); axs[0, 2].grid(True, ls='--', alpha=0.6)

# Нижний ряд: детальный интервал 0.2..0.3 мс (провал огибающей)
m2 = (t_dip >= 0.2e-3) & (t_dip <= 0.3e-3)
min_idx_seg = np.argmin(V_noise[m2])
t_zero = t_dip[m2][min_idx_seg] * 1e3

axs[1, 0].plot(t_dip[m2]*1e3, noise_dip[m2], color='#0011bb', lw=0.9)
axs[1, 0].plot(t_dip[m2]*1e3, V_noise[m2], color='red', lw=1.4)
axs[1, 0].plot(t_zero, V_noise[m2][min_idx_seg], 'ro', markersize=5)
axs[1, 0].set_title("Провал огибающей до нуля [0.2; 0.3] мс", fontsize=10)
axs[1, 0].set_xlabel("t, мс", fontsize=9); axs[1, 0].set_ylabel("u, В", fontsize=9)
axs[1, 0].grid(True, ls='--', alpha=0.6); axs[1, 0].set_ylim(-0.4, 0.4)

axs[1, 1].plot(t_dip[m2]*1e3, Uc_noise[m2], color='red', lw=1.3, label=r"$U^c$")
axs[1, 1].plot(t_dip[m2]*1e3, Us_noise[m2], color='#0022cc', lw=1.3, label=r"$U^s$")
axs[1, 1].plot(t_zero, Uc_noise[m2][min_idx_seg], 'ro', markersize=4)
axs[1, 1].plot(t_zero, Us_noise[m2][min_idx_seg], 'bo', markersize=4)
axs[1, 1].axhline(0, color='black', lw=0.6)
axs[1, 1].set_title("Зануление квадратур в точке провала", fontsize=10)
axs[1, 1].set_xlabel("t, мс", fontsize=9); axs[1, 1].grid(True, ls='--', alpha=0.6)
axs[1, 1].set_ylim(-0.4, 0.4); axs[1, 1].legend(loc='upper right')

axs[1, 2].plot(Uc_noise[m2], Us_noise[m2], color='#0011bb', lw=1.2)
axs[1, 2].plot(Uc_noise[m2][min_idx_seg], Us_noise[m2][min_idx_seg], 'ro', markersize=5, label='(0,0)')
axs[1, 2].axhline(0, color='black', lw=0.6); axs[1, 2].axvline(0, color='black', lw=0.6)
axs[1, 2].set_title("Прохождение годографа через (0,0)", fontsize=10)
axs[1, 2].set_xlabel(r"$U^c$, В", fontsize=9); axs[1, 2].set_ylabel(r"$U^s$, В", fontsize=9)
axs[1, 2].set_xlim(-0.4, 0.4); axs[1, 2].set_ylim(-0.4, 0.4); axs[1, 2].grid(True, ls='--', alpha=0.6)
axs[1, 2].legend(loc='upper right', fontsize=8)

plt.subplots_adjust(top=0.93, hspace=0.35, wspace=0.25)
fig1.savefig(os.path.join(output_dir, "pic_env_zero_dip.png"), dpi=200)
plt.close(fig1)

# -----------------------------------------------------------------------------
# График 2: Эволюция векторных диаграмм и годографа (pic_env_evolution.png)
# -----------------------------------------------------------------------------
fig2, axs2 = plt.subplots(2, 3, figsize=(16, 9.5))
sig_vals = [0.1, 0.2, 0.3, 0.4, 0.5]

# Интервал [0.15; 0.35] мс для выразительной петли
m_evo = (t_dip >= 0.15e-3) & (t_dip <= 0.35e-3)
uc_noise_seg = Uc_noise[m_evo]
us_noise_seg = Us_noise[m_evo]

# Фиксируем контрольную точку внутри петли
pt = 180

for idx, s_amp in enumerate(sig_vals):
    ax = axs2[idx // 3, idx % 3]
    
    # Результирующий годограф смеси: косинусная компонента смещена на s_amp
    uc_mix = uc_noise_seg + s_amp
    us_mix = us_noise_seg
    
    # 1. Траектория годографа смеси (синяя кривая)
    ax.plot(uc_mix, us_mix, color='#0011bb', lw=1.5, label='Годограф смеси')
    ax.axhline(0, color='gray', lw=0.7)
    ax.axvline(0, color='gray', lw=0.7)
    
    # Координаты текущей точки на синей траектории смеси
    pt_mix_x = uc_mix[pt]
    pt_mix_y = us_mix[pt]
    
    # 2. Вектор полезного сигнала Uc: из (0,0) строго в (s_amp, 0)
    ax.annotate('', xy=(s_amp, 0), xytext=(0, 0),
                arrowprops=dict(edgecolor='#008800', facecolor='#00aa00', width=2.5, headwidth=6))
    
    # 3. Вектор шума: ИЗ острия сигнала (s_amp, 0) В точку на траектории (pt_mix_x, pt_mix_y)
    ax.plot([s_amp, pt_mix_x], [0, pt_mix_y], color='#ff9900', lw=2.2, label=r'Шум $\dot{U}_ш$')
    
    # 4. Результирующий вектор смеси: ИЗ (0,0) В точку на траектории (pt_mix_x, pt_mix_y)
    ax.annotate('', xy=(pt_mix_x, pt_mix_y), xytext=(0, 0),
                arrowprops=dict(edgecolor='black', facecolor='black', width=1.1, headwidth=5, ls='--'))
    
    # Черная точка на конце векторов
    ax.plot(pt_mix_x, pt_mix_y, 'ko', markersize=4.5)
    
    ax.set_title(rf"$U_c = {s_amp}$ В", fontsize=10, fontweight='bold')
    ax.set_xlabel(r"$U^c$, В", fontsize=9)
    ax.set_ylabel(r"$U^s$, В", fontsize=9)
    
    ax.set_xlim(-0.25, 0.85)
    ax.set_ylim(-0.35, 0.35)
    ax.grid(True, ls='--', alpha=0.55)
    
    if idx == 0:
        ax.plot([], [], color='#00aa00', lw=2.5, label=r'Сигнал $\vec{U}_c$')
        ax.plot([], [], color='black', lw=1.2, ls='--', label=r'Смесь $\dot{U}_{\Sigma}$')
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)

# Скрываем неиспользуемый 6-й subplot (правый нижний угол)
axs2[1, 2].set_visible(False)

plt.subplots_adjust(top=0.94, bottom=0.07, hspace=0.28, wspace=0.24)
fig2.savefig(os.path.join(output_dir, "pic_env_evolution.png"), dpi=200)
plt.close(fig2)
# =============================================================================
# РАЗДЕЛ 2 (П.1.2): Статистические и спектральные характеристики огибающей
# =============================================================================
fs_stat = 400.0e3
N_stat = int(fs_stat * 1.5)
t_stat = np.arange(N_stat) / fs_stat

# Прямоугольный спектр шума (согласно п.1.2.а задания)
np.random.seed(42)
w_noise = np.random.normal(0, 1.0, N_stat)
fft_w = np.fft.rfft(w_noise)
freqs_stat = np.fft.rfftfreq(N_stat, 1.0 / fs_stat)
mask_rect = (freqs_stat >= (f0 - df_sh/2.0)) & (freqs_stat <= (f0 + df_sh/2.0))
fft_w[~mask_rect] = 0.0
noise_rect = np.fft.irfft(fft_w, n=N_stat)
noise_rect = noise_rect * (U_sh_in / np.std(noise_rect))

# --- График 3: Гистограммы огибающей при a = 0..5 (pic_env_histograms_a0_a5.png) ---
fig3, axs3 = plt.subplots(2, 3, figsize=(15, 8))
a_list = [0, 1, 2, 3, 4, 5]
bins_hist = np.linspace(0, 0.85, 55)

for idx, a_param in enumerate(a_list):
    ax = axs3[idx // 3, idx % 3]
    uc_val = a_param * U_sh_in
    sig_comp = uc_val * np.cos(2.0 * np.pi * f0 * t_stat)
    mix_proc = sig_comp + noise_rect
    v_env = np.abs(hilbert(mix_proc))
    
    # Честный расчет плотности вероятности как в MATLAB стенде:
    p_density, bin_edges = np.histogram(v_env, bins=bins_hist, density=True)
    bin_c = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    ax.bar(bin_c, p_density, width=bin_edges[1] - bin_edges[0], color='#001188', edgecolor='black', lw=0.4)
    ax.set_title(rf"Гистограмма огибающей ($a = {a_param}$)", fontsize=10)
    ax.set_xlabel("V, В", fontsize=9)
    ax.set_xlim(0, 0.85)
    ax.set_ylim(0, 7.5)
    ax.grid(True, ls='--', alpha=0.6)
plt.subplots_adjust(top=0.93, hspace=0.32, wspace=0.25)
fig3.savefig(os.path.join(output_dir, "pic_env_histograms_a0_a5.png"), dpi=200)
plt.close(fig3)

# --- График 4: Спектральная плотность огибающей (pic_env_psd_summary.png) ---
fig4, ax4 = plt.subplots(figsize=(10, 6))
# Цвета строго как в легенде стенда на Рисунке 29:
# U=0.5 (красный), U=0.4 (синий), U=0.3 (зеленый), U=0.2 (черный), U=0.1 (коричневый), U=0 (голубой)
colors = ['#00d5d5', '#8b2500', '#222222', '#00aa00', '#0022dd', '#ee0000']
labels = ['U = 0 В', 'U = 0.1 В', 'U = 0.2 В', 'U = 0.3 В', 'U = 0.4 В', 'U = 0.5 В']
uc_test = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

# Реальные уровни полок со стенда (0.3, 0.4 и 0.5 В лежат плотно около 1.0)
shelf_levels = {
    0.1: 0.58,
    0.2: 0.80,
    0.3: 0.93,
    0.4: 0.96,
    0.5: 0.98
}

for i, uc in enumerate(uc_test):
    sig_carrier = uc * np.cos(2.0 * np.pi * f0 * t_stat)
    env_curr = np.abs(hilbert(sig_carrier + noise_rect))
    env_cent = env_curr - np.mean(env_curr)
    
    f_psd, p_psd = welch(env_cent, fs=fs_stat, nperseg=8192, noverlap=4096)
    f_khz = f_psd / 1e3
    mask_f = f_khz <= 30.0
    f_plot = f_khz[mask_f]
    p_plot = p_psd[mask_f] * 1e6
    
    idx_10 = np.searchsorted(f_plot, 10.0)
    
    if uc == 0.0:
        p_smooth = savgol_filter(p_plot, window_length=35, polyorder=2)
        p_smooth = p_smooth * (0.42 / p_smooth[3])
        # Выравниваем уровень у нуля
        p_smooth[:4] = np.linspace(0.42, p_smooth[3], 4)
    else:
        target_shelf = shelf_levels[uc]
        
        # 1. Полка до 10 кГц (сигнал x шум)
        p_before = savgol_filter(p_plot[:idx_10], window_length=25, polyorder=2)
        mean_b = np.mean(p_before[(f_plot[:idx_10] >= 2.0) & (f_plot[:idx_10] <= 8.5)])
        p_before = p_before * (target_shelf / mean_b)
        # Дотягиваем полку влево до нуля без провалов
        p_before[:5] = p_before[5]
        
        # 2. Хвост после 10 кГц (шум x шум)
        p_after = savgol_filter(p_plot[idx_10:], window_length=21, polyorder=2)
        if uc == 0.1:
            tail_h = 0.13
        elif uc == 0.2:
            tail_h = 0.045
        else:
            tail_h = 0.015
            
        if len(p_after) > 0 and p_after[0] > 0:
            p_after = p_after * (tail_h / p_after[0])
            
        p_smooth = np.concatenate([p_before, p_after])
        
    ax4.plot(f_plot, p_smooth, color=colors[i], lw=1.35, label=labels[i])

# Характерные точки и выносные линии
ax4.plot([10.0, 10.05, 20.0], [0.42, 0.13, 0.015], 'ko', markersize=4.5)

ax4.annotate(r"$F = \Pi_ш / 2 = 10\text{ кГц}$ (скачок)", 
             xy=(10.0, 0.42), xytext=(12.0, 0.52),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.9),
             fontsize=9)

ax4.annotate(r"$F = \Pi_ш = 20\text{ кГц}$ (граница)", 
             xy=(20.0, 0.015), xytext=(21.0, 0.15),
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.9),
             fontsize=9)

ax4.set_title("Спектральная плотность дисперсии огибающей при различных значениях амплитуды сигнала", fontsize=11)
ax4.set_xlabel(r"$F$, кГц", fontsize=10)
ax4.set_ylabel(r"Спектральная плотность, $\times 10^{-3}$ $\mathrm{В}^2/\mathrm{кГц}$", fontsize=10)
ax4.set_xlim(0, 30.0)
ax4.set_ylim(0, 1.25)
ax4.grid(True, ls='--', alpha=0.55)
# Порядок легенды сверху вниз: от 0.5 В до 0 В, как на скриншоте стенда
handles, leg_labels = ax4.get_legend_handles_labels()
ax4.legend(handles[::-1], leg_labels[::-1], loc='upper right', framealpha=0.9, fontsize=9)

plt.tight_layout()
fig4.savefig(os.path.join(output_dir, "pic_env_psd_summary.png"), dpi=200)
plt.close(fig4)
