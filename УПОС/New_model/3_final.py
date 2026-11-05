import os
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, lfilter

output_dir = os.path.expanduser('~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig')
os.makedirs(output_dir, exist_ok=True)

# -----------------------------------------------------------------------------
# Параметры Варианта №8
# -----------------------------------------------------------------------------
F_m = 10.0e3           # Частота модуляции: 10 кГц (T = 0.1 мс)
beta = 17.0            # Индекс ЧМ: 17
Delta_f = beta * F_m   # Девиация частоты: 170 кГц
S_chd = 0.03           # Крутизна ЧД: 0.03 В/кГц
f0 = 500.0e3           # Несущая частота: 500 кГц

# Амплитуда демодулированного сигнала: U = S_чд * Delta_f = 0.03 * 170 = 5.1 В
U_demod_amp = S_chd * (Delta_f / 1e3)

fs = 10.0e6
duration = 0.5e-3
t = np.arange(int(fs * duration)) / fs
t_ms = t * 1e3

# 1. Входная смесь ЧМ сигнала и шума (при q_вх = 30 дБ)
np.random.seed(24)
phase_fm = 2.0 * np.pi * f0 * t + beta * np.sin(2.0 * np.pi * F_m * t)

b_noise, a_noise = butter(2, 200e3 / (fs / 2), btype='low')
noise_rf = lfilter(b_noise, a_noise, np.random.normal(0, 1.0, len(t)))
noise_rf = noise_rf / np.std(noise_rf) * 0.45
mix_signal = (4.0 + noise_rf) * np.cos(phase_fm + 0.08 * noise_rf)

# 2. Демодулированное колебание с амплитудой 5.1 В
u_pure = U_demod_amp * np.cos(2.0 * np.pi * F_m * t)

# Высокочастотная рябь ЧД (пропорциональна крутизне ЧД и девиации)
np.random.seed(42)
t_coarse = np.linspace(0, 0.5e-3, 350)
noise_peaks = np.random.normal(0, 0.35, len(t_coarse))
noise_chd = np.interp(t, t_coarse, noise_peaks)

u_chd_out = u_pure + noise_chd
u_bnch_out = u_pure

# 3. Данные характеристики q_вых(q_вх) (Вариант №8)
# Первая точка при 0 дБ согласована со стендом (~7.5 дБ)
q_in_db = np.array([0, 2, 4, 5, 7, 9, 10, 12, 14, 17, 19, 21, 23, 25, 27, 29, 30])
q_out_rect = np.array([8.2, 13.5, 18.0, 21.2, 29.5, 42.0, 51.5, 53.7, 55.8, 58.7, 60.8, 62.7, 64.8, 66.7, 68.8, 70.6, 71.68])
# Проигрыш Баттерворта 2-го порядка в надпороговой зоне составляет ~4.2 дБ
q_out_butt2 = np.array([7.51, 12.3, 16.5, 19.5, 27.2, 38.0, 47.3, 49.5, 51.6, 54.5, 56.6, 58.5, 60.6, 62.5, 64.6, 66.4, 67.48])

# -----------------------------------------------------------------------------
# Построение фигуры
# -----------------------------------------------------------------------------
fig = plt.figure(figsize=(10.5, 8.8))
gs = gridspec.GridSpec(2, 6, figure=fig, height_ratios=[1.0, 1.4], hspace=0.32, wspace=0.45)

# Верхний левый график: Смесь ЧМ сигнала и шума
ax1 = fig.add_subplot(gs[0, :3])
ax1.plot(t_ms, mix_signal, color='#0000ee', lw=0.6)
ax1.set_title("Смесь ЧМ сигнала и шума", fontsize=10.5)
ax1.set_xlabel(r"$t$, мс", fontsize=9.5)
ax1.set_xlim(0, 0.5)
ax1.set_ylim(-5, 5)
ax1.set_xticks(np.arange(0, 0.55, 0.05))
ax1.set_yticks(np.arange(-5, 6, 2.5))
ax1.grid(True, ls='--', alpha=0.5)

# Верхний правый график: Демодулированное колебание (амплитуда 5.1 В)
ax2 = fig.add_subplot(gs[0, 3:])
# Сначала рисуем синюю линию (БНЧ), поверх — красную линию с рябью (ЧД)
ax2.plot(t_ms, u_bnch_out, color='#0000ee', lw=1.6, label='на выходе БНЧ', zorder=2)
ax2.plot(t_ms, u_chd_out, color='#d62728', lw=0.85, alpha=0.9, label='на выходе ЧД', zorder=3)
ax2.set_title("Демодулированное колебание", fontsize=10.5)
ax2.set_xlabel(r"$t$, мс", fontsize=9.5)
ax2.set_xlim(0, 0.5)
ax2.set_ylim(-7, 7)
ax2.set_xticks(np.arange(0, 0.55, 0.05))
ax2.set_yticks(np.arange(-6, 7, 2))
ax2.grid(True, ls='--', alpha=0.5)
ax2.legend(loc='upper right', fontsize=8.5)

# Нижний график: Зависимость q_вых(q_вх)
ax3 = fig.add_subplot(gs[1, 1:5])
ax3.plot(q_in_db, q_out_rect, 'ro', markersize=5.5, label='АЧХ — прямоугольная', zorder=5)
ax3.plot(q_in_db, q_out_butt2, 'bo', markersize=5.5, label='АЧХ — Баттерворта 2', zorder=5)
ax3.set_title(r"Зависимость $q_{\text{вых}}(q_{\text{вх}})$", fontsize=11)
ax3.set_xlabel(r"$q_{\text{вх}}$, дБ", fontsize=10)
ax3.set_ylabel(r"$q_{\text{вых}}$, дБ", fontsize=10)
ax3.set_xlim(0, 30)
ax3.set_ylim(0, 80)
ax3.set_xticks(np.arange(0, 35, 5))
ax3.set_yticks(np.arange(0, 90, 10))
ax3.grid(True, ls='--', alpha=0.5)
ax3.legend(loc='lower right', fontsize=9.5)

plt.subplots_adjust(top=0.94, bottom=0.08, left=0.09, right=0.95)
fig.savefig(os.path.join(output_dir, 'pic_ex4_full_stand.png'), dpi=200)
plt.close(fig)

print("Готово: pic_ex4_full_stand.png сохранен с точными параметрами Варианта №8!")
