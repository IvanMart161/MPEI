# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.patches import Rectangle
# from scipy.signal import correlate
#
# output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
# os.makedirs(output_dir, exist_ok=True)
#
# # -------------------------------------------------------------
# # 1. Параметры 8-го варианта
# # -------------------------------------------------------------
# f0 = 10.0e6           # Центральная частота 10 МГц
# df_fss = 0.55e6       # Полоса ФСС 0.55 МГц
# G_in = 3.0e-11        # В^2/Гц (спектральная плотность дисперсии шума)
# R_val = 2.8           # кОм
#
# # Расчет шумовой полосы и СКО (действующего напряжения)
# pi_noise = (np.pi / 3.0) * df_fss
# sigma_2 = G_in * pi_noise
# sigma_u = np.sqrt(sigma_2)  # ~4.16 мВ
#
# # -------------------------------------------------------------
# # 2. Моделирование сигналов
# # -------------------------------------------------------------
# fs = 100.0e6
# N = 1_015_000
# t = np.arange(N) / fs
#
# freqs = np.fft.rfftfreq(N, 1.0 / fs)
# x = 2.0 * (freqs - f0) / df_fss
# H_fss = 1.0 / np.sqrt(1.0 + x**6)  # 3-звенный ФСС
#
# np.random.seed(42)
# white_noise = np.random.normal(0, np.sqrt(G_in * fs / 2.0), N)
# noise_fss = np.fft.irfft(np.fft.rfft(white_noise) * H_fss, n=N)
# noise_fss = noise_fss * (sigma_u / np.std(noise_fss))  # нормировка по СКО
#
# # -------------------------------------------------------------
# # 3. Подготовка данных для 6 графиков
# # -------------------------------------------------------------
# # 1) АЧХ ФСС
# f_ach = np.linspace(9.0e6, 11.5e6, 1000)
# x_ach = 2.0 * (f_ach - f0) / df_fss
# ach_curve = 1.0 / np.sqrt(1.0 + x_ach**6)
#
# # 2) Реализация во времени (0..35 мкс)
# idx_35us = np.where(t <= 35.0e-6)[0]
# t_us = t[idx_35us] * 1e6
# u_mv = noise_fss[idx_35us] * 1e3
#
# # 3) Гистограмма
# hist_vals, hist_bins = np.histogram(noise_fss * 1e3, bins=45, range=(-15, 15), density=False)
# hist_prob = hist_vals / len(noise_fss)
# bin_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])
#
# # 4) Спектральная плотность (метод Уэлча)
# seg_len = 8192
# n_segs = N // seg_len
# win = np.hanning(seg_len)
# scale_psd = 2.0 / (fs * np.sum(win**2))
# reshaped = noise_fss[:n_segs * seg_len].reshape(n_segs, seg_len)
# psd_arr = np.mean(np.abs(np.fft.rfft(reshaped * win, axis=1))**2 * scale_psd, axis=0)
# f_psd = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e6
# psd_scaled = psd_arr * 1e6  # В^2/МГц
# psd_scaled = psd_scaled * ((G_in * 1e6) / np.max(psd_scaled)) * 1.05
#
# # 5) АКФ шума (tau от -4 до 4 мкс)
# M_corr = 8000
# sig_cut = noise_fss[:M_corr]
# sig_cut = sig_cut - np.mean(sig_cut)
# akf_raw = correlate(sig_cut, sig_cut, mode='same')
# lags = (np.arange(len(akf_raw)) - len(akf_raw)//2) / fs
# akf_curve = akf_raw * (sigma_2 / np.max(akf_raw))
#
# mask_4us = (lags >= -4.0e-6) & (lags <= 4.0e-6)
# lags_4us = lags[mask_4us] * 1e6
# akf_4us = akf_curve[mask_4us]
#
# mask_1us = (lags >= -1.0e-6) & (lags <= 1.0e-6)
# lags_1us = lags[mask_1us] * 1e6
# akf_1us = akf_curve[mask_1us]
#
# # -------------------------------------------------------------
# # 4. Отрисовка стенда
# # -------------------------------------------------------------
# bg_color = '#b0c8c8'
# fig = plt.figure(figsize=(16, 8.5), facecolor=bg_color)
#
# # Левая панель с параметрами
# ax_left = fig.add_axes([0.0, 0.0, 0.17, 1.0], facecolor=bg_color)
# ax_left.axis('off')
#
# def draw_box(ax, x, y, w, h, text, title=None, unit=None):
#     rect = Rectangle((x, y), w, h, facecolor='white', edgecolor='gray', lw=1)
#     ax.add_patch(rect)
#     ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=9, fontweight='bold')
#     if title:
#         ax.text(x + w/2, y + h + 0.012, title, ha='center', va='bottom', fontsize=8)
#     if unit:
#         ax.text(x + w + 0.015, y + h/2, unit, ha='left', va='center', fontsize=8)
#
# draw_box(ax_left, 0.28, 0.90, 0.44, 0.04, f"{R_val}", "Rг=Rн", "кОм")
# draw_box(ax_left, 0.08, 0.79, 0.38, 0.04, "2.450", "L1", "мкГн")
# draw_box(ax_left, 0.54, 0.79, 0.38, 0.04, "95.833", "C1", "пФ")
# draw_box(ax_left, 0.08, 0.70, 0.38, 0.04, "1.225", "L2", "мкГн")
# draw_box(ax_left, 0.54, 0.70, 0.38, 0.04, "191.667", "C2", "пФ")
# draw_box(ax_left, 0.28, 0.60, 0.44, 0.04, "5.684", "Cсв", "пФ")
#
# ax_left.plot([0, 1], [0.56, 0.56], color='#7a9292', lw=1.5)
# ax_left.text(0.5, 0.52, "Спектральная плотность\nдисперсии шума", ha='center', va='center', fontsize=8)
# draw_box(ax_left, 0.22, 0.44, 0.56, 0.04, "3E-011", unit="В^2/Гц")
#
# ax_left.plot([0, 1], [0.41, 0.41], color='#7a9292', lw=1.5)
# ax_left.text(0.5, 0.37, "Количество отсчётов", ha='center', va='center', fontsize=8)
# draw_box(ax_left, 0.15, 0.32, 0.70, 0.035, "1 015 000")
#
# ax_left.text(0.5, 0.27, "Эффективное напряжение шума\nна выходе ФСС", ha='center', va='center', fontsize=8)
# draw_box(ax_left, 0.25, 0.20, 0.50, 0.04, f"{sigma_u*1e3:.2f}", unit="мВ")
#
# draw_box(ax_left, 0.15, 0.11, 0.70, 0.04, "АЧХ")
# draw_box(ax_left, 0.15, 0.06, 0.70, 0.04, "Пуск")
# draw_box(ax_left, 0.15, 0.01, 0.70, 0.04, "Закончить")
# ax_left.plot([1.0, 1.0], [0, 1], color='#607878', lw=1.5)
#
# # -------------------------------------------------------------
# # Схема фильтра (верхняя центральная область)
# # -------------------------------------------------------------
# ax_sch = fig.add_axes([0.18, 0.68, 0.58, 0.29], facecolor='white')
# ax_sch.set_xlim(-0.5, 11.5)
# ax_sch.set_ylim(-1.5, 3.5)
# ax_sch.axis('off')
#
# # Шины
# ax_sch.plot([0, 11], [0, 0], color='green', lw=2)
# ax_sch.plot([0, 11], [2, 2], color='green', lw=2)
# for g_x in [0, 4, 8]:
#     ax_sch.plot([g_x, g_x], [0, -0.4], color='green', lw=1.5)
#     ax_sch.plot([g_x-0.3, g_x+0.3], [-0.4, -0.4], color='green', lw=1.5)
#
# # Источник
# ax_sch.plot([0, 0], [0, 2], color='green', lw=2)
# circle = plt.Circle((0, 1), 0.45, fill=False, edgecolor='blue', lw=2)
# ax_sch.add_patch(circle)
# ax_sch.text(0, 1, "~", ha='center', va='center', fontsize=18, color='blue')
#
# # Резисторы
# ax_sch.add_patch(Rectangle((0.8, 1.8), 0.8, 0.4, facecolor='white', edgecolor='blue', lw=1.5))
# ax_sch.text(1.2, 2.4, r"$R_\Gamma$", color='green', fontsize=11, fontweight='bold', ha='center')
# ax_sch.add_patch(Rectangle((10.8, 0.6), 0.4, 0.8, facecolor='white', edgecolor='blue', lw=1.5))
# ax_sch.text(10.3, 1.0, r"$R_н$", color='green', fontsize=11, fontweight='bold', ha='right')
#
# # 3 параллельных контура
# x_stages = [3.2, 6.4, 9.6]
# labels_c = [r"$C_1$", r"$C_2$", r"$C_1$"]
# labels_l = [r"$L_1$", r"$L_2$", r"$L_1$"]
#
# for i, xs in enumerate(x_stages):
#     ax_sch.plot([xs, xs], [0, 2], color='green', lw=1.5)
#     ax_sch.plot([xs-0.6, xs+0.6], [1.7, 1.7], color='green', lw=1.5)
#     ax_sch.plot([xs-0.6, xs+0.6], [0.3, 0.3], color='green', lw=1.5)
#
#     # Конденсатор
#     ax_sch.plot([xs-0.6, xs-0.6], [1.7, 1.1], color='green', lw=1.5)
#     ax_sch.plot([xs-0.6, xs-0.6], [0.3, 0.9], color='green', lw=1.5)
#     ax_sch.plot([xs-0.85, xs-0.35], [1.1, 1.1], color='blue', lw=2)
#     ax_sch.plot([xs-0.85, xs-0.35], [0.9, 0.9], color='blue', lw=2)
#     ax_sch.text(xs-0.95, 1.0, labels_c[i], color='green', fontsize=11, fontweight='bold', ha='right')
#
#     # Катушка индуктивности
#     ax_sch.plot([xs+0.6, xs+0.6], [1.7, 0.3], color='blue', lw=2.5, ls='--')
#     ax_sch.text(xs+0.95, 1.0, labels_l[i], color='green', fontsize=11, fontweight='bold', ha='left')
#
#     # Красные узлы
#     ax_sch.plot([xs], [2], 'ro', markersize=5)
#     ax_sch.plot([xs], [0], 'ro', markersize=5)
#
# # Емкости связи C_св
# for c_x in [4.8, 8.0]:
#     ax_sch.plot([c_x-0.3, c_x-0.3], [1.8, 2.2], color='blue', lw=2)
#     ax_sch.plot([c_x+0.3, c_x+0.3], [1.8, 2.2], color='blue', lw=2)
#     ax_sch.text(c_x, 2.4, r"$C_{\mathrm{св}}$", color='green', fontsize=11, fontweight='bold', ha='center')
#
# # -------------------------------------------------------------
# # 6 измерительных графиков
# # -------------------------------------------------------------
# # 1. АЧХ ФСС (справа вверху)
# ax1 = fig.add_axes([0.79, 0.68, 0.18, 0.27])
# ax1.plot(f_ach / 1e6, ach_curve, color='#0011dd', lw=1.5)
# ax1.set_title("АЧХ ФСС", fontsize=9)
# ax1.set_xlabel("f, МГц", fontsize=8)
# ax1.set_xlim(9.0, 11.5)
# ax1.set_ylim(0, 1.05)
# ax1.grid(True, ls='--', color='gray')
#
# # 2. Шум на выходе ФСС (широкий по центру)
# ax2 = fig.add_axes([0.22, 0.37, 0.52, 0.22])
# ax2.plot(t_us, u_mv, color='#0011dd', lw=0.4)
# ax2.set_title("Шум на выходе ФСС", fontsize=9)
# ax2.set_xlabel("t, мкс", fontsize=8)
# ax2.set_ylabel("u, мВ", fontsize=8)
# ax2.set_xlim(0, 35)
# ax2.set_ylim(-18, 18)
# ax2.grid(True, ls='--', color='gray')
#
# # 3. Гистограмма (справа в центре)
# ax3 = fig.add_axes([0.79, 0.37, 0.18, 0.22])
# ax3.bar(bin_centers, hist_prob, width=0.45, color='black', edgecolor='black')
# ax3.set_title("Гистограмма", fontsize=9)
# ax3.set_xlabel("u, мВ", fontsize=8)
# ax3.set_xlim(-15, 15)
# ax3.set_ylim(0, 0.18)
# ax3.grid(True, ls='--', color='gray')
#
# # 4. СП шума на выходе ФСС (слева снизу)
# ax4 = fig.add_axes([0.22, 0.07, 0.23, 0.22])
# ax4.plot(f_psd, psd_scaled, color='#cc0000', lw=1.2)
# ax4.set_title("Спектральная плотность\nшума на выходе ФСС", fontsize=9)
# ax4.set_xlabel("f, МГц", fontsize=8)
# ax4.set_ylabel(r"Спектральная плотность ($\mathrm{В}^2/\mathrm{МГц}$)", fontsize=8)
# ax4.set_xlim(8.0, 12.0)
# ax4.set_ylim(0, 3.8e-5)
# ax4.text(-0.04, 1.03, r"$\times 10^{-5}$", transform=ax4.transAxes, fontsize=8)
# ax4.grid(True, ls='--', color='gray')
#
# # 5. АКФ шума (по центру снизу)
# ax5 = fig.add_axes([0.50, 0.07, 0.23, 0.22])
# ax5.plot(lags_4us, akf_4us, color='#0011dd', lw=0.6)
# ax5.set_title("АКФ шума на выходе ФСС", fontsize=9)
# ax5.set_xlabel(r"$\tau$, мкс", fontsize=8)
# ax5.set_xlim(-4, 4)
# ax5.set_ylim(-1.8e-5, 2.0e-5)
# ax5.text(-0.04, 1.03, r"$\times 10^{-5}$", transform=ax5.transAxes, fontsize=8)
# ax5.grid(True, ls='--', color='gray')
#
# # 6. Центральная часть АКФ (справа снизу)
# ax6 = fig.add_axes([0.79, 0.07, 0.18, 0.22])
# ax6.plot(lags_1us, akf_1us, color='#0011dd', lw=0.9)
# ax6.set_title("Центральная часть АКФ шума\nна выходе ФСС", fontsize=9)
# ax6.set_xlabel(r"$\tau$, мкс", fontsize=8)
# ax6.set_xlim(-1, 1)
# ax6.set_ylim(-1.8e-5, 2.0e-5)
# ax6.text(-0.04, 1.03, r"$\times 10^{-5}$", transform=ax6.transAxes, fontsize=8)
# ax6.grid(True, ls='--', color='gray')
#
# save_path = os.path.join(output_dir, "pic_fss_stand.png")
# plt.savefig(save_path, dpi=200)
# plt.close()
# print(f"Стенд ФСС успешно сгенерирован: {save_path}")
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------------------
# 1. Параметры 8-го варианта ФСС
# -------------------------------------------------------------
f0 = 10.0e6           # Центральная частота 10 МГц
df_fss = 0.55e6       # Полоса ФСС 0.55 МГц
G_in = 3.0e-11        # В^2/Гц

# СКО шума на выходе ФСС
pi_noise = (np.pi / 3.0) * df_fss
sigma_2 = G_in * pi_noise
sigma_u = np.sqrt(sigma_2)  # ~4.16 мВ

# -------------------------------------------------------------
# 2. Моделирование процесса
# -------------------------------------------------------------
fs = 100.0e6
N = 1_015_000
t = np.arange(N) / fs

freqs = np.fft.rfftfreq(N, 1.0 / fs)
x = 2.0 * (freqs - f0) / df_fss
H_fss = 1.0 / np.sqrt(1.0 + x**6)

np.random.seed(42)
white_noise = np.random.normal(0, np.sqrt(G_in * fs / 2.0), N)
noise_fss = np.fft.irfft(np.fft.rfft(white_noise) * H_fss, n=N)
noise_fss = noise_fss * (sigma_u / np.std(noise_fss))

# -------------------------------------------------------------
# 3. Подготовка данных для графиков
# -------------------------------------------------------------
# 1) АЧХ ФСС
f_ach = np.linspace(8.5e6, 11.5e6, 1000)
x_ach = 2.0 * (f_ach - f0) / df_fss
ach_curve = 1.0 / np.sqrt(1.0 + x_ach**6)

# 2) Реализация во времени (0..35 мкс)[cite: 17]
idx_35us = np.where(t <= 35.0e-6)[0]
t_us = t[idx_35us] * 1e6
u_mv = noise_fss[idx_35us] * 1e3

# 3) Гистограмма распределения[cite: 17]
hist_vals, hist_bins = np.histogram(noise_fss * 1e3, bins=45, range=(-15, 15), density=False)
hist_prob = hist_vals / len(noise_fss)
bin_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])

# 4) СПМ шума на выходе ФСС[cite: 17]
seg_len = 8192
n_segs = N // seg_len
win = np.hanning(seg_len)
scale_psd = 2.0 / (fs * np.sum(win**2))
reshaped = noise_fss[:n_segs * seg_len].reshape(n_segs, seg_len)
psd_arr = np.mean(np.abs(np.fft.rfft(reshaped * win, axis=1))**2 * scale_psd, axis=0)
f_psd = np.fft.rfftfreq(seg_len, 1.0 / fs) / 1e6
psd_scaled = psd_arr * 1e6
psd_scaled = psd_scaled * ((G_in * 1e6) / np.max(psd_scaled)) * 1.05

# 5) АКФ шума (полная и центральная)[cite: 17]
M_corr = 8000
sig_cut = noise_fss[:M_corr] - np.mean(noise_fss[:M_corr])
akf_raw = correlate(sig_cut, sig_cut, mode='same')
lags = (np.arange(len(akf_raw)) - len(akf_raw)//2) / fs
akf_curve = akf_raw * (sigma_2 / np.max(akf_raw))

mask_4us = (lags >= -4.0e-6) & (lags <= 4.0e-6)
lags_4us = lags[mask_4us] * 1e6
akf_4us = akf_curve[mask_4us]

mask_1us = (lags >= -1.0e-6) & (lags <= 1.0e-6)
lags_1us = lags[mask_1us] * 1e6
akf_1us = akf_curve[mask_1us]

# -------------------------------------------------------------
# 4. Чистая сетка 2x3 без серого интерфейса
# -------------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))

# [0, 0] АЧХ ФСС
axs[0, 0].plot(f_ach / 1e6, ach_curve, color='#0011dd', lw=1.5)
axs[0, 0].set_title("АЧХ трехзвенного ФСС", fontsize=9)
axs[0, 0].set_xlabel("f, МГц", fontsize=8)
axs[0, 0].set_xlim(8.5, 11.5)
axs[0, 0].set_ylim(0, 1.05)
axs[0, 0].grid(True, ls='--', color='lightgray')

# [0, 1] Шум во времени
axs[0, 1].plot(t_us, u_mv, color='#0011dd', lw=0.4)
axs[0, 1].set_title("Шум на выходе ФСС", fontsize=9)
axs[0, 1].set_xlabel("t, мкс", fontsize=8)
axs[0, 1].set_ylabel("u, мВ", fontsize=8)
axs[0, 1].set_xlim(0, 35)
axs[0, 1].set_ylim(-16, 16)
axs[0, 1].grid(True, ls='--', color='lightgray')

# [0, 2] Гистограмма
axs[0, 2].bar(bin_centers, hist_prob, width=0.45, color='black', edgecolor='black')
axs[0, 2].set_title("Гистограмма распределения шума", fontsize=9)
axs[0, 2].set_xlabel("u, мВ", fontsize=8)
axs[0, 2].set_xlim(-15, 15)
axs[0, 2].set_ylim(0, 0.18)
axs[0, 2].grid(True, ls='--', color='lightgray')

# [1, 0] СП шума
axs[1, 0].plot(f_psd, psd_scaled, color='#cc0000', lw=1.2)
axs[1, 0].set_title("Спектральная плотность шума на выходе ФСС", fontsize=9)
axs[1, 0].set_xlabel("f, МГц", fontsize=8)
axs[1, 0].set_ylabel(r"СП ($\mathrm{В}^2/\mathrm{МГц}$)", fontsize=8)
axs[1, 0].set_xlim(8.0, 12.0)
axs[1, 0].set_ylim(0, 3.8e-5)
axs[1, 0].text(-0.02, 1.03, r"$\times 10^{-5}$", transform=axs[1, 0].transAxes, fontsize=8)
axs[1, 0].grid(True, ls='--', color='lightgray')

# [1, 1] АКФ шума
axs[1, 1].plot(lags_4us, akf_4us, color='#0011dd', lw=0.6)
axs[1, 1].set_title("АКФ шума на выходе ФСС", fontsize=9)
axs[1, 1].set_xlabel(r"$\tau$, мкс", fontsize=8)
axs[1, 1].set_xlim(-4, 4)
axs[1, 1].set_ylim(-1.8e-5, 2.0e-5)
axs[1, 1].text(-0.02, 1.03, r"$\times 10^{-5}$", transform=axs[1, 1].transAxes, fontsize=8)
axs[1, 1].grid(True, ls='--', color='lightgray')

# [1, 2] Центральная часть АКФ
axs[1, 2].plot(lags_1us, akf_1us, color='#0011dd', lw=0.9)
axs[1, 2].set_title("Центральная часть АКФ шума", fontsize=9)
axs[1, 2].set_xlabel(r"$\tau$, мкс", fontsize=8)
axs[1, 2].set_xlim(-1, 1)
axs[1, 2].set_ylim(-1.8e-5, 2.0e-5)
axs[1, 2].text(-0.02, 1.03, r"$\times 10^{-5}$", transform=axs[1, 2].transAxes, fontsize=8)
axs[1, 2].grid(True, ls='--', color='lightgray')

plt.tight_layout()
save_path = os.path.join(output_dir, "pic_fss_results.png")
plt.savefig(save_path, dpi=200)
plt.close()
print(f"График сохранен: {save_path}")
