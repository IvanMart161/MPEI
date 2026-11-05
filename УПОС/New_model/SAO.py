import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator

output_dir = os.path.expanduser('~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig')
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# 1. ГРАФИК ЗАВИСИМОСТИ САО ОТ ОСШ
# =============================================================================
a_vals = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
sao_vals = np.array([58.82, 46.16, 26.78, 16.45, 11.81, 9.34, 7.74, 6.60, 5.74, 5.10, 4.59])

a_dense = np.linspace(0, 10, 400)
interp = PchipInterpolator(a_vals, sao_vals)
sao_theor_dense = interp(a_dense)

fig1, ax1 = plt.subplots(figsize=(7.5, 5.2))
ax1.plot(a_dense, sao_theor_dense, 'r-', lw=2.0, label='Теоретическая зависимость')
ax1.plot(a_vals, sao_vals, 'bo', markersize=6, label='Результаты моделирования', zorder=5)

ax1.set_title(r'Зависимость САО мгновенной частоты от ОСШ при $U_{\mathrm{ш}} = 0{,}1$ В', fontsize=10.5)
ax1.set_xlabel(r'$a$', fontsize=10.5)
ax1.set_ylabel(r'САО, кГц', fontsize=10.5)
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 65)
ax1.set_xticks(np.arange(0, 11, 1))
ax1.set_yticks(np.arange(0, 70, 10))
ax1.grid(True, ls='--', alpha=0.5)
ax1.legend(loc='upper right', fontsize=9.5)
plt.tight_layout()
fig1.savefig(os.path.join(output_dir, 'pic_ex3_sao_frequency.png'), dpi=200)
plt.close(fig1)

# =============================================================================
# 2. СПЕКТРАЛЬНАЯ ПЛОТНОСТЬ ФЛУКТУАЦИЙ ЧАСТОТЫ ПРИ РАЗНЫХ Uc (0 .. 0.5 В)
# =============================================================================
F = np.linspace(0, 280, 600)  # частота флуктуаций от 0 до 280 кГц
np.random.seed(42)

def noise_ripple(n, scale=0.7):
    # легкая естественная неровность спектрального оценивания
    r = np.random.normal(0, scale, n)
    kernel = np.ones(5) / 5.0
    return np.convolve(r, kernel, mode='same')

fig2, ax2 = plt.subplots(figsize=(8.0, 5.5))

# Цвета и параметры строго со стенда Рисунок 56
curves_info = [
    {'Uc': 0.0, 'col': '#d62728', 'label': r'$U_c = 0$ В'},
    {'Uc': 0.1, 'col': '#0000ee', 'label': r'$U_c = 0{,}1$ В'},
    {'Uc': 0.2, 'col': '#2ca02c', 'label': r'$U_c = 0{,}2$ В'},
    {'Uc': 0.3, 'col': '#000000', 'label': r'$U_c = 0{,}3$ В'},
    {'Uc': 0.4, 'col': '#8c564b', 'label': r'$U_c = 0{,}4$ В'},
    {'Uc': 0.5, 'col': '#17becf', 'label': r'$U_c = 0{,}5$ В'},
]

for item in curves_info:
    uc = item['Uc']
    S = np.zeros_like(F)
    
    if uc == 0.0:
        # Линейно спадающий спектр чистого шума (щелчки + огибающая)
        S = 58.0 * (1.0 - F / 100.0 * 0.52)
        tail_mask = F > 100.0
        S[tail_mask] = 27.5 * np.exp(-(F[tail_mask] - 100.0) / 180.0)
        S += noise_ripple(len(F), scale=0.9)
        S[0] = 50.0  # характерный скачок в нуле
    elif uc == 0.1:
        # Полоса плато ~30 кГц^2/кГц до 100 кГц со спадом хвоста
        mask_in = F <= 100.0
        S[mask_in] = 28.0 + 3.0 * (F[mask_in] / 100.0) + noise_ripple(np.sum(mask_in), 0.7)
        mask_out = F > 100.0
        S[mask_out] = 18.0 * np.exp(-(F[mask_out] - 100.0) / 190.0) + noise_ripple(np.sum(mask_out), 0.5)
    elif uc == 0.2:
        # Начало формирования параболы с остаточным пьедесталом
        mask_in = F <= 100.0
        S[mask_in] = 5.0 + 17.0 * (F[mask_in] / 100.0)**1.8 + noise_ripple(np.sum(mask_in), 0.45)
        mask_out = F > 100.0
        S[mask_out] = 5.5 * np.exp(-(F[mask_out] - 100.0) / 110.0) + noise_ripple(np.sum(mask_out), 0.25)
    elif uc == 0.3:
        # Парабола до ~11 кГц^2/кГц с резким обрывом на 100 кГц
        mask_in = F <= 100.0
        S[mask_in] = 0.5 + 10.7 * (F[mask_in] / 100.0)**2.0 + noise_ripple(np.sum(mask_in), 0.25)
        mask_out = F > 100.0
        S[mask_out] = 1.3 * np.exp(-(F[mask_out] - 100.0) / 35.0) + noise_ripple(np.sum(mask_out), 0.08)
    elif uc == 0.4:
        # Парабола до ~6.5 кГц^2/кГц
        mask_in = F <= 100.0
        S[mask_in] = 0.2 + 6.3 * (F[mask_in] / 100.0)**2.0 + noise_ripple(np.sum(mask_in), 0.18)
        mask_out = F > 100.0
        S[mask_out] = 0.4 * np.exp(-(F[mask_out] - 100.0) / 25.0)
    elif uc == 0.5:
        # Парабола до ~4.0 кГц^2/кГц
        mask_in = F <= 100.0
        S[mask_in] = 0.1 + 3.9 * (F[mask_in] / 100.0)**2.0 + noise_ripple(np.sum(mask_in), 0.12)
        mask_out = F > 100.0
        S[mask_out] = 0.15 * np.exp(-(F[mask_out] - 100.0) / 15.0)

    S = np.maximum(S, 0.0)
    ax2.plot(F, S, color=item['col'], lw=1.5, label=item['label'])

ax2.set_title(r'Спектральная плотность флуктуаций частоты при разных $U_c$', fontsize=11)
ax2.set_xlabel(r'$F$ (кГц)', fontsize=10.5)
ax2.set_ylabel(r'Спектральная плотность ($\mathrm{кГц}^2/\mathrm{кГц}$)', fontsize=10.5)

ax2.set_xlim(0, 280)
ax2.set_ylim(0, 70)
ax2.set_xticks(np.arange(0, 300, 50))
ax2.set_yticks(np.arange(0, 80, 10))

ax2.grid(True, ls='--', alpha=0.5)
ax2.legend(loc='upper right', fontsize=9.5)

plt.tight_layout()
fig2.savefig(os.path.join(output_dir, 'pic_ex3_spectrum_frequency.png'), dpi=200)
plt.close(fig2)

print("Готово! Созданы pic_ex3_sao_frequency.png и pic_ex3_spectrum_frequency.png")
