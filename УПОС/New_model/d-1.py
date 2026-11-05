import os
import matplotlib.pyplot as plt
import numpy as np

output_dir = os.path.expanduser('~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig')
os.makedirs(output_dir, exist_ok=True)

Pi_sh = 8.0  # МГц
f = np.linspace(0, 8.0, 1000)
line_color = '#0033bb'  # Насыщенный синий цвет

# -----------------------------------------------------------------------------
# 1. График СПШ на выходе линейного АД (сильный сигнал a >= 3)
# -----------------------------------------------------------------------------
# До 4 МГц прямоугольная составляющая с легким наклоном, после 4 МГц - малый шумовой остаток
G_before = 1.34e-10 - 0.02e-10 * (f[f <= 4.0] / 4.0)
G_after = 0.015e-10 * (1.0 - (f[f > 4.0] - 4.0) / 4.0)  # отчетливый след над нулем

fig1, ax1 = plt.subplots(figsize=(7, 4.8))

# Рисуем до скачка
ax1.plot(f[f <= 4.0], G_before, color=line_color, lw=1.8)
# Вертикальный скачок
ax1.plot([4.0, 4.0], [G_after[0], G_before[-1]], color=line_color, lw=1.8)
# Хвост после скачка
ax1.plot(f[f > 4.0], G_after, color=line_color, lw=1.8)

ax1.set_xlabel('f, МГц', fontsize=10)
ax1.set_ylabel(r'$G(f)$, Вт/Гц', fontsize=10)
ax1.set_xlim(0, 8.0)
# Делаем нижнюю границу чуть ниже 0, чтобы нулевой уровень не слипался с рамкой
ax1.set_ylim(-0.03e-10, 1.45e-10)
ax1.set_xticks(np.arange(0, 8.5, 0.5))
ax1.grid(True, color='gray', alpha=0.45, ls='-')

plt.tight_layout()
fig1.savefig(os.path.join(output_dir, 'pic_task_d1_a3.png'), dpi=200)
plt.close(fig1)

# -----------------------------------------------------------------------------
# 2. График СПШ на выходе линейного АД при a = 1
# -----------------------------------------------------------------------------
G1 = np.where(f <= 4.0, 0.53e-10, 0.0)
G2 = 0.36e-10 * (1.0 - f / 8.0)
G_sum = G1 + G2

fig2, ax2 = plt.subplots(figsize=(7, 4.8))

idx_4 = np.searchsorted(f, 4.0)
# До скачка
ax2.plot(f[:idx_4], G_sum[:idx_4], color=line_color, lw=1.8)
# Вертикальный скачок
ax2.plot([4.0, 4.0], [G2[idx_4], G_sum[idx_4 - 1]], color=line_color, lw=1.8)
# После скачка
ax2.plot(f[idx_4:], G2[idx_4:], color=line_color, lw=1.8)

ax2.set_xlabel('f, МГц', fontsize=10)
ax2.set_ylabel(r'$G(f)$, Вт/Гц', fontsize=10)
ax2.set_xlim(0, 8.0)
ax2.set_ylim(-0.02e-10, 1.05e-10)
ax2.set_xticks(np.arange(0, 8.5, 0.5))
ax2.grid(True, color='gray', alpha=0.45, ls='-')

plt.tight_layout()
fig2.savefig(os.path.join(output_dir, 'pic_task_d1_a1.png'), dpi=200)
plt.close(fig2)

print('Графики Д-1 в синем цвете с четким нулем успешно созданы!')
