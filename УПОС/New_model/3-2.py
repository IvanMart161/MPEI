import os
import matplotlib.pyplot as plt
import numpy as np

output_dir = os.path.expanduser('~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig')
os.makedirs(output_dir, exist_ok=True)

# Частотная сетка в кГц
f_khz = np.linspace(0, 15.0, 1000)
f_hz = f_khz * 1e3

# Расчет СПШ по формуле G(f) = 4.17e-20 * f^2
G_val = np.where(f_khz <= 12.0, 4.17e-20 * (f_hz**2), 0.0)

fig, ax = plt.subplots(figsize=(7, 4.8))

# График параболической спектральной плотности
ax.plot(f_khz[f_khz <= 12.0], G_val[f_khz <= 12.0], color='#cc0000', lw=1.8)
# Вертикальный обрыв на 12 кГц
ax.plot([12.0, 12.0], [0.0, 4.17e-20 * (12e3**2)], color='#cc0000', lw=1.8)
# Нулевая линия после 12 кГц
ax.plot(f_khz[f_khz >= 12.0], G_val[f_khz >= 12.0], color='#cc0000', lw=1.8)

ax.set_xlabel(r'$f$, кГц', fontsize=10)
ax.set_ylabel(r'$G(f)$, $\mathrm{В}^2/\mathrm{Гц}$', fontsize=10)
ax.set_xlim(0, 15.0)
ax.set_ylim(0, 6.5e-12)

# Сетка и метки как на эталонном графике
ax.set_xticks(np.arange(0, 16.5, 1.5))
ax.grid(True, color='gray', alpha=0.5, ls='-')

plt.tight_layout()
fig.savefig(os.path.join(output_dir, 'pic_task_3_2.png'), dpi=200)
plt.close(fig)

print("График задачи 3.2 сохранен в exports_fig/pic_task_3_2.png")
