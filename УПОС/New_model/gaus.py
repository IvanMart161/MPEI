import os
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.path.expanduser("~/Documents/MPEI/УПОС/РГР1_Иван/exports_fig")
os.makedirs(output_dir, exist_ok=True)

# Исходные данные задачи 2-1
G0 = 5e-19           # Вт/Гц
Kp_lin = 10**(6/10)  # 6 дБ в разах (~3.981)
Pi_sh = 1e6          # 1 МГц
R_n = 50.0           # Ом

# Расчет параметров шума
G_out = G0 * Kp_lin
P_sh = G_out * Pi_sh
U_sh = np.sqrt(P_sh * R_n)  # ~9.98 мкВ
sigma2 = U_sh**2

# Массивы напряжений для графиков
u = np.linspace(-4 * U_sh, 4 * U_sh, 1000)
U_env = np.linspace(0, 4 * U_sh, 1000)

# Плотности вероятностей
w_u = (1.0 / (U_sh * np.sqrt(2 * np.pi))) * np.exp(-u**2 / (2 * sigma2))
w_U = (U_env / sigma2) * np.exp(-U_env**2 / (2 * sigma2))

# Отрисовка графиков
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# 1. Нормальное распределение мгновенных значений
ax1.plot(u * 1e5, w_u * 1e-4, color='#0022aa', lw=2)
ax1.set_title("Нормальное распределение (Гаусс)", fontsize=11)
ax1.set_xlabel(r"Мгновенное напряжение $u$, $\times 10^{-5}$ В", fontsize=10)
ax1.set_ylabel(r"Плотность вероятности $w(u)$, $\times 10^{4}$ 1/В", fontsize=10)
ax1.set_xlim(-4, 4)
ax1.set_ylim(0, 4.2)
ax1.grid(True, ls='--', color='lightgray')

# 2. Распределение Рэлея огибающей
ax2.plot(U_env * 1e5, w_U * 1e-4, color='#cc1111', lw=2)
ax2.set_title("Распределение Рэлея", fontsize=11)
ax2.set_xlabel(r"Огибающая $U$, $\times 10^{-5}$ В", fontsize=10)
ax2.set_ylabel(r"Плотность вероятности $w(U)$, $\times 10^{4}$ 1/В", fontsize=10)
ax2.set_xlim(0, 4)
ax2.set_ylim(0, 6.5)
ax2.grid(True, ls='--', color='lightgray')

plt.tight_layout()
save_path = os.path.join(output_dir, "pic_task_2_1.png")
plt.savefig(save_path, dpi=200)
plt.close()
print(f"График сохранен: {save_path}")
