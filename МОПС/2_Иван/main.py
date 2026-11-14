import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# Исходные константы (Лекция 9, слайд 10)
# ==========================================
alpha = 1.0                     # с^-1
T = 0.01                        # 10 мс
f0 = 1602e6                     # Гц
omega0 = 2.0 * np.pi * f0       # рад/с
c = 3e8                         # м/с
k_dop = omega0 / c              # коэффициент пересчета (рад/с) / (м/с)


def calc_N0_tilde(q_lin):
    return (2.0 / (q_lin * T**2)) * (1.0 + 1.0 / (2.0 * q_lin * T))


def calc_S_xi(sigma_a):
    return 2.0 * (sigma_a**2) * alpha * (k_dop**2)


def steady_state_params(q_dBHz, sigma_a):
    q_lin = 10.0**(0.1 * q_dBHz)
    N0_t = calc_N0_tilde(q_lin)
    S_xi = calc_S_xi(sigma_a)

    # Дисперсия ошибки фильтрации частоты D11 и коэффициенты фильтра
    D11 = (alpha * N0_t / 2.0) * (
        np.sqrt(1.0 + 2.0 * np.sqrt(S_xi / (alpha**2 * N0_t))) - 1.0
    )
    sigma_omega = np.sqrt(D11)

    K1 = 2.0 * D11 / N0_t
    K2 = (K1**2) / 2.0

    # Эквивалентная шумовая полоса ЧАП (аналитический интеграл)
    num = (K1**2) * (alpha + K1) + (alpha * K1 + K2)**2
    den = 4.0 * (alpha + K1) * (alpha * K1 + K2)
    dF_fll = num / den

    return sigma_omega, dF_fll, D11


# ==========================================
# Пункты 1 и 2: Аналитические зависимости
# ==========================================
q_grid_dBHz = np.linspace(14.0, 50.0, 200)
sigma_omega_q, dF_q, _ = steady_state_params(q_grid_dBHz, sigma_a=10.0)

sigma_a_grid = np.linspace(1.0, 30.0, 200)
sigma_omega_a, dF_a, _ = steady_state_params(q_dBHz=30.0, sigma_a=sigma_a_grid)

# ==========================================
# Пункты 3–8: Дискретное моделирование (Калман)
# ==========================================
q_sim_dBHz = 30.0
sigma_a_sim = 10.0
q_sim_lin = 10.0**(0.1 * q_sim_dBHz)

N0_t_sim = calc_N0_tilde(q_sim_lin)
S_xi_sim = calc_S_xi(sigma_a_sim)

sigma_xi_d = np.sqrt(S_xi_sim / (2.0 * T))     # СКО дискретного формирующего шума
sigma_n_d = np.sqrt(N0_t_sim / (2.0 * T))      # СКО дискретного шума наблюдений

T_sim = 1000.0                                 # 1000 с
N_steps = int(T_sim / T)
t_arr = np.arange(N_steps) * T

# Матрицы дискретной системы
F = np.array([
    [1.0, T],
    [0.0, 1.0 - alpha * T],
])
G = np.array([
    [0.0],
    [alpha * T],
])
H = np.array([[1.0, 0.0]])
Q_d = (G @ G.T) * (sigma_xi_d**2)
R_d = sigma_n_d**2

# Начальные условия
rng = np.random.default_rng(42)
x_true = np.array([100.0, 100.0])                # Omega_0 = 100 рад/с, v_0 = 0
x_est = np.array([0.0, 0.0])                   # Оценка в начале = 0
D = np.array([
    [34.0**2, 0.0],
    [0.0, 340.0**2],
])

# Массивы для сохранения траекторий
omega_true_hist = np.zeros(N_steps)
omega_est_hist = np.zeros(N_steps)
D11_hist = np.zeros(N_steps)

xi_noise = rng.normal(0.0, sigma_xi_d, size=N_steps)
n_noise = rng.normal(0.0, sigma_n_d, size=N_steps)

for k in range(N_steps):
    omega_true_hist[k] = x_true[0]
    omega_est_hist[k] = x_est[0]
    D11_hist[k] = D[0, 0]

    x_true = F @ x_true + np.array([0.0, alpha * T * xi_noise[k]])
    y_meas = x_true[0] + n_noise[k]

    x_pred = F @ x_est
    D_pred = F @ D @ F.T + Q_d

    S_inn = D_pred[0, 0] + R_d
    K = D_pred[:, 0] / S_inn

    innov = y_meas - x_pred[0]
    x_est = x_pred + K * innov
    D = (np.eye(2) - np.outer(K, H)) @ D_pred

err_omega = omega_est_hist - omega_true_hist

# ==========================================
# Сохранение отдельных фигур: 1.png .. 5.png
# ==========================================

# 1.png: Истинная доплеровская частота (0..1000 с)
fig1, ax1 = plt.subplots(figsize=(9, 4.5))
ax1.plot(t_arr, omega_true_hist, linewidth=0.8)
ax1.set_title(r'Истинная доплеровская частота $\Omega_k(t_k)$, $t \in [0, 1000]\text{ с}$', fontsize=11)
ax1.set_xlabel(r'$t$, с')
ax1.set_ylabel(r'$\Omega$, рад/с')
ax1.grid(True)
fig1.tight_layout()
fig1.savefig("1.png", dpi=300)

# 2.png: СКО фильтрации (0..1 с)
idx_1s = int(1.0 / T)
fig2, ax2 = plt.subplots(figsize=(9, 4.5))
ax2.plot(t_arr[:idx_1s], np.sqrt(D11_hist[:idx_1s]), linewidth=2, color='tab:orange')
ax2.set_title(r'СКО фильтрации частоты $\sqrt{D_{11}}(t_k)$ на переходном участке, $t \in [0, 1]\text{ с}$', fontsize=11)
ax2.set_xlabel(r'$t$, с')
ax2.set_ylabel(r'$\sqrt{D_{11}}$, рад/с')
ax2.grid(True)
fig2.tight_layout()
fig2.savefig("2.png", dpi=300)

# 3.png: Мгновенная ошибка фильтрации (0..10 с)
idx_10s = int(10.0 / T)
fig3, ax3 = plt.subplots(figsize=(9, 4.5))
ax3.plot(t_arr[:idx_10s], err_omega[:idx_10s], linewidth=1.0, color='tab:green')
ax3.set_title(r'Мгновенная ошибка фильтрации $\epsilon_\Omega(t_k) = \hat{\Omega}_k - \Omega_k$, $t \in [0, 10]\text{ с}$', fontsize=11)
ax3.set_xlabel(r'$t$, с')
ax3.set_ylabel(r'$\epsilon_\Omega$, рад/с')
ax3.grid(True)
fig3.tight_layout()
fig3.savefig("3.png", dpi=300)

# 4.png: Зависимости от q_{c/n0} (2 графика)
fig4, axes4 = plt.subplots(1, 2, figsize=(12, 4.8))
axes4[0].plot(q_grid_dBHz, sigma_omega_q, linewidth=2)
axes4[0].set_title(r'СКО фильтрации частоты $\sigma_\Omega(q_{c/n_0})$ при $\sigma_a = 10\text{ м/с}^2$', fontsize=10)
axes4[0].set_xlabel(r'$q_{c/n_0}$, дБГц')
axes4[0].set_ylabel(r'$\sigma_\Omega$, рад/с')
axes4[0].grid(True)

axes4[1].plot(q_grid_dBHz, dF_q, linewidth=2, color='tab:orange')
axes4[1].set_title(r'Оптимальная полоса ЧАП $\Delta F_{\text{ЧАП}}(q_{c/n_0})$ при $\sigma_a = 10\text{ м/с}^2$', fontsize=10)
axes4[1].set_xlabel(r'$q_{c/n_0}$, дБГц')
axes4[1].set_ylabel(r'$\Delta F_{\text{ЧАП}}$, Гц')
axes4[1].grid(True)
fig4.tight_layout(pad=2.0)
fig4.savefig("4.png", dpi=300)

# 5.png: Зависимости от sigma_a (2 графика)
fig5, axes5 = plt.subplots(1, 2, figsize=(12, 4.8))
axes5[0].plot(sigma_a_grid, sigma_omega_a, linewidth=2, color='tab:green')
axes5[0].set_title(r'СКО фильтрации частоты $\sigma_\Omega(\sigma_a)$ при $q_{c/n_0} = 30\text{ дБГц}$', fontsize=10)
axes5[0].set_xlabel(r'$\sigma_a$, $\text{м/с}^2$')
axes5[0].set_ylabel(r'$\sigma_\Omega$, рад/с')
axes5[0].grid(True)

axes5[1].plot(sigma_a_grid, dF_a, linewidth=2, color='tab:red')
axes5[1].set_title(r'Оптимальная полоса ЧАП $\Delta F_{\text{ЧАП}}(\sigma_a)$ при $q_{c/n_0} = 30\text{ дБГц}$', fontsize=10)
axes5[1].set_xlabel(r'$\sigma_a$, $\text{м/с}^2$')
axes5[1].set_ylabel(r'$\Delta F_{\text{ЧАП}}$, Гц')
axes5[1].grid(True)
fig5.tight_layout(pad=2.0)
fig5.savefig("5.png", dpi=300)

# ==========================================
# Пункт 8: Сравнение в установившемся режиме
# ==========================================
steady_start_idx = int(5.0 / T)
D11_discrete_steady = D11_hist[-1]
_, _, D11_analytic_steady = steady_state_params(q_sim_dBHz, sigma_a_sim)
D11_empirical = np.var(err_omega[steady_start_idx:])

print("=" * 60)
print("СРАВНЕНИЕ В УСТАНОВИВШЕМСЯ РЕЖИМЕ (t > 5 с):")
print(
    f"Теоретическая D11 (дискретный Калман): {D11_discrete_steady:.2f} (рад/с)^2  "
    f"(СКО = {np.sqrt(D11_discrete_steady):.2f} рад/с)"
)
print(
    f"Теоретическая D11 (непрерывная модель): {D11_analytic_steady:.2f} (рад/с)^2  "
    f"(СКО = {np.sqrt(D11_analytic_steady):.2f} рад/с)"
)
print(
    f"Выборочная дисперсия по реализации:    {D11_empirical:.2f} (рад/с)^2  "
    f"(СКО = {np.sqrt(D11_empirical):.2f} рад/с)"
)
print("=" * 60)

plt.show()
