import numpy as np
import matplotlib.pyplot as plt

Fd = 33e6               
f0 = 8e6                
T = 1e-3                
N = int(T * Fd)         

k = np.arange(N)
t = k / Fd

code_str = (
    "1111111000001111011111000101110011001000001001010011101101000"
    "111100111110011011000101010010001110001101101010111000100110001000"
    "100000000100001000110000100111001010101100001101111010011011100100"
    "010100001010110100111111011001001001011011111100100110101001100110"
    "000000110001100101000110100101111111010001011000111010110010110011"
    "110001111101110100000110101101101110110000010110101111101010101000"
    "000101001010111100101110111000000111111010010011110101110101000111"
    "001000011001110000101111011011001101000011101111000011"
)

dk_symbols = np.array([1.0 if b == '0' else -1.0 for b in code_str])
chip_indices = np.clip(np.floor(k * 511.0 / N).astype(int), 0, 510)
g_ref = dk_symbols[chip_indices]

y = np.loadtxt("Yvar8.txt")
if len(y) > N:
    y = y[:N]

num_tau = 1022
tau_ms = np.linspace(0.0, 1.0, num_tau, endpoint=False)  
f_dop = np.linspace(-2000.0, 2000.0, 9)                  
G_conj = np.conj(np.fft.fft(g_ref))
X2_matrix = np.zeros((len(f_dop), num_tau))

sample_tau_indices = np.round(np.linspace(0, N - 1, num_tau, endpoint=False)).astype(int)

for i, fd in enumerate(f_dop):
    s_demod = y * np.exp(-1j * 2.0 * np.pi * (f0 + fd) * t)
    
    corr = np.fft.ifft(np.fft.fft(s_demod) * G_conj)
    
    X2_full = np.abs(corr)**2
    X2_matrix[i, :] = X2_full[sample_tau_indices]

max_val = np.max(X2_matrix)
fd_idx, tau_idx = np.unravel_index(np.argmax(X2_matrix), X2_matrix.shape)
est_fd = f_dop[fd_idx]
est_tau = tau_ms[tau_idx]

noise_var = np.median(X2_matrix) / np.log(2)
P_fa = 1e-4
h_threshold = -noise_var * np.log(P_fa)

print("\n" + "=" * 50)
print(f"Максимум корреляционного отклика: {max_val:.2e}")
print(f"Порог обнаружения (P_fa = 1e-4):  {h_threshold:.2e}")
print(f"Оценка доплеровской частоты f_д:  {est_fd:.1f} Гц")
print(f"Оценка задержки сигнала tau_з:    {est_tau:.4f} мс")
print("=" * 50 + "\n")

TAU, FD = np.meshgrid(tau_ms, f_dop)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(TAU, FD, X2_matrix, cmap='coolwarm', edgecolor='none', alpha=0.9)

ax.set_xlabel(r'$\tau_з$, мс', fontsize=11, labelpad=10)
ax.set_ylabel(r'$f_д$, Гц', fontsize=11, labelpad=10)
ax.set_zlabel(r'$X^2$', fontsize=11, labelpad=10)
ax.set_title(f'Поле поиска сигнала ГЛОНАСС L1OF (Вариант 8)\n$f_д = {est_fd:.0f}$ Гц, $\\tau_з = {est_tau:.4f}$ мс', fontsize=12)

ax.view_init(elev=35, azim=-125)
fig.colorbar(surf, shrink=0.5, aspect=10)

plt.tight_layout()
plt.savefig("search_surface_var8.png", dpi=300)
plt.show()
