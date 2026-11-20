import scipy.io as sio, numpy as np, matplotlib.pyplot as plt, scipy.signal as signal
import psd

# ================= Задание 1: синхронизация x и y =================
mat = sio.loadmat('test1_25p0.mat')
x = mat['pdin'].ravel()
y = mat['pdout'].ravel()
Fs = 245.76e6
df = float(mat['rx_freq'][0, 0]) - float(mat['tx_freq'][0, 0])

y = y - np.mean(y)                                           # 1. убрать постоянную составляющую
y = y * np.exp(2j*np.pi*df/Fs*np.arange(len(y)))             # 2. убрать частотный сдвиг
delay = np.argmax(np.abs(np.correlate(y, x, 'valid')))       # 3. временное выравнивание
y = y[delay:delay + len(x)]
y = y * (np.vdot(y, x) / np.vdot(y, y))                      # 4. комплексный коэффициент g

scale = np.sqrt(np.mean(np.abs(x)**2))                       # нормируем мощность x к 1,
x, y = x/scale, y/scale                                      # чтобы σ² и µ были "безразмерными"

# ================= Задание 2 =================
d = y - x                       # опорный сигнал: ошибка тракта (линейные + нелинейные искажения)
D = 7                           # фильтр на 2D+1 = 15 коэффициентов, центр модели в D

# Матрица U: строка n = [x[n+D], x[n+D-1], ..., x[n-D]]  (вектор состояния u_n)
xp = np.concatenate([np.zeros(D), x, np.zeros(D)])
U = np.array([xp[n:n + 2*D + 1][::-1] for n in range(len(x))])

nmse = lambda e: 10*np.log10(np.vdot(e, e).real / np.vdot(d, d).real)   # ошибка относительно d, дБ

# ---------- LS (7.26) и LS с регуляризацией Тихонова (8.6) ----------
UH = U.conj().T
w_ls = np.linalg.solve(UH @ U, UH @ d)                       # w = (U^H U)^-1 U^H d

sigma2 = 1e-2 * len(x)          # σ² (для суммы по N отсчётов; это 1e-2 на один отсчёт)
w_reg = np.linalg.solve(UH @ U + sigma2*np.eye(2*D + 1), UH @ d)   # w = (U^H U + σ² I)^-1 U^H d

# ---------- LMS с утечкой (8.15): w = (1-λ) w + µ u* e ----------
mu = 0.003                      # шаг адаптации
lam = mu * 1e-2                 # λ = µσ²
w = np.zeros(2*D + 1, complex)
for epoch in range(20):          # данных мало, поэтому проходим по ним несколько раз
    for n in range(len(x)):
        e_n = d[n] - U[n] @ w                   # ошибка с ТЕКУЩИМИ коэффициентами
        w = (1 - lam)*w + mu*np.conj(U[n])*e_n  # шаг стохастического градиента
w_lms = w

# Проверка: коэффициенты зафиксированы, ошибка считается по всему блоку (как в лекции про переобучение)
for name, w in [('LS', w_ls), ('LS + регуляризация', w_reg), ('LMS с утечкой', w_lms)]:
    print(f'{name:20s} NMSE = {nmse(d - U@w):6.2f} дБ,  ||w|| = {np.linalg.norm(w):.3f}')

# ================= Графики =================
w = w_reg                         # итоговый фильтр
e = d - U @ w                     # остаток после компенсации

spec = lambda s: psd.psd_welch(s*scale, 2048, signal.windows.blackmanharris(2048), 1024, Fs)
H = lambda w: 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(w, 2048))))

plt.figure('PSD')
for s in (x, y, d, e):
    P, f = spec(s)
    plt.plot(f/1e6, P)
plt.plot(f/1e6, H(w))                                        # АЧХ корректора |W(f)|, дБ
plt.legend(['x', 'y', 'd = y - x', 'остаток e', '|W(f)|'])
plt.xlabel('freq, MHz'); plt.ylabel('PSD, dB/Hz'); plt.grid(True)

plt.figure('Коэффициенты')
plt.subplot(211); plt.title('Без регуляризации'); plt.plot(w_ls.real); plt.plot(w_ls.imag)
plt.grid();
plt.subplot(212); plt.title('С регуляризацией');  plt.plot(w_reg.real); plt.plot(w_reg.imag)
plt.grid();
plt.show()
