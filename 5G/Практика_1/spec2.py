import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as signal

ORD = 64 
Fs = 10

h = signal.firwin(ORD+1, 0.4, window = 'hamming')

w = np.linspace(-np.pi, np.pi, 1000)

f, H = signal.freqz(h, 1, worN = w)

f = f * Fs / (2 * np.pi)

# plt.plot(f, np.abs(H))

Phi = np.unwrap(np.angle(H))

plt.figure()
plt.plot(f, 10*np.log10(np.abs(H)**2))
plt.grid()
plt.show()

plt.figure()
plt.plot(f, Phi)
plt.grid()
plt.show()
