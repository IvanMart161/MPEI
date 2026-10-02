import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as signal

def psd(x, nfft = 1024, fs = 1.0): 
    f, S = signal.welch(s, fs = fs, window ='blackmanharris', nperseg = nfft, nfft = nfft, noverlap = int(nfft/2), \
             detrend = False, return_onesided = False)

    f = np.fft.fftshift(f)
    S = 10*np.log10(np.fft.fftshift(S))
    return f, S

N = 65536*8

t = np.arange(N)
noise = 0.1 * np.random.randn(N)

s = np.cos(2 * np.pi * 0.2 * t) + 0.1 * np.cos(2 * np.pi * 0.26 *t) + noise

f, S = psd(s, fs = 50)

plt.figure()
plt.plot(f, S)
plt.show()
