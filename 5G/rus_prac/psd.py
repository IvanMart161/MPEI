# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 11:17:07 2021

@author: Main_PC
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal


"""
Функция расчета СПМ методом коррелограммы 
Внимание! Функция не оптимизирована! Вычислительно не использует FFT
"""    
def psd_correlogram(x, D, w, nfft, fs):
    N = len(x)

    rxx = np.zeros(2*D+1, dtype = 'complex')


    # расчет значений автокорреляционной функции для положительных m
    for m in range(D+1):
        ind = np.arange(0, N-m-1, 1, dtype = 'int') 
        rxx[m+D] = np.sum(x[ind+m]*np.conj(x[ind]))/(N-m)
        
    
    # расчет значений автокорреляционной функции для отрицательных m    
    for m in np.arange(-D, 0, 1):
        ind = np.arange(0, N+m-1, 1, dtype = 'int')  
        rxx[m+D] = np.sum(np.conj(x[ind-m])*x[ind])/(N+m)    

    # Оконное взвешивание
    rxx = rxx*w
    
    # расчет СПМ
    X = np.zeros(nfft, dtype = 'complex')
    frq = fs * np.arange(-nfft/2, nfft/2+1, 1) / nfft
    m = np.arange(-D,D+1,1)
    X = np.zeros(nfft+1, dtype = 'complex')
    
    for m in range(2*D+1):
        X = X + rxx[m]*np.exp(-1j*2*np.pi*m*frq / fs) / fs
        
    X = 10.0*np.log10(np.abs(X))
    return X, frq
    


    

"""
Функция расчета модифицированной периодограммы случайного процесса
"""
def psd_periodogram(x, w, fs, log_flag = True):
    n = len(x)
    wn =np.sum(w**2)

    X = (np.abs(np.fft.fft(x*w))**2) / (wn * fs)
    X = np.fft.fftshift(X)
    
    if(log_flag):
        X = 10.0 * np.log10(X)

    frq = fs * np.arange(-n/2, n/2, 1) / n
    
    return (X, frq)



"""
Функция расчета периодограммы Бартлетта случайного процесса
"""
def psd_bartlett(x, nfft, fs):
    n = len(x)
    w = np.ones((nfft))
    p = 0
    X = 0
    c = 0
    while(p+nfft <= n):
        (X0, frq)  = psd_periodogram(x[p:p+nfft], w, fs, log_flag = False)
        X = X + X0
        p = p + nfft
        c = c + 1
    
    if(p < n):
        tmp = np.zeros(nfft)
        tmp[0:n-p] = x[p:n]
        (X0, frq)  = psd_periodogram(tmp, w, fs, log_flag = False)
        X = X + X0
        c = c + 1
    
    X = 10.0 * np.log10(X / c)    
    return (X, frq)



"""
Функция расчета периодограммы Уэлча случайного процесса
"""
def psd_welch(x, nfft, w, shift, fs):
    n = len(x)
    p = 0
    X = 0
    c = 0
    while(p+nfft <= n):
        (X0, frq)  = psd_periodogram(x[p:p+nfft], w, fs, log_flag = False)
        X = X + X0
        p = p + shift
        c = c + 1
    
    if(p < n):
        tmp = np.zeros(nfft, dtype=complex)
        tmp[0:n-p] = x[p:n]
        (X0, frq)  = psd_periodogram(tmp, w, fs, log_flag = False)
        X = X + X0
        c = c + 1
    
    X = 10.0 * np.log10(X / c)    
    return (X, frq)