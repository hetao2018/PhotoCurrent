import numpy as np
from Two_color_laser import *
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftshift


def W(E):
    '''电离率 经验公式'''
    Ee = 0.0000001 + abs(E)             #  防止分母为零
    Ei, Eh = 15.6, 13.6       ##Eh and ⁢𝐸i are the ionization potentials of hydrogen and the atom in question.
    W = 4 * (Ei / Eh) ** 2.5 * (1 / Ee) * np.exp(-2 / 3 * (Ei / Eh) ** 1.5 * (1 / Ee))
    return W

def Ne(E):
    '''等离子体浓度'''
    Ng = 2.4 * 10 ** 19 * (0.5292 * 10 ** (-8)) ** 3   #initial density of neutral particles
    Ne = np.zeros(len(t),float)
    Ne[0] = Ng * W(E)[0]
    # for i in range(len(time)-1):

    for i in range(len(t) - 1):
        Ne[i + 1] = Ne[i] + (Ng - Ne[i]) * W(E[i + 1]) * (t[i+1] - t[i])
    return Ne

def Je(E):
    dJ_dt = Ne(E) * E
    J = np.zeros(len(t),float)
    J[0] = Ne(E)[0] * E[0]
    for i in range(len(t)-1):
        J[i+1] = J[i] + dJ_dt[i] * (t[i+1]-t[i])
    return J

def Fourier(E):
    dJ_dt = Ne(E) * E
    num = 10 * len(t)
    dJ_dt = np.pad(dJ_dt, (num, num), 'constant')
    FFt = fft(dJ_dt)
    return fftshift(FFt)

plt.xlim(90000,120000)
co_rotating = co_P(t) + co_S(t)
counter_rotating = counter_P(t) + counter_S(t)
plt.plot(abs(Fourier(counter_P(t))))
plt.show()

