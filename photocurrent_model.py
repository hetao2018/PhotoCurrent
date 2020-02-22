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
    # num = 10 * len(t)
    # dJ_dt = np.pad(dJ_dt, (num, num), 'constant')
    FFt = fft(E)
    return fftshift(FFt)

co_rotating,counter_rotating = np.zeros(t.size,dtype=complex),np.zeros(t.size,dtype=complex)
co_rotating.real = co_P(t)
co_rotating.imag = co_S(t)
counter_rotating.real = counter_P(t)
counter_rotating.imag = counter_S(t)

dJ_dt_co_P = Ne(co_rotating) * co_P(t)
dJ_dt_co_S = Ne(co_rotating) * co_S(t)
dJ_dt_counter_P = Ne(counter_rotating) * counter_P(t)
dJ_dt_counter_S = Ne(counter_rotating) * counter_S(t)


plt.plot(Fourier(dJ_dt_co_S))
plt.show()



'''
co_theta = np.angle(co_rotating)
co_abs = abs(co_rotating)
counter_theta = np.angle(counter_rotating)
counter_abs = abs(counter_rotating)
Size = 12
plt.figure(figsize=(10,5))
plt.subplot(121,polar=True)
plt.plot(co_theta,co_abs,linewidth=1.5,color="y",label="Co-rotating")
plt.rgrids(np.arange(0.01,1.1*max(co_abs),0.03),angle=180)
plt.legend(loc=(0.7,1),prop={'size': Size})
plt.subplot(122,polar=True)
plt.plot(counter_theta,counter_abs,linewidth=1.5,color="c",label="Counter-rotating")
plt.legend(loc=(0.7,1),prop={'size': Size})
plt.rgrids(np.arange(0.01,1.1*max(counter_abs),0.03),angle=180)
plt.tight_layout(pad=1.0, w_pad=4)
plt.show()
'''


