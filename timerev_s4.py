"""
@author: Y R Lin @CCU PHY
"""
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import math
from sympy import *

h_ = 1.054  # (h/2pi)
Nstates = 2  # Hilbert space size, choose N=2 for spin-1/2
a = destroy(Nstates)
a_dagger = a.dag()
sm = sigmaz()
Ec = 385*10**-3*np.pi*2
Ej = 8.9*np.pi*2  # Josephson energy
r = 0.956*2*np.pi*10**-3  # decoherence rate    #1/T2
Gamma = 1.686*2*np.pi*10**-3  # relaxation rate   #1/T1
w_q = np.sqrt(8*Ej*Ec)-Ec  # transition frequency = 4.8GHz
I_q = qeye(Nstates)  # I
ground = basis(Nstates, 0)
excited = basis(Nstates, 1)
# theroitical k   atom-field coupling constant
k = 1.6*0.4*np.sqrt(50)/h_*(Ej/2/Ec)**0.25*10000
t_pulse = 2.63

alpha = Ec
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2  # Hamiltonaian

# 定義exp_pulse


def step_func(t):
    return (-1*np.tanh(1000*(t-t_pulse))+1)/2


def exp_pulse(t, args):
    t0 = args['gate time']  # pulse time
    w_d = args['drive frequency']
    V = args['amp']
    tau = args['characteristic time']

    return V*np.exp((t-t0)/tau) * step_func(t)


def omega(t, args):
    return np.sqrt((exp_pulse(t, args))**2/100)*k/4/1j


def kappa(t, args):
    return 0.5*np.sqrt(r**2/4-16*(np.abs(omega(t, args))**2))


def sigmam(t, args):
    return -2*1j*omega(t, args)*(r/(r**2+8*omega(t, args)**2))*(1-np.exp(-3*r*t/4)*(np.cosh(kappa(t, args)*t)+(kappa(t, args)/r+3*r/16/kappa(t, args))*np.sinh(kappa(t, args)*t)))


def exp_out(t, args):
    t0 = args['gate time']  # pulse time
    w_d = args['drive frequency']
    V = args['amp']
    tau = args['characteristic time']
    P = (V*np.exp((t-t0)/tau))**2 / 2/50
    return V*np.exp((t-t0)/tau) * step_func(t) - np.sqrt(Gamma)*sigmam(t, args)


# exp_pulse 參數
args1 = {}
t_pulse = 2.63
V = 9.8
tau = 0.6
t0 = t_pulse
args1['gate time'] = t_pulse
args1['amp'] = V
args1['drive frequency'] = w_q
args1['characteristic time'] = tau

t_list = np.linspace(0, 2*t_pulse, 50000)
# exp_pulse 參數
args2 = {}
t_pulse = 2.63
V = 15.76
tau = 0.23
t0 = t_pulse
args2['gate time'] = t_pulse
args2['amp'] = V
args2['drive frequency'] = w_q
args2['characteristic time'] = tau

# exp_pulse 參數
args3 = {}
t_pulse = 2.63
V = 18.32
tau = 0.17
t0 = t_pulse
args3['gate time'] = t_pulse
args3['amp'] = V
args3['drive frequency'] = w_q
args3['characteristic time'] = tau

# exp_pulse 參數
args4 = {}
t_pulse = 2.63
V = 36.84
tau = 0.04
t0 = t_pulse
args4['gate time'] = t_pulse
args4['amp'] = V
args4['drive frequency'] = w_q
args4['characteristic time'] = tau

# 畫出exp_pulse
t_list = np.linspace(0, 2*t_pulse, 50000)
plt.plot(t_list, exp_pulse(t_list, args1), label='tau = 600ns')
plt.plot(t_list, exp_pulse(t_list, args2), label='tau = 230ns')
plt.plot(t_list, exp_pulse(t_list, args3), label='tau = 170ns')
plt.plot(t_list, exp_pulse(t_list, args4), label='tau = 40ns')
plt.legend(loc='upper right')
plt.title('Voffres to t (Fig.2a)')
plt.ylabel('Voffres (nV)')
plt.xlabel('Time (micros)')
plt.show()


t_list = np.linspace(0, 3*t_pulse, 50000)
plt.plot(t_list, np.abs(sigmam(t_list, args1)), label='tau = 600ns')
plt.plot(t_list, np.abs(exp_out(t_list, args2)), label='tau = 230ns')
plt.plot(t_list, exp_out(t_list, args3), label='tau = 170ns')
plt.plot(t_list, exp_out(t_list, args4), label='tau = 40ns')
plt.legend(loc='upper right')
plt.title('Vres to t (Fig.2a)')
plt.ylabel('Vres (nV)')
plt.xlabel('Time (micros)')
plt.show()


v_list = np.linspace(0, 50, 1000)
t = np.linspace(0, t0, 1000)
k = 0
c = 0
for v in v_list:
    c = 0
    for i in t:
        c = c + (v*np.exp((i-t0)/tau))**2*t0/1000
    N = c/w_q/h_/10  # 光子數
    if N > 0.0899 and N < 0.0901:
        k = v
        print(k)  # control N 計算v的值
