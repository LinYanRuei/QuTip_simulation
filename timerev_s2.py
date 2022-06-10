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
print(a)

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


# exp_pulse 參數
args1 = {}
t_pulse = 2.63
V = 1.58
tau = 0.6
t0 = t_pulse
args1['gate time'] = t_pulse
args1['amp'] = V
args1['drive frequency'] = w_q
args1['characteristic time'] = tau


w_list = np.linspace(4.83, 4.87, 1000)
ref = []
ref = np.abs(1-Gamma/(r + 1j*(w_q-2*np.pi*w_list)))


plt.plot(w_list, ref)
plt.title('Magnitude of reflection to prob frequency')
plt.ylabel('|r|')
plt.xlabel('Wprob/2pi  (GHz)')
plt.show()

P_list = np.linspace(-165, -115, 1000)
k = 1.6*0.4*np.sqrt(50)/h_*(Ej/2/Ec)**0.25*10000  # theroitical k
ref2 = []
ref2 = np.abs(1-(Gamma**2/(Gamma*r+k*k*10**(P_list/10))))
Pin = -144
Omega = k*np.sqrt(10**(Pin/10))
print(Omega)
plt.plot(P_list, ref2)
plt.title('Magnitude of reflection to resonant power P')
plt.ylabel('|r|')
plt.xlabel('Pin  (dBm)')
plt.show()
