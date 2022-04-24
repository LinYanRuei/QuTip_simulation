import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import math
from sympy import *
t = Symbol('t')

# rabi freq = sqrt(4*g**2+(drive_freq-W_t)**2)
Nstates = 2  # Hilbert space size, choose N=2 for spin-1/2
a = destroy(Nstates)
a_dagger = a.dag()
sm = sigmaz()
w_q = 5*2*np.pi  # cavity frequency = 5GHz
wa = 1.0 * 2 * np.pi  # atom frequency
g = 0.05 * 2 * np.pi  # coupling strength
W_t = 1.9*2*np.pi  # transition frequency
kappa = 0.005          # cavity dissipation rate
gamma = 0.05           # atom dissipation rate
I_q = qeye(Nstates)  # I
ground = basis(Nstates, 0)
excited = basis(Nstates, 1)


alpha = -300*0.001*2*np.pi
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2   # Hamiltonaian
print(a)


def step_func(t):  # step function
    return (-1*np.tanh(1000*t)+1)/2


def exp_pulse(t, args):      # 製作脈衝
    t_g = args['gate time']  # pulse time
    w_d = args['drive frequency']
    amp = args['amp']
    return np.exp(amp*t) * np.cos(w_d * t)*step_func(t)


H_d = [a+a.dag(), exp_pulse]
H = [H_0, H_d]
ini = ground


# 設定脈衝的參數 plot 1
args1 = {}
t_pulse = 30
args1['gate time'] = t_pulse
args1['amp'] = 0.32
args1['drive frequency'] = w_q
t_list = np.linspace(-t_pulse, t_pulse, 10000)

# 畫出expnential pulse
plt.plot(t_list, exp_pulse(t_list, args1))
plt.title('Gaussian Pulse')
plt.ylabel('Amplitude')
plt.xlabel('Time')
plt.show()

result = mesolve(H, ini, t_list, [], args=args1)

# 觀測量 plot 1
state0 = ground*ground.dag()
state1 = excited*excited.dag()
result0 = expect(state0, result.states)
result1 = expect(state1, result.states)

# Plot 1
plt.plot(t_list, result0, label='Ground')
plt.plot(t_list, result1, label='Excited')
plt.title('Population Evolution')
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend(loc='right')
plt.show()
