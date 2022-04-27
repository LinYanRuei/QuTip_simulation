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


def step_func(t):
    return (-1*np.tanh(1000*(t-t_pulse))+1)/2


def exp_pulse(t, args):      # 製作脈衝
    t_g = args['gate time']  # pulse time
    w_d = args['drive frequency']
    amp = args['amp']
    sigma = t_g/4.0

    return amp*np.exp(t/sigma) * np.cos(w_d * t)*step_func(t)


H_d = [a+a.dag(), exp_pulse]
H = [H_0, H_d]
ini = ground


# 設定脈衝的參數
args3 = {}
t_pulse = 35
args3['gate time'] = t_pulse
args3['amp'] = 1
args3['drive frequency'] = w_q

# 繪圖
t_list = np.linspace(0, 2*t_pulse, 10000)
plt.plot(t_list, step_func(t_list))
plt.show()
plt.plot(t_list,  exp_pulse(t_list, args3))
plt.title('Exponential Pulse')
plt.ylabel('Amplitude')
plt.xlabel('Time')
plt.show()

result = mesolve(H, ini, t_list, [], args=args3)

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


# 設定脈衝的參數 plot 2
result_population = []  # 紀錄結果的集合
ini = ground

amp_range = np.linspace(0, 0.1*np.pi, 100)  # 脈衝強度範圍
t_list1 = np.linspace(0, t_pulse, 100)
args2 = {}
args2['gate time'] = t_pulse
args2['drive frequency'] = w_q

state2 = excited*excited.dag()
for item in amp_range:
    args2['amp'] = item
    result = mesolve(H, ini, t_list1, [], args=args2)
    result_population.append(expect(state2, result.states[-1]))
    if expect(state2, result.states[-1]) > 0.99:
        print(item)


# Plot 2
plt.plot(amp_range, result_population, label='Excited')
plt.title('Rabi experiment')
plt.xlabel('Amplititude')
plt.ylabel('Population')
plt.legend(loc='upper right')
plt.show()
