"""
@author: Y R Lin @CCU PHY
"""
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import math
from sympy import *
t = Symbol('t')

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

alpha = -300*0.001*2*np.pi
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2   # Hamiltonaian


def step_func(t):
    return (-1*np.tanh(1000*(t-t_pulse))+1)/2


def exp_pulse(t, args):
    t0 = args['gate time']  # pulse time
    w_d = args['drive frequency']
    V = args['amp']
    tau = args['characteristic time']

    return V*np.exp((t-t0)/tau) * (np.cos(w_d * t))*step_func(t)


H_d = [a+a.dag(), exp_pulse]
H = [H_0, H_d]
ini = ground


# 設定脈衝的參數
# exp_pulse 參數
args3 = {}
t_pulse = 30
V = 2.12
tau = 1.4
t0 = t_pulse
args3['gate time'] = t_pulse
args3['amp'] = V
args3['drive frequency'] = w_q
args3['characteristic time'] = tau
### N = 0.09
# 繪圖
t_list = np.linspace(0, 2*t_pulse, 10000)
plt.plot(t_list,  exp_pulse(t_list, args3))
plt.title('Exponential Pulse')
plt.ylabel('Amplitude (nV)')
plt.xlabel('Time (ns)')
plt.show()

result = mesolve(H, ini, t_list, [], args=args3)

# 觀測量 plot 1
state0 = ground*ground.dag()
state1 = excited*excited.dag()
result0 = expect(state0, result.states)
result1 = expect(state1, result.states)
#################
c = 0
t = np.linspace(0, t_pulse, 1000)
for i in t:
    c = c + (V*np.exp((i-t0)/tau))**2*t_pulse/1000  # * (np.cos(w_q * i))

N = c/w_q/h_/10  # 光子數
print('c = ', c)  # c = energy*R(100歐姆)
print('N=', N)
######################
# Plot 1
plt.plot(t_list, result0, label='Ground')
plt.plot(t_list, result1, label='Excited')
plt.title('Population Evolution (Exponential)')
plt.xlabel('Time (ns)')
plt.ylabel('Population')
plt.legend(loc='right')
plt.show()


# 設定脈衝的參數 plot 2
result_population = []  # 紀錄結果的集合
ini = ground

amp_range = np.linspace(0, 2*V, 100)  # 脈衝強度範圍
t_list1 = np.linspace(0, t_pulse, 100)
args2 = {}
args2['gate time'] = t_pulse
args2['amp'] = V
args2['drive frequency'] = w_q
args2['characteristic time'] = tau
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
plt.xlabel('Amplititude (nV)')
plt.ylabel('Population')
plt.legend(loc='upper right')
plt.show()
