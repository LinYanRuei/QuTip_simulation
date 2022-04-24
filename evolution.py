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
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2   # Hamiltonaian from Qiskit


def drive_pulse(t, args):      # 製作高斯脈衝
    t_g = args['gate time']  # pulse time
    amp = args['amp']
    w_d = args['drive frequency']
    tau = t_g/2.0
    sigma = t_g/4.0

    return amp*np.exp(-(t-tau)**2/(sigma**2)) * np.cos(w_d * t)


# 設定高斯脈衝的參數 plot 1
args1 = {}
t_pulse = 30
amp = 0.238
args1['gate time'] = t_pulse
args1['amp'] = amp
args1['drive frequency'] = w_q
t_list = np.linspace(0, t_pulse, 10000)

# 畫出Gaussian pulse
t_list = np.linspace(0, 30, 50000)
plt.plot(t_list, drive_pulse(t_list, args1))
plt.title('Gaussian Pulse')
plt.ylabel('Amplitude')
plt.xlabel('Time')
plt.show()

H_d = [a+a.dag(), drive_pulse]
H = [H_0, H_d]
ini = ground
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


# Plot 2
# 針對不同的脈衝強度做計算

result_population = []  # 紀錄結果的集合
ini = ground

# 設定高斯脈衝的參數 plot 2
t_pulse = 30
amp_range = np.linspace(0, 0.5*np.pi, 100)  # 脈衝強度範圍
t_list1 = np.linspace(0, t_pulse, 100)
args2 = {}
args2['gate time'] = t_pulse
args2['drive frequency'] = w_q

state2 = excited*excited.dag()
for item in amp_range:
    args2['amp'] = item
    result = mesolve(H, ini, t_list1, [], args=args2)
    result_population.append(expect(state2, result.states[-1]))
    if expect(state2, result.states[-1]) > 0.999:
        print(item)


# Plot 2
plt.plot(amp_range, result_population, label='Excited')
plt.title('Rabi experiment')
plt.xlabel('Amplititude')
plt.ylabel('Population')
plt.legend(loc='upper right')
plt.show()
