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
#t_pulse = 2.63

alpha = -300*0.001*2*np.pi
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2  # Hamiltonaian from Qiskit


def drive_pulse(t, args):      # 製作高斯脈衝
    t_g = args['gate time']  # pulse time
    amp = args['amp']
    w_d = args['drive frequency']
    tau = 15
    sigma = 4/2.355

    return amp*np.exp(-(t-tau)**2/2*(sigma**2)) * np.cos(w_d * t)


# 設定高斯脈衝的參數 plot 1
args1 = {}
t_pulse = 30
amp = 8
tau = 15
sigma = 4/2.355
args1['gate time'] = t_pulse
args1['amp'] = amp
args1['drive frequency'] = w_q
t_list = np.linspace(0, t_pulse, 10000)


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
plt.title('Population Evolution (Gaussian)')
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend(loc='lower right')
plt.show()

t = np.linspace(0, t_pulse, 1000)
k = 0
c = 0

for i in t:
    c = c + (amp*np.exp(-(i-tau)**2/2*(sigma**2)))**2*t_pulse/1000
N = c/w_q/h_/10  # 光子數
print('N=', N)
print('c = ', c)

# Plot 2
# 針對不同的脈衝強度做計算

result_population = []  # 紀錄結果的集合
ini = ground

# 設定高斯脈衝的參數 plot 2
amp_range = np.linspace(0, 2*amp, 100)  # 脈衝強度範圍
t_list1 = np.linspace(0, 2*t_pulse, 100)
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


z = 0
for i in t_list:
    z = amp*np.exp(-(i-tau)**2/2*(sigma**2))
    if z > 1.04 and z < 1.06:
        print('i', i)
