"""
@author: Y R Lin @CCU PHY
"""
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import math
from sympy import *
t = Symbol('t')

# rabi freq = sqrt(4*g**2+(drive_freq-W_t)**2)
h_ = 1.054  # (h/2pi)
Nstates = 2  # Hilbert space size, choose N=2 for spin-1/2
a = destroy(Nstates)
a_dagger = a.dag()
sm = sigmaz()
w_q = 5*2*np.pi  # atom frequency = 5GHz
I_q = qeye(Nstates)  # I
ground = basis(Nstates, 0)
excited = basis(Nstates, 1)
ini = ground
Gamma = 1.686*2*np.pi*10**-3  # relaxation rate   #1/T1
alpha = -300*0.001*2*np.pi
Ej = (w_q+alpha)**2/8/alpha
k = 1.6*0.4*np.sqrt(50)/h_*(Ej/2/alpha)**0.25  # theroitical k

alpha = -300*0.001*2*np.pi
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2   # Hamiltonaian


def drive_pulse(t, args):      # 製作高斯脈衝
    t_g = args['gate time']  # pulse time
    amp = args['amp']
    w_d = args['drive frequency']
    tau = 15
    sigma = 4/2.355

    return amp*np.exp(-(t-tau)**2/2*(sigma**2))


def drive_pulse2(t, args):      # 製作高斯脈衝
    t_g = args['gate time']  # pulse time
    amp = args['amp']
    w_d = args['drive frequency']
    tau = t_g/2.0
    sigma = t_g/4.0

    return amp*np.exp(-(t-tau)**2/(sigma**2))


H_d = [a+a.dag(), drive_pulse]
H = [H_0, H_d]


# 設定高斯脈衝的參數 plot 1
args1 = {}
t_pulse = 30
steps = 10000
amps = 2.12
tau = 15
sigma = 4/2.355
args1['gate time'] = t_pulse
args1['amp'] = amps
args1['drive frequency'] = w_q
t_list = np.linspace(0, t_pulse, steps)

# 畫出Gaussian pulse
plt.plot(t_list, drive_pulse(t_list, args1))
plt.title('Gaussian Pulse')
plt.ylabel('Amplitude')
plt.xlabel('Time')
plt.show()

c = 0
for x in t_list:
    c = c+(amps*np.exp(-(x-tau)**2/(sigma**2))
           * np.cos(w_q * x))**2*t_pulse/steps
    # 2* c == pulse's energy
d = c/100
N = d/w_q/h_   # number of photon
print('N=', N)
print('amplitude = ', amps)

e = 0
rabi_max = []
amp_range = np.linspace(0, 40, 1000)


for item in amp_range:
    f = 0
    for x in t_list:
        if f < k*np.sqrt((item*np.exp(-(x-tau)**2/(sigma**2)))
                         ** 2/100):
            f = k*np.sqrt((item*np.exp(-(x-tau)**2/(sigma**2)))
                          ** 2/100)
        e = e + k*np.sqrt((item*np.exp(-(x-tau)**2/(sigma**2)))
                          ** 2/100)*t_pulse/steps
    rabi_max.append(f)

plt.plot(amp_range, rabi_max)
plt.title('Rabi experiment')
plt.ylabel('Mazimum of Rabi freq')
plt.xlabel('amplitude')
plt.show()
print(f)
print(e/30)
plt.plot(t_list, k*np.sqrt((drive_pulse(t_list, args1))**2/100))
plt.title('Rabi frequency to time')
plt.ylabel('Rabi frequency (GHz)')
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
