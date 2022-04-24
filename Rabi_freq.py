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
ini = ground


alpha = -300*0.001*2*np.pi
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2   # Hamiltonaian


def drive_pulse(t, args):      # 製作高斯脈衝
    t_g = args['gate time']  # pulse time
    amp = args['amp']
    w_d = args['drive frequency']
    tau = t_g/2.0
    sigma = t_g/4.0

    return amp*np.exp(-(t-tau)**2/(sigma**2)) * np.cos(w_d * t)


H_d = [a+a.dag(), drive_pulse]
H = [H_0, H_d]


# 設定高斯脈衝的參數 plot 1
args1 = {}
t_pulse = 30
amps = 35
steps = 10000
t_g = t_pulse
tau = t_g/2.0
sigma = t_g/4.0
args1['gate time'] = t_pulse
args1['amp'] = amps
args1['drive frequency'] = w_q
t_list = np.linspace(0, t_pulse, steps)

# 畫出Gaussian pulse
t_list = np.linspace(0, 30, 50000)
plt.plot(t_list, drive_pulse(t_list, args1))
plt.title('Gaussian Pulse')
plt.ylabel('Amplitude')
plt.xlabel('Time')
# plt.show()

c = 0
for x in t_list:
    c = c+amps*np.exp(-(x-tau)**2/(sigma**2)) * amps * \
        np.exp(-(x-tau)**2/(sigma**2)) * t_pulse/steps
    # c == pulse's energy

n = c/wa  # number of photon
rabi_freq = sqrt(4*(n+1)*alpha**2+(w_q-wa)**2)

print('c=', c)
print('n=', n)
print('amp=', amps, 'rabi freq=', rabi_freq/69.37385873)


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
# plt.show()
