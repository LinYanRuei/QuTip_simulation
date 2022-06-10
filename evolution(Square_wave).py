"""
@author: Y R Lin @CCU PHY
"""
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from scipy import integrate, special


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


w_q = 5*2*np.pi
alpha = -300*0.001*2*np.pi
H_0 = w_q*a*a.dag()+alpha*a*a*a.dag()*a.dag()/2   # Hamiltonaian


def square_pulse(t, args):
    t_g = args['gate time']  # pulse time
    w_d = args['drive frequency']
    amp = args['amps']
    return amp*((1*np.tanh(10000*(t))+1)/2-(1*np.tanh(10000*(t-t_g))+1)/2)*np.cos(w_d*t)


H_d = [a+a.dag(), square_pulse]
H = [H_0, H_d]
ini = ground
args_test = {}
t_pulse = 1.4
amps = 2.14
args_test['gate time'] = t_pulse
args_test['drive frequency'] = w_q
args_test['amps'] = amps

t_list = np.linspace(-10, 40, 10000)
plt.plot(t_list, square_pulse(t_list, args_test))
plt.title('Square wave Pulse')
plt.ylabel('Amplitude')
plt.xlabel('Time')
plt.show()
#######
# Plot 1 H
result = mesolve(H, ini, t_list, [], args=args_test)

# 觀測量 plot 1
state0 = ground*ground.dag()
state1 = excited*excited.dag()
result0 = expect(state0, result.states)
result1 = expect(state1, result.states)

# Plot 1
plt.plot(t_list, result0, label='Ground')
plt.plot(t_list, result1, label='Excited')
plt.title('Population Evolution (square wave)')
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend(loc='right')
plt.show()
###############
t = np.linspace(0, t_pulse, 1000)
k = 0
c = 0
for i in t:
    c = c + (amps*((1*np.tanh(10000*(i))+1)/2 -
             (1*np.tanh(10000*(i-t_pulse))+1)/2))**2*t_pulse/1000
N = c/w_q/h_/10  # 光子數
print('N=', N)
print('c = ', c)

# Plot2
result_population = []  # 紀錄結果的集合
ini = ground

# 設定脈衝的參數 plot 2
amp_range = np.linspace(0, np.pi, 100)  # 脈衝強度範圍
t_list1 = t_list
args2 = {}
args2['gate time'] = t_pulse
args2['drive frequency'] = w_q
state2 = excited*excited.dag()
for item in amp_range:
    args2['amps'] = item
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
