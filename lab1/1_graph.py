
import matplotlib.pyplot as plt

# 1 is seq, other - par, 5k * 5k

# d = '<', '<'
# D = '1', '8'

x = [1, 2, 3, 4]
time = [0.32, 0.22, 0.19, 0.18]

eff = [time[0] / (time[i] * x[i]) for i in range (len(x))]

plt.plot(x, eff)
plt.xlabel('cores')
plt.ylabel('time in s')
plt.title('efficiency')
plt.show()

speed = [time[0] / (time[i]) for i in range (len(x))]

plt.plot(x, speed)
plt.xlabel('cores')
plt.ylabel('time in s')
plt.title('speed')
plt.show()