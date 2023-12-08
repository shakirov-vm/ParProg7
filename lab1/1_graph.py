
import matplotlib.pyplot as plt

# 1 is seq, other - par, 15k * 15k

# d = '<', '<'
# D = '1', '8'

x = [1, 2, 4]
time = [2.961830, 1.476772, 0.743388]

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