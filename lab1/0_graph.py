
import matplotlib.pyplot as plt

# d = '<', '<'
# D = '1', '8'

x = [1, 2, 3, 4]
time = [2.77, 1.40, 0.94, 0.74]

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