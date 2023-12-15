
import matplotlib.pyplot as plt

# max_iter = 10000

x = [1, 2, 3, 4, 5, 6]
time = [840, 615, 495, 455, 410, 380]

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