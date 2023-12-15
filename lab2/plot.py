
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

fig = plt.figure(figsize = [12, 8])

df = pd.read_csv("result_para.txt", sep = " ")
x = np.array(df.iloc[:, 0].tolist())

b_max = 10
for i in range(0, b_max):
    y = np.array(df.iloc[:, 1 + i].tolist())
    plt.plot(x, y, label = f'b = {i / 10}', linewidth = 2)

plt.legend(loc = "best")
plt.grid()

plt.xlabel("x")
plt.ylabel("y")

plt.show()