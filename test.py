import numpy as np
import matplotlib.pyplot as plt

labels = ['C_r 03', 'C_r 05', 'C_r 0.1', 'C_r 0.2', 'C_r 0.5', 'C_r 1', 'C_r 2', 'Unconfined']

my_list = [np.random.uniform(10, 30, 5) for _ in labels]
my_mean = [values.mean() for values in my_list]

plt.plot(np.arange(len(my_mean)) + 1, my_mean, color='r')
print (my_list, my_mean)
plt.boxplot(my_list, labels=labels)
plt.show()