import pandas as pd
import matplotlib.pyplot as plt

cDF = pd.read_csv('c_times.csv')
fDF = pd.read_csv('fortran_times.csv')

fig, axs = plt.subplots(1,2,figsize=(12,4), sharey=False)

cDF.plot.line(x = 'n', y = 'Time_ij', ax = axs[0], ylabel='Time (seconds)')
cDF.plot.line(x = 'n', y = 'Time_ji', ax = axs[0])
axs[0].set_title('C',fontsize=15)

fDF.plot.line(x = 'n', y = 'Time_ij', ax = axs[1], ylabel='Time (seconds)')
fDF.plot.line(x = 'n', y = 'Time_ji', ax = axs[1])
axs[1].set_title('Fortran',fontsize=15)

# plt.show()
fig.savefig('graph.jpeg')
