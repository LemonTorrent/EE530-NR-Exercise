import matplotlib.pyplot as plt
import numpy as np
import time

tolerances = [0.1, 0.01, 0.05, 0.001, 0.0001]
times = [0.0069, 0.0057, 0.0141, 0.0044, 0.0044] # From provided files
iterations = [2, 2, 2, 3, 4] # From provided files

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Tolerance ($\epsilon$) (log scale)')
ax1.set_ylabel('Processing Time (seconds)', color=color)
ax1.plot(tolerances, times, marker='o', color=color, label="Time")
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xscale('log')

ax2 = ax1.twinx()  
color = 'tab:blue'
ax2.set_ylabel('Iterations', color=color)  
ax2.plot(tolerances, iterations, marker='s', color=color, linestyle='--', label="Iterations")
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0, 4)
ax2.set_yticks([1, 2, 3, 4])

fig.tight_layout()  
plt.title('Newton-Raphson Performance vs Tolerance')
plt.savefig('nr_performance.png')
print("Plot saved to nr_performance.png")