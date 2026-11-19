import matplotlib.pyplot as plt

tolerances = [0.1, 0.01, 0.05, 0.001, 0.0001]
times = [0.0069, 0.0057, 0.0141, 0.0044, 0.0044] # From provided files
iterations = [2, 2, 2, 3, 4] # From provided files

data = sorted(zip(tolerances, times, iterations))
sorted_tolerances, sorted_times, sorted_iterations = zip(*data)

fig, ax = plt.subplots()
ax.plot(sorted_tolerances, sorted_times, marker='o', color='tab:red')
ax.set_xlabel('Tolerance ($\\epsilon$) (log scale)')
ax.set_ylabel('Processing Time (seconds)')
ax.set_xscale('log')
ax.set_title('Newton-Raphson Processing Time vs Tolerance')
fig.tight_layout()
fig.savefig('nr_time_performance.png')
plt.close(fig)
print("Plot saved to nr_time_performance.png")

fig, ax = plt.subplots()
ax.plot(sorted_tolerances, sorted_iterations, marker='s', color='tab:blue')
ax.set_xlabel('Tolerance ($\\epsilon$) (log scale)')
ax.set_ylabel('Iterations')
ax.set_xscale('log')
ax.set_ylim(0, 4)
ax.set_yticks([1, 2, 3, 4])
ax.set_title('Newton-Raphson Iterations vs Tolerance')
fig.tight_layout()
fig.savefig('nr_iterations_performance.png')
plt.close(fig)
print("Plot saved to nr_iterations_performance.png")