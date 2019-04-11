import numpy as np
import matplotlib.pyplot as plt

n = np.array([64, 256, 1024, 4096])
p = np.array([1,2,4,6,9,12,18, 36])

time_vs_np = np.array([
[0.0118445, 0.218464, 4.04543, 74.0699],
[0.00705552, 0.273146, 2.36198, 42.8971],
[0.00373796, 0.0639896, 1.147880, 20.7168],
[0.00314088, 0.0449188, 0.785950, 14.3255],
[0.00225658, 0.0295443, 0.532675, 9.55536],
[0.00713626, 0.0226034, 0.399351, 7.19222],
[0.00598702, 0.0352874, 0.264528, 4.82325],
[0.00442977, 0.0104089, 0.141798, 2.31187]
])

time_vs_p = np.array([
[0.0116838, 0.218034, 4.051550, 78.9013],
[0.00657731, 0.120411, 2.23435, 42.0498],
[0.0038715, 0.0673242, 1.2312600, 22.3115],
[0.00298387, 0.0470334, 0.838947, 15.3274],
[0.00233850, 0.0335307, 0.585164, 10.4454],
[0.00222564, 0.0266532, 0.454868, 8.1042],
[0.00300185, 0.0341911, 0.484426, 8.83917],
[0.00364528, 0.026544, 0.329556, 5.81617] ]) 

speedup_np = time_vs_np[0, :] / time_vs_np
speedup_p = time_vs_p[0, :] / time_vs_p

efficiency_np = speedup_np / p[:, None]
efficiency_p = speedup_p / p[:, None]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, (ax3_1, ax3_2) = plt.subplots(ncols=2)
fig4, (ax4_1, ax4_2) = plt.subplots(ncols=2)

ax1.loglog(n, time_vs_np[0, :], label='#p = 1', linestyle='-', marker='x')
ax1.loglog(n, time_vs_np[1, :], label='#p = 2', linestyle='-', marker='x') 
ax1.loglog(n, time_vs_np[2, :], label='#p = 4', linestyle='-', marker='x')
ax1.loglog(n, time_vs_np[3, :], label='#p = 6', linestyle='-', marker='x')
ax1.loglog(n, time_vs_np[4, :], label='#p = 9', linestyle='-', marker='x')
ax1.loglog(n, time_vs_np[5, :], label='#p = 12', linestyle='-', marker='x')
ax1.loglog(n, time_vs_np[6, :], label='#p = 18', linestyle='-', marker='x')
ax1.loglog(n, time_vs_np[7, :], label='#p = 36', linestyle='-', marker='x')

ax1.legend()
ax1.set_title('Time vs. #processes for #threads=1')
ax1.set_xlabel('n')
ax1.set_ylabel('Time [s]')
ax1.grid()

ax2.loglog(n, time_vs_p[0, :], label='#t = 1', linestyle='-', marker='x')
ax2.loglog(n, time_vs_p[1, :], label='#t = 2', linestyle='-', marker='x') 
ax2.loglog(n, time_vs_p[2, :], label='#t = 4', linestyle='-', marker='x')
ax2.loglog(n, time_vs_p[3, :], label='#t = 6', linestyle='-', marker='x')
ax2.loglog(n, time_vs_p[4, :], label='#t = 9', linestyle='-', marker='x')
ax2.loglog(n, time_vs_p[5, :], label='#t = 12', linestyle='-', marker='x') 
ax2.loglog(n, time_vs_p[6, :], label='#t = 18', linestyle='-', marker='x') 
ax2.loglog(n, time_vs_p[7, :], label='#t = 36', linestyle='-', marker='x')  


ax2.legend()
ax2.set_xlabel('n')
ax2.set_ylabel('Time [s]')
ax2.set_title('Time vs. #threads for #processes=1')
ax2.grid() 

n = n**2


ax3_1.plot(n, speedup_np[0, :], label='#p = 1', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[1, :], label='#p = 2', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[2, :], label='#p = 4', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[3, :], label='#p = 6', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[4, :], label='#p = 9', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[5, :], label='#p = 12', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[6, :], label='#p = 18', linestyle='-', marker='x')
ax3_1.plot(n, speedup_np[7, :], label='#p = 36', linestyle='-', marker='x')

ax3_1.legend()
ax3_1.grid()
ax3_1.set_xlabel(r'$n^2$')
ax3_1.set_ylabel('Speedup')
ax3_1.set_title(r'Speedup vs. $n^2$, #t=1')


ax3_2.plot(n, efficiency_np[0, :], label='#p = 1', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[1, :], label='#p = 2', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[2, :], label='#p = 4', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[3, :], label='#p = 6', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[4, :], label='#p = 9', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[5, :], label='#p = 12', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[6, :], label='#p = 18', linestyle='-', marker='x')
ax3_2.plot(n, efficiency_np[7, :], label='#p = 36', linestyle='-', marker='x')

ax3_2.legend()
ax3_2.grid()
ax3_2.set_xlabel(r'$n^2$')
ax3_2.set_ylabel('Efficiency')
ax3_2.set_title(r'Efficiency vs. $n^2$, #t=1')

ax4_1.plot(n, speedup_p[0, :], label='#t = 1', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[1, :], label='#t = 2', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[2, :], label='#t = 4', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[3, :], label='#t = 6', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[4, :], label='#t = 9', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[5, :], label='#t = 12', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[6, :], label='#t = 18', linestyle='-', marker='x')
ax4_1.plot(n, speedup_p[7, :], label='#t = 36', linestyle='-', marker='x')

ax4_1.legend()
ax4_1.grid()
ax4_1.set_xlabel(r'$n^2$')
ax4_1.set_ylabel('Speedup')
ax4_1.set_title(r'Speedup vs. $n^2$, #p=1')


ax4_2.plot(n, efficiency_p[0, :], label='#t = 1', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[1, :], label='#t = 2', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[2, :], label='#t = 4', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[3, :], label='#t = 6', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[4, :], label='#t = 9', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[5, :], label='#t = 12', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[6, :], label='#t = 18', linestyle='-', marker='x')
ax4_2.plot(n, efficiency_p[7, :], label='#t = 36', linestyle='-', marker='x')

ax4_2.legend()
ax4_2.grid()
ax4_2.set_xlabel(r'$n^2$')
ax4_2.set_ylabel('Efficiency')
ax4_2.set_title(r'Efficiency vs. $n^2$, #p=1')


plt.show()


