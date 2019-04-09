import numpy as np
import matplotlib.pyplot as plt

n = np.array([64, 256, 1024, 4096])**2

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
[0.00364528, 0.026544, 0.329556, 5.81617]
])



speedup_np = time_vs_np[0, :] / time_vs_np
speedup_p = time_vs_p[0, :] / time_vs_p

fig1, (ax11, ax12) = plt.subplots(ncols=2)
fig2, (ax21, ax22) = plt.subplots(ncols=2)

ax11.loglog(n, time_vs_np[0, :], label='np = 1')
ax11.loglog(n, time_vs_np[1, :], label='np = 2')
ax11.loglog(n, time_vs_np[2, :], label='np = 4')
ax11.loglog(n, time_vs_np[3, :], label='np = 6')
ax11.loglog(n, time_vs_np[4, :], label='np = 9')
ax11.loglog(n, time_vs_np[5, :], label='np = 12')
ax11.loglog(n, time_vs_np[6, :], label='np = 18')
ax11.loglog(n, time_vs_np[7, :], label='np = 36')

ax11.legend()
ax11.set_xlabel('n^2')
ax11.set_ylabel('Time [s]')
ax11.set_title('Time vs. number of processes')
ax11.grid()

ax12.plot(n, speedup_np[0, :], label='np = 1')
ax12.plot(n, speedup_np[1, :], label='np = 2')
ax12.plot(n, speedup_np[2, :], label='np = 4')
ax12.plot(n, speedup_np[3, :], label='np = 6')
ax12.plot(n, speedup_np[4, :], label='np = 9')
ax12.plot(n, speedup_np[5, :], label='np = 12')
ax12.plot(n, speedup_np[6, :], label='np = 18')
ax12.plot(n, speedup_np[7, :], label='np = 36')

ax12.legend()
ax12.set_xlabel('n^2')
ax12.set_ylabel('Speedup')
ax12.set_title('Speedup vs. number of processes')
ax12.grid()


ax21.plot(n, time_vs_p[0, :], label='p = 1')
ax21.plot(n, time_vs_p[1, :], label='p = 2')
ax21.plot(n, time_vs_p[2, :], label='p = 4')
ax21.plot(n, time_vs_p[3, :], label='p = 6')
ax21.plot(n, time_vs_p[4, :], label='p = 9')
ax21.plot(n, time_vs_p[5, :], label='p = 12')
ax21.plot(n, time_vs_p[6, :], label='p = 18')
ax21.plot(n, time_vs_p[7, :], label='p = 36')

ax21.legend()
ax21.set_xlabel('n^2')
ax21.set_ylabel('Time [s]')
ax21.set_title('Time vs. number of threads')


ax22.plot(n, speedup_p[0, :], label='p = 1')
ax22.plot(n, speedup_p[1, :], label='p = 2')
ax22.plot(n, speedup_p[2, :], label='p = 4')
ax22.plot(n, speedup_p[3, :], label='p = 6')
ax22.plot(n, speedup_p[4, :], label='p = 9')
ax22.plot(n, speedup_p[5, :], label='p = 12')
ax22.plot(n, speedup_p[6, :], label='p = 18')
ax22.plot(n, speedup_p[7, :], label='p = 36')

ax22.legend()
ax22.set_xlabel('n^2')
ax22.set_ylabel('Speedup')
ax22.set_title('Speedup vs. number of threads')

plt.show()


