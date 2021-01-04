from matplotlib import pyplot as plt
import pandas
pandas.set_option("mode.chained_assignment", None)
from pandas import read_csv
import numpy as np
import os.path as osp

data = read_csv("datalog_demo4.csv", low_memory=False)
data = data.rename(columns=lambda x: x.strip())
data.loc[data.index.size-1].iloc[1:] = [0]*(data.columns.size - 1)

# print data

plt.figure(figsize=[3,2.5])
# plt.figure()
plt.plot(data["timestamp"]*1e-6, data["rf_0"], 'g')
plt.plot(data["timestamp"]*1e-6, data["rf_1"], 'r')
plt.plot(data["timestamp"]*1e-6, data["rf_2"], 'b')
plt.gca().set_xlim([2.0,12.0])
plt.gca().set_ylim([-5,12.0])
plt.savefig(osp.join(osp.expanduser('~'), 'Downloads/roller_to_arm_force.png'))

plt.figure(figsize=[3,2.5])
# plt.figure()
plt.plot(data["timestamp"]*1e-6, data["gf_0"], 'r')
plt.plot(data["timestamp"]*1e-6, -data["gf_1"], 'g')
plt.plot(data["timestamp"]*1e-6, data["gf_2"], 'b')
plt.gca().set_xlim([2.0,12.0])
plt.gca().set_ylim([-5,12.0])
plt.savefig(osp.join(osp.expanduser('~'), 'Downloads/ground_to_roller_force.png'))

# plt.figure()
# plt.plot(data["timestamp"]*1e-6, data["dist_0"], 'g')
# plt.plot(data["timestamp"]*1e-6, data["dist_1"], 'r')


# plt.gca().set_xlim([0.9,2.1])
# plt.gca().set_ylim([-0.6,0.6])
# plt.savefig(osp.join(osp.expanduser('~'), 'Downloads/mesh_omega.png'))

# FOR PLOTTING DEMO 1 ANGULAR VELOCITY
# plt.figure(figsize=[3,2.5])
# # plt.figure()
# plt.plot(data["timestamp"]*1e-6/9, data["omega_0"], 'r')
# # plt.plot(data["timestamp"]*1e-6/9, data["omega_1"], 'g')
# plt.plot(data["timestamp"]*1e-6/9, data["omega_2"], 'b')
# plt.gca().set_xlim([0.9,2.1])
# plt.gca().set_ylim([-0.6,0.6])
# plt.savefig(osp.join(osp.expanduser('~'), 'Downloads/mesh_omega.png'))

# fig, axes = plt.subplots(5,1, sharex=True)
# axes[0].set_ylabel("Force (N)")
# plot torque
# axes[1].plot(data["timestamp"]*1e-6, data["torque_0"]/30)
# axes[1].set_ylabel("Torque (Nm)")
# plot speed
# print (data["timestamp"].iloc[5]-data["timestamp"].iloc[4])
# speed = np.diff(data["speed_0"]*1e6, prepend=0)/(data["timestamp"].iloc[5]-data["timestamp"].iloc[4])
# axes[2].plot(data["timestamp"]*1e-6, speed)
# axes[2].set_ylabel("Shaft speed (RPM)")
# plot pressure
# axes[3].plot(data["timestamp"]*1e-6, data["pressure_0"]*0.0063125)
# axes[3].set_ylabel("Pressure (bar)")
# plot temperature
# axes[4].plot(data["timestamp"]*1e-6, data["front_temp_0"]/10)
# axes[4].set_ylabel("Front temp (deg C)")
# axes[4].set_xlabel("Time (sec)")
plt.show()