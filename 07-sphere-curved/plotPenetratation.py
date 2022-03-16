# plotPenetratation.py

import numpy as np
import pandas
pandas.set_option("mode.chained_assignment", None)
from pandas import read_csv
from matplotlib import pyplot as plt

# data = read_csv("log_curvature.csv", low_memory=False)
# data = data.rename(columns=lambda x: x.strip()) # this removes white spaces padding the column names

# pene0 = data[data["pene_0"] > 0]["pene_0"].iloc[0]
# print pene0

# plt.figure()
# plt.plot(data['timestamp']*1e-6, 1e3*(data["pene_0"] - pene0))
# plt.show()

data_nocomp = read_csv("log_curvature_nocompensation.csv", low_memory=False)
data_nocomp = data_nocomp.rename(columns=lambda x: x.strip()) # this removes white spaces padding the column names

pene0_nocomp = data_nocomp[data_nocomp["pene_0"] > 0]["pene_0"].iloc[0]

data_fullcomp = read_csv("log_curvature_fullcompensation.csv", low_memory=False)
data_fullcomp = data_fullcomp.rename(columns=lambda x: x.strip()) # this removes white spaces padding the column names

pene0_fullcomp = data_fullcomp[data_fullcomp["pene_0"] > 0]["pene_0"].iloc[0]

plt.figure()
plt.plot(data_nocomp['timestamp']*1e-6, 1e3*(data_nocomp["pene_0"] - pene0_nocomp), 'r')
plt.plot(data_fullcomp['timestamp']*1e-6, 1e3*(data_fullcomp["pene_0"] - pene0_fullcomp), 'g')
plt.show()