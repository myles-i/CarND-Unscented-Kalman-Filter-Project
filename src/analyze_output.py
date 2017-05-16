import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import pdb as pdb

my_data = pd.read_csv('../build/output.txt', sep = '\t')
my_data.NIS.plot()
plt.show()



input_data = pd.read_csv("../data/obj_pose-laser-radar-synthetic-input.txt","\t")

	