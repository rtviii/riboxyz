import os
import pickle
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from alphashape import alphashape


RCSB_ID      = '6Z6K'
cluster_path = "/home/rtviii/dev/riboxyz/mesh_generation/{}_cluster.npy".format(RCSB_ID)
cluster      = np.load(cluster_path)

import plotly.graph_objects as go


# pts = np.loadtxt(np.DataSource().open('https://raw.githubusercontent.com/plotly/datasets/master/mesh_dataset.txt'))
x, y, z = cluster.T
_ = go.Mesh3d(x=x, y=y, z=z,
                   alphahull=3,
                   opacity=0.4,
                   color='cyan')


print(len(_.x))
print(len(_.y))
print(len(_.z))
fig = go.Figure(data=[_])
fig.show()