import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval as make_tuple
from mpl_toolkits import mplot3d

def make_3D_plot(filename : str):
  f_desc = open(filename, 'r')
  x_data = []
  y_data = []
  z_data = []
  for line in f_desc.readlines():
    (x,y,z) = line.replace('(','').replace(')','').split(',')
    #if (float(z) != 0):
    #  print(f'({x},{y},{z})')
    #  continue
    x_data.append(float(x))
    y_data.append(float(y))
    z_data.append(float(z))
  f_desc.close()

  fig = plt.figure()
  ax = plt.axes(projection="3d")
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  ax.scatter3D(x_data, y_data, z_data)

  plt.show()

make_3D_plot('./output')