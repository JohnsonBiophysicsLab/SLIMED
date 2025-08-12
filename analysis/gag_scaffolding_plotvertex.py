#!/usr/bin/env python
# coding: utf-8

# In[1]:


#@python 3.8.10
#@author moon ying
#plot continuum membrane results with matplotlib 3d plot and trimesh package
#dependencies are listed as follows

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import trimesh as tr
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import seaborn as sns
sns.set_context("talk")

import pyvista as pv
import vtk
from pyvistaqt import BackgroundPlotter

# use one of...
from pykdtree.kdtree import KDTree
# from scipy.spatial import cKDTree as KDTree


# In[2]:


infile = "vertexfinal.csv"
trangulationfile = "face.csv"
outfile = "392"


# In[3]:


vertices = pd.read_csv(infile, header = None)
vertices.columns = ["x","y","z"]
vertices


# In[4]:


df_faces = pd.read_csv(trangulationfile, header = None)
faces = df_faces.to_numpy()
#faces.columns = ["x","y","z"]
faces


# In[5]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')

#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],
                       triangles=faces, cmap=cm.coolwarm,
                       edgecolor = "k")
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')



ax.view_init(elev = 10, azim = -45)
plt.savefig(outfile + "gagsonly_10deg_with_axis.png", dpi = 300)


# In[6]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"], triangles = faces,                                cmap = "Blues", alpha = 0.3, edgecolor="k", label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors


ax.set_xlim3d(-40.0, 40.0)
ax.set_ylim3d(-40.0, 40.0)
ax.set_zlim3d(-40.0, 40.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 20, azim = 45)

#plt.savefig(outfile + "membrane_10deg_uniform_scale.png", dpi = 300)


# In[7]:


# VTK requires mesh faces info stores as [num_vertices, vertex1, vertex2, ...]
df_faces.insert(0, 'num', 3)
vtk_faces = df_faces.to_numpy()
vtk_faces


# In[8]:


subdiv_times = 4
subfilter = "loop"

arr_x = vertices["x"].to_numpy()
arr_y = vertices["y"].to_numpy()
arr_z = vertices["z"].to_numpy()

npdat = np.transpose(np.array([arr_x, arr_y, arr_z]))

# Convert the numpy array into a pyvista mesh object
mesh = pv.PolyData(npdat, vtk_faces)

# Subdivide with pyvista module
# currently `subfilter` field supports three schemes: linear, loop, butterfly
meshsub = mesh.subdivide(subdiv_times, subfilter=subfilter)
meshsub


# In[9]:


vertices_subdiv = meshsub.points.T
vertices_subdiv


# In[10]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})



surf_subdiv = ax.plot_trisurf(vertices_subdiv[0],                                 vertices_subdiv[1],
                                vertices_subdiv[2], \
                               cmap = "Blues", alpha = 1.0, edgecolor="none", label = "Limit Surface")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_subdiv._facecolors2d=surf_subdiv._facecolors
surf_subdiv._edgecolors2d=surf_subdiv._edgecolors

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"], triangles = faces,                                cmap = "Blues", alpha = 0.24, linewidth = 0.5, edgecolor="k", label = "Control Mesh")
surf_vertices._facecolors2d=surf_vertices._edgecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors


ax.set_xlim3d(-50.0, 50.0)
ax.set_ylim3d(-50.0, 50.0)
ax.set_zlim3d(-50.0, 50.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 20, azim = 45)

plt.savefig(outfile + "ls_cm_20deg_uniform_scale.png", dpi = 900)


# In[11]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})



surf_subdiv = ax.plot_trisurf(vertices_subdiv[0],                                 vertices_subdiv[1],
                                vertices_subdiv[2], \
                               cmap = "Blues", alpha = 1.0, edgecolor="none", label = "Limit Surface")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_subdiv._facecolors2d=surf_subdiv._facecolors
surf_subdiv._edgecolors2d=surf_subdiv._edgecolors

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"], triangles = faces,                                cmap = "Blues", alpha = 0.24, linewidth = 0.5, edgecolor="k", label = "Control Mesh")
surf_vertices._facecolors2d=surf_vertices._edgecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors


ax.set_xlim3d(-50.0, 50.0)
ax.set_ylim3d(-50.0, 50.0)
ax.set_zlim3d(-50.0, 50.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 90, azim = 0)

plt.savefig(outfile + "ls_cm_90deg_uniform_scale.png", dpi = 900)


# In[12]:


# Plot the PolyData with PyVista
plotter = pv.plot(meshsub)

# Get the Matplotlib axes from the PyVista plotter
mpl_axes = plotter.show_axes()

# Customize the PyVista plot if needed
# plotter.add_mesh(...) or other PyVista plotting functions

# Create a Matplotlib plot using the same axes
fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': '3d'}, num='My Overlayed Plot')
ax.axis('off')  # Turn off axes if needed

# Your Matplotlib plotting code goes here
# For example:
# ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 0, 0], color='red', linewidth=2)
# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"], triangles = faces,                                cmap = "Blues", alpha = 0.3, edgecolor="k", label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors


# Show the Matplotlib plot
plt.show()


# Energy and Force
# ==

# In[ ]:


df_ef = pd.read_csv("EnergyForce.csv", index_col = False, header = None, skiprows = 1)
df_ef.columns = ["E_curv", "E_area", "E_reg", "E_tot", "F_mean"]
df_ef


# In[ ]:


plt.plot(df_ef["E_tot"])


# In[ ]:


plt.yscale("log")
plt.plot(df_ef["E_reg"])


# In[ ]:


plt.yscale("log")
plt.plot(df_ef["E_tot"])


# In[ ]:


gags_r


# In[ ]:





# In[ ]:




