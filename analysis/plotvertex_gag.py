#!/usr/bin/env python
# coding: utf-8

# In[12]:


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


# In[13]:


infile = "vertexfinal.csv"
outfile = "392"


# In[14]:


vertices = pd.read_csv(infile, header = None)
vertices.columns = ["x","y","z"]
vertices


# In[15]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')

#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"], cmap=cm.coolwarm)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.view_init(elev = 10, azim = -45)
plt.savefig(outfile + "gagsonly_10deg_with_axis.png", dpi = 300)


# In[16]:


gags = pd.read_csv("../COM.csv", header = None, skiprows=1)
gags.columns = ["x","y","z"]
gags


# In[17]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],                                 cmap="Blues", alpha = 0.7, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')

ax.set_xlim3d(-100.0, 100.0)
ax.set_ylim3d(-100.0, 100.0)
ax.set_zlim3d(-100.0, 100.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 10, azim = 45)

plt.savefig(outfile + "gags_membrane_10deg_uniform_scale.png", dpi = 300)


# In[18]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],                                 cmap="Blues", alpha = 0.95, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_xlim3d(-100.0, 100.0)
ax.set_ylim3d(-100.0, 100.0)
ax.set_zlim3d(-100.0, 100.0)
#ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 10, azim = 45)

plt.savefig(outfile + "gags_membrane_10deg_uniform_scale_with_axis.png", dpi = 300)


# In[19]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],                                 cmap="Blues", alpha = 0.7, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')

ax.set_xlim3d(-100.0, 100.0)
ax.set_ylim3d(-100.0, 100.0)
ax.set_zlim3d(-100.0, 100.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 0, azim = 45)

plt.savefig(outfile + "gags_membrane_0deg_uniform_zscale.png", dpi = 300)


# In[20]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],                                 cmap="Blues", alpha = 0.7, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')

ax.set_xlim3d(-100.0, 100.0)
ax.set_ylim3d(-100.0, 100.0)
ax.set_zlim3d(-20.0, 20.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 0, azim = 45)

plt.savefig(outfile + "gags_membrane_0deg_5x_zscale.png", dpi = 300)


# In[21]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],                                 cmap="Blues", alpha = 0.7, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')

ax.set_xlim3d(-100.0, 100.0)
ax.set_ylim3d(-100.0, 100.0)
ax.set_zlim3d(-20.0, 20.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 10, azim = 45)

plt.savefig(outfile + "gags_membrane_10deg_5x_zscale.png", dpi = 300)


# In[22]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"],                                 cmap="Blues", alpha = 0.7, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
surf_vertices._facecolors2d=surf_vertices._facecolors
surf_vertices._edgecolors2d=surf_vertices._edgecolors

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')

ax.set_xlim3d(-100.0, 100.0)
ax.set_ylim3d(-100.0, 100.0)
ax.set_zlim3d(-100.0, 100.0)
ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 90, azim = 45)

plt.savefig(outfile + "gags_membrane_topdown_uniform_zscale.png", dpi = 300)


# In[ ]:





# In[23]:


faces = pd.read_csv("face.csv", header = None)-1
faces.columns = ["x","y","z"]
faces


# In[24]:


mesh = tr.Trimesh(vertices.to_numpy(), faces.to_numpy())
mesh.show()


# In[25]:


mesh = tr.Trimesh(-vertices.to_numpy(), faces.to_numpy())
mesh.show()


# In[26]:


#Move membrane above the gag complex (-17 to +10)
#Add scaling constant k (0.01, +0.005 per iter)


# Energy and Force
# ==

# In[27]:


df_ef = pd.read_csv("EnergyForce.csv", index_col = False, header = None, skiprows = 1)
df_ef.columns = ["E_curv", "E_area", "E_reg", "E_tot", "F_mean"]
df_ef


# In[ ]:


plt.yscale("log")
plt.plot(df_ef["E_reg"])


# In[ ]:


plt.yscale("log")
plt.plot(df_ef["E_tot"])


# In[ ]:


fig = plt.figure(figsize=(24, 12))
ax = plt.axes(projection='3d')


#ax.scatter3D(vertices["x"], vertices["y"], vertices["z"], c = vertices["z"], s = 16.0);

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
#surf_gags = ax.plot_trisurf(gags["x"], gags["y"], gags["z"], cmap="Oranges")
#surf_vertices = ax.plot_trisurf(vertices["x"], vertices["y"], vertices["z"], \
#                                cmap="Blues", alpha = 0.95, label = "membrane")
# https://stackoverflow.com/questions/54994600/pyplot-legend-poly3dcollection-object-has-no-attribute-edgecolors2d
#surf_vertices._facecolors2d=surf_vertices._facecolors
#surf_vertices._edgecolors2d=surf_vertices._edgecolors

gags = pd.read_csv("COM.csv")
gags.columns = ["x","y","z"]

scatter_gags = ax.scatter(gags["x"], gags["y"], gags["z"], s = 3, c = "#000000", label = "gag lattice")

gags_r = pd.read_csv("Spline0.csv", header = None)
gags_r.columns = ["x","y","z"]

scatter_gags_r = ax.scatter(gags_r["x"], gags_r["y"], gags_r["z"], s = 3, c = "#FF0000", label = "gag lattice 0")

gags_s = pd.read_csv("Spline200.csv", header = None)
gags_s.columns = ["x","y","z"]

scatter_gags_s = ax.scatter(gags_s["x"], gags_s["y"], gags_s["z"], s = 3, c = "#00FF00", label = "gag lattice 150")



ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#ax.set_xlim3d(-100.0, 100.0)
#ax.set_ylim3d(-100.0, 100.0)
#ax.set_zlim3d(-100.0, 100.0)
#ax.set_axis_off()
ax.legend(prop={'size': 12})
ax.view_init(elev = 10, azim = 45)

plt.savefig(outfile + "gags_membrane_10deg_uniform_scale_with_axis.png", dpi = 300)


# In[ ]:


gags_r


# In[ ]:





# In[ ]:




