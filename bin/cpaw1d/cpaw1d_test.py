import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.gridspec as gridspec
import types
from pegasus_read import vtk_read as vtk


num=3000
filename = "id0/cpaw1d."+"%04d"%num+".all.vtk"
data = vtk(filename)
X = data['x']
X = X[:,0,0]
X_cc = (X[1:]+X[:-1])/2.0
B = data['cell_centered_B']
By = B[:,0,0,1]

width = 512.11743/72.2
font = 9
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rc('font', family = 'serif')
mpl.rcParams['xtick.labelsize']=font-1
mpl.rcParams['ytick.labelsize']=font-1
fig=plt.figure(figsize=(1,1))
fig.set_figheight((np.sqrt(5.0)-1.0)/2.0 * width/2)
fig.set_figwidth(width/2)
gs1 = gridspec.GridSpec(1,1,height_ratios = [1])
gs1.update(wspace=0.0, hspace=0.)

ax1 = plt.subplot(gs1[0])
ax1.plot(X_cc,By,'-',color='b',linewidth=0.5,label=r"$B_y$")
ax1.scatter(X_cc,By,color='b',s=0.5,linewidth=0.5)
ax1.legend(loc = "upper right",fontsize=font,frameon=False)
plt.xlabel(r"$x/d_{\rm i}$",fontsize=font)
plt.ylabel(r"$B/B_0$",fontsize=font)

#plt.yscale('log')
#plt.xscale('log')

plt.savefig("plots_cpaw1d/plot."+"%04d"%num+".pdf",bbox_inches='tight')
