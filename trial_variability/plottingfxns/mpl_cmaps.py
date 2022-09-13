import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat

a = np.flipud(plt.cm.RdBu_r(range(256))[:,:3]) # flip so blue is positive and red is negative
b = plt.cm.Blues(range(256))[:,:3]
r = plt.cm.Reds(range(256))[:,:3]
c = plt.cm.cividis(range(256))[:,:3]

#cmap_low = np.array([0,86,115])
#cmap_high = np.array([216,68,0])
cmap_def = np.array([[0,86,115, 255],[255,255,255,255],[216,68,0,255]]) # go from color 1 to white to col 2
from matplotlib.colors import ListedColormap
N = 256
vals1 = np.ones((N, 4)) # go from color 1 to white
vals1[:, 0] = np.linspace(cmap_def[0,0]/256, 1, N)
vals1[:, 1] = np.linspace(cmap_def[0,1]/256, 1, N)
vals1[:, 2] = np.linspace(cmap_def[0,2]/256, 1, N)
vals2 = np.ones((N, 4)) # go from white to color 2
vals2[:, 0] = np.linspace(1, cmap_def[2,0]/256, N)
vals2[:, 1] = np.linspace(1, cmap_def[2,1]/256, N)
vals2[:, 2] = np.linspace(1, cmap_def[2,2]/256, N)
cmap_array1 = np.concatenate((vals1,vals2),axis=0)

cmap = ListedColormap(cmap_array1)
np.save(file='data/colors/ejc_custom1.npy',arr=cmap_array1)

cmap_def = np.array([[98.,88.,159., 255.],[255.,255.,255.,255.],[147.,173.,144.,255.]]) # go from color 1 to white to col 2
cmap_def[2,:] = 1*cmap_def[2,:]
from matplotlib.colors import ListedColormap
N = 256
vals1 = np.ones((N, 4)) # go from color 1 to white
vals1[:, 0] = np.linspace(cmap_def[0,0]/256, 1, N)
vals1[:, 1] = np.linspace(cmap_def[0,1]/256, 1, N)
vals1[:, 2] = np.linspace(cmap_def[0,2]/256, 1, N)
vals2 = np.ones((N, 4)) # go from white to color 2
vals2[:, 0] = np.linspace(1, cmap_def[2,0]/256, N)
vals2[:, 1] = np.linspace(1, cmap_def[2,1]/256, N)
vals2[:, 2] = np.linspace(1, cmap_def[2,2]/256, N)
cmap_array2 = np.concatenate((vals1,vals2),axis=0)

cmap = ListedColormap(cmap_array2)
np.save(file='data/colors/ejc_custom2.npy',arr=cmap_array2)

# make reds colormap with grey for values == 0
r0 = r.copy()
r0[0,:] = [0.5,0.5,0.5] 
np.save(file='data/colors/Reds_gray0.npy',arr=r0)

# flip colormap
r0_flip = np.flipud(r0)
np.save(file='data/colors/Reds_gray0_flip.npy',arr=r0_flip)


savemat(file_name='data/colors/mpl_cmaps.mat',mdict={'RdBu_r':a,'Blues':b,'cividis':c,'Reds_gray0':r0,'Reds_gray0_flip':r0_flip,'ejc_custom1':cmap_array1[:,:3],'ejc_custom2':cmap_array2[:,:3]})
