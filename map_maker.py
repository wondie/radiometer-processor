import cmocean

import matplotlib.pyplot as plt
import numpy as np
print (cmocean.cm.cmapnames)

cmap = cmocean.cm.topo
fig, axes = plt.subplots(1, 2, figsize=(8,4))
A = np.random.randint(-50, 201, (10,10))
salinity = 'D:/MSU/RA/acidification/acidification/salinity/salinity.tif'
mappable = axes[0].pcolormesh(A, vmin=-200, vmax=200, cmap=cmap)
axes[0].set_title('No values<-50, but still\nshow possibility in colorbar')
fig.colorbar(mappable, ax=axes[0])

newcmap = cmocean.tools.crop(cmap, -50, 200, 0)
mappable = axes[1].pcolormesh(A, vmin=-50, vmax=200, cmap=newcmap)
axes[1].set_title('Colorbar only shows color\nrange used by data')
fig.colorbar(mappable, ax=axes[1])
fig.show()
#
# ax0.coastlines()
# ax0.add_feature(cfeature.BORDERS)
# ax0.add_feature(cfeature.STATES)
# ax0.set_title('Mean 500mb Heights (dam)', fontsize=12)
# fig = plt.subplot(ax0)

fig.savefig('500mbHGT.png')
