import cmocean
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy import feature
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import georaster

extent = (-93.99, -80.3, 15.507,  31.955)
tran = ccrs.PlateCarree()
# proj = tran
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])

fig, axes = plt.subplots(2, 3, figsize=(45,37), subplot_kw={'projection':ccrs.AlbersEqualArea(central_lon, central_lat)})
row1 = axes[0]
row2 = axes[1]

def read_tif_to_array(tif_path):
    map = georaster.SingleBandRaster(tif_path)
    map.r[map.r < 0.0] = np.nan
    return map

sal = read_tif_to_array('D:/MSU/RA/acidification/acidification/salinity/salinity.tif')
ph = read_tif_to_array('D:/MSU/RA/acidification/acidification/ph/pH.tif')
ast = read_tif_to_array('D:/MSU/RA/acidification/acidification/ast/omega_Ar.tif')
css = read_tif_to_array('D:/MSU/RA/acidification/acidification/css/omega_Ca.tif')
spcp = read_tif_to_array('D:/MSU/RA/acidification/acidification/spcp/spCO2_micro.tif')
ta = read_tif_to_array('D:/MSU/RA/acidification/acidification/ta/TA.tif')

def add_plot(data, ax, title, colorbar, code, tick=None):
    tran = ccrs.PlateCarree()
    # plt.figure(figsize=(8,9))
    ax.set_extent(extent)
    norm = simple_norm(data.r, 'linear')
    colorbar = cmocean.tools.crop_by_percent(
        colorbar, 2, which='both', N=None
    )

    # ax.add_feature(feature.LAND, edgecolor='black', zorder=11)
    ax.add_feature(feature.LAND, zorder=11)
    mappable = ax.pcolormesh(data.r, cmap=colorbar, norm=norm)
    ax.imshow(data.r, extent=data.extent, transform=tran,
              interpolation='nearest', zorder=10, cmap=colorbar, norm=norm)
    cbbox = inset_axes(ax, '99%', '16%', loc=4)
    cbbox.tick_params(
        axis='both', left=False, top=False, right=False, bottom=False,
        labelleft=False, labeltop=False, labelright=False, labelbottom=False
    )


    cbbox.set_facecolor([1, 1, 1, 1])
    cbaxes = inset_axes(cbbox, '94%', '30%', loc=10)
    cb = fig.colorbar(mappable, cax=cbaxes, ax=ax,orientation='horizontal')
    cb.ax.set_title(title, fontsize=50, y=1.2)
    cb.ax.tick_params(labelsize=40, length=10, width=4, which='major')
    ax.text(0.07, 0.94, code, ha='center', va='center',fontsize=50, transform=ax.transAxes,
      zorder=14)
    if tick is not None:
        tick_locator = ticker.MaxNLocator(nbins=tick)
        cb.locator = tick_locator
        cb.update_ticks()

haline = cmocean.cm.haline_r
thermal = cmocean.cm.thermal
balance = cmocean.cm.balance
matter = cmocean.cm.matter
add_plot(sal, axes[0,0], 'Salinity', haline, '(a)')
add_plot(ph, axes[0,1], 'pH', thermal, '(b)')
add_plot(spcp, axes[0,2], r'Surface Pressure of $CO_{2}$ in µatm', balance, '(c)')
add_plot(ta, axes[1,0], 'Total Alkalinity in µmol/kg',matter, '(d)')
add_plot(css, axes[1,1], 'Calcite Saturation State',matter, '(e)')
add_plot(ast, axes[1,2], 'Aragonite Saturation State',matter, '(f)')
fig.tight_layout(pad=0.2)
# plt.show()
fig.savefig('acid1.png', dpi=150)


add_plot(sal, axes[0,0], 'Salinity', thermal, '(a)')
add_plot(ph, axes[0,1], 'pH', thermal, '(b)')
add_plot(spcp, axes[0,2], r'Surface Pressure of $CO_{2}$ in µatm', balance, '(c)')
add_plot(ta, axes[1,0], 'Total Alkalinity in µmol/kg',haline, '(d)')
add_plot(css, axes[1,1], 'Calcite Saturation State',haline, '(e)')
add_plot(ast, axes[1,2], 'Aragonite Saturation State',haline, '(f)')
# fig.tight_layout(pad=0.2)
# plt.show()
fig.savefig('acid2.png', dpi=150)




