import cmocean
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy import feature
from cartopy.feature import NaturalEarthFeature, AdaptiveScaler
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import georaster
from cartopy.io.shapereader import Reader

plt.rcParams["figure.autolayout"] = True
extent = (-93.99, -80.3, 15.507, 31.955)
tran = ccrs.PlateCarree()
# proj = tran
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])

fig, axes = plt.subplots(2, 3, figsize=(45, 37), subplot_kw={
    'projection': ccrs.AlbersEqualArea(central_lon, central_lat)})
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


def add_plot(data, ax, title, colorbar, code, min_cut=None,max_cut=None,min_p=None, max_p=None, tick=None):
    tran = ccrs.PlateCarree()
    # plt.figure(figsize=(8,9))
    ax.set_extent(extent)
    colorbar = cmocean.tools.crop_by_percent(
        colorbar, 2, which='both', N=None
    )
    # ax.add_feature(feature.LAND, edgecolor='black', zorder=11)
    auto_scaler = AdaptiveScaler('110m', (('50m', 50), ('10m', 15)))
    land = NaturalEarthFeature('physical', 'land', auto_scaler,
                               edgecolor='face', facecolor=np.array((220, 220, 220)) / 256., zorder=-1)
    ax.add_feature(land, zorder=11)

    bathym = feature.NaturalEarthFeature(
        name='bathymetry_K_200', scale='10m', category='physical')
    ax.add_feature(
        bathym, facecolor='none', edgecolor='black',
        linestyle='solid', linewidth=1, zorder=11)

    # shape_feature = feature.ShapelyFeature(Reader(
    #     'D:/MSU/RA/acidification/bathymetry/GoM_200m.shp').geometries(),
    #                                ccrs.PlateCarree(), facecolor='none')
    # ax.add_feature(shape_feature, facecolor='none',
    #                edgecolor='black', linestyle='solid', linewidth=1, zorder=11)

    cbbox = inset_axes(ax, '99%', '16%', loc=4)
    cbbox.tick_params(
        axis='both', left=False, top=False, right=False, bottom=False,
        labelleft=False, labeltop=False, labelright=False, labelbottom=False
    )

    cbbox.set_facecolor([1, 1, 1, 1])
    cbaxes = inset_axes(cbbox, '94%', '30%', loc=10)
    # if max_cut is not None:
    #     norm = simple_norm(data.r, 'linear', max_cut=max_cut)
    #     print('1norm')
    # elif max_p is not None:
    #     print ('2norm')
    #     norm = simple_norm(data.r, 'linear', min_percent=min_p, max_percent=max_p)
    # else:
    #     print ('3norm')
    #     norm = simple_norm(data.r, 'linear')
    norm = simple_norm(data.r, 'linear', min_percent=min_p, max_percent=max_p, min_cut=min_cut, max_cut=max_cut)
    mappable = ax.pcolormesh(data.r, cmap=colorbar, norm=norm)
    ax.imshow(data.r, extent=data.extent, transform=tran,
              interpolation='nearest', zorder=10, cmap=colorbar, norm=norm)

    cb = fig.colorbar(mappable, cax=cbaxes, ax=ax, orientation='horizontal')
    cb.ax.set_title(title, fontsize=40, y=1.2)
    cb.ax.tick_params(labelsize=30, length=10, width=4, which='major')
    ax.text(0.07, 0.94, code, ha='center', va='center',
            fontsize=50, transform=ax.transAxes, zorder=14)
    if tick is not None:
        tick_locator = ticker.MaxNLocator(nbins=tick)
        cb.locator = tick_locator
        cb.update_ticks()


haline = cmocean.cm.haline_r
thermal = cmocean.cm.thermal_r
balance = cmocean.cm.balance
matter = cmocean.cm.matter
dense = cmocean.cm.dense
speed = cmocean.cm.speed



# add_plot(sal, axes[0, 0], 'Salinity', haline, '(a)', max_cut=40)
# add_plot(ta, axes[0, 1], 'Total Alkalinity in µmol/kg', haline, '(b)')
# add_plot(spcp, axes[0, 2], r'Surface Pressure of $CO_{2}$ in µatm', haline, '(c)')
# add_plot(ph, axes[1, 0], 'pH', haline, '(d)')
# add_plot(css, axes[1, 1], 'Calcite Saturation State', haline, '(e)')
# add_plot(ast, axes[1, 2], 'Aragonit Saturation State', haline, '(f)')
# fig.tight_layout(pad=0.2)
# plt.show()
# fig.savefig('acid1.png', dpi=150)

add_plot(sal, axes[0, 0], 'Salinity', haline, '(a)', min_cut=30, max_cut=40)
add_plot(ta, axes[0, 1], 'Total Alkalinity in µmol/kg', dense, '(b)', max_cut=2500, tick=5)
add_plot(spcp, axes[0, 2], r'Partial Pressure of $CO_{2}$ in µatm', thermal, '(c)')
add_plot(ph, axes[1, 0], 'pH', speed, '(d)')
add_plot(css, axes[1, 1], 'Calcite Saturation State', matter, '(e)')
add_plot(ast, axes[1, 2], 'Aragonite Saturation State', matter, '(f)')
fig.tight_layout(pad=0.5)
# fig.colorbar(im, ax=axs, shrink=0.6)
# plt.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
# fig.set_tight_layout(False)
# plt.show()
fig.savefig('acid2022_c.png', dpi=150, bbox_inches='tight')
#
# add_plot(sal, axes[0, 0], 'Salinity', haline, '(a)')
# add_plot(ta, axes[0, 1], 'Total Alkalinity in µmol/kg', matter, '(b)')
# add_plot(spcp, axes[0, 2], r'Surface Pressure of $CO_{2}$ in µatm', matter, '(c)')
# add_plot(ph, axes[1, 0], 'pH', thermal, '(d)')
# add_plot(css, axes[1, 1], 'Calcite Saturation State', thermal, '(e)')
# add_plot(ast, axes[1, 2], 'Aragonite Saturation State', thermal, '(f)')
#
# # plt.show()
# fig.savefig('acid3.png', dpi=150)
#
# add_plot(sal, axes[0, 0], 'Salinity', haline, '(a)')
# add_plot(ta, axes[0, 1], 'Total Alkalinity in µmol/kg', thermal, '(b)')
# add_plot(spcp, axes[0, 2], r'Surface Pressure of $CO_{2}$ in µatm', thermal, '(c)')
# add_plot(ph, axes[1, 0], 'pH', matter, '(d)')
# add_plot(css, axes[1, 1], 'Calcite Saturation State', matter, '(e)')
# add_plot(ast, axes[1, 2], 'Aragonite Saturation State', matter, '(f)')
#
# # plt.show()
# fig.savefig('acid4.png', dpi=150)
#
# add_plot(sal, axes[0, 0], 'Salinity', thermal, '(a)')
# add_plot(ta, axes[0, 1], 'Total Alkalinity in µmol/kg', matter, '(c)')
# add_plot(spcp, axes[0, 2], r'Surface Pressure of $CO_{2}$ in µatm', balance, '(c)')
# add_plot(ph, axes[1, 0], 'pH', haline, '(d)')
# add_plot(css, axes[1, 1], 'Calcite Saturation State', haline, '(e)')
# add_plot(ast, axes[1, 2], 'Aragonite Saturation State', haline, '(f)')
#
# # plt.show()
# fig.savefig('acid5.png', dpi=150, bbox_inches='tight')
#
# add_plot(sal, axes[0, 0], 'Salinity', haline, '(a)', min_p=60)
# add_plot(ta, axes[0, 1], 'Total Alkalinity in µmol/kg', haline, '(b)', max_cut=2500, tick=5)
# add_plot(spcp, axes[0, 2], r'Surface Pressure of $CO_{2}$ in µatm', haline, '(c)')
# add_plot(ph, axes[1, 0], 'pH', haline, '(d)')
# add_plot(css, axes[1, 1], 'Calcite Saturation State', haline, '(e)')
# add_plot(ast, axes[1, 2], 'Aragonite Saturation State', haline, '(f)')
# # plt.show()
# fig.savefig('acid6.png', dpi=150)
