import warnings
from collections import OrderedDict
from datetime import datetime

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
import pandas as pd
from cartopy.io.shapereader import Reader
# plt.rcParams['text.usetex'] = True

plt.rcParams["figure.autolayout"] = True
extent = (-93.99, -80.3, 15.507, 31.955)


def read_tif_to_array(tif_path):
    map = georaster.SingleBandRaster(tif_path)
    map.r[map.r < 0.0] = np.nan
    return map



def add_plot(data, ax, title, colorbar, code, min_cut=None,max_cut=None,min_p=None, max_p=None, tick=None):
    tran = ccrs.PlateCarree()
    # proj = tran
    central_lon = np.mean(extent[:2])
    central_lat = np.mean(extent[2:])

    fig, axes = plt.subplots(2, 3, figsize=(45, 37), subplot_kw={
        'projection': ccrs.AlbersEqualArea(central_lon, central_lat)})
    # row1 = axes[0]
    # row2 = axes[1]
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
    norm = simple_norm(data.r, 'linear', min_percent=min_p, max_percent=max_p,
                       min_cut=min_cut, max_cut=max_cut)
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

def save_plot(plt, size, file_name=None):
    style = plot_size(size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # print("Failed")
        if file_name is None:
            file_name = 'boxplot'
        plt.savefig(r'D:\MSU\codes\radiometer_processor\data\{}.png'.format(file_name), dpi = style['dpi'])
        plt.show()
def create_scatter_plot(actual, estimated, validation, model, training_y, training_x, r2, equation_latex, size='paper'):

    style = plot_size(size)

    fig, axes = plt.subplot_mosaic(
        [['left', 'right']], constrained_layout=False, figsize=style['fig_size'],
        # gridspec_kw=dict(top = 0.95, bottom = 0, right = 1.04, left = 0.0,
        # hspace = 1, wspace = 0)
    )

    val_ax = axes['right']
    add_validation_plot(actual, estimated, val_ax, validation, style)
    alg_ax = axes['left']
    add_algorithm_plot(training_y, training_x,model,  alg_ax, r2, equation_latex, style)
    # fig.tight_layout(h_pad=0.8, w_pad=0.32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # print("Failed")
        plt.savefig(r'D:\MSU\codes\radiometer_processor\data\result_{}.png'.format(r2),
                dpi = style['dpi'])
        fig.show()


def plot_size(size):
    style = {'fig_size': (8, 4), 'point_size': 10, 'plot_font_size': 10,
             'label_font_size': 8, 'dpi': 200}
    if size == 'poster':
        style = {'fig_size': (20, 10), 'point_size': 50, 'plot_font_size': 24,
                 'label_font_size': 30, 'dpi': 500}
    if size == 'ppt':
        style = {'fig_size': (10, 5), 'point_size': 30, 'plot_font_size': 10,
                 'label_font_size': 15, 'dpi': 300}
    if size == 'paper':
        style = {'fig_size': (8, 4), 'point_size': 10, 'plot_font_size': 8,
                 'label_font_size': 10, 'dpi': 200}
    return style


def add_validation_plot(actual, estimated, axis, validation, style):
    colors = ['m']
    axis.scatter(actual, estimated, c=colors, s=style['point_size'])
    p1 = max(max(estimated), max(actual)) + 10
    p2 = 0
    axis.axline((p2, p2), (p1, p1))
    axis.set_xlabel('Measured SPM(mg/L)', fontsize=style['label_font_size'])
    axis.set_ylabel('Estimated SPM(mg/L)', fontsize=style['label_font_size'])
    # plt.title('VALIDATION', fontsize=11)
    y = 75
    for val_type, val in validation.items():
        # print(val_type)
        axis.text(1, y, '{}: {}'.format(val_type, f'{val:.2f}'),
                    fontsize=style['plot_font_size'], zorder=34)
        y = y - 5  # to keep adding from top of y axis to bottom(0), x is constant

    # axis.text(80, 80, '(b)', fontsize=style['plot_font_size'], zorder=34)
    axis.text(-0.15, 1.0, '(b)', transform=axis.transAxes, size=12, weight='bold')

    axis.axis('square')

def add_algorithm_plot(training_y, training_x, model, axis, r2, equation_latex, style):
    colors = ['b']
    axis.scatter(training_x,training_y,  c=colors, s=style['point_size'])

    z = np.polyfit(training_x, training_y, 1)
    p = np.poly1d(z)

    #add trendline to plot
    axis.plot(training_x, p(training_x))

    axis.set_ylabel('Measured SPM(mg/L)', fontsize=style['label_font_size'])
    # axis.set_xlabel('Reflectance(Rrs)', fontsize=10)
    axis.set_xlabel('$%s$' % equation_latex, fontsize=style['label_font_size'])
    # axis.text(0.6, 0.6, r'$\mathcal{A}\mathrm{sin}(2 \omega t)$',
    #          fontsize=20)
    axis.set_box_aspect(1)
    # plt.title('VALIDATION', fontsize=11)
    equation = 'y={:.2f}x+{:.2f}'.format(np.ravel(model.coef_)[0],np.ravel(model.intercept_)[0])

    if np.ravel(model.intercept_)[0] < 0:
        equation = 'y={:.2f}x{:.2f}'.format(np.ravel(model.coef_)[0], np.ravel(model.intercept_)[0])

    # axis.text(max(training_x)+5, max(training_y), '(a)', fontsize=style['plot_font_size'], zorder=34)
    axis.text(-0.15, 1.0, '(a)', transform=axis.transAxes, size=12, weight='bold')
    axis.text(min(training_x), max(training_y), equation,
              fontsize=style['plot_font_size'], zorder=34)
    axis.text(min(training_x), max(training_y)-5, r'$R^{2}$'+'={:.2f}'.format(r2),
              fontsize=style['plot_font_size'], zorder=34)


def read_xy_data(x_path, y_path):
    # read data
    try:
        x_df = pd.read_excel(x_path)
    except Exception as ex:
        x_df = pd.read_csv(x_path, 'rb', encoding='utf8')
    x_df = x_df.reset_index()  # make sure indexes pair with number of rows
    # read data
    try:
        y_df = pd.read_excel(y_path)
    except Exception as ex:
        y_df = pd.read_csv(x_path, 'rb', encoding='utf8')
    y_df = y_df.where((pd.notnull(y_df)), None)
    y_df = y_df.reset_index()  # make sure indexes pair with number of rows
    # organize by date and average SPM data - x data
    return x_df, y_df

def read_data(data_path):
    # read data.

    try:
        df = pd.read_excel(data_path)
    except Exception as ex:
        df = pd.read_csv(data_path, 'rb', encoding='utf8')
    df = df.reset_index()  # make sure indexes pair with number of rows
    df = df.where((pd.notnull(df)), None)

    # read data
    return df
def group_data_by_date(spm_df, discharge_df, spm_date_field, uas_spm_field, discharge_date_field,
                       discharge_fields, equal_dates=False, insitu_spm_field=None):
    y_daily_values_cont = []
    insitu_spm_values = OrderedDict()
    # if len(y_fields) > 1:
    for discharge_field in discharge_fields:
        discharge_daily_values = OrderedDict()
        for index, row in discharge_df.iterrows():
            curr_date = row[discharge_date_field].to_pydatetime()
            curr_date_str = curr_date.strftime("%m-%d-%Y")
            # print (row[y_field])
            if row[discharge_field] is not None:
                discharge_daily_values[curr_date_str] = float(row[discharge_field])
            else:
                discharge_daily_values[curr_date_str] = 0

        y_daily_values_cont.append(discharge_daily_values)


    uas_spm_daily_data = OrderedDict()
    uas_spm_daily_values = OrderedDict()

    for index, row in spm_df.iterrows():
        curr_date = row[spm_date_field].to_pydatetime()
        curr_date_str = curr_date.strftime("%m-%d-%Y")
        if insitu_spm_field is not None:
            insitu_spm_values[curr_date_str] = row[insitu_spm_field]
        if curr_date_str not in uas_spm_daily_data.keys():
            uas_spm_daily_data[curr_date_str] = []

        uas_spm_daily_data[curr_date_str].append(row[uas_spm_field])
    # average SPM data - x data
    date_objs = []
    if equal_dates:
        daily_x_dates = uas_spm_daily_data.keys()
        for curr_y_date, y_data in y_daily_values_cont[0].items():
            date_objs.append(datetime.strptime(curr_y_date, '%m-%d-%Y'))
            if curr_y_date in daily_x_dates:
                # print (uas_spm_daily_data[curr_y_date])
                uas_spm_daily_values[curr_y_date] = np.mean(uas_spm_daily_data[curr_y_date])
            else:
                uas_spm_daily_values[curr_y_date] = None
    else:
        for curr_date, x_data in uas_spm_daily_data.items():
            # print (x_data)
            uas_spm_daily_values[curr_date] = np.mean(x_data)
            date_objs.append(datetime.strptime(curr_date, '%m-%d-%Y'))
    # organize discharge data by date
    # print(insitu_spm_values)
    return date_objs, uas_spm_daily_values, y_daily_values_cont, insitu_spm_values

# spm_data, spm_uas_data, spm_uas_mean, spm_mean = setup_boxplot_data(
#     None, 'Date', '%m-%d-%Y'
# )
def run_plot():
    central_lon = np.mean(extent[:2])
    central_lat = np.mean(extent[2:])
    fig, axes = plt.subplots(2, 3, figsize=(45, 37), subplot_kw={
        'projection': ccrs.AlbersEqualArea(central_lon, central_lat)})
    sal = read_tif_to_array('D:/MSU/RA/acidification/acidification/salinity/salinity.tif')
    ph = read_tif_to_array('D:/MSU/RA/acidification/acidification/ph/pH.tif')
    ast = read_tif_to_array('D:/MSU/RA/acidification/acidification/ast/omega_Ar.tif')
    css = read_tif_to_array('D:/MSU/RA/acidification/acidification/css/omega_Ca.tif')
    spcp = read_tif_to_array('D:/MSU/RA/acidification/acidification/spcp/spCO2_micro.tif')
    ta = read_tif_to_array('D:/MSU/RA/acidification/acidification/ta/TA.tif')

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
    fig.savefig('D:/MSU/RA/acidification/acidif'
                'ication/acid2022_c.png', dpi=150, bbox_inches='tight')
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

# run_plot()