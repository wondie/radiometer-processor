import csv
from collections import OrderedDict
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
from datetime import datetime, timedelta
import pandas as pd
from plot import save_plot
# TODO change stats label date to monthly to change the plot x values to monthly
TIME_SERIES = r'D:\MSU\dissertation\SPM_Multi-spectral\data\discharge\monthly\sites_SPM_monthly.xlsx'
MONTHS = ['2018_03', '2018_05', '2018_06', '2018_07', '2018_12', '2019_06',
          '2019_07', '2021_07']
stats = {'2018-03': {'label': '03-01-2018', 'mean': 48.49413536058343,
                     'iqr': 17.808972000000004, 'cilo': 47.60989,
                     'cihi': 47.60989, 'whishi': 70.75, 'whislo': 18.75,
                     'fliers': [], 'q1': 39.423992, 'med': 47.609985,
                     'q3': 57.232964},
         '2018-05': {'label': '05-01-2018', 'mean': 50.66261584940906,
                     'iqr': 17.96060175000001, 'cilo': 50.934883,
                     'cihi': 50.934883, 'whishi': 70.75, 'whislo': 18.75,
                     'fliers': [], 'q1': 42.03115124999999, 'med': 50.9356,
                     'q3': 59.991753},
         '2018-06': {'label': '06-01-2018', 'mean': 53.106575937526294,
                     'iqr': 15.364933, 'cilo': 53.151974, 'cihi': 53.151974,
                     'whishi': 70.74999, 'whislo': 22.800732, 'fliers': [],
                     'q1': 45.848103, 'med': 53.15095, 'q3': 61.213036},
         '2018-07': {'label': '07-01-2018', 'mean': 48.4147312334088,
                     'iqr': 17.568292, 'cilo': 47.791264, 'cihi': 47.791264,
                     'whishi': 70.75, 'whislo': 18.750008, 'fliers': [],
                     'q1': 39.791958, 'med': 47.790653, 'q3': 57.36025},
         '2018-12': {'label': '12-01-2018', 'mean': 54.76202082100937,
                     'iqr': 15.137205999999999, 'cilo': 54.438995,
                     'cihi': 54.438995, 'whishi': 70.75, 'whislo': 24.823877,
                     'fliers': [], 'q1': 47.52966, 'med': 54.439846,
                     'q3': 62.666866},
         '2019-06': {'label': '06-01-2019', 'mean': 39.861843893526085,
                     'iqr': 17.317583, 'cilo': 37.521755, 'cihi': 37.521755,
                     'whishi': 70.75, 'whislo': 18.75, 'fliers': [],
                     'q1': 30.319347, 'med': 37.52079, 'q3': 47.63693},
         '2019-07': {'label': '07-01-2019', 'mean': 37.61548303090227,
                     'iqr': 15.914034999999998, 'cilo': 34.929714,
                     'cihi': 34.929714, 'whishi': 68.45436, 'whislo': 18.75,
                     'fliers': [], 'q1': 28.669275, 'med': 34.9295,
                     'q3': 44.58331},
         '2021-07': {'label': '07-01-2021', 'mean': 43.84053343948744,
                     'iqr': 19.364843, 'cilo': 43.622818, 'cihi': 43.622818,
                     'whishi': 70.75, 'whislo': 18.75, 'fliers': [],
                     'q1': 33.70563, 'med': 43.62234, 'q3': 53.070473}
         }


def read_pixel_value_from_file(file_path, date_str):
    pixel_values = []

    # with open(file_path, 'r') as f:
    #     pixel_values = f.read().split(',')
    with open(file_path, 'r', newline='') as plentydata:
        reader = csv.reader(plentydata, delimiter=',',
                            quoting=csv.QUOTE_NONNUMERIC)
        # data = list(reader)
        # print('Read', file_path)
        for row in reader:
            pixel_values.append(row)
        print('finished reading ', file_path)
    stats = cbook.boxplot_stats(pixel_values, labels=[date_str], bootstrap=1)

    # print(stats)
    return stats


def txt_file_to_pixel_values():
    pixel_values = OrderedDict()
    for month in MONTHS:
        date_str = month.split('_')
        curr_date_str = '{}-{}'.format(date_str[1], date_str[0])
        file_path = 'E:/SPMOutput/{}.txt'.format(curr_date_str)
        stats = read_pixel_value_from_file(file_path, curr_date_str)
        if len(stats) > 0:
            pixel_values[curr_date_str] = stats[0]

    return pixel_values


def raster_to_pixel_values():
    pixel_values = OrderedDict()
    for month in MONTHS:
        raster = "E:/SPMOutput/SPM_{}.tif".format(month)
        ds = gdal.Open(raster)
        date_str = month.split('_')
        curr_date_str = '{}-{}'.format(date_str[1], date_str[0])
        # print(ds.GetRasterBand(1).ReadAsArray())
        # print(curr_date_str)

        pixel_values[curr_date_str] = [element for v in
                                       ds.GetRasterBand(1).ReadAsArray() for
                                       element in v if element > 0]
        # print('Finish ', curr_date_str)
        with open('E:/SPMOutput/{}.txt'.format(curr_date_str), 'w') as f:
            f.write(','.join([str(x) for x in pixel_values[curr_date_str]]))

    return pixel_values


def setup_boxplot_data_daily(path=None, date='Date', format="%Y-%m"):
    if path is None:
        path = TIME_SERIES
    # global box_plot_mean, uas_spm_values, index, insitu_spm_values
    df = pd.read_excel(path)
    df = df.replace(np.nan, None)

    df = df.reset_index()  # make sure indexes pair with number of rows
    insitu_spm_data = OrderedDict()
    insitu_spm_mean = []
    uas_spm_data = OrderedDict()
    uas_spm_mean = []
    for index, row in df.iterrows():
        curr_date = row[date].to_pydatetime()
        curr_date_str = curr_date.strftime(format)
        if curr_date_str not in insitu_spm_data.keys():
            # uas_spm_values.append(row['UAS SPM'])
            insitu_spm_data[curr_date_str] = []
            uas_spm_data[curr_date_str] = []
        if row['UAS SPM'] is not None:
        # if row['UAS SPM'] is not None:
            uas_spm_data[curr_date_str].append(row['UAS SPM'])
        if row['Insitu SPM'] is not None:
            insitu_spm_data[curr_date_str].append(row['Insitu SPM'])
    # labels = insitu_spm_data.keys()
    # print(uas_spm_data)
    insitu_spm_values = insitu_spm_data.values()
    uas_spm_values = uas_spm_data.values()
    for ls in insitu_spm_values:
        insitu_spm_mean.append(np.mean(ls))
    for sv in uas_spm_values:
        uas_spm_mean.append(np.mean(sv))
    return insitu_spm_data, uas_spm_data, uas_spm_mean, insitu_spm_mean

def setup_boxplot_data_monthly(path=None, date='Date', month_format="%Y-%m", convert_to_date=True, date_format='%m-%d-%Y'):
    if path is None:
        path = TIME_SERIES
    # global box_plot_mean, uas_spm_values, index, insitu_spm_values
    df = pd.read_excel(path)
    df = df.replace(np.nan, None)

    df = df.reset_index()  # make sure indexes pair with number of rows
    insitu_spm_data = OrderedDict()
    insitu_spm_mean = OrderedDict()
    uas_spm_data = OrderedDict()
    uas_spm_mean = OrderedDict()
    for index, row in df.iterrows():
        curr_date = row[date].to_pydatetime()

        curr_date_str = curr_date.strftime(month_format)
        print (curr_date)


        if curr_date_str not in insitu_spm_data.keys():
            # uas_spm_values.append(row['UAS SPM'])
            # Group months by mission. The starting month will be the mission month
            close_date = curr_date - timedelta(days=5)
            curr_date_str = close_date.strftime(month_format)
            if curr_date_str not in insitu_spm_data.keys():
                insitu_spm_data[curr_date_str] = []
                uas_spm_data[curr_date_str] = []
        if convert_to_date:
            del insitu_spm_data[curr_date_str]
            del uas_spm_data[curr_date_str]
            curr_date_str = datetime.strptime(curr_date_str, month_format).\
                date().strftime(date_format)
            if curr_date_str not in insitu_spm_data.keys():
                insitu_spm_data[curr_date_str] = []
                uas_spm_data[curr_date_str] = []
        if row['UAS SPM'] is not None:
            # if row['UAS SPM'] is not None:
            uas_spm_data[curr_date_str].append(row['UAS SPM'])
        if row['Insitu SPM'] is not None:
            insitu_spm_data[curr_date_str].append(row['Insitu SPM'])
    # labels = insitu_spm_data.keys()
    # print(uas_spm_data)
    # insitu_spm_values = insitu_spm_data.values()
    # uas_spm_values = uas_spm_data.values()
    for dt, ls in insitu_spm_data.items():
        if dt not in insitu_spm_mean.keys():
            insitu_spm_mean[dt] = []
        insitu_spm_mean[dt] = np.mean(ls)
    for dt, ls in uas_spm_data.items():
        if dt not in uas_spm_mean.keys():
            uas_spm_mean[dt] = []
        uas_spm_mean[dt] = np.mean(ls)

    return insitu_spm_data, uas_spm_data, uas_spm_mean, insitu_spm_mean

insitu_spm_data, spm_uas_data, spm_uas_mean, spm_mean = setup_boxplot_data_monthly(
    r'D:\MSU\dissertation\SPM_Multi-spectral\data\sites_SPM.xlsx', 'Date', '%m-%Y'
)
print(insitu_spm_data)
print(spm_uas_data)
print(spm_uas_mean)
print(spm_mean)
def add_box_plots(plt, axes, insitu_spm_values, spm_uas_mean, spm_mean, format):
    plt.figure(figsize=(10, 5))

    axes.plot(np.arange(len(spm_uas_mean)) + 1, spm_uas_mean, marker='o',
              color='#bf5700', label='UAS')
    axes.plot(np.arange(len(spm_mean)) + 1, spm_mean, marker='o',
              color='#2177e2', label='In situ')
    # plt.legend()
    axes.legend(loc='upper right')
    boxprops = dict(color="#2177e2", linewidth=1.5)
    medianprops = dict(color="#2177e2", linewidth=1.5)
    flierprops = dict(marker='o', markerfacecolor='none', markersize=7,
                      linestyle='none', markeredgecolor='#2177e2')
    # print (box_plot_values)
    # if format is not None:
    #     temp_insitu_spm_values = OrderedDict()
    #     for f
    axes.boxplot(
        insitu_spm_values, labels=['', '', '', '', '', '', '', ''],
        boxprops=boxprops, medianprops=medianprops,
        capprops=dict(color="#2177e2"),
        whiskerprops=dict(color="#2177e2"),
        flierprops=flierprops
    )
    boxprops2 = dict(color='#bf5700', linewidth=1.5)
    medianprops2 = dict(color='#bf5700', linewidth=1.5)
    flierprops2 = dict(marker='o', markerfacecolor='none', markersize=7,
                       linestyle='none', markeredgecolor='#bf5700')
    # raster_arrays = raster_to_arrays()
    # pixel_values = txt_file_to_pixel_values()
    pixel_values = stats
    pixel_values_value = pixel_values.values()
    # plt.boxplot(pixel_values_value, labels=labels,  boxprops=boxprops2, medianprops=medianprops2,
    #             capprops=dict(color="red"),
    #             whiskerprops=dict(color="red") )
    axes.bxp(
        pixel_values_value, widths=(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
        boxprops=boxprops2, medianprops=medianprops2,
        capprops=dict(color="#bf5700"),
        whiskerprops=dict(color="#bf5700"),
        flierprops=flierprops2
    )
    # axes.set_title('Default')
    # y_axis = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    # y_values = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    # plt.yticks(y_axis, y_values)
    # plt.ylabel('SPM (mg/L)')


if __name__ == '__main__':
    insitu_spm_values, uas_spm_values, spm_uas_mean, spm_mean = setup_boxplot_data_daily()
    # print (uas_spm_values)'
    fig, axes = plt.subplots(1, 1)
    add_box_plots(plt, axes, insitu_spm_values.values(), spm_uas_mean, spm_mean, format='%Y-%m')
    axes.grid(True)
    # plt.xlabel('Sampling Periods')
    axes.set(xlabel="Sampling Periods", ylabel='SPM (mg/L)')
    plt.tight_layout()

    # plt.show()
    save_plot(fig, 'poster')
