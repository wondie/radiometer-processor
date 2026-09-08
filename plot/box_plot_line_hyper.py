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
TIME_SERIES = r'G:\Other computers\My Laptop\dissertation\SPM_Hyperspectral\data\sites_SPM_hyper.xlsx'
MONTHS = ['2018_03', '2018_05', '2018_06', '2018_07', '2018_12', '2019_06',
          '2019_07']
stats = {'2018-03': {'label': '03-01-2018', 'mean': 48.92021942138672,
                     'iqr': 16.231426239013672, 'cilo': 46.18420247191051,
                     'cihi': 46.18786021119496, 'whishi': 79.80335235595703, 'whislo': 17.424707412719727,
                     'fliers': [], 'q1': 39.22479248046875, 'med': 46.186031341552734,
                     'q3': 55.45621871948242},
         '2018-05': {'label': '05-01-2018', 'mean': 47.81859588623047,
                     'iqr': 17.5078125, 'cilo': 43.98210815210446,
                     'cihi': 43.98552032689945, 'whishi': 80.37709045410156, 'whislo': 17.356876373291016,
                     'fliers': [], 'q1': 36.60756301879883, 'med': 43.98381423950195,
                     'q3': 54.11537551879883},
         '2018-06': {'label': '06-01-2018', 'mean': 47.594200134277344,
                     'iqr': 16.250728607177734, 'cilo': 42.06770646397376,
                     'cihi': 42.07066023524499, 'whishi': 76.36460876464844, 'whislo': 16.938583374023438,
                     'fliers': [], 'q1': 35.737789154052734, 'med': 42.069183349609375,
                     'q3': 51.98851776123047},
         '2018-07': {'label': '07-01-2018', 'mean': 67.10711669921875,
                     'iqr': 28.97795867919922, 'cilo': 60.02605834751214,
                     'cihi': 60.03306946010505, 'whishi': 119.656005859375, 'whislo': 16.895000457763672,
                     'fliers': [], 'q1': 47.21111297607422, 'med': 60.029563903808594,
                     'q3': 76.18907165527344},
         '2018-12': {'label': '12-01-2018', 'mean': 86.19490051269531,
                     'iqr': 44.246673583984375, 'cilo': 71.68429910347876,
                     'cihi': 71.69581259085717, 'whishi': 164.64617919921875, 'whislo': 16.895000457763672,
                     'fliers': [], 'q1': 54.02949523925781, 'med': 71.69005584716797,
                     'q3': 98.27616882324219},
         '2019-06': {'label': '06-01-2019', 'mean': 45.17489242553711,
                     'iqr': 15.114555358886719, 'cilo': 40.98647488248668,
                     'cihi': 40.98969088899769, 'whishi': 72.73553466796875, 'whislo': 17.165565490722656,
                     'fliers': [], 'q1': 34.94915008544922, 'med': 40.98808288574219,
                     'q3': 50.06370544433594},
         '2019-07': {'label': '07-01-2019', 'mean': 38.72906494140625,
                     'iqr': 11.808418273925781, 'cilo': 35.71976508579002,
                     'cihi': 35.72148857631935, 'whishi': 60.74414825439453, 'whislo': 17.133949279785156,
                     'fliers': [], 'q1': 31.223102569580078, 'med': 35.72062683105469,
                     'q3': 43.03152084350586}
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

def add_box_plots(plt, axes, insitu_spm_values, uas_spm_values, spm_uas_mean, spm_mean, format):
    insitu_spm_values = list(insitu_spm_values)
    uas_spm_values = list(uas_spm_values)
    positions = np.arange(len(insitu_spm_values)) + 1
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
        insitu_spm_values, positions=positions - 0.15, widths=0.25,
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
    pixel_values_value = uas_spm_values
    # plt.boxplot(pixel_values_value, labels=labels,  boxprops=boxprops2, medianprops=medianprops2,
    #             capprops=dict(color="red"),
    #             whiskerprops=dict(color="red") )
    axes.boxplot(
        pixel_values_value, positions=positions + 0.15, widths=0.25,
        boxprops=boxprops2, medianprops=medianprops2,
        capprops=dict(color="#bf5700"),
        whiskerprops=dict(color="#bf5700"),
        flierprops=flierprops2
    )
    axes.set_xticks(positions)
    axes.set_xticklabels([''] * len(positions))
    # axes.set_title('Default')
    # y_axis = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    # y_values = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    # plt.yticks(y_axis, y_values)
    # plt.ylabel('SPM (mg/L)')


if __name__ == '__main__':
    # Use monthly grouping (mm-YYYY format)
    insitu_spm_data, uas_spm_data, spm_uas_mean, spm_mean = setup_boxplot_data_daily(format='%m-%Y')
    
    # Filter out 2021 data
    filtered_insitu = OrderedDict()
    filtered_uas = OrderedDict()
    filtered_insitu_mean = []
    filtered_uas_mean = []
    
    for (period, insitu_vals), (_, uas_vals), insitu_m, uas_m in zip(
        insitu_spm_data.items(), uas_spm_data.items(), spm_mean, spm_uas_mean):
        if not period.endswith('-2021'):
            filtered_insitu[period] = insitu_vals
            filtered_uas[period] = uas_vals
            filtered_insitu_mean.append(insitu_m)
            filtered_uas_mean.append(uas_m)
    
    insitu_spm_data = filtered_insitu
    uas_spm_data = filtered_uas
    spm_mean = filtered_insitu_mean
    spm_uas_mean = filtered_uas_mean
    
    # print (uas_spm_values)'
    fig, axes = plt.subplots(1, 1, figsize=(10, 5))
    
    insitu_spm_values = list(insitu_spm_data.values())
    uas_spm_values = list(uas_spm_data.values())
    positions = np.arange(len(insitu_spm_values)) + 1
    
    print(f"Number of periods: {len(positions)}")
    print(f"Insitu values: {len(insitu_spm_values)}")
    print(f"UAS values: {len(uas_spm_values)}")
    
    # Check data content in detail
    print("\n" + "="*60)
    print("DATA CHECK")
    print("="*60)
    for i, (period_name, insitu, uas) in enumerate(zip(insitu_spm_data.keys(), insitu_spm_values, uas_spm_values)):
        print(f"\nPeriod {i+1} ({period_name}):")
        print(f"  Insitu: {len(insitu)} values -> {insitu}")
        print(f"  UAS:    {len(uas)} values -> {uas}")
        if len(insitu) <= 1:
            print(f"  WARNING: Insitu has {len(insitu)} value(s) - boxplot needs multiple values!")
        if len(uas) <= 1:
            print(f"  WARNING: UAS has {len(uas)} value(s) - boxplot needs multiple values!")
    print("="*60 + "\n")
    
    # Filter out empty lists
    insitu_spm_values_filtered = [v if len(v) > 0 else [0] for v in insitu_spm_values]
    uas_spm_values_filtered = [v if len(v) > 0 else [0] for v in uas_spm_values]
    
    # Plot boxplot ONLY for Insitu SPM (centered at positions)
    # Styled to match Chapter 1's boxplot_line.svg reference: peachpuff box,
    # brown outline/whiskers, blue "UAS" / green "In situ" mean lines,
    # horizontal-only gridlines, non-rotated period labels.
    boxprops = dict(color="#915f0e", linewidth=1.5, facecolor='#ffdab9', alpha=0.5)
    medianprops = dict(color="#915f0e", linewidth=1.5)
    flierprops = dict(marker='o', markerfacecolor='none', markersize=7,
                      linestyle='none', markeredgecolor='#915f0e')
    bp1 = axes.boxplot(
        insitu_spm_values_filtered, positions=positions, widths=0.5,
        patch_artist=True,
        boxprops=boxprops, medianprops=medianprops,
        capprops=dict(color="#915f0e", linewidth=1.5),
        whiskerprops=dict(color="#915f0e", linewidth=1.5),
        flierprops=flierprops,
        zorder=1
    )

    # Plot line graphs for mean values ON TOP
    axes.plot(positions, spm_uas_mean, marker='o', color='#2c42b3', label='UAS',
              linewidth=1.5, markersize=6, zorder=2)
    axes.plot(positions, spm_mean, marker='o', color='#2cb330', label='In situ',
              linewidth=1.5, markersize=6, zorder=2)
    axes.legend(loc='upper right')

    # Set x-axis labels with sampling periods (horizontal, not rotated)
    axes.set_xticks(positions)
    axes.set_xticklabels(list(insitu_spm_data.keys()))

    axes.grid(axis='y', color='#b0b0b0', linewidth=0.8, zorder=0)
    axes.set(xlabel="Sampling Periods", ylabel='SPM (mg/L)')
    plt.tight_layout()

    # Save PNG for quick preview and SVG for paper/publication use (vector,
    # matching the reference image's own format).
    output_path = r'G:\Other computers\My Laptop\codes\radiometer_processor\data\output\boxplot_insitu_uas.png'
    svg_output_path = r'G:\Other computers\My Laptop\codes\radiometer_processor\data\output\boxplot_insitu_uas.svg'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(svg_output_path, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")
    print(f"Paper/vector version saved to: {svg_output_path}")
    
    plt.show()
    plt.close()
