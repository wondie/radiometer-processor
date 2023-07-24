from collections import OrderedDict
from datetime import datetime, date

import matplotlib.pyplot as plt
import matplotlib.dates
import pandas as pd
import seaborn
import numpy as np

from box_plot_line import setup_boxplot_data_monthly, add_box_plots, stats

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
plt.rcParams.update(params)
from plot import save_plot, read_xy_data, group_data_by_date

def plot_data(spm_path, discharge_path, insitu_spm_field, uas_spm_field,
              spm_date_field, y1_field, y2_field, y3_field, y4_field,
              spm_label,discharge1_label, discharge2_label, discharge_date_field,
              residence_time=None):
    spm_df, discharge_df = read_xy_data(spm_path, discharge_path)

    date_objs, spm_daily_values, discharge_daily_values_cont, uas_spm_values = \
        group_data_by_date(spm_df, discharge_df, spm_date_field,
        uas_spm_field, discharge_date_field,
        [y1_field, y2_field, y3_field, y4_field], False, insitu_spm_field)
    new_discharge_daily_values_cont = []
    insitu_spm_data, spm_uas_data, spm_uas_mean, insitu_spm_mean = \
        setup_boxplot_data_monthly(
            spm_path, 'Date', '%m-%Y'
        )
    spm_date_objs = []
    spm_date_dict = OrderedDict()

    for cur_date, spm_dat in insitu_spm_data.items():
        curr_date_obj = datetime.strptime(cur_date, '%m-%d-%Y')
        spm_date_objs.append(curr_date_obj)
        spm_date_dict[cur_date] = curr_date_obj

    for i, discharge_daily_values in enumerate(discharge_daily_values_cont):
        discharge_plot_data = OrderedDict()
        for curr_date, discharge_data in discharge_daily_values.items():
            curr_date_obj = datetime.strptime(curr_date, '%m-%d-%Y')
            # Add the residence day from the curr_date of discharge to find matching
            # SPM data recorded

            if i == 3: # Index 2 is Bonnet Caree
                if discharge_data is None:
                    discharge_plot_data[curr_date_obj] = 0
                else:
                    discharge_plot_data[curr_date_obj] = discharge_data
            else:
                discharge_plot_data[curr_date_obj] = discharge_data

        new_discharge_daily_values_cont.append(discharge_plot_data)

    fig, ax = plt.subplots()
    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines.right.set_position(("axes", 1.06))

    # dates = matplotlib.dates.date2num(date_objs)
    # ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b %Y'))
    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())

    left = date(2018, 2, 1)
    right = date(2021,8, 1)
    print ('Pearson correlation coefficients: ', np.corrcoef(np.array(list(insitu_spm_mean.values())), np.array(list(spm_uas_mean.values()))))
    p, = ax.plot(spm_date_objs, insitu_spm_mean.values(),
              color="#5c0f7c",  markersize=14, marker='.',
                 linestyle = '-',label=insitu_spm_field)

    p0, = ax.plot(spm_date_objs, spm_uas_mean.values(),
                  color="#c284dc",  markersize=14,marker='.',
                  linestyle = '-',label=uas_spm_field)

    boxprops = dict(color="#2177e2", linewidth=1.5)
    medianprops = dict(color="#2177e2", linewidth=1.5)
    flierprops = dict(marker='o', markerfacecolor='none', markersize=7,
                      linestyle='none', markeredgecolor='#2177e2')


# print(np.array( [ np.random.normal( i, 1, 10 ) for i in range(3) ] ))
    # ax.set_xticklabels(list(insitu_spm_data.keys()), rotation=45 )
    boxprops2 = dict(color='#bf5700', linewidth=1.5)
    medianprops2 = dict(color='#bf5700', linewidth=1.5)
    flierprops2 = dict(marker='o', markerfacecolor='none', markersize=7,
                       linestyle='none', markeredgecolor='#bf5700')
    # raster_arrays = raster_to_arrays()
    # pixel_values = txt_file_to_pixel_values()
    # pixel_values = stats
    # pixel_values_value = pixel_values.values()
    # final_pixel_values = []
    # # Replace string date with date obj
    # for val in pixel_values_value:
    #     val['label'] = spm_date_dict[val['label']]
    #     final_pixel_values.append(val)

    # plt.boxplot(pixel_values_value, labels=labels,  boxprops=boxprops2, medianprops=medianprops2,
    #             capprops=dict(color="red"),
    #             whiskerprops=dict(color="red") )
    # print(final_pixel_values)
    # ax.bxp(
    #     final_pixel_values, widths=(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
    #     boxprops=boxprops2, medianprops=medianprops2,
    #     capprops=dict(color="#bf5700"),
    #     whiskerprops=dict(color="#bf5700"),
    #     flierprops=flierprops2
    # )
    # ax.set(xlabel="Sampling Periods", ylabel='SPM (mg/L)')
    # add_box_plots(ax, insitu_spm_data.values())
    p1, = twin1.plot(new_discharge_daily_values_cont[0].keys(),
                     new_discharge_daily_values_cont[0].values(), markersize=0,
                          linestyle = '-', color="#d4a373", label=y1_field)

    p2, = twin1.plot(new_discharge_daily_values_cont[1].keys(),
                     new_discharge_daily_values_cont[1].values(), markersize=0,
                          linestyle = '-',color="#bc4b51", label=y2_field)

    p3, = twin2.plot(new_discharge_daily_values_cont[2].keys(),
                     new_discharge_daily_values_cont[2].values(), markersize=0,
                          linestyle = '-',color="#0ead69", label=y3_field)

    values_y3_data = []
    for j in new_discharge_daily_values_cont[3].values():
        values_y3_data.append(float(j))
    p4 = twin2.bar(list(new_discharge_daily_values_cont[3].keys()),
                   values_y3_data, color= '#7cb518', label=y4_field)
    # twin2.set_ylim(twin2.get_ylim())

    plt.gcf().set_size_inches(20, 12)
    plt.gca().set_xbound(left, right)
    ax.set_xlabel(spm_date_field)
    ax.set_ylabel(spm_label)

    ax.set_ylim(10, max(spm_daily_values.values())*1.5)
    twin1.set_ylim(0, max(new_discharge_daily_values_cont[1].values())*1.5)
    twin2.set_ylim(0, max(values_y3_data)*1.5)

    twin1.invert_yaxis()
    twin2.invert_yaxis()

    twin1.set_ylabel(discharge1_label)
    twin2.set_ylabel(discharge2_label)
    ax.yaxis.label.set_color(p.get_color())
    # ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())

    ax.tick_params(axis='y', colors=p.get_color())
    twin1.tick_params(axis='y', colors=p2.get_color())
    twin2.tick_params(axis='y', colors=p3.get_color())

    ax.tick_params(axis='y', colors=p.get_color())
    twin1.tick_params(axis='y', colors=p2.get_color())
    twin2.tick_params(axis='y', colors=p3.get_color())
    # print (insitu_spm_data.values())
    spm_pos = pd.DatetimeIndex(discharge_daily_values_cont[0].keys())
    # print ('Pos ', spm_pos)
    for dt, val in discharge_daily_values_cont[0].items():
        if dt not in insitu_spm_data.keys():
            insitu_spm_data[dt] = []

    spm_indexes = [spm_pos.get_loc(d) for d in insitu_spm_data.keys()]
    print(spm_indexes)
    # ax.boxplot(
    #     list(insitu_spm_data.values()), positions=spm_indexes,
    #     boxprops=boxprops, medianprops=medianprops,
    #     capprops=dict(color="#2177e2"),
    #     whiskerprops=dict(color="#2177e2"),
    #     flierprops=flierprops
    # )
    # seaborn.boxplot(x=list(insitu_spm_data.keys()), y=list(insitu_spm_data.values()))
    # ax.set_xticklabels(list(discharge_daily_values_cont[0].keys()), rotation=45 )

    plt.gcf().autofmt_xdate()

    plt.tight_layout()
    ax.legend(handles=[p, p0, p1, p2, p3, p4], loc=4)
    save_plot(plt, 'poster', 'SPM and Discharge')

SPM_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\sites_SPM.xlsx'
discharge_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\discharge_daily_combined.xlsx'

plot_data(SPM_path, discharge_path, 'Insitu SPM', 'UAS SPM', 'Date', 'Jourdan River', 'Wolf River','Pearl River','Bonnet Carre Spillway', 'SPM (mg/L)', 'Gauge Height','Discharge (cpfs)','Date', 4)

