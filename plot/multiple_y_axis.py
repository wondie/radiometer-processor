from collections import OrderedDict
from datetime import datetime
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates
import numpy as np

from box_plot_line import setup_boxplot_data, add_box_plots, stats

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
plt.rcParams.update(params)
from plot import save_plot, read_xy_data, group_data_by_date

def plot_data(y0_path, y_path, y01_field, y0_field, y0_date_field, y1_field, y2_field,
              y3_field,y4_field, y0_label,y1_label, y2_label, y_date_field, residence_time=None):
    y0_df, y_df = read_xy_data(y0_path, y_path)

    date_objs, y0_daily_values, y_daily_values_cont, y01_values = group_data_by_date(y0_df, y_df, y0_date_field,
                                                        y0_field, y_date_field,
                                                        [y1_field, y2_field, y3_field, y4_field], False,y01_field=y01_field)
    new_y_daily_values_cont = []
    spm_data, spm_uas_data, spm_uas_mean, spm_mean = setup_boxplot_data(
        'Date', '%m-%d-%Y'
    )
    spm_date_objs = []
    spm_date_dict = OrderedDict()
    for cur_date, spm_dat in spm_data.items():
        curr_date_obj = datetime.strptime(cur_date, '%m-%d-%Y')
        spm_date_objs.append(curr_date_obj)
        spm_date_dict[cur_date] = curr_date_obj

    for i, y_daily_values in enumerate(y_daily_values_cont):
        y_plot_data = OrderedDict()
        for curr_date, y_data in y_daily_values.items():
            curr_date_obj = datetime.strptime(curr_date, '%m-%d-%Y')
            # Add the residence day from the curr_date of discharge to find matching
            # SPM data recorded

            if i == 3: # Index 2 is Bonnet Caree
                if y_data is None:
                    y_plot_data[curr_date_obj] = 0
                else:
                    y_plot_data[curr_date_obj] = y_data
            else:
                y_plot_data[curr_date_obj] = y_data

        new_y_daily_values_cont.append(y_plot_data)

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

    left = dt.date(2018, 3, 1)
    right = dt.date(2021,8, 1)


    p, = ax.plot(spm_date_objs, spm_data.values(),
              color="#5c0f7c",  markersize=14, marker='.',
                 linestyle = '-',label=y01_field)

    p0, = ax.plot(spm_date_objs, spm_uas_data.values(),
                  color="#c284dc",  markersize=14,marker='.',
                  linestyle = '-',label=y0_field)

    boxprops = dict(color="#2177e2", linewidth=1.5)
    medianprops = dict(color="#2177e2", linewidth=1.5)
    flierprops = dict(marker='o', markerfacecolor='none', markersize=7,
                      linestyle='none', markeredgecolor='#2177e2')

    # print (box_plot_values)
    # ax.boxplot(
    #     spm_data.values(), labels=['', '', '', '', '', '', '', ''],
    #     boxprops=boxprops, medianprops=medianprops,
    #     capprops=dict(color="#2177e2"),
    #     whiskerprops=dict(color="#2177e2"),
    #     flierprops=flierprops
    # )
    #
    # boxprops2 = dict(color='#bf5700', linewidth=1.5)
    # medianprops2 = dict(color='#bf5700', linewidth=1.5)
    # flierprops2 = dict(marker='o', markerfacecolor='none', markersize=7,
    #                    linestyle='none', markeredgecolor='#bf5700')
    # raster_arrays = raster_to_arrays()
    # pixel_values = txt_file_to_pixel_values()
    pixel_values = stats
    pixel_values_value = pixel_values.values()
    final_pixel_values = []
    # Replace string date with date obj
    for val in pixel_values_value:
        val['label'] = spm_date_dict[val['label']]
        final_pixel_values.append(val)
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
    # add_box_plots(ax, spm_data.values())
    p1, = twin1.plot(new_y_daily_values_cont[0].keys(),
                     new_y_daily_values_cont[0].values(), markersize=0,
                          linestyle = '-', color="#d4a373", label=y1_field)

    p2, = twin1.plot(new_y_daily_values_cont[1].keys(),
                     new_y_daily_values_cont[1].values(), markersize=0,
                          linestyle = '-',color="#bc4b51", label=y2_field)

    p3, = twin2.plot(new_y_daily_values_cont[2].keys(),
                     new_y_daily_values_cont[2].values(), markersize=0,
                          linestyle = '-',color="#0ead69", label=y3_field)

    values_y3_data = []
    for j in new_y_daily_values_cont[3].values():
        values_y3_data.append(float(j))
    p4 = twin2.bar(list(new_y_daily_values_cont[3].keys()),
                   values_y3_data, color= '#7cb518', label=y4_field)
    # twin2.set_ylim(twin2.get_ylim())

    plt.gcf().set_size_inches(20, 12)
    plt.gca().set_xbound(left, right)
    ax.set_xlabel(y0_date_field)
    ax.set_ylabel(y0_label)

    ax.set_ylim(10, max(y0_daily_values.values())*1.5)
    twin1.set_ylim(0, max(new_y_daily_values_cont[1].values())*1.5)
    twin2.set_ylim(1,  max(values_y3_data)*1.5)

    twin1.invert_yaxis()
    twin2.invert_yaxis()

    twin1.set_ylabel(y1_label)
    twin2.set_ylabel(y2_label)
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

    plt.gcf().autofmt_xdate()
    plt.tight_layout()
    ax.legend(handles=[p, p0, p1, p2, p3, p4], loc=4)
    save_plot(plt, 'poster', 'SPM and Discharge')

SPM_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\sites_SPM.xlsx'
y_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\discharge_daily_combined.xlsx'

plot_data(SPM_path, y_path, 'Insitu SPM', 'UAS SPM', 'Date', 'Jourdan River', 'Wolf River','Pearl River','Bonnet Carre Spillway', 'SPM (mg/L)', 'Gauge Height','Discharge (cpfs)','Date', 4)
