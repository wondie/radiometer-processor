import csv
from collections import OrderedDict
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from plot import save_plot, read_xy_data, group_data_by_date


def prepare_SPM_discharge_corr_data(x_path, y_path, x_field_name, y_field_name,
                                    x_date_field, y_date_field,
                                    residence_time=None):
    x_df, y_df = read_xy_data(x_path, y_path)
    # print (x_date_field,
    #        x_field_name, y_date_field,
    #        [y_field_name])
    date_objs, x_daily_values, y_daily_values_cont, insitu_spm_table = group_data_by_date(x_df, y_df, x_date_field,
                                                 x_field_name, y_date_field,
                                                 [y_field_name])
    y_daily_values = y_daily_values_cont[0]
    # create final x and y data with available dates for SPM
    # with the consideration of residence time
    x_plot_data = []
    y_plot_data = []
    for curr_date, x_data in x_daily_values.items():
        curr_date_obj = datetime.strptime(curr_date, '%m-%d-%Y')
        # subtract the residence day from the curr_date of SPM /x to find matching
        # discharge data that passed through the gauging station
        corrected_date = curr_date_obj - timedelta(days=residence_time)
        corrected_date_str = corrected_date.strftime("%m-%d-%Y")
        if corrected_date_str in y_daily_values.keys():
            if y_daily_values[corrected_date_str] is not None and y_daily_values[corrected_date_str] != 0:
                x_plot_data.append(x_data)
                y_plot_data.append(y_daily_values[corrected_date_str])

    return x_plot_data, y_plot_data




def assignIDs(list):
    '''Take a list of strings, and for each unique value assign a number.
    Returns a map for "unique-val"->id.
    '''
    sortedList = sorted(list)

    # taken from
    # http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order/480227#480227
    seen = set()
    seen_add = seen.add
    uniqueList = [x for x in sortedList if x not in seen and not seen_add(x)]

    return dict(zip(uniqueList, range(len(uniqueList))))


def create_scatter_plot(x, y, color, x_field_name, y_field_name):
    # xMap = assignIDs(x)
    # xAsInts = np.array([xMap[i] for i in x])
    x_np_arr = np.array(x)
    pearR = np.corrcoef(x, y)[1, 0]
    # least squares from:
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, np.array(y), rcond=None)[0]

    plt.scatter(x, y, color=color)
    plt.plot(x_np_arr, x_np_arr * m + c, color=color,
             label="r = {}".format(round(pearR, 2)))
    # plt.xticks(xMap.values(),xMap.keys())
    plt.legend(loc=1)
    plt.xlabel(x_field_name)
    plt.ylabel(y_field_name)
    plt.title('Fit of {} and {}'.format(x_field_name, y_field_name))
    # plt.show()
    save_plot(plt, 'poster', '{} {}'.format(x_field_name, y_field_name))


SPM_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\sites_SPM_daily.xlsx'
discharge_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\discharge_daily_combined.xlsx'
pearl_river_path = r'D:\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\Pearl.xlsx'
x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
                                       'Jourdan River', 'Date', 'Date', 16)
create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Jourdan River')

x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
                                       'Wolf River', 'Date', 'Date', 13)
create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Wolf River')
x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
                                       'Pearl River', 'Date', 'Date', 39)
create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Pearl River')

x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
                                       'Bonnet Carre Spillway', 'Date', 'Date',
                                       17)
create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Bonnet Carre Spillway')




x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
                                       'Jourdan River', 'Date', 'Date', 16)
create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Jourdan River')

x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
                                       'Wolf River', 'Date', 'Date', 13)
create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Wolf River')
x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
                                       'Pearl River', 'Date', 'Date', 39)
create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Pearl River')

x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
                                       'Bonnet Carre Spillway', 'Date', 'Date', 17)
create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Bonnet Carre Spillway')