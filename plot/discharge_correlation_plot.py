import csv
import os.path
from collections import OrderedDict
from datetime import datetime, timedelta

from fontTools.misc.cython import returns
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
import pandas as pd
from plot import save_plot, read_xy_data, group_data_by_date


def prepare_SPM_discharge_corr_data(x_path, y_path, x_field_name, y_field_name,
                                    x_date_field, y_date_field,
                                    residence_time=None):
    x_df, y_df = read_xy_data(x_path, y_path)
    # print (x_date_field,
    #        x_field_name, y_date_fie-ld,
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


def create_scatter_plot(x, y, color, x_field_name, y_field_name, transt_d, cor='nan'):
    # xMap = assignIDs(x)
    # xAsInts = np.array([xMap[i] for i in x])
    x_np_arr = np.array(x)
    pearR = np.corrcoef(x, y)[1, 0]
    # least squares from:
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
    A = np.vstack([x, np.ones(len(x))]).T
    try:
        m, c = np.linalg.lstsq(A, np.array(y), rcond=None)[0]
    except Exception as ex:
        return
    plt.scatter(x, y, color=color)
    plt.plot(x_np_arr, x_np_arr * m + c, color=color,
             label="r = {}".format(round(pearR, 2)))
    # plt.xticks(xMap.values(),xMap.keys())
    plt.legend(loc=1)
    plt.xlabel(x_field_name)
    plt.ylabel(y_field_name)
    plt.title('Fit of {} and {}'.format(x_field_name, y_field_name))
    # plt.show()
    save_plot(plt, 'poster', '{} {} {} {}'.format(x_field_name, y_field_name, transt_d, cor[0]))


def transit_time_correlation_plot(correlation, param_type, x_label_name, y_label_name, riv):
    for k, v in list(correlation.items()):

        if len(v) > 0:
            if v[0] < 0:
                del correlation[k]
        else:
            del correlation[k]

    x = np.array(list(correlation.keys()))
    y = list(chain.from_iterable(list(correlation.values())))

    plt.rcParams['figure.figsize'] = [6, 3] # [width, height]
    fig, ax = plt.subplots()
    ax.set_xticks(list(range(0, 110, 10)))
    plt.plot(x, y)
    plt.xlabel("{} (Days)".format(x_label_name)) # "Transit time (Days)"
    plt.ylabel("{} (r)".format(y_label_name)) # "correlation coefficient (r)"
    plt.title("Transit Time Vs Correlation, {}".format(riv))

    plt.ylim([0, 1])
    save_plot(plt, 'poster', 'transit_time_vs_correlation {}_{}'.format(param_type, riv))
    return y

def transit_time_correlation_boxplots(correlation, param_type):
    for c in (correlation.keys()):
        correlation_list = correlation[c]
        del_index = []
        for r in range(len(correlation_list)):
            if correlation_list[r] < 0:
                del_index.append(r)
        for ele in sorted(del_index, reverse=True):
            del correlation[c][ele]

    plt.ylim([0, 1])
    plt.boxplot(correlation.values(), labels=correlation.keys(), showfliers=False)

    save_plot(plt, 'poster', 'transit_time_vs_correlation_{}_boxplots'.format(param_type))

SPM_path = r'C:\Users\andex\OneDrive\Documents\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\sites_SPM_daily.xlsx'
discharge_path = r'C:\Users\andex\OneDrive\Documents\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\discharge_daily_combined.xlsx'
pearl_river_path = r'C:\Users\andex\OneDrive\Documents\MSU\dissertation\SPM_Multi-spectral\data\discharge\daily\Pearl.xlsx'

def loop_correlation(river, SPM_type, transit_time, correlations, x_y):
    x, y = prepare_SPM_discharge_corr_data(
        SPM_path, discharge_path, '{} SPM'.format(SPM_type), river, 'Date',
        'Date', transit_time
    )

    if transit_time not in correlations.keys():
        correlations[transit_time] = []
    corr = np.corrcoef(x, y)[0][1]

    if corr > 0:
        correlations[transit_time].append(corr)
        x_y[transit_time] = [x, y]


def create_transit_time_plots(river, SPM_type):
    correlations = {}
    x_y = {}
    for i in range(1, 100):
        loop_correlation(river, SPM_type, i, correlations, x_y)

    if len(correlations.keys()) > 0:

        best_corr =  max(correlations.values())
        transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
        print (river, transit_time, best_corr)
        create_scatter_plot(x_y[transit_time][0], x_y[transit_time][1],
                            'blue', SPM_type, river, transit_time)
        correlation = transit_time_correlation_plot(
            correlations, SPM_type, 'Transit Time', 'Correlation Coefficient', river)
        river_correlations[river] = correlation

river_correlations = {}
create_transit_time_plots('Wolf River', 'Insitu')
create_transit_time_plots('Bonnet Carre Spillway', 'Insitu')
create_transit_time_plots('Pearl River', 'Insitu')
create_transit_time_plots('Jourdan River', 'Insitu')
transit_time_correlation_boxplots(river_correlations, 'Insitu')

create_transit_time_plots('Wolf River', 'UAS')
create_transit_time_plots('Bonnet Carre Spillway', 'UAS')
create_transit_time_plots('Pearl River', 'UAS')
create_transit_time_plots('Jourdan River', 'UAS')
transit_time_correlation_boxplots(river_correlations, 'UAS')
