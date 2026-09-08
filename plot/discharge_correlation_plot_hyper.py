import csv
import os.path
from collections import OrderedDict
from datetime import datetime, timedelta
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
import pandas as pd
from plot import save_plot, read_xy_data, group_data_by_date, plot_size


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
        # print(v)
        if len(v) > 0:
            if v[0] < 0:
                del correlation[k]
        else:
            del correlation[k]
    print (correlation.values())
    x = np.array(list(correlation.keys()))
    y = list(chain.from_iterable(list(correlation.values())))
    # print (y_intermid)
    # y = [y[0] for y in y_intermid]
    #print(x)
    plt.rcParams['figure.figsize'] = [6, 3] # [width, height]
    fig, ax = plt.subplots()
    ax.set_xticks(list(range(0, 110, 10)))
    plt.plot(x, y)
    plt.xlabel("{} (Days)".format(x_label_name)) # "Transit time (Days)"
    plt.ylabel("{} (r)".format(y_label_name)) # "correlation coefficient (r)"
    plt.title("Transit Time Vs Correlation, {}".format(riv))
    #plt.show()
    plt.ylim([0, 1])
    save_plot(plt, 'poster', 'transit_time_vs_correlation {}_{}'.format(param_type, riv))
    return y

def transit_time_correlation_all_rivers(correlations_dict, param_type, x_label_name, y_label_name):
    plt.figure(figsize=(10, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for idx, (river_name, correlation) in enumerate(correlations_dict.items()):
        filtered_corr = {}
        for k, v in list(correlation.items()):
            if len(v) > 0 and v[0] >= 0:
                filtered_corr[k] = v
        
        if len(filtered_corr) > 0:
            x = np.array(list(filtered_corr.keys()))
            y = list(chain.from_iterable(list(filtered_corr.values())))
            plt.plot(x, y, label=river_name, color=colors[idx % len(colors)], linewidth=2)
    
    plt.xlabel("{} (Days)".format(x_label_name))
    plt.ylabel("{} (r)".format(y_label_name))
    plt.title("Transit Time Vs Correlation - {} SPM (All Rivers)".format(param_type))
    plt.ylim([0, 1])
    plt.xlim([0, 100])
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    save_plot(plt, 'poster', 'transit_time_vs_correlation_{}_rivers'.format(param_type))

def transit_time_correlation_boxplots(correlation, param_type):
    for c in (correlation.keys()):
        correlation_list = correlation[c]
        del_index = []
        for r in range(len(correlation_list)):
            if correlation_list[r] < 0:
                del_index.append(r)
        for ele in sorted(del_index, reverse=True):
            del correlation[c][ele]
    plt.figure(figsize=(20, 10))
    plt.rcParams.update({
        'font.size': 28, 
        'font.weight': 'bold',
        'xtick.labelsize': 28, 
        'ytick.labelsize': 28,
        'axes.labelweight': 'bold',
        'axes.linewidth': 2.5
    })
    plt.ylim([0, 1])
    box = plt.boxplot(correlation.values(), labels=correlation.keys(), showfliers=False,
                      widths=0.6,
                      patch_artist=False,
                      boxprops=dict(linewidth=2.5),
                      whiskerprops=dict(linewidth=2.5),
                      capprops=dict(linewidth=2.5),
                      medianprops=dict(linewidth=3.0, color='#ff7f0e'))
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    plt.tight_layout(pad=2.0)
    plt.savefig(r'G:\Other computers\My Laptop\codes\radiometer_processor\data\output\transit_time_vs_correlation_{}_boxplots.png'.format(param_type), dpi=150, bbox_inches='tight')
    plt.close()


def transit_time_correlation_boxplots_combined(insitu_correlations, uas_correlations,
                                                river_order=('Bonnet Carre Spillway', 'Jourdan River',
                                                             'Pearl River', 'Wolf River')):
    """Combine the Insitu (panel A) and UAS (panel B) transit-time correlation
    boxplots into a single two-panel figure, for the case where the two aren't
    significantly different and don't need separate figures.
    """
    def clean(correlations):
        cleaned = OrderedDict()
        for river in river_order:
            values = [v for v in correlations.get(river, []) if v is not None and v >= 0]
            cleaned[river] = values
        return cleaned

    labels = [r.replace('Bonnet Carre Spillway', 'Bonnet Carre\nSpillway') for r in river_order]

    fig, axes = plt.subplots(1, 2, figsize=(15, 3.5), sharey=True)
    boxprops = dict(linewidth=1.5)
    whiskerprops = dict(linewidth=1.5)
    capprops = dict(linewidth=1.5)
    medianprops = dict(linewidth=2.0, color='#ff7f0e')

    for ax, correlations, panel_label in zip(axes, (clean(insitu_correlations), clean(uas_correlations)), ('A)', 'B)')):
        ax.boxplot(
            list(correlations.values()), labels=labels, showfliers=False,
            widths=0.6, patch_artist=False,
            boxprops=boxprops, whiskerprops=whiskerprops,
            capprops=capprops, medianprops=medianprops,
        )
        ax.set_ylim(0, 1)
        ax.set_xlabel('Rivers')
        ax.tick_params(axis='x', labelsize=10)
        ax.text(0.95, 0.95, panel_label, transform=ax.transAxes,
                ha='right', va='top', fontweight='bold')

    axes[0].set_ylabel('Correlation (r)')

    plt.tight_layout()
    fig.subplots_adjust(wspace=0.15)
    output_path = r'G:\Other computers\My Laptop\codes\radiometer_processor\data\output\transit_time_vs_correlation_combined_boxplots.svg'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    return output_path


root = r'G:\Other computers\My Laptop\dissertation\SPM_Hyperspectral'
SPM_path = r'{}\data\discharge\daily\sites_SPM_daily.xlsx'.format(root)
discharge_path = r'{}\data\discharge\daily\discharge_daily_combined.xlsx'.format(root)
pearl_river_path = r'{}\data\discharge\daily\Pearl.xlsx'.format(root)
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
#                                        'Jourdan River', 'Date', 'Date', 16)
# create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Jourdan River')
#
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
#                                        'Wolf River', 'Date', 'Date', 13)
# create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Wolf River')
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
#                                        'Pearl River', 'Date', 'Date', 39)
# create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Pearl River')
#
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'Insitu SPM',
#                                        'Bonnet Carre Spillway', 'Date', 'Date',
#                                        17)
# create_scatter_plot(x, y, 'blue', 'Insitu SPM', 'Bonnet Carre Spillway')

# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
#                                        'Jourdan River', 'Date', 'Date', 16)
# create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Jourdan River')
#
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
#                                        'Wolf River', 'Date', 'Date', 13)
# create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Wolf River')
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
#                                        'Pearl River', 'Date', 'Date', 39)
# create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Pearl River')
#
# x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, 'UAS SPM',
#                                        'Bonnet Carre Spillway', 'Date', 'Date', 17)
# create_scatter_plot(x, y, 'blue', 'UAS SPM', 'Bonnet Carre Spillway')

def loop_correlation(river, SPM_type, transit_time, correlations):
    # file_path = r'D:\MSU\codes\radiometer_processor\data\{} {} {}.png'.format(SPM_type, river, transit_time)
    # if os.path.isfile(file_path):
    #     return
    x, y = prepare_SPM_discharge_corr_data(SPM_path, discharge_path, '{} SPM'.format(SPM_type),
                                           river, 'Date', 'Date', transit_time)

    if transit_time not in correlations.keys():
        correlations[transit_time] = []
    corr = np.corrcoef(x, y)[0][1]
    print(corr)
    if corr > 0:
        create_scatter_plot(x, y, 'blue', SPM_type, river, transit_time)
        correlations[transit_time].append(corr)


correlations = {}
for i in range(1, 100):
    loop_correlation('Jourdan River', 'Insitu', i, correlations)

print ('Insitu')
param_type = 'Insitu'

print ('Corr ', correlations)
correlations_jourdan = correlations.copy()
if len(correlations.keys()) > 0:

    best_corr =  max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Jourdan River', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Jourdan River')
    correlation_uas_jourdan = correlation #list(chain.from_iterable(list(correlations.values())))


correlations = {}
for i in range(1, 100):
    loop_correlation('Wolf River', 'Insitu', i, correlations)

print (correlations)
correlations_wolf = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr =  max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Wolf River', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Wolf River')
    correlation_uas_wolf = correlation #list(chain.from_iterable(list(correlations.values())))


correlations = {}
for i in range(1, 100):
    loop_correlation('Pearl River', 'Insitu', i, correlations)
correlations_pearl = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr = max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Pearl River', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Pearl River')
    correlation_uas_pearl = correlation #list(chain.from_iterable(list(correlations.values())))

correlations = {}
for i in range(1, 100):
    loop_correlation('Bonnet Carre Spillway', 'Insitu', i, correlations)
correlations_spillway = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr = max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Bonnet Carre Spillway', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Bonnet Carre Spillway')
    correlation_uas_spillway = correlation #list(chain.from_iterable(list(correlations.values())))

# Captured (not plotted standalone) so it can be combined with UAS below into
# a single two-panel A)/B) figure instead of two separate figures.
insitu_river_correlations = {'Jourdan River': correlation_uas_jourdan, 'Wolf River': correlation_uas_wolf,
                             'Pearl River': correlation_uas_pearl, 'Bonnet Carre Spillway': correlation_uas_spillway}


correlations = {}
for i in range(1, 100):
    loop_correlation('Jourdan River', 'UAS', i, correlations)

print ('UAS')
param_type = 'UAS'
correlations_jourdan = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr = max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Jourdan River', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Jourdan River')
    correlation_uas_jourdan = correlation #list(chain.from_iterable(list(correlations.values())))

correlations = {}
for i in range(1, 100):
    loop_correlation('Wolf River', 'UAS', i, correlations)

correlations_wolf = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr = max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Wolf River', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Wolf River')
    correlation_uas_wolf = correlation #list(chain.from_iterable(list(correlations.values())))

correlations = {}
for i in range(1, 100):
    loop_correlation('Pearl River', 'UAS', i, correlations)

correlations_pearl = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr = max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Pearl River', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Pearl River')
    correlation_uas_pearl = correlation #list(chain.from_iterable(list(correlations.values())))

correlations = {}
for i in range(1, 100):
    loop_correlation('Bonnet Carre Spillway', 'UAS', i, correlations)

correlations_spillway = correlations.copy()
if len(correlations.keys()) > 0:
    best_corr = max(correlations.values())
    transit_time = list(correlations.keys())[list(correlations.values()).index(best_corr)]
    print ('Bonnet Carre Spillway', transit_time, best_corr)
    correlation = transit_time_correlation_plot(correlations, param_type, 'Transit Time', 'Correlation Coefficient', 'Bonnet Carre Spillway')
    correlation_uas_spillway = correlation #list(chain.from_iterable(list(correlations.values())))

uas_river_correlations = {'Jourdan River': correlation_uas_jourdan, 'Wolf River': correlation_uas_wolf,
                          'Pearl River': correlation_uas_pearl, 'Bonnet Carre Spillway': correlation_uas_spillway}

combined_output_path = transit_time_correlation_boxplots_combined(insitu_river_correlations, uas_river_correlations)
print('Combined Insitu/UAS boxplot saved to:', combined_output_path)