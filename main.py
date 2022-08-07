###Copyright (C) Wondimagegn Tesfaye Beshah - All Rights Reserved
# Unauthorized copying of this file and software as a whole,
# via any medium is strictly prohibited
# Proprietary and confidential
# Written by Wondimagegn Tesfaye Beshah <wondim81@gmail.com>, September 2021
###
import glob
import json
import os
# import re
from functools import reduce
from itertools import groupby
# import threading
from collections import OrderedDict
import math
# from datetime import datetime
from os import remove, path
from os.path import expanduser
# import random

from natsort import os_sorted
# from sklearn.metrics import mean_squared_error
import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (
    QMainWindow, QApplication,
    QFileDialog, QTableWidgetItem,
    QComboBox, QDialogButtonBox, QHeaderView, QLineEdit, QWidget, QHBoxLayout,
    QToolButton, QTableWidget, QVBoxLayout, QSizePolicy, QSpacerItem, QPlainTextEdit, QMessageBox, QPushButton,
    QTreeWidgetItem)
from ui.processor import Ui_RrsProcessor
data_path = 'D:/MSU/codes/radiometer_processor/sampledata/processed/'
PREFIX = 'water'
def prepare_irradiance(path):
    """
    [Use calculate reflectance] Loop through site irradiance observation asc files and create a csv
    file containing all the observation and the average. You will need to see
    :param path: The path containing all raw radiometer data
    :type path: String
    :return:
    :rtype:
    """
    irr_files = glob.glob('{}*{}*.asc'.format(path, '[iI][rR][rR]'))
    irr_files_sorted = sorted(irr_files)
    data0 = None
    irr_cols = []
    for file_path in irr_files_sorted:
        suff = file_path[file_path.rindex('_') + 1:].replace('.asc', '')
        irr_col = 'raw_irr_{}'.format(suff)
        irr_cols.append(irr_col)
        headers = ['wavelength', 'reference', irr_col]
        if suff == '000':
            data0 = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                               names=headers)
        else:
            data123 = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                                names=headers)
            if data0 is not None:
                data0[irr_col] = data123[irr_col]

    cols = []
    for col in irr_cols:
        new_col = col.replace('raw_', '')
        data0[new_col] = data0[col]/(10**6)
        cols.append(new_col)

    data0['Avg_Irr'] = data0[cols].mean(1)
    output_irr = irr_files_sorted[0].replace('_000', '').replace('.asc', '.xlsx')
    create_xlsx_chart(data0, output_irr)


def detect_outlier(data_1):
    if isinstance(data_1, OrderedDict) or isinstance(data_1, dict):
        data = list(data_1.values())
    else:
        data = data_1

    outliers = OrderedDict()
    threshold = 2.0e-05
    # print (data)
    mean_1 = np.mean(data)
    std_1 = np.std(data)

    for key, y in data_1.items():
        z_score = (y - mean_1) / std_1
        if np.abs(z_score) > threshold:
            outliers[key] = y
    return outliers

def calculate_reflectance(path):
    """
    Loop through site water observation asc files and create a csv
    file containing all the observation and the average. You will
    need to see the chart to choose the best spectra.
    :param path: The path containing all raw radiometer data
    :type path: String
    :return:
    :rtype:
    """
    water_files_sorted = glob.glob('{}/*.asc'.format(path))
    water_files_sorted = os_sorted(water_files_sorted)
    # print ('w',water_files_sorted)
    grouped_files = [(j,list(i)) for j, i in groupby(water_files_sorted,
                  lambda a: '{}/{}_{}'.format(
                      os.path.dirname(a),
                      os.path.basename(a).split('_')[0],
                      os.path.basename(a).split('_')[1])
                                                     )]
    # print ('g',grouped_files)
    processed_files_cols = OrderedDict()
    processed_files_df = OrderedDict()
    # print (path)

    for output_path, water_files in grouped_files:
        df = None
        for file_path in water_files:
            # print ('f',file_path)
            suff = file_path[file_path.rindex('_') + 1:].replace('.asc', '')
            water_col = 'raw_{}_{}'.format(PREFIX, suff)

            headers = ['wavelength', 'reference', water_col]

            if suff == '000':

                df = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                                   names=headers)
            elif suff == '003':
                # water_cols[file_path].append(water_col)
                data12 = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                                     names=headers)
                if df is not None:
                    df['sky'] = data12[water_col]

            else:# suffix is 001, or 002
                # water_cols[file_path.replace('_{}'.format(suff), '')].append(water_col)
                data12 = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                                    names=headers)
                if df is not None:
                    df[water_col] = data12[water_col]

        cols_2 = []
        for col in df.columns:
            if col.startswith('raw'):
                new_col = col.replace('raw_', '')

                df[new_col] =  (df[col] - (0.02*df['sky']))/((100/99)*df['reference']*math.pi)

                cols_2.append(new_col)

        if df is None:
            print ('Unable to find water files in ', path)
            return

        processed_files_cols[output_path] = [c for c in list(df.columns)
                                             if c.startswith(PREFIX)]
        processed_files_df[output_path] = df
        # create_xlsx_chart(df, output_path)
    return processed_files_cols, processed_files_df

def process_reflectance(root_path):
    """
    Loop though each site folder containing radiometer raw data.
    :param root_path: The root folder
    :type root_path: String
    :return:
    :rtype:
    """
    dirs = os.listdir(root_path)
    for dr in dirs:
        path = os.path.join(root_path, dr)
        calculate_reflectance(path)

def create_xlsx_chart(df, output_xlsx, cols=None):
    """
    Creates an Excel file with Sheet1 containing data and Chart sheet.
    :param df: The dataframe containing the data
    :type df: pd.Dataframe
    :param output_xlsx: The output Excel file path
    :type output_xlsx: String
    :param chart_start: The start column of the chart
    :type chart_start: Integer
    :return:
    :rtype:
    """
    writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1')
    if cols is not None:
        create_cols_chart(df, writer, cols, len(df.columns) + 1)
    else:
        chart_start = len(df.columns) - 3

        create_chart(df, writer, chart_start, len(df.columns) + 1)
    writer.save()
    # print ('Written: ', output_xlsx)

def create_cols_chart(df, writer, cols, x_label=1):
    """
    Creates a chart from DataFrame.
    :param df: The dataframe containing the data.
    :type df: pd>DataFrame
    :param writer: The Excel writer
    :type writer: pd.ExcelWriter
    :param cols: List of columns to be on the Y axis
    :type start: List
    :param x_label: The label column
    :type x_label: Integer
    :return:
    :rtype:
    """
    workbook = writer.book
    chart = workbook.add_chart({'type': 'line'})
    for col in cols:
        i = df.columns.get_loc(col)
        chart.add_series({
            'name': ['Sheet1', 0, i],
            'values': ['Sheet1', 1, i, len(df.index), i],
            'categories': ['Sheet1', 1, x_label, len(df.index), x_label]
        })
    chart.set_x_axis({'name': 'wavelength'})
    chart.set_y_axis({'name': 'Rrs'})
    chart_sheet = workbook.add_chartsheet()
    # Configure the chart.
    chart_sheet.set_chart(chart)

def create_chart(df, writer, start, end, x_label=1):
    """
    Creates a chart from DataFrame.
    :param df: The dataframe containing the data.
    :type df: pd>DataFrame
    :param writer: The Excel writer
    :type writer: pd.ExcelWriter
    :param start: The start column of the chart
    :type start: Integer
    :param end: The end column of the chart
    :type end: Integer
    :param x_label: The label column
    :type x_label: Integer
    :return:
    :rtype:
    """
    workbook = writer.book
    chart = workbook.add_chart({'type': 'line'})
    for i in range(start, end):
        chart.add_series({
            'name': ['Sheet1', 0, i],
            'values': ['Sheet1', 1, i, len(df.index), i],
            'categories': ['Sheet1', 1, x_label, len(df.index), x_label]
        })
    chart.set_x_axis({'name': 'wavelength'})
    chart.set_y_axis({'name': 'Rrs'})
    chart_sheet = workbook.add_chartsheet()
    # Configure the chart.
    chart_sheet.set_chart(chart)

# After creating average using good observation if there is a bad data,
# copy all to one folder in this case Site. The are required to be xlsx file contaning WMS and Water
# all_rrs = 'D:/MSU/codes/radiometer_processor/sampledata/2020/sites'
# all_rrs = r'D:\MSU\SPM Testing\WMS_2020_2021\WMS_2020_2021'
def collect_sites_rrs_to_single_sheet(path, exclusions=[]):
    """
    Creates site RRS from different sites into single sheet
    picking the average from each site files.
    :param path: The path containing all the files of different sites xlsx.
    :type path: String
    :param exclusions: List of files to exclude
    :type exclusions: List
    :return:
    :rtype:
    """
    water_files = glob.glob('{}/{}'.format(
        path, '[wW][mM][sS]*[wW][aA][tT][eE][rR]*.xlsx'
    ))
    output_xlsx = os.path.join(path, 'WMS_Rrs.xlsx')
    # print (water_files)
    df_cols = OrderedDict()
    writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')
    for i, f in enumerate(water_files):
        file_name = os.path.basename(f)
        col_names = file_name.split('_')
        if col_names[1] in exclusions:
            continue
        col_name = '{}_{}'.format(col_names[0], col_names[1])
        # print (col_name)
        df = pd.read_excel(f, sheet_name='Sheet1')
        if i == 0:
            df_cols['wavelength'] = (df['wavelength'])
        df_cols[col_name] = df['avg_Rrs']

    new_df = pd.concat(df_cols.values(), axis=1, keys=df_cols.keys())
    new_df.to_excel(writer, sheet_name='Sheet1')
    create_chart(new_df, writer, 2, len(new_df.columns) + 1)
    writer.save()

def collect_rrs_to_single_sheet(processed_df):
    """
    Creates site RRS from different sites into single sheet
    picking the average from each site files.
    :param water_files: A list containing all the files of different sites xlsx.
    :type path: List
    :return:
    :rtype:
    """
    # water_files = glob.glob('{}/{}'.format(
    #     path, '[wW][mM][sS]*[wW][aA][tT][eE][rR]*.xlsx'
    # ))
    if (len(processed_df.keys())) == 0:
        return
    first_path = list(processed_df.keys())[0]
    path = os.path.dirname(first_path)
    file_name = os.path.basename(first_path)
    site = file_name.split('_')[0]
    output_xlsx = os.path.join(path, '{}_Rrs.xlsx'.format(site))
    # print (water_files)
    df_cols = OrderedDict()
    writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')
    for i, (f, df)  in enumerate(processed_df.items()):
        file_name = os.path.basename(f)
        col_names = file_name.split('_')

        col_name = '{}_{}'.format(col_names[0], col_names[1].replace('.xlsx', ''))

        if i == 0:
            df_cols['wavelength'] = (df['wavelength'])
        df_cols[col_name] = df['avg_Rrs']

    new_df = pd.concat(df_cols.values(), axis=1, keys=df_cols.keys())
    new_df.to_excel(writer, sheet_name='Sheet1')
    create_chart(new_df, writer, 2, len(new_df.columns) + 1)
    writer.save()
    return output_xlsx

# collect_sites_rrs_to_single_sheet(all_rrs, exclusions=['57', '58', '63'])
# process_reflectance('D:/MSU/codes/radiometer_processor/sampledata/2020/WMS')
# rrs = 'D:/MSU/codes/radiometer_processor/sampledata/2020/sites/WMS_Rrs.xlsx'

def interpolate(xlsx_file_path, min, max):
    """
    Interpolate the site Rrs data into one step increment
    :param xlsx_file_path: The file bath of all sites combined Rrs xlsx file
    :type xlsx_file_path: String
    :param min: The minimum wavelength value
    :type min: Integer
    :param max: The maximum wavelength value
    :type max: Integer
    :return:
    :rtype:
    """
    df = pd.read_excel(xlsx_file_path, sheet_name='Sheet1', index_col=None)
    for i in range(min, max+1):
        df = df.append({'wavelength': i}, ignore_index=True)
    df.sort_values(by=['wavelength'], inplace=True)

    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df = df.interpolate()
    y = df[df['wavelength'].isin(range(min+1, max+1))]

    output_xlsx = xlsx_file_path.replace('.xlsx', '_interpolate.xlsx')

    writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')
    y.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(y, writer, 1, len(y.columns), 0)
    writer.save()

#interpolate(rrs, 277, 1094)
def merge_two_spreadsheet(xlsx_path1, xlsx_path2, combined_path):
    """
    Merges two spreadsheets.
    :param xlsx_path1: xlsx one path
    :type xlsx_path1: String
    :param xlsx_path2: xlsx two path
    :type xlsx_path2: String
    :param combined_path: output xlsx path
    :type combined_path: String
    :return:
    :rtype:
    """

    df1 = pd.read_excel(xlsx_path1, sheet_name='Sheet1', index_col=None)
    df2 = pd.read_excel(xlsx_path2, sheet_name='Sheet1', index_col=None)

    df3 = pd.merge(df1, df2, on="wavelength")
    writer = pd.ExcelWriter(combined_path, engine='xlsxwriter')
    df3.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(df3, writer, 1, len(df3.columns), 0)
    writer.save()


def merge_spreadsheet(list_of_files, combined_path):
    """
    Merges two spreadsheets.
    :param xlsx_path1: xlsx one path
    :type xlsx_path1: String
    :param xlsx_path2: xlsx two path
    :type xlsx_path2: String
    :param combined_path: output xlsx path
    :type combined_path: String
    :return:
    :rtype:
    """
    # print (list_of_files)
    # df0 = pd.read_excel(list_of_files[0], sheet_name='Sheet1', index_col=None, engine='openpyxl')
    # df2 = pd.read_excel(list_of_files[1], sheet_name='Sheet1', index_col=None, engine='openpyxl')
    # df1 = pd.merge(df0, df2, on="wavelength")
    dfs = []
    for  xlsx_path in list_of_files:
        df = pd.read_excel(xlsx_path)
        # xl = pd.ExcelFile(xlsx_path)
        # df = xl.parse("Sheet1")
        # # import openpyxl
        #
        # wb = openpyxl.load_workbook("/Users/Scott/Desktop/Workbook1.xlsx")
        dfs.append(df)
    # print (dfs)
    df_merged = reduce(lambda left, right: pd.merge(left, right, on=['wavelength'],
                                                    how='inner'), dfs)

    #     if i== 0:
    #         continue
    #     df2 = pd.read_excel(xlsx_path, sheet_name='Sheet1', index_col=None, engine='openpyxl')
    #     # df2 = pd.read_excel(xlsx_path2, sheet_name='Sheet1', index_col=None, engine='openpyxl')
    #     print ('Merge')
    #     df1 = pd.merge(df1, df2, on="wavelength")
        # df3 = df1

    writer = pd.ExcelWriter(combined_path, engine='xlsxwriter')
    # df3.to_excel(writer, sheet_name='Sheet1', index=False)
    df_merged.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(df_merged, writer, 1, len(df_merged.columns), 0)
    writer.save()


# merge_spreadsheet(
#     'D:/MSU/codes/radiometer_processor/sampledata/2020/Radiometer_rrs.xlsx',
#     'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_interpolate.xlsx',
#   'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_combined.xlsx'
# )


def convert_hyperspectral_to_multispectral(srf_path, rrs_path, output_rrs):
    """
    Applies a spatial response function of multi-spectral sensors
    on radiometer reflectance values to come up with multi-spectral
    reflectance values.
    :param srf_csv: srf_csv is a csv file containing the following columns:
    wavelength, response value for each subsequent bands from the start
    to end of the wavelengths.
    :type srf_csv: String
    :param rrs_csv: contains a csv file with reflectance value for each sites.
    The columns are wavelength(should match srf_csv wavelength),
    sites reflectance values for each wavelength.
    :type rrs_csv: String
    :param output_rrs: The output rrs path
    :type output_rrs: String
    :return:
    :rtype:
    """
    if srf_path.endswith('.csv'):
        srf = pd.read_csv(srf_path)
    elif srf_path.endswith(('.xls', '.xlsx')):
        srf = pd.read_excel(srf_path)
    else:
        return
    if rrs_path.endswith('.csv'):
        rrs = pd.read_csv(rrs_path)
    elif rrs_path.endswith(('.xls', '.xlsx')):
        rrs = pd.read_excel(rrs_path)
    else:
        return
    try:
        new_df = pd.merge(rrs, srf, on='wavelength')
    except Exception as ex:
        print (ex)
        return
    final_cols = ['sites']
    final_cols.extend([c for c in srf.columns if c != 'wavelength'])

    final_dic = {}
    final_dic['sites'] = [c for c in rrs.columns if c != 'wavelength']

    for col_rrs in rrs:
        if col_rrs == 'wavelength':
            continue

        for col_srf in srf:
            if col_srf == 'wavelength':
                continue

            p1 = new_df[col_srf] * new_df[col_rrs]
            if col_srf in final_dic.keys():
                # add the first calculated band
                final_dic[col_srf].append(sum(p1) / sum(new_df[col_srf]))
            else:
                final_dic[col_srf] = [sum(p1) / sum(new_df[col_srf])]

    final_df = pd.DataFrame(final_dic, columns=final_cols)
    writer = pd.ExcelWriter(output_rrs, engine='xlsxwriter')
    final_df.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(final_df, writer, 1, len(final_df.columns), 0)
    writer.save()

# spectral_response_function_file = 'D:/MSU/codes/radiometer_processor/sampledata/Spectral_response_function_micasense_277_1094.csv'
# site_rrs_path = 'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_interpolate.csv'
# output = 'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense.xlsx'
# convert_hyperspectral_to_multispectral(spectral_response_function_file, site_rrs_path, output)

spectral_response_function_file = r'D:\MSU\SPM_MODIS\Spectral_response_function_modis.csv'
site_rrs_path = r'D:\MSU\SPM_MODIS\WMS_Rrs_combined.csv'
output = r'D:\MSU\SPM_MODIS\WMS_Rrs_modis.xlsx'
# convert_hyperspectral_to_multispectral(spectral_response_function_file, site_rrs_path, output)


def merge_two_column(xlsx_path1, xlsx_path2, common_column, output_path):
    """
    Merge two excel sheets by common column values.
    :param xlsx_path1: Xlsx path 1
    :type xlsx_path1: String
    :param xlsx_path2: Xlsx path 2
    :type xlsx_path2: String
    :param common_column: The common column name that exists in both sheets
    :type common_column: String
    :param output_path: The output path
    :type output_path: String
    :return:
    :rtype:
    """
    path1 = pd.read_excel(xlsx_path1, 'Sheet1', engine='openpyxl')
    path2 = pd.read_excel(xlsx_path2, 'Sheet1', engine='openpyxl')
    # print (path1)
    new_df = pd.merge(path1, path2, on=common_column)
    # print (new_df)
    # final_df = pd.DataFrame(final_dic, columns=final_cols)
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    new_df.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(new_df, writer, 1, len(new_df.columns), 0)
    writer.save()
    # print (new_df)

spm = 'D:/MSU/codes/radiometer_processor/sampledata/2020/spm.xlsx'
# micasense_path = 'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense.xlsx'
# combined_path = 'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense_SPM.xlsx'
# merge_two_column(spm, micasense_path, 'sites', combined_path)

# micasense_path = 'D:/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense.xlsx'
combined_path = r'D:\MSU\SPM_MODIS\WMS_Rrs_modis_SPM.xlsx'
# merge_two_column(spm, output, 'sites', combined_path)


class RadiometerProcessor(QMainWindow, Ui_RrsProcessor):
    log = pyqtSignal(str)
    configFolder = pyqtSignal(str)

    def __init__(self):
        """
        The user interface dialog that loads settings and
        gives the users options to configure and run the data capture.
        """
        QMainWindow.__init__(self, None)
        self.setupUi(self)
        self._raw_data_folder = None
        self.raw_data_folder_btn.clicked.connect(
            lambda: self.file_dialog(self.raw_data_folder_le)
        )
        self.combine_folder_btn.clicked.connect(
            lambda: self.file_dialog(self.combine_folder_le)
        )
        self.combine_output_btn.clicked.connect(
             self.save_file_dialog
        )
        self.srf_btn.clicked.connect(self.open_a_file_dialog)
        self.site_rrs_btn.clicked.connect(self.open_a_file_dialog)
        self.converted_rrs_btn.clicked.connect(self.save_file_dialog)
        self.convert_rrs_btn.clicked.connect(self.convert_rrs_to_multispectral)
        self.srf_path = None
        self.site_rrs_path = None
        self.converted_rrs_path = None
        self.last_used_folder = None
        self.init_gui()
        self.combine_output_path = None
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.toolbar_cont.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)

        self.processed_df = None
        self.processed_cols = None

    def convert_rrs_to_multispectral(self):
        convert_hyperspectral_to_multispectral(
            self.srf_path,self.site_rrs_path, self.converted_rrs_path
        )

    def plot(self, file_name):
        ''' plot some random stuff '''

        self.figure.clear()
        # create an axis
        ax = self.figure.add_subplot(111)
        # plot data
        self.processed_df[file_name].plot(
            kind='line', x ='wavelength', y=self.processed_cols[file_name], ax=ax, title=os.path.basename(file_name))
        # refresh canvas
        self.canvas.draw()

    def generate(self):
        root = self.raw_tree.invisibleRootItem()
        child_count = root.childCount()

        for i in range(child_count):
            item = root.child(i)

            file_name = item.text(0)  # text at first (0) column
            state = item.checkState(0)
            if state == Qt.Unchecked:
                del self.processed_cols[file_name]
                del self.processed_df[file_name]
                continue

        for file_name, cols in self.processed_cols.items():
            # self.processed_df[file_name]
            self.processed_df[file_name]['avg_Rrs'] = self.processed_df[file_name][cols].mean(1)
            cols.append('avg_Rrs')
            # create_xlsx_chart(self.processed_df[file_name], file_name, cols)
        output_site_rrs = collect_rrs_to_single_sheet(self.processed_df)
        interpolate(output_site_rrs, 277, 1094)
            # print (self.processed_cols)
                # print (url)
            # df['avg_Rrs'] = df[cols_2].mean(1)

    def init_gui(self):
        """
        Initializes the widgets with signals
        :return:
        """
        self.log_boxes = {}

        self.setWindowFlags(
            self.windowFlags() |
            Qt.WindowSystemMenuHint |
            Qt.WindowMaximizeButtonHint | Qt.WindowMinimizeButtonHint
        )
        icon_path = self.resource_path('favicon.ico')
        self.setWindowIcon(QIcon(icon_path))

        config_path = self.resource_path('config.json')
        if path.isfile(config_path):
            self.config_path = config_path.replace('\\', '/')
            with open(self.config_path) as f:
                self.config = json.load(f)
        self.last_used_folder = self.config['raw_data_folder']

        self.raw_tree.itemClicked.connect(self.load_plot)
        self.raw_tree.itemChanged.connect(self.update_plot)
        self.generate_btn.clicked.connect(self.generate)

        self.combine_btn.clicked.connect(self.combine)

    def update_plot(self, item):
        state = item.checkState(0)
        text = item.text(0)
        if os.path.isfile(text):
            pass
        else:
            parent_item = item.parent()
            if parent_item is None:
                return
            file_name = parent_item.text(0)
            if state == Qt.Unchecked:
                self.processed_cols[file_name].remove(text)

            else:
                self.processed_cols[file_name].append(text)

            self.plot(parent_item.text(0))

    def load_plot(self, text):
        if isinstance(text, str):
            file_path = text
        else:
            file_path = text.text(0)

        if file_path in self.processed_df.keys():
            self.plot(file_path)

    def fill_item(self, item, value):
        def new_item(self, parent, text, val=None):
            # print (text)
            child = QTreeWidgetItem([text])
            child.setFlags(child.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsSelectable)
            child.setCheckState(0, Qt.Checked)
            self.fill_item(child, val)
            parent.addChild(child)
            child.setExpanded(True)

        if value is None:
            return

        elif isinstance(value, dict) or isinstance(value, OrderedDict):
            for key, val in value.items():
                new_item(self, item, key, val)
        elif isinstance(value, (list, tuple)):
            for val in value:
                text = (str(val) if not isinstance(val, (dict, list, tuple, OrderedDict))
                        else '[%s]' % type(val).__name__)
                new_item(self, item, text, val)

    def fill_widget(self, widget, value):
        if value is not None:
            widget.clear()
        self.fill_item(widget.invisibleRootItem(), value)


    @staticmethod
    def resource_path(relative_path):
        """
         Get absolute path to resource, works for dev and for PyInstaller
        :param relative_path:
        :return:
        """
        try:
            # PyInstaller creates a temp folder and stores path in _MEIPASS
            base_path = sys._MEIPASS
        except Exception:
            base_path = path.abspath(".")

        return path.join(base_path, relative_path)

    def show_message(self, title, msg):
        """
        Show messages as a popup.
        :param title: The type of the message. "Error", "Information"
        :type title: String
        :param msg: The message to be displayed.
        :return:
        """
        msgBox = QMessageBox()
        if title == 'Error':
            icon = QMessageBox.Critical
        elif title == "Information":
            icon = QMessageBox.Information
        else:
            icon = QMessageBox.Information
        msgBox.setIcon(icon)
        msgBox.setText(msg)
        msgBox.setWindowTitle('Radiometer Processor {}'.format(title))
        msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        icon_path = self.resource_path('favicon.ico')
        msgBox.setWindowIcon(QIcon(icon_path))
        msgBox.exec_()

    def save_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(
            self, "Write File Name", self.last_used_folder,
            "Spreadsheet Files (*.csv *.xls *.xlsx)", options=options
        )
        if file_name:
            # print(file_name)
            if not file_name.endswith('.xlsx'):

                file_name = '{}.xlsx'.format(file_name)
            sender = self.sender().objectName()
            if sender == 'combine_output_btn':
                self.combine_output_le.setText(file_name)
                self.combine_output_path = file_name
            elif sender == 'converted_rrs_btn':
                self.converted_rrs_le.setText(file_name)
                self.converted_rrs_path = file_name
            self.last_used_folder = os.path.dirname(file_name).replace('\\', '/')



    def open_a_file_dialog(self):
        # print (self.sender().objectName())
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select a File", self.last_used_folder,
            "Spreadsheet Files (*.csv *.xls *.xlsx)", options=options)

        if file_name:
            # print(file_name)
            if not file_name.endswith(('.xls','.xlsx', '.csv')):
                return
            sender = self.sender().objectName()
            if sender == 'srf_btn':
                self.srf_le.setText(file_name)
                self.srf_path = file_name
                self.last_used_folder = os.path.dirname(file_name).replace('\\', '/')

            if sender == 'site_rrs_btn':
                self.site_rrs_le.setText(file_name)
                self.site_rrs_path = file_name
                self.last_used_folder = os.path.dirname(file_name).replace('\\', '/')

    def file_dialog(self, line_edit):
        """
        Displays a file dialog for a user to specify a GPX Folder.
        :param line_edit: The line edit in which the folder is going to be set.
        :type line_edit: QLineEdit
        """
        title = QApplication.translate(
            "SensorsBridge",
            "Select a Folder"
        )
        if self.last_used_folder is None:
            last_path = expanduser("~")
        else:
            last_path =  self.last_used_folder

        path = QFileDialog.getExistingDirectory(
            self.parent(),
            title,
            last_path,
            QFileDialog.ShowDirsOnly
        )
        path = path.replace('\\', '/')
        if len(path) == 0:
            return
        if len(path) > 0:
            line_edit.setText(path)

        if line_edit == self.raw_data_folder_le:
            self._raw_data_folder = path
            self.last_used_folder = path
            # raw_path = self.config_folder_le.text()
            self.processed_cols, self.processed_df = calculate_reflectance(path)
            # print (processed)
            self.fill_widget(self.raw_tree, self.processed_cols)
            if len(self.processed_cols.keys()) >0:
                self.load_plot(list(self.processed_cols.keys())[0])
        if line_edit == self.combine_folder_le:

            combine_path = self.combine_folder_le.text()
            self.last_used_folder = combine_path
            combine_files = OrderedDict()
            for root, dirs, files in os.walk(combine_path):
                for i, fi in enumerate(files):
                    file_path = os.path.join(root, fi)
                    if fi.endswith('interpolate.xlsx') and os.path.isfile(file_path) and not fi.startswith('~$'):
                        combine_files[file_path] = fi

            self.fill_widget(self.combine_tree, combine_files)
            # if len(self.processed_cols.keys()) > 0:
            #     self.load_plot(list(self.processed_cols.keys())[0])


    def combine(self):
        root = self.combine_tree.invisibleRootItem()
        child_count = root.childCount()
        combine_files = []
        for i in range(child_count):
            item = root.child(i)

            file_name = item.text(0)  # text at first (0) column
            state = item.checkState(0)
            if state == Qt.Checked:
                combine_files.append(file_name)

        merge_spreadsheet(combine_files, self.combine_output_path)
        # for file_name, cols in self.processed_cols.items():
        #     # self.processed_df[file_name]
        #     self.processed_df[file_name]['avg_Rrs'] = self.processed_df[file_name][cols].mean(1)
        #     cols.append('avg_Rrs')
        #     # create_xlsx_chart(self.processed_df[file_name], file_name, cols)
        # output_site_rrs = collect_rrs_to_single_sheet(self.processed_df)


if __name__ == "__main__":
    import sys
    # Prevent the multiprocessing window from executiving this several times.

    processes = []
    threads = []


    def main():
        """
        The main method that opens the application dialog and
        connects run and terminate process signals.
        :return:
        """
        app = QApplication(sys.argv)
        main = RadiometerProcessor()
        main.show()
        sys.exit(app.exec_())

    main()
#
# raw_path = r'D:\MSU\SPM Testing\WMS_2020_2021\WMS_1_12_06022021'
# calculate_reflectance(raw_path)