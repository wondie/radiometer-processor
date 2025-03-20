import glob
import os
from collections import OrderedDict
import math
import pandas as pd
import time

path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/processed/'

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
        irr_col = 'Raw_Irr_{}'.format(suff)
        irr_cols.append(irr_col)
        headers = ['Wavelength', 'Reference', irr_col]
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
        new_col = col.replace('Raw_', '')
        data0[new_col] = data0[col]/(10**6)
        cols.append(new_col)

    data0['Avg_Irr'] = data0[cols].mean(1)
    output_irr = irr_files_sorted[0].replace('_000', '').replace('.asc', '.xlsx')
    create_xlsx_chart(data0, output_irr)

def calculate_reflectance(path):
    """
    Loop through site water observation asc files and create a csv
    file containing all the observation and the average. You will need to see
    :param path: The path containing all raw radiometer data
    :type path: String
    :return:
    :rtype:
    """
    water_files = glob.glob('{}/*{}*.asc'.format(path, '[wW][aA][tT][eE][rR]'))
    water_files_sorted = sorted(water_files)
    df = None
    water_cols = []
    for file_path in water_files_sorted:
        suff = file_path[file_path.rindex('_') + 1:].replace('.asc', '')
        water_col = 'Raw_water_{}'.format(suff)

        headers = ['Wavelength', 'Reference', water_col]
        if suff == '000':
            water_cols.append(water_col)
            df = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                               names=headers)
        elif suff == '003':
            data12 = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                                 names=headers)
            if df is not None:
                df['SKY'] = data12[water_col]

        else:
            water_cols.append(water_col)
            data12 = pd.read_fwf(file_path, skiprows=14, skipfooter=1,
                                names=headers)
            if df is not None:
                df[water_col] = data12[water_col]

    cols = []

    for col in water_cols:
        new_col = col.replace('Raw_', '')

        df[new_col] =  (df[col] - (0.02*df['SKY']))/((100/99)*df['Reference']*math.pi)

        cols.append(new_col)

    if df is None:
        print ('Unable to find water files in ', path)
        return
    df['Avg_Rrs'] = df[cols].mean(1)

    output_xlsx = water_files_sorted[0].replace('_000', '').replace('.asc', '.xlsx')
    create_xlsx_chart(df, output_xlsx)

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

def create_xlsx_chart(df, output_xlsx, chart_start=None):
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
    if chart_start is None:
        chart_start = len(df.columns) - 3
    df.to_excel(writer, sheet_name='Sheet1')
    create_chart(df, writer, chart_start, len(df.columns) + 1)
    # writer.close()
    print ('Written: ', output_xlsx)


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
    chart.set_x_axis({'name': 'Wavelength'})
    chart.set_y_axis({'name': 'Rrs'})
    chart_sheet = workbook.add_chartsheet()
    # Configure the chart.
    chart_sheet.set_chart(chart)

# After creating average using good observation if there is a bad data,
# copy all to one folder in this case Site. They are required to be xlsx file containing WMS and Water
all_rrs = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/Sites'
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
    print (water_files)
    df_cols = OrderedDict()
    writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')
    for i, f in enumerate(water_files):
        file_name = os.path.basename(f)
        col_names = file_name.split('_')
        if col_names[1] in exclusions:
            continue
        col_name = '{}_{}'.format(col_names[0], col_names[1])
        print (col_name)
        df = pd.read_excel(f, sheet_name='Sheet1')
        if i == 0:
            df_cols['Wavelength'] = (df['Wavelength'])
        df_cols[col_name] = df['Avg_Rrs']

    new_df = pd.concat(df_cols.values(), axis=1, keys=df_cols.keys())
    new_df.to_excel(writer, sheet_name='Sheet1')
    create_chart(new_df, writer, 2, len(new_df.columns) + 1)
    # writer.close()

# collect_sites_rrs_to_single_sheet(all_rrs, exclusions=['57', '58', '63'])
# process_reflectance('C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS')
rrs = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/Sites/WMS_Rrs.xlsx'

def interpolate(xlsx_file_path, min, max):
    """
    Interpolate the site Rrs data into one-step increment
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
        df = df.append({'Wavelength': i}, ignore_index=True)
    df.sort_values(by=['Wavelength'], inplace=True)

    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df = df.interpolate()
    y = df[df['Wavelength'].isin(range(min+1, max+1))]

    output_xlsx = xlsx_file_path.replace('.xlsx', '_interpolate.xlsx')

    writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')
    y.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(y, writer, 1, len(y.columns), 0)
    # writer.close()

#interpolate(rrs, 277, 1094)
def merge_spreadsheet(xlsx_path1, xlsx_path2, combined_path):
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

    df3 = pd.merge(df1, df2, on="Wavelength")
    writer = pd.ExcelWriter(combined_path, engine='xlsxwriter')
    df3.to_excel(writer, sheet_name='Sheet1', index=False)
    create_chart(df3, writer, 1, len(df3.columns), 0)
    # writer.save()

# merge_spreadsheet(
#     'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/Radiometer_rrs.xlsx',
#     'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_interpolate.xlsx',
#   'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_combined.xlsx'
# )


def convert_hyperspectral_to_multispectral(srf_csv, rrs_csv, output_rrs):
    """
    Applies a spatial response function of multi-spectral sensors
    on radiometer reflectance values to come up with multi-spectral
    reflectance values.
    :param srf_csv: srf_csv is a csv file containing the following columns:
    wavelength, response value for each subsequent bands from the start
    to end of the wavelengths.
    :type srf_csv: String
    :param rrs_csv: contains a csv file with reflectance value for each site.
    The columns are wavelength(should match srf_csv wavelength),
    sites reflectance values for each wavelength.
    :type rrs_csv: String
    :param output_rrs: The output rrs path
    :type output_rrs: String
    :return:
    :rtype:
    """
    srf = pd.read_csv(srf_csv)
    rrs = pd.read_csv(rrs_csv)

    new_df = pd.merge(rrs, srf, on='Wavelength')
    final_cols = ['Sites']
    final_cols.extend([c for c in srf.columns if c != 'Wavelength'])

    final_dic = {}
    final_dic['Sites'] = [c for c in rrs.columns if c != 'Wavelength']

    for col_rrs in rrs:
        if col_rrs == 'Wavelength':
            continue

        for col_srf in srf:
            if col_srf == 'Wavelength':
                continue

            p1 = new_df[col_srf] * new_df[col_rrs]
            if col_srf in final_dic.keys():
                # add the first calculated band
                final_dic[col_srf].append(sum(p1) / sum(new_df[col_srf]))
            else:
                final_dic[col_srf] = [sum(p1) / sum(new_df[col_srf])]

    final_df = pd.DataFrame(final_dic, columns=final_cols)
    with pd.ExcelWriter(output_rrs, engine='xlsxwriter') as writer:
        final_df.to_excel(writer, sheet_name='Sheet1', index=False)
        create_chart(final_df, writer, 1, len(final_df.columns), 0)
        writer.close()

# spectral_response_function_file = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/Spectral_response_function_micasense_277_1094.csv'
# site_rrs_path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_interpolate.csv'
# output = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense.xlsx'
# convert_hyperspectral_to_multispectral(spectral_response_function_file, site_rrs_path, output)

spectral_response_function_file = r'C:\Users\andex\OneDrive\Documents\MSU\SPM_MODIS\Spectral_response_function_modis.csv'
site_rrs_path = r'C:\Users\andex\OneDrive\Documents\MSU\SPM_MODIS\WMS_Rrs_combined.csv'
output = r'C:\Users\andex\OneDrive\Documents\MSU\SPM_MODIS\WMS_Rrs_modis.xlsx'
# convert_hyperspectral_to_multispectral(spectral_response_function_file, site_rrs_path, output)


def merge_two_column(xlsx_path1, xlsx_path2, common_column, output_path):
    """
    Merge two Excel sheets by common column values.
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
    # print ('path1', path1)
    # print ('path2', path2)
    new_df = pd.merge(path1, path2, on=common_column)
    # print ('new_df', new_df)
    # final_df = pd.DataFrame(final_dic, columns=final_cols)
    print (output_path)
    # try:
    #     with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
    #         print ('err1')
    #         new_df.to_excel(writer, sheet_name='Sheet1', index=False)
    #         print ('err2')
    #         create_chart(new_df, writer, 1, len(new_df.columns), 0)
    # except Exception as err:

    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        print ('err1')
        time.sleep(15)
        new_df.to_excel(writer, sheet_name='Sheet1', index=False)
        print ('err2')
        create_chart(new_df, writer, 1, len(new_df.columns), 0)

    print ('new_df')

spm = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/spm.xlsx'
# micasense_path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense.xlsx'
# combined_path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense_SPM.xlsx'
# merge_two_column(spm, micasense_path, 'Sites', combined_path)

# micasense_path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/sampledata/2020/WMS_Rrs_micasense.xlsx'
combined_path = r'C:\Users\andex\OneDrive\Documents\MSU\SPM_MODIS\WMS_Rrs_modis_SPM.xlsx'
merge_two_column(spm, output, 'Sites', combined_path)
