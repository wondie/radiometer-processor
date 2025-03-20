
import os.path
import time
from collections import OrderedDict
import pandas as pd
import numpy as np

from sklearn.metrics import mean_squared_error, mean_absolute_error, mean_absolute_percentage_error
# from sklearn import model_selection
# from matplotlib import     pyplot as plt
from sympy.printing.latex import LatexPrinter

start = time.time()
import re


from itertools import combinations_with_replacement
from scipy import optimize, stats



from mgsub import mgsub
# Importing required libraries
# from sklearn.datasets import load_breast_cancer

from sklearn.model_selection import KFold, cross_val_score
from sklearn.linear_model import LinearRegression
# from sklearn.metrics import accuracy_score
from plot import plot
from sympy.core.function import _coeff_isneg as neg
from sympy import sympify, latex
from sympy import Add
# TODO algorithm have different r squared based on the number of placeholder bands used, eg. if a,b,c vs a,b,c,d,e with the same function see fun1
# TODO Simplified function ends up missing multiplication between two bands
BANDS = {'red':'R_{rs}668', 'rededge':'R_{rs}717', 'NIR': 'R_{rs}842'}

TESTING_OUTPUT_PATH = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/data/testing_output.xlsx'
def minimize_1(coef, df, bands, return_df=False):
    # print (bands)
    a, b, c = bands
    x = ((df[c] * df[a] ** coef[0]) / (df[b] **
                                       coef[1])) + (df[b] / df[a]) ** coef[2]
    # formula = eval(x)
    # if  (df['SPM'].corr(x)) ** 2 >= 0.7:
    #     x = (df[c] * df[a] ** round(coef[0], 4) / (df[b] **round(coef[1], 4))) + (df[b] / df[a]) ** round(coef[2], 3)
    #     # print (a, b, c)
    #     # print (bands, '(({0} * {1} ** {2}) / ({3} **{4})) + ({3} / {1}) ** {5}'.format(c,a,coef[0],b,coef[1],coef[2]))
    #     print(df['SPM'].corr(x)** 2)
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2
def minimize_1_1(coef, df, bands, return_df=False):
    # print (bands)
    a, b, c, d, e = bands
    x = ((df[a] * df[b] ** coef[0]) / (df[c] **
                                       coef[1])) + (df[d] / df[e]) ** coef[2]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_1_2(coef, df, bands, return_df=False):
    # print (bands)
    a, b, c = bands
    # (a**2)**5.478/((b**1.133)*c)**5.049
    x = (df[a] ** coef[0])** coef[1]/((df[b] ** coef[2]) * df[c] **coef[3])
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_2(coef, df, bands, return_df=False):
    a, b, c, d, e, f, g, h = bands
    x = (((df[a] / df[b]) ** coef[0] / (df[c] / df[d]) ** coef[1]) ** coef[2]) * (
                (df[e] / df[f]) ** coef[3] / (df[g] / df[h]) ** coef[4])
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_2_2(coef, df, bands, return_df=False):
    a, b, c = bands
    x = (((df[a] / df[b]) ** coef[0] / (df[c] / df[a]) ** coef[1]) ** coef[2]) * (
            (df[b] / df[c]) ** coef[3] / (df[a] / df[b]) ** coef[4])
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_3(coef, df, bands, return_df=False):
    a, b, c = bands
    # 'a/(b-c)'
    x = df[a] ** coef[0] / (df[b] - df[c]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_4(coef, df, bands, return_df=False):
    a, b, c = bands
    # 'a/(b+c)'
    x = df[a] ** coef[0] / (df[b] + df[c]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_5(coef, df, bands, return_df=False):
    a, b, c = bands
    # 'a/(b*c)'
    x = df[a] ** coef[0] / (df[b] * df[c]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_6(coef, df, bands, return_df=False):
    a, b, c, d = bands
    # (a-b)**2/(c-d)**2
    x = (df[a] - df[b]) ** coef[0] / (df[c] - df[d]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_7(coef, df, bands, return_df=False):
    a, b, c, d = bands
    # (a*b)**2/(c*d)**2
    x = (df[a] * df[b]) ** coef[0] / (df[c] * df[d]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_8(coef, df, bands, return_df=False):
    a, b, c, d = bands
    # (a*b)**2/(c*d)**2
    x = (df[a] * df[b]) ** coef[0] / (df[c] * df[d]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_9(coef, df, bands, return_df=False):
    a, b, c, d = bands
    # (a-b)**2/(c+d)**2
    x = (df[a] - df[b]) ** coef[0] / (df[c] + df[d]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_10(coef, df, bands, return_df=False):
    a, b, c, d = bands
    # (a+b)**2/(c-d)**2
    x = (df[a] + df[b]) ** coef[0] / (df[c] - df[d]) ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_11(coef, df, bands, return_df=False):
    a, b = bands
    # 'a**2-b**2'
    x = df[a] ** coef[0] - df[b] ** coef[1]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_12(coef, df, bands, return_df=False):
    a, b, c, d, e, f = bands
    # (((a*b)^3.1/(c*d)^3.1)^3.1)+(e/f)^3.1
    x = (((df[a] * df[b]) ** coef[0] / (df[c] * df[d]) ** coef[1]) ** coef[2]) + (df[e] / df[f]) ** coef[3]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_13(coef, df, bands, return_df=False):
    a, b, c, d, e, f = bands
    # '(((a*b^3.1)/(c*d^3.2)))+(e/f)^3.1'
    x = (((df[a] * df[b] ** coef[0]) / (df[c] * df[d] ** coef[1]))) + (df[e] / df[f]) ** coef[2]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_14(coef, df, bands, return_df=False):
    a, b, c, d, e, f = bands
    # (((a*b^3.1)/(c*d^3.2)))-(e/f)^3.1
    x = (((df[a] * df[b] ** coef[0]) / (df[c] * df[d] ** coef[1]))) - (df[e] / df[f]) ** coef[2]
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_15(coef, df, bands, return_df=False):
    a, b, c, d, e, f, g, h = bands
    # (((a/b)^3.1/(c/d)^3.1) * ((e/f)^3.1/(g/h)^3.1))
    x = (((df[a] / df[b]) ** coef[0] / (df[c] / df[d]) ** coef[1]) * (
                (df[e] / df[f]) ** coef[2] / (df[g] / df[h]) ** coef[3]))
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_16(coef, df, bands, return_df=False):
    # print(bands)
    a, b, c, d = bands
    # 'a**2-b**2'
    x = (coef[0] ** np.log(df[a]) /coef[1] ** np.log(df[b])) - (coef[2] ** np.log(df[c]) /coef[3] ** np.log(df[d]))
    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_17(coef, df, bands, return_df=False):

    a, b, c, d = bands
    # 'a**2-b**2'
    # x =  (coef[0] ** np.log(df[a]) - coef[1] ** np.log(df[b]))
    x = ((df[a]  ** coef[0]) / (df[b]** coef[1])) - (df[c] ** coef[2] / df[d] ** coef[3])

# we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    r2 = (df['SPM'].corr(x)) ** 2
    # if r2 >= 0.7:
    #     print ('((df["{0}"]  ** {1}) / (df["{2}"]** {3})) - (df["{4}"] ** {5}/ df["{6}"] ** {7})'.format(
    #         a, coef[0], b, coef[1], c, coef[2], d, coef[3]
    #     ))
    negative_r2 = -1 * r2

    return negative_r2

def minimize_18(coef, df, bands, return_df=False):
    # print(bands)
    a, b= bands
    # 'a**2-b**2'
    # x =  (coef[0] ** np.log(df[a]) - coef[1] ** np.log(df[b]))
    x = (df[a]  ** coef[0] / df[b]** coef[1]) - 1
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2

def minimize_19(coef, df, bands, return_df=False):

    a, b, c, d = bands
    # 'a**2-b**2'
    # x =  (coef[0] ** np.log(df[a]) - coef[1] ** np.log(df[b]))
    x = ((df[a]  ** coef[0]) / (df[b]** coef[1])) + (df[c] ** coef[2] / df[d] ** coef[3])

    # we want to maximize, so we have to multiply by -1
    if return_df:
        return x
    r2 = (df['SPM'].corr(x)) ** 2
    # if r2 >= 0.7:
    #     print ('((df["{0}"]  ** {1}) / (df["{2}"]** {3})) - (df["{4}"] ** {5}/ df["{6}"] ** {7})'.format(
    #         a, coef[0], b, coef[1], c, coef[2], d, coef[3]
    #     ))
    negative_r2 = -1 * r2

    return negative_r2

def compute_algorithm(df, return_df=False):
    x = (df["+"
            "red"]  ** 6.887001053122411 / df["rededge"]** 7.550940051329842) - (df["rededge"] ** 1.5186888746782579/ df["NIR"] ** 2.2402455268389008)
    if return_df:
        return x
    return -1 * (df['SPM'].corr(x)) ** 2


def solver(df, fun, bands, coefficients, minimum_rsquared):
    """
    Execute the minimize function
    :param df:
    :param fun:
    :param bands:
    :param coefficients:
    :param minimum_rsquared:
    :return:
    """
    # read your dataframe from somewhere. e.g. csv
    bounds = None
    if len(coefficients) > 0:
        bounds = [(-100, 100) for i in range(0, len(coefficients))]
    # print (fun.__name__, bounds, coefficients)
    res = optimize.minimize(
        fun, coefficients, args=(df, bands,),
        bounds=bounds, options={"disp": False},
        method='L-BFGS-B'
    )

    if res.success:
        if (-1 * res.fun) >= minimum_rsquared and (-1 * res.fun) < 1:
            # that's the highest rsqared we can get
            return {
                'rsquared': -1 * res.fun,
                'coefficients': list(res.x),
                # 'bands_f': bands_formatted,
                'bands': bands
            }
        # else:
        #     return {
        #         'rsquared': -1 * res.fun,
        #         'coefficients': list(res.x),
        #         # 'bands_f': bands_formatted,
        #         'bands': bands
        #     }
    return None
        # print("Sorry, the optimization was not successful. Try with another initial"
        #       " guess or optimization method")


# solver('C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/data/mic.xlsx')
# def mean_absolute_percentage_error(y_true, y_pred):
#     y_true, y_pred = np.array(y_true), np.array(y_pred)
#     return np.mean(np.abs((y_true - y_pred) / y_true)) * 100
def root_mean_squared_error_percent(y_true, y_pred):
    '''
    Compute Root Mean Square Percentage Error between two arrays.
    '''
    actual, estimated = np.array(y_true), np.array(y_pred)
    loss = np.sqrt(np.mean(np.square(((actual - estimated) / actual)), axis=0))

    return loss


def get_model_errors(actual, predicted):
    rmse = mean_squared_error(actual, predicted, squared=False)
    rmsep = root_mean_squared_error_percent(actual, predicted)*100

    mae = mean_absolute_error(actual, predicted)
    maep = mean_absolute_percentage_error(actual, predicted)*100

    return rmse, rmsep, mae, maep


def validate_model(model, coef, bands, test_df, fun, i, r2, training_y, training_x, equation_latex):
    """
    Validate a model using a testing data.
    @param model: The LinearRegression model to be tested.
    @param coef: The coefficients of the model
    @param bands: The bands used in the model
    @param test_df: Testing data frame
    @param fun: The function containing the source question of the model
    @param i: Index of the testing
    @param r2: R squared of the model
    @return:
    """
    x_out = fun(coef, test_df, bands, True)
    test_df['x_test_{}'.format(i)] = x_out
    predicted = list((model.coef_ * test_df['x_test_{}'.format(i)]) + model.intercept_)
    actual = list(test_df['SPM'])
    rmse, rmsep, mae, mape = get_model_errors(actual, predicted)

    df = pd.DataFrame({'actual': actual, 'predicted':predicted})

    writer = pd.ExcelWriter(TESTING_OUTPUT_PATH, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', header=True)
    writer.save()

    print('Equation: ', equation_latex, 'R2: ', f'{r2:.3f}', 'RMSE: ', f'{rmse:.3f}','RMSE %: ',
          f'{rmsep:.3f}', 'MAE: ', f'{mae:.3f}', 'MAEP: ',
          f'{mape:.3f}', fun.__name__)

    validation = OrderedDict([('RMSE(mg/L)',rmse), ('RMSE(%)',rmsep),
                              ('MAE(mg/L)', mae), ('MAE(%)',mape)])



    plot.create_scatter_plot(actual, predicted, validation, model, training_y,
                             training_x, r2, equation_latex, 'ppt')

    validation['Equation'] = equation_latex
    return validation
def negative_to_last(expr):
    if neg(expr.args[0]):

        args = list(expr.args)
        args.append(args.pop(0))
        # print (args, LatexPrinter({'order':'none'})._print_Add(Add(*args,evaluate=False)))
        return LatexPrinter({'order':'none'})._print_Add(Add(*args,evaluate=False))
    else:
        return latex(expr, mul_symbol='times', order='none')

def run_regression(X, y):
    # y = y.values.reshape(-1, 1)
    x = X.values.reshape(-1, 1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(X,y)
    # print(r_value, p_value, std_err)
    model = LinearRegression()
    model.fit(x, y)
    r2 = model.score(x, y)
    return model, r2, r_value, p_value, std_err

def kfold_validation(X, y):
    #
    # # Loading the dataset
    # data = load_breast_cancer(as_frame=True)
    # df = data.frame
    # X = df.iloc[:, :-1]
    # y = df.iloc[:, -1]

    # Implementing cross validation

    k = 10
    kf = KFold(n_splits=k, random_state=None)
    model = LinearRegression()

    acc_score = []

    for train_index, test_index in kf.split(X):
        # print (train_index, test_index)
        X_train, X_test = X[train_index], X[test_index]

        y_train, y_test = y[train_index], y[test_index]
        # print (X_train.values.reshape(-1, 1))
        try:
            model.fit(X_train.values.reshape(-1, 1), y_train)
        except ValueError:
            return None, None
        # pred_values = model.predict(X_test.values.reshape(-1, 1))
        # print(pred_values)
        acc = model.score(X_train.values.reshape(-1, 1), y_train)
        acc_score.append(acc)

    r2 = sum(acc_score) / k
    # if avg_acc_score > 0.75:
        # pyplot.errorbar(folds, means, yerr=[mins, maxs], fmt='o')
        # # plot the ideal case in a separate color
        # pyplot.plot(folds, [ideal for _ in range(len(folds))], color='r')
        # # show the plot
        # pyplot.show()
    # print('accuracy of each fold - {}'.format(acc_score))
    # print('Avg accuracy : {}'.format(avg_acc_score))
    return model, r2


def permituate_equation_by_band(band_list, function_string, fun, df, df_t, output_path, minimum_rsquared=0.71):
    """
    Uses an equation and tries different combinations of columns in the input_list to
    come up with list of equation above the R-squared provided.
    :param headers: Is the list containing all bands
    :type headers: List
    :param function_string: is an equation to be tested.
    :type function_string: String
    :param fun: The equation functional form as python function.
    :type fun: Function Object
    :param df: The dataframe containing the Xs and Y of the data.
    :param minimum_rsquared: the minimum r-squared used to filter the equations.
        If it is above the minimum r-squared, the equation is saved in a result file for review.
    :return:
    """
    coefficients_str = re.findall(r"\d+[.]\d+|\d+", function_string)
    coefficients = [float(a) for a in coefficients_str]
    all_place_holders = re.findall(r"[a-z]", function_string)
    place_holders = list(set(all_place_holders))
    place_count = len(place_holders)
    place_holders.sort()

    all_permutation = list(combinations_with_replacement(band_list, place_count))
    # print (input_list, place_count, all_permutation)
    # return
    for i, band_permutation in enumerate(all_permutation):
        for a in range(-10, 10):
            new_coefficients = np.random.uniform(a,a+1, len(coefficients))
            result = solver(df, fun, band_permutation, new_coefficients, minimum_rsquared)
            if result is not None:
                res_coef = [str(round(co, 4)) for co in result['coefficients']]
                eq_1 = mgsub(function_string, coefficients_str, res_coef)

                eq_2 = mgsub(eq_1, place_holders, list(result['bands']))
                print (eq_2)
                equation_latex = negative_to_last(sympify(eq_2))
                # print (equation_latex)
                for key in BANDS.keys():
                    equation_latex = re.sub(r"\b%s\b" % key, BANDS[key], equation_latex)

                X = fun(result['coefficients'], df, result['bands'], True)
                df['x_{}'.format(i)] = X

                y = df['SPM']
                # print (df)
                # model, r2 = kfold_validation(X, y)
                model, r2, r_value, p_value, std_err = run_regression(X, y)
                if model is None:
                    continue
                validation = validate_model(model, result['coefficients'], result['bands'], df_t,
                               fun, i, r2, y, X, equation_latex)
                # print(res.pvalue)
                validation['P Value'] = p_value
                validation.update({'R-squared':r2})
                validation.move_to_end('R-squared', last=False)
                validation.update({'Statistical Results':'Value'})
                validation.move_to_end('Statistical Results', last = False)

                validation_data = pd.DataFrame(validation, index=[0]).T
                output_dir = os.path.dirname(output_path)
                validation_path = os.path.join(output_dir, 'result_{}.xlsx'.format(r2))
                writer = pd.ExcelWriter(validation_path, engine='xlsxwriter')
                validation_data.to_excel(writer, sheet_name='Sheet1', header=None)
                writer.save()

# X, y = make_classification(n_samples=1000, n_features=20, n_informative=15, n_redundant=5, random_state=1)
# print('X', X)
mic_headers = ['red', 'rededge', 'NIR']
### IMPORTANT - coeffients can be any number but one equation shouldn't have duplicate coefficents
equations = {
    'c *(a ** 1 /b **2) + (b/ a) ** 3': minimize_1, # a, b, c, a, b
    # '(a *( b ** 1 / c **2)) + (d/ e) ** 3': minimize_1_1, # a, b, ,c, d, e
    # '(a**2)**5.478/((b**1.133)*c)**5.049': minimize_1_2,
    # '(((a/b)**3.1/(c/a)**3.2)**0.3)*((b/c)**3.3/(a/b)**3.4)': minimize_2_2,
    # 'a**2/(b-c)**3.1': minimize_3,
    # 'a**2/(b+c)**3.1': minimize_4,
    # 'a**2/(b*c)**3.1': minimize_5,
    # '(a-b)**2/(c-d)**4': minimize_6,
    # '(a*b)**2/(c*d)**4': minimize_7,
    # '(a+b)**2/(c+d)**4': minimize_8,
    # '(a-b)**2/(c+d)**4': minimize_9,
    # '(a+b)**2/(c-d)**4': minimize_10,
    # 'a**2-b**4': minimize_11,
    # '(((a*b)^3.1/(c*d)^3.2)^3.3)+(e/f)^3.4': minimize_12,
    # '(((a*b^3.1)/(c*d^3.2)))+(e/f)^3.3': minimize_13,
    # '(((a*b^3.1)/(c*d^3.2)))-(e/f)^3.3': minimize_14,
    # '(((a/b)^3.1/(c/d)^3.2) * ((e/f)^3.3/(g/h)^3.4))': minimize_15,
    # '3.1^(a)/3.2^(b) - 3.3^(e)/3.4^(f)': minimize_16,
    # '(a**3.1/b**3.2)-(c**3.3/d**3.4)': minimize_17,
    #
    # '(a^3.1/b^3.2)': minimize_18,
    # '(a**3.1/b**3.2)+(c**3.3/d**3.4)': minimize_19,
    # '(a/b)': minimize_19
}

file_path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/data/mic.xlsx'
testing_path = 'C:/Users/andex/OneDrive/Documents/MSU/dissertation/SPM Testing/Mic4_3_final_testing.xlsx'
output_path = 'C:/Users/andex/OneDrive/Documents/MSU/codes/radiometer_processor/data/output.xlsx'

algorithm_data_df = pd.read_excel(file_path)
testing_data_df = pd.read_excel(testing_path)
# print (compute_algorithm(algorithm_data_df))
for eq, fun in equations.items():
    permituate_equation_by_band(mic_headers, eq, fun, algorithm_data_df, testing_data_df, output_path)
# kfold_test(df)
end = time.time()
time_taken = end - start
print('Time taken in Sec: {}'.format(time_taken))
