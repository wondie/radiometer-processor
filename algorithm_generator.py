import math
import time
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error, mean_absolute_percentage_error
from sklearn import model_selection
from matplotlib import pyplot as plt
start = time.time()
import re

from itertools import permutations, combinations, combinations_with_replacement
from scipy import optimize
import pandas as pd
from mgsub import mgsub
# Importing required libraries
from sklearn.datasets import load_breast_cancer

from sklearn.model_selection import KFold, cross_val_score
from sklearn.linear_model import LinearRegression
from sklearn.metrics import accuracy_score


def minimize_1(coef, df, bands, return_df=False):
    a, b, c = bands
    x = ((df[c] * df[a] ** coef[0]) / (df[b] **
                                       coef[1])) + (df[b] / df[a]) ** coef[2]
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
    # print(bands)
    a, b = bands
    # 'a**2-b**2'
    x =  (coef[0] ** np.log(df[a]) - coef[1] ** np.log(df[b]))
    # we want to maximize, so we have to multiply by -1
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
    bnds = [(-100, 100) for i in range(0, len(coefficients))]
    res = optimize.minimize(
        fun, coefficients, args=(df, bands,),
        bounds=bnds, options={"disp": False},
        method='L-BFGS-B'
    )

    if res.success:
        if (-1 * res.fun) >= minimum_rsquared:
            # that's the highest rsqared we can get
            bands_formatted = ['df.{}'.format(b) for b in bands]
            return {
                'rsquared': -1 * res.fun,
                'coefficients': list(res.x),
                'bands_f': bands_formatted,
                'bands': bands

            }
    else:
        pass
        # print("Sorry, the optimization was not successful. Try with another initial"
        #       " guess or optimization method")


# solver('D:/MSU/codes/radiometer_processor/data/mic.xlsx')

def get_model_errors(actual, predicted):
    # actual = [0, 1, 2, 0, 3]
    # predicted = [0.1, 1.3, 2.1, 0.5, 3.1]
    mse = mean_squared_error(actual, predicted)
    msep = mean_absolute_percentage_error(actual, predicted)
    mae = mean_absolute_error(actual, predicted)
    rmse = math.sqrt(mse)
    return rmse, msep, mae


def validate_model(model, coef, bands, test_df, fun, i, r2):
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
    rmse, msep, mae = get_model_errors(actual, predicted)
    plt.scatter(actual, predicted, c='crimson')

    p1 = max(max(predicted), max(actual))
    p2 = min(min(predicted), min(actual))
    plt.axline((p2, p2), (p1, p1))
    plt.xlabel('Measured', fontsize=15)
    plt.ylabel('Predicted', fontsize=15)
    plt.axis('square')
    print('R2: ', f'{r2:.3f}', 'RMSE: ', f'{rmse:.3f}','RMSE %: ',
          f'{msep:.3f}', 'MAE: ', f'{mae:.3f}')
    plt.show()

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
        model.fit(X_train.values.reshape(-1, 1), y_train)

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

def permituate_equation_by_band(input_list, function_string, fun, df, df_t, minimum_rsquared=0.70):
    """
    Uses an equation and tries different combinations of columns in the input_list to
    come up with list of equation above the R-squared provided.
    :param input_list: is the file path containing the first column as
                       a Y axis and the subsequent columns as X axis.
    :type input_list: List
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
    place_holders = re.findall(r"[a-z]", function_string)
    place_count = len(place_holders)
    all_permutation = list(combinations_with_replacement(input_list, place_count))

    for i, permutation in enumerate(all_permutation):
        result = solver(df, fun, permutation, coefficients, minimum_rsquared)
        if type(result) == dict:
            res_coef = [str(co) for co in result['coefficients']]
            eq_1 = mgsub(function_string, coefficients_str, res_coef)
            eq_2 = mgsub(eq_1, place_holders, list(result['bands']))
            print(result['rsquared'], eq_2, fun.__name__)
            # x_out = pd.eval(eq_2)
            x_out = fun(result['coefficients'], df, result['bands'], True)
            df['x_{}'.format(i)] = x_out
            # print(type(df['x_{}'.format(i)]))
            X = df['x_{}'.format(i)]
            # print (X)
            y = df['SPM']
            # print (type(y))
            model, r2 = kfold_validation(X, y)
            # if (result['rsquared'] > 0.65):
            validate_model(model, result['coefficients'], result['bands'], df_t, fun, i, r2)


# X, y = make_classification(n_samples=1000, n_features=20, n_informative=15, n_redundant=5, random_state=1)
# print('X', X)
mic_headers = ['green', 'red', 'rededge', 'NIR']

equations = {
    '(a**2)**5.478/((b**1.133)*c)**5.049': minimize_1,
    '(((a/b)**3.1/(c/d)**3.1)**0.3)*((e/f)**3.1/(g/h)**3.1)': minimize_2,
    'a**2/(b-c)**3.1': minimize_3,
    'a**2/(b+c)**3.1': minimize_4,
    'a**2/(b*c)**3.1': minimize_5,
    '(a-b)**2/(c-d)**2': minimize_6,
    '(a*b)**2/(c*d)**2': minimize_7,
    '(a+b)**2/(c+d)**2': minimize_8,
    '(a-b)**2/(c+d)**2': minimize_9,
    '(a+b)**2/(c-d)**2': minimize_10,
    'a**2-b**2': minimize_11,
    '(((a*b)^3.1/(c*d)^3.1)^3.1)+(e/f)^3.1': minimize_12,
    '(((a*b^3.1)/(c*d^3.2)))+(e/f)^3.1': minimize_13,
    '(((a*b^3.1)/(c*d^3.2)))-(e/f)^3.1': minimize_14,
    '(((a/b)^3.1/(c/d)^3.1) * ((e/f)^3.1/(g/h)^3.1))': minimize_15,
    '3.3^(a)/3.3^(b) - 3.3^(e)/3.3^(f)': minimize_16,
    '3.3^(a)/3.3^(b)': minimize_17
}

file_path = 'D:/MSU/codes/radiometer_processor/data/mic.xlsx'
testing_path = 'D:/MSU/dissertation/SPM Testing/Mic4_3_final_testing.xlsx'
algorithm_data_df = pd.read_excel(file_path)
testing_data_df = pd.read_excel(testing_path)

for eq, fun in equations.items():
    permituate_equation_by_band(mic_headers, eq, fun, algorithm_data_df, testing_data_df)
# kfold_test(df)
end = time.time()
time_taken = end - start
print('Time taken in Sec: {}'.format(time_taken))
