import pandas as pd
from pandas import DataFrame,Series
import math
import copy
from sklearn.preprocessing import MinMaxScaler
from scipy import signal
from numpy import ndarray
from pathlib import Path as P
from typing import Any
from scipy.special import voigt_profile
import numpy as np
from functools import partial
from scipy.optimize import minimize
import nest_asyncio
import logger_config
from logging import Logger

log: Logger = logger_config.get_logger(__name__)
nest_asyncio.apply()

# 定义要拟合的函数列表
class peak_funcs:
    
    @staticmethod
    def exp(x, a, b, c):
        return a * np.exp(b * (x-c))
    
    @staticmethod
    def gauss(x, A, mu, sigma):
        return A * np.exp(-(x - mu)**2 / (2 * sigma**2))
    
    @staticmethod
    def lorentz(x, A, mu, gamma):
        return A / (1 + ((x - mu) / gamma)**2)
    
    @staticmethod
    def voigt(x, A, mu, sigma, gamma):
        return A * voigt_profile(x - mu, sigma, gamma)

models = [getattr(peak_funcs, i) for i in dir(peak_funcs) if not i.startswith('__')]

# 定义损失函数
def loss(params, x, y, func):
    y_pred = func(x, *params)
    loss_v = np.sum((y - y_pred) + np.maximum(0.001, 10*(y_pred - y)))
    return loss_v

# 读取excel
def get_data_from_excel(file:P)-> dict[str, ndarray]:
    filename: str= file.stem
    df: DataFrame = pd.read_excel(file, header=None)
    df.dropna(axis=0, how='all', inplace=True)
    df.dropna(axis=1, how='all', inplace=True)
    mask: DataFrame = df.applymap(lambda x: 'nm' == str(x).strip())    ## 在df中寻找值为字符串'nm'索引
    mask_series: DataFrame | Series = mask.stack()   ## 二维数据打平转化为series(row,col,value)
    indices = mask_series[lambda x: x].index # type: ignore
    if len(indices) > 1:
        raise Exception(f'{filename}: 存在多个值为nm的单元格{indices}')
    elif len(indices) == 0:
        raise Exception(f'{filename}: 不存在值为nm的单元格')
    else:
        arr: ndarray = df.loc[indices[0][0]:,indices[0][1]:].to_numpy()
        name: str = f'{filename}-{str(arr[0][1])}'
        arr[0,1] = name
        return {filename: arr}

# 数据前处理
def pre_process(data:ndarray)-> tuple[ndarray, MinMaxScaler]:
    scaler = MinMaxScaler()
    arr_normalized = scaler.fit_transform(data.reshape(-1,1)).reshape(-1)
    arr_normalized = signal.savgol_filter(arr_normalized, window_length=10, polyorder=2)
    return arr_normalized, scaler

## 寻找峰值
def get_peaks(data:ndarray, threshold=10)-> list[int]:

    peaks_normal: ndarray
    _property:dict
    peaks_normal, _property = signal.find_peaks(data, prominence=0.002, distance=10)
    peaks_cwt: ndarray = signal.find_peaks_cwt(data, np.arange(1, 10), min_length=4, min_snr=1)
    ## 合并去重,过滤低值
    peaks_merged: list[Any] = sorted(list(set(peaks_normal.tolist() + peaks_cwt.tolist())))
    peaks=[i for i in peaks_merged if data[i] > 0.05]
    ## 筛选主峰
    diffs = np.diff(peaks)
    separators = np.where(diffs >= threshold)[0] + 1
    subarrays= np.split(peaks, separators)
    peaks=[]
    ## 密集区域稀疏化
    for sub in subarrays:
        if len(sub) == 1:
            sub = sub[0]
        else:
            value_in_peaks_normal =np.array([i for i in sub if i in peaks_normal])
            if len(value_in_peaks_normal) == 0:
                sub = sub.mean()
            else:
                index = np.argmin(value_in_peaks_normal - sub.mean())
                sub= value_in_peaks_normal[index]
        peaks.append(sub)
    print('peaks:',peaks)
    return peaks

# 迭代寻找峰值主函数
def iter_peaks(x_data, y_data, iter_num:int|None = None, results:list[dict] = []) -> list[dict]:
    """
    find the best fitting model for each peak.

    Args:
        x_data: The x-axis data points.
        y_data: The y-axis data points.
        iter_num: 最大迭代次数 (optional).
        results: 输出的结果 (optional).

    Returns:
        A list of fitting results, where each result contains:
            - name: The name of the model used for fitting.
            - params: The optimal parameters found for the model.
    """
    try:
        # 识别峰位
        peak_indexs = get_peaks(y_data)
        iter_num = iter_num if iter_num else len(peak_indexs)

        # 计算最高峰位的相关信息
        scale = len(y_data)
        max_peak_index= np.argmax(y_data[peak_indexs])
        max_intensity = y_data[peak_indexs[max_peak_index]]
        center = peak_indexs[max_peak_index] / scale
        _width_scipy=signal.peak_widths(y_data, [peak_indexs[max_peak_index]], rel_height=0.5)[0][0] / scale
        width = _width_scipy if _width_scipy > 0.02 else 0.02

        # 设置不同模型拟合函数和初猜值
        tasks = []
        for model in models:
            initial_func_guess=[]
            if model.__name__ in ['gauss','lorentz']:
                initial_func_guess = [max_intensity,center,width]
            elif model.__name__ == 'voigt':
                initial_func_guess = [max_intensity/4, center, width-0.01, width/2-0.01]
            elif model.__name__ == 'exp':
                initial_func_guess = [1.0, -10.0, -0.01]
            params = {
                'fun': partial(loss, func=model),
                'x0': initial_func_guess,
                'args': (x_data, y_data)
            }
            tasks.append({'name': model.__name__, 'params': copy.deepcopy(params)})

        # 并行加速运行拟合函数，并行失败，待研究
        ## task_results = Parallel(n_jobs=-1)(delayed(minimize)(**task['params']) for task in tasks)
        task_results=[minimize(**task['params']) for task in tasks]

        # 过滤拟合失败的结果
        task_results_filtered= [result for result in task_results if not math.isnan(result.fun)]

        # 选择拟合最好的模型
        optimal_fit_info = min(task_results_filtered, key=lambda x: x.fun)
        optimal_index = task_results.index(optimal_fit_info)
        optimal_params= optimal_fit_info.x
        model_func = models[optimal_index]

        # 保存当前拟合的最优模型参数
        results.append({
            'name': model_func.__name__,
            'params': optimal_params,
        })

        # 初始数据减去拟合函数的值，生成新的待拟合数据
        y_fit= model_func(x_data, *optimal_params)
        y_new = y_data - y_fit

        # 递归拟合上一步的残差, 直至iter_num == 0
        iter_num -= 1
        if iter_num != 0:
            return iter_peaks(x_data, y_new, iter_num, results)
        else:
            return results
    except Exception as e:
        print(f'peak process error in the last {iter_num} iteration: {e}')
        return results