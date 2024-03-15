import plotly.express as px
import plotly.graph_objects as go
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import numpy as np
import pickle
import requests
from scipy.stats import pearsonr, spearmanr

def render(df_data, col1='eg',col2='2d', label=''):
    # 创建数据框
    df_data = df_data.dropna(axis=0)
    df_data = df_data.reset_index(drop=True)
    x_data = df_data[col1]
    y_data = df_data[col2]
    names = df_data[label if (label != '') else col1]

    # 计算指标
    corr_p, _p_value = pearsonr(x_data, y_data)
    corr_sp, _p_value = spearmanr(x_data, y_data)
    rmse = mean_squared_error(x_data, y_data, squared=False)
    mae = mean_absolute_error(x_data, y_data)
    r2 = r2_score(y_data, x_data)

    # 画布
    fig = px.scatter(
        x=x_data,
        y=y_data,
        text=names,
        labels={'x': col1, 'y': col2},
        hover_name=names,
    )
    fig.update_layout(
        width=1200,  # 设置画布宽度
        height=1200,  # 设置画布高度
        xaxis=dict(scaleanchor="y", scaleratio=1),  # 设置x轴的比例与y轴相同
        yaxis=dict(scaleanchor="x", scaleratio=1),  # 设置y轴的比例与x轴相同
        shapes=[
            go.layout.Shape(
                type="line",
                x0=min(x_data),
                y0=min(x_data),
                x1=max(x_data),
                y1=max(x_data),
                line=dict(color="gray"),
            )
        ],
        annotations=[
            go.layout.Annotation(
                text=f"{col1.upper()}_{col2.upper()}<br>RMSE: {rmse:.2f}<br>MAE: {mae:.2f}<br>R2: {r2:.2f}<br>Corr_P: {corr_p:.2f}<br>Corr_SP: {corr_sp:.2f}",
                xref="paper",
                yref="paper",
                x=0,
                y=1,
                showarrow=False,
                font=dict(size=12),
            )
        ]
    )
    fig.show()

def getDatafromModel(inputData):
    url = "http://192.168.2.233:5050/unimol"
    headers = {
        "accept": "application/json, text/plain, */*",
        "accept-language": "zh-CN,zh;q=0.9,en;q=0.8",
        "content-type": "application/json",
        "sec-ch-ua": "\"Chromium\";v=\"118\", \"Google Chrome\";    v=\"118\", \"Not=A?Brand\";v=\"99\"",
        "sec-ch-ua-mobile": "?0",
        "sec-ch-ua-platform": "\"Windows\"",
        "sec-fetch-dest": "empty",
        "sec-fetch-mode": "cors",
        "sec-fetch-site": "same-site",
    }
    payload = {
        "atoms": inputData['atoms'],
        "coordinates": inputData['coordinates'],
        "results": None,
        "models": ["repr"],#"LUMO","Eg"
        "names": inputData['names'], # dict_mol['name_list'][:1000],
        "smiles": None,
        "molBlocks": None, #dict_mol['molBlock_list'][:1000],
    }
    response_init_struct = requests.post(url, headers=headers,  json=payload)
    return response_init_struct.json()