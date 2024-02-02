from typing import Any
from fastapi import APIRouter, Request, Body
from dataCheck.spectrum.utils import get_peaks,pre_process
from dataCheck.spectrum.type import RES, UVdata, GetData
import numpy as np
import asyncio
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix

semaphore = asyncio.Semaphore(1)

router = APIRouter(
    prefix="/spectrum", # absolute path : localhost:port/dataCheck/spectrum
    tags=["spectrum"],
)

@router.get("/get_fig_name_list")
async def get_fig_name_list(req: Request) -> RES[list[str]]: 
    uv_data: dict[str,UVdata] = req.app.state.uv_data
    nameList: list[str] = [k for k in uv_data]
    return RES[list[str]](data = nameList)

@router.post("/get_uv_data")
async def get_uv_data(req: Request, data: GetData= Body()) -> RES[UVdata]: 
    uv_res = UVdata(name='',raw_arr=[],peaks_arr=[])
    uv_data: dict[str,UVdata] = req.app.state.uv_data
    if data.name in uv_data:
        uv_res: UVdata = uv_data[data.name]
        return RES[UVdata](data = uv_res)
    return RES[UVdata](error = f'parameter error: no data named {data.name}')
    

@router.post("/del_uv_data")
async def del_uv_data(req: Request, data: GetData = Body()) -> RES[str]: 
    uv_data: dict[str,UVdata] = req.app.state.uv_data
    if data.name in uv_data:
        async with semaphore:
            del uv_data[data.name]
            return RES[str](data = f'successful delete {data.name}')
    return RES[str](data = f'no data named {data.name}')

@router.post("/put_uv_data")
async def put_uv_data(req: Request, data: UVdata = Body()) -> RES[str]: 
    ## print(data,query(data))
    uv_data: dict[str,UVdata] = req.app.state.uv_data
    if len(data.raw_arr) != 401:
        return RES[str](error = 'parameter error: raw_arr length must be 401')
    try:
        y_data,_scalar = pre_process(data.raw_arr[:,1][::-1]) # type: ignore
        peaks_indices = get_peaks(y_data)
        peaks_arr= np.zeros(401)
        peaks_arr[peaks_indices] = y_data[peaks_indices]
        data.peaks_arr = peaks_arr.tolist()
        async with semaphore:
            uv_data[data.name]=data
        return RES[str](data=f'success updated {data.name}')
    except Exception as e:
        return RES[str](error = f'fail to update with error: {e}')

@router.post("/check_uv_data")
async def check_uv_data(req: Request, data: UVdata = Body()) -> RES[list[float]]:
    uv_data: dict[str,UVdata] = req.app.state.uv_data 
    ## 先put uvdata
    _result: RES[str] = await put_uv_data(req, data)
    ## 对比返回结果
    all_arr: list[Any] = [uv_data[i].peaks_arr for i in uv_data]
    matrix = np.concatenate(all_arr).reshape(-1,401)
    data_matrix_sparse = csr_matrix(matrix)
  
    # 稀疏数组
    target_array_sparse = csr_matrix(data.peaks_arr)  # (1, 401)
  
    # 将稀疏矩阵转换为密集格式
    data_matrix_dense = data_matrix_sparse.toarray()
    target_array_dense = target_array_sparse.toarray()
  
    # 计算余弦相似度
    similarities = 1 - cdist(target_array_dense, data_matrix_dense, metric='cosine')
  
    return RES[list[float]](data=similarities[0])