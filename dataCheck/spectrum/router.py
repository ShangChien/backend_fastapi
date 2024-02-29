from typing import Any
from fastapi import APIRouter, Request, Body
from dataCheck.spectrum.utils import get_peaks,pre_process,cosine_similarity
from dataCheck.spectrum.type import RES, UVData, GetData, Result
import numpy as np
import asyncio
import logging

logger = logging.getLogger(__name__)

semaphore = asyncio.Semaphore(1)

router = APIRouter(
    prefix="/spectrum", # absolute path : localhost:port/dataCheck/spectrum
    tags=["spectrum"],
)

@router.get("/get_names")
async def get_fig_name_list(req: Request) -> RES[list[str]]: 
    uv_data: dict[str,dict] = req.app.state.uv_data
    nameList: list[str] = [k for k in uv_data]
    return RES[list[str]](data = nameList)

@router.post("/get")
async def get_uv_data(req: Request, data: GetData= Body()) -> RES[UVData]: 
    uv_data: dict[str,dict] = req.app.state.uv_data
    if data.name in uv_data:
        uv_res= UVData(**uv_data[data.name])
        return RES[UVData](data = uv_res)
    return RES[UVData](error = f'parameter error: no data named {data.name}', success=False)
    

@router.post("/del")
async def del_uv_data(req: Request, data: GetData = Body()) -> RES[str]: 
    uv_data: dict[str,dict] = req.app.state.uv_data
    if data.name in uv_data:
        async with semaphore:
            del uv_data[data.name]
            return RES[str](data = f'successful delete {data.name}')
    return RES[str](data = f'no data named {data.name}', success=False)

async def put_actions(data: UVData, req: Request) -> RES[list[float]]:
    uv_data: dict[str,dict] = req.app.state.uv_data
    if len(data.raw_arr) != 401:
        return RES(error = 'parameter error: raw_arr length must be 401', success=False)
    try:
        y_data,_scalar = pre_process(np.array(data.raw_arr)[::-1]) # type: ignore
        peaks_indices = get_peaks(y_data)
        peaks_arr= np.zeros(401)
        peaks_arr[peaks_indices] = np.array(y_data)[peaks_indices]
        async with semaphore:
            uv_data[data.name] = {'name':data.name, 'raw_arr':y_data, 'peaks_arr':peaks_arr.tolist()}
        return RES(data=peaks_arr.tolist())
    except Exception as e:
        return RES(error = f'fail to update with error: {e}', success=False)


@router.post("/put")
async def put_uv_data(req: Request, data: UVData = Body()) -> RES[list[float]]: 
    # print(data)
    return await put_actions(req=req, data=data)

@router.post("/check")
async def check_uv_data(req: Request, data: UVData = Body()) -> RES[list[Result]]:
    uv_data: dict[str,dict] = req.app.state.uv_data 

    ## 先put uvdata
    put_op: RES[list[float]] = await put_actions(req=req, data=data)
    if not put_op.success or not put_op.data:
        return RES(error=put_op.error, success=False)
    ## 获取所有数据拼接为矩阵

    names=[key for key in uv_data]
    all_arr: list[Any] = [uv_data[i]['peaks_arr'] for i in uv_data]
    matrix = np.concatenate(all_arr).reshape(-1,401)
    similarities: list[float] = cosine_similarity(put_op.data, matrix)
    res: list[Result]=[Result(name=k,similarity=v) for k,v in zip(names,similarities)]
    return RES(data=res)