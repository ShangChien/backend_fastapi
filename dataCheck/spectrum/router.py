from fastapi import APIRouter, Request, Body
from dataCheck.spectrum.utils import get_peaks,pre_process
from pydantic import BaseModel
from typing import Any, TypeVar, Generic, Optional
import numpy as np
import asyncio
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix

semaphore = asyncio.Semaphore(1)

router = APIRouter(
    prefix="/spectrum", # absolute path : localhost:port/dataCheck/spectrum
    tags=["spectrum"],
)

class getUVdata(BaseModel):
  name: str|None = None
  index: int|None = None

class UVdata(BaseModel):
  name:str
  raw_arr:list[Any]
  peaks_arr:list[int|float|None|Any]|None = []

T = TypeVar('T')
class RES(BaseModel, Generic[T]):
    success: bool = True
    data: Optional[T] = None
    error: str= ''

@router.get("/get_fig_name_list")
async def get_fig_name_list(req: Request) -> RES[list[str]]: 
  nameList: list[str] = [uv.name for uv in req.app.state.uv_data]
  return RES[list[str]](data = nameList)

@router.post("/get_uv_data")
async def get_uv_data(req: Request, data: getUVdata = Body()) -> RES[UVdata]: 
  uv_res = UVdata(name='',raw_arr=[],peaks_arr=[])
  if data.name:
    uv_res: UVdata = [uv for uv in req.app.state.uv_data if uv.name == data.name][0]
  elif data.index:
    uv_res: UVdata = req.app.state.uv_data[data.index]
  else:
    return RES[UVdata](error = f'parameter error{data.model_dump_json()}')
  return RES[UVdata](data = uv_res)

@router.post("/del_uv_data")
async def del_uv_data(req: Request, data: getUVdata = Body()) -> RES[str]: 
  ##
  if data.name:
    target = [uv for uv in req.app.state.uv_data if uv.name == data.name]
    if len(target) > 0:
      async with semaphore:
        req.app.state.uv_data.pop(target[0])
        return RES[str](data = f'successful delete {data.name}')
    return RES[str](data = f'no data named {data.name}')
  return RES[str](data = f'field name is empty, data:{data.model_dump_json()}')

@router.post("/put_uv_data")
async def put_uv_data(req: Request, data: UVdata = Body()) -> RES[str]: 
  ## print(data,query(data))
  action:str=''
  target: list[UVdata] = [uv for uv in req.app.state.uv_data if uv.name == data.name]
  length = len(target)
  y_data,_scalar = pre_process(data.raw_arr[1:,1][::-1]) # type: ignore
  data.peaks_arr = get_peaks(y_data)
  if length == 0:
    action='add'
    async with semaphore:
      req.app.state.uv_data.append(data)
  else:
    action='update'
    async with semaphore:
      target[0].raw_arr = data.raw_arr
      target[0].peaks_arr = data.peaks_arr # type: ignore
  return RES[str](data=f'operation {action} for {data.name}')

@router.post("/check_uv_data")
async def check_uv_data(req: Request, data: UVdata = Body()) -> RES[list[float]]: 
  ## 先put uvdata
  _result: RES[str] = await put_uv_data(req, data)
  ## 对比返回结果
  matrix = np.concatenate([i['peaks_arr'] for i in req.app.state.uv_data]).reshape(-1,401)
  data_matrix_sparse = csr_matrix(matrix)

  # 假设results[1]['peaks_arr']是你要比较的稀疏数组
  target_array_sparse = csr_matrix(data.peaks_arr)  # (1, 401)

  # 将稀疏矩阵转换为密集格式
  data_matrix_dense = data_matrix_sparse.toarray()
  target_array_dense = target_array_sparse.toarray()

  # 计算余弦相似度
  similarities = 1 - cdist(target_array_dense, data_matrix_dense, metric='cosine')

  return RES[list[float]](data=similarities[0])