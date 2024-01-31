from typing import Any, TypeVar, Generic, Optional
from pydantic import BaseModel


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
