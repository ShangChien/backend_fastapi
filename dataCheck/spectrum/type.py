from typing import Any, TypeVar, Generic, TypeGuard
from pydantic import BaseModel
from numpy import ndarray

class GetData(BaseModel):
    name: str|None = None
    index: int|None = None

class UVData(BaseModel):
    name: str
    raw_arr: list[int|float]
    peaks_arr: list[int|float] | None = []
    
class Result(BaseModel):
    name: str
    similarity: float

T = TypeVar('T')
class RES(BaseModel, Generic[T]):
    success: bool = True
    data: T | None = None
    error: str = ''

def is_ndarray(val: Any) -> TypeGuard[ndarray]:
    return isinstance(val, ndarray)