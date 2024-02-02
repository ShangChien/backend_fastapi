from typing import Any, TypeVar, Generic, Optional, TypeGuard, Sequence as Seq
from pydantic import BaseModel
from numpy import ndarray

class GetData(BaseModel):
    name: str|None = None
    index: int|None = None

class UVdata(BaseModel):
    name: str
    raw_arr: Seq[int|float]
    peaks_arr: Seq[int|float] | None = []

T = TypeVar('T')
class RES(BaseModel, Generic[T]):
    success: bool = True
    data: Optional[T] = None
    error: str = ''

def is_ndarray(val: Any) -> TypeGuard[ndarray]:
    return isinstance(val, ndarray)