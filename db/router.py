from fastapi import APIRouter, Body
from .dataSQL import searchItem,db_operate,query,get_mol_by_id,molReturn
from typing import Any

class itemhasSmiles(searchItem):
    smiles:str

router = APIRouter(
    prefix="/db",
    tags=["dataBase"],
)

@router.post("/search")
async def search(data: list[searchItem] = Body()) -> dict[str, list[ dict[ str, int | str ] | None ]] | Any:
    result: list[dict[str, int | str] | None] = db_operate(query(data))
    return {'data':result}

@router.post("/molDetail", status_code=200)
async def getMolById(data: dict[str,int] = Body()) -> dict[str,molReturn| Any] :
    result = get_mol_by_id(data['data'])
    return {'data':result}