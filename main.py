import uvicorn
import asyncio

from fastapi import FastAPI, APIRouter, Body, Response, status
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from pydantic import BaseModel
from typing import Any, Dict, Annotated
from joblib import Parallel,delayed
from rdkit.Chem import rdDepictor # type: ignore

from dependencies import enumData,enum_atoms_smiles,layout_cdxml,smi2cdxml,create_queue,destroy_queue
from db.dataSQL import searchItem,db_operate,query,get_mol_by_id,molReturn
from unimol.api import dataUnimol,inference,RES
import logger_config
from logging import Logger

log: Logger = logger_config.get_logger(__name__)
rdDepictor.SetPreferCoordGen(True)

router = APIRouter()

semaphore = asyncio.Semaphore(1)

class layout(BaseModel):
  cols:int=10
  rows:int=20
  scale:float=1.0
  singlePage:bool=False

class smi2cdx(BaseModel):
  smiles:list[str]=[]
  layout:layout

class itemhasSmiles(searchItem):
  smiles:str

@router.post("/enum")
async def enum_molecule(enumData: enumData = Body()) -> dict[str,Any]:
	smis=await enum_atoms_smiles(enumData)
	return {
    'data':smis
  }

@router.post("/cdxml")
async def convert2cdxml(data: smi2cdx = Body()) -> dict[str, list[Any]]:
  cdxmls=Parallel(n_jobs=10)(delayed(smi2cdxml)(smi) for smi in data.smiles)
  slides=layout_cdxml(cdxmls = cdxmls, # type: ignore
                        cols = data.layout.cols,
                        rows = data.layout.rows,
                       scale = data.layout.scale,
                 single_page = data.layout.singlePage)
  return {'data':slides}

@router.post("/search")
async def search(data: list[searchItem] = Body()) -> dict[str, list[ dict[ str, int | str ] | None ]] | Any: 
  ## print(data,query(data))
  result: list[Dict[str, int | str] | None] = db_operate(query(data))
  return {'data':result}

@router.post("/molDetail", status_code=200)
async def getMolById(data: dict[str,int] = Body()) -> dict[str,molReturn| Any] :
  result = get_mol_by_id(data['data'])
  return {'data':result}

@router.post("/unimol", status_code=200)
async def predictor(*, data: Annotated[dataUnimol, Body()], response: Response) -> RES[dataUnimol]:
    res_result: RES[dataUnimol] = await asyncio.to_thread(inference, data)
    response.status_code = status.HTTP_200_OK if res_result.success else status.HTTP_400_BAD_REQUEST
    return res_result

@asynccontextmanager
async def lifespan(app: FastAPI):
    # on_startup event
    await create_queue(app)

    yield
    # on_shutdown event
    await destroy_queue(app)

app = FastAPI(lifespan=lifespan)
app.add_middleware(
  CORSMiddleware,
  allow_origins=["*"],
  allow_credentials=True,
  allow_methods=["*"],
  allow_headers=["*"],
)
app.include_router(router)


if __name__ == "__main__":
  uvicorn.run(app="main:app",
              host="0.0.0.0",
              port=5050, 
              reload=True)
