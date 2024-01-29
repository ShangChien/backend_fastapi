from fastapi import APIRouter, Body
from .utils import enumData,enum_atoms_smiles,layout_cdxml,smi2cdxml,smi2cdx
from joblib import Parallel, delayed
from typing import Any


router = APIRouter(
    prefix="/enum_smiles",
    tags=["enum_smiles"],
)

@router.post("/enum_by_atoms")
async def enum_molecule(enumData: enumData = Body()) -> dict[str,Any]:
	smis: list[str]=await enum_atoms_smiles(enumData)
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