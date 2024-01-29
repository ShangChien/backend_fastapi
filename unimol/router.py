import asyncio
from fastapi import  APIRouter, Body, Response, status
from unimol.api import dataUnimol,inference,RES
from typing import  Annotated

router = APIRouter(
    prefix="/ml",
    tags=["machine_learning"],
)

@router.post("/unimol", status_code=200)
async def predictor(*, data: Annotated[dataUnimol, Body()], response: Response) -> RES[dataUnimol]:
    res_result: RES[dataUnimol] = await asyncio.to_thread(inference, data)
    response.status_code = status.HTTP_200_OK if res_result.success else status.HTTP_400_BAD_REQUEST
    return res_result