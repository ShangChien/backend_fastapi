from fastapi import APIRouter
from .spectrum.router import router as spectrum_router

router = APIRouter(
    prefix="/data_check",
    tags=["data_check"],
)

router.include_router(spectrum_router)
