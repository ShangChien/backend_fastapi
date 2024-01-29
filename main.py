import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
from dataCheck.router import router as data_check_router
from enumSmiles.router import router as enum_router
from db.router import router as db_router
from unimol.router import router as unimol_router

from dependencies import create_queue,destroy_queue
import logger_config
from logging import Logger

log: Logger = logger_config.get_logger(__name__)

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
app.include_router(unimol_router)
app.include_router(data_check_router)
app.include_router(enum_router)
app.include_router(db_router)


if __name__ == "__main__":
  uvicorn.run(app="main:app",
              host="0.0.0.0",
              port=5050, 
              reload=True)
