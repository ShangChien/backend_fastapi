from fastapi import APIRouter, Body, Request
from queue_hook import taskQueue, Task
from fastapi.responses import JSONResponse
from typing import TypeVar,Optional,Generic,Any
from pydantic import BaseModel
from datetime import datetime, timedelta

T = TypeVar('T')
class RES(JSONResponse, Generic[T]):
    success: bool = True
    content: Optional[T] = None
    error: str= ''

class queueInfo(BaseModel):
    running: list[Any] = []
    waiting: list[Any] = []
    canceled: list[Any] = []
    completed: list[Any] = []

router = APIRouter(
    prefix="/queue",
    tags=["queue"],
)

@router.post("/add" )
async def addTask(req: Request, data: dict[str,int] = Body()) -> RES[str]:
    queue:taskQueue = req.app.state.queue
    task=Task(data=data)
    await queue.add(task)
    return RES[str](content = f'task {task.id} added')

@router.post("/cancel" )
async def cancelTask(req: Request, task_id:str  = Body()) -> RES[str]:
    queue:taskQueue = req.app.state.queue
    res:str =''
    if task_id in queue.waiting:
        task = queue.waiting[task_id]
        await queue.cancel(task)
        res = f'task {task.id} canceled from waiting'
    else:
        res = f'task {task_id} is not support cancel'
    return RES[str](content = res)

@router.post("/modify" )
async def modifyTask(req: Request, data: Task = Body()) -> RES[str]:
    queue:taskQueue = req.app.state.queue
    res:str =''
    if data.id in queue.waiting:
        task = queue.waiting[data.id]
        task.prior = data.prior
        await queue.modify(task)
        res = f'task {task.id} has modified: canceled from waiting, then add a new'
    else:
        res = f'task {data.id} is not support modified'
    return RES[str](content = res)

pre_time = datetime.now()
status_res:Any = None
@router.get("/status" )
async def getStatus(req: Request) -> RES[queueInfo] :
    cur_time = datetime.now()
    global pre_time
    global status_res
    if cur_time - pre_time < timedelta(seconds=5):
        return status_res
    pre_time = cur_time
    queue:taskQueue = req.app.state.queue
    queue_info=queueInfo(
        running = [queue.running[k] for k in queue.running],
        waiting = [queue.waiting[k] for k in queue.waiting ],
        canceled = [queue.canceled[k] for k in queue.canceled],
        completed = [queue.completed[k] for k in queue.completed]
    )
    status_res = RES[queueInfo](content=queue_info.model_dump_json())
    return status_res