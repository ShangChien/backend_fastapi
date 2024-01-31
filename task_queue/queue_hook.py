import asyncio
import pickle
from datetime import datetime
from asyncio import PriorityQueue
from pydantic import BaseModel, Field
from uuid import uuid4
from dataclasses import dataclass,field
from fastapi import FastAPI
from typing import Any,Callable
import logger_config
from logging import Logger



log: Logger = logger_config.get_logger(__name__)

class Task(BaseModel):
    data: Any = None
    func: str|None = None
    id: str = Field(default_factory = lambda: uuid4().hex)
    time: float = Field(default_factory= lambda: datetime.now().timestamp())
    prior: list[int] = [5,5,5,5]
    result: Any = None
    
    def priorValue(self) -> int:
        val=int(''.join([str(i) for i in self.prior]))
        return val
    
    def __lt__(self, other) -> bool:
        return self.priorValue() < other.priorValue()
    def __le__(self, other) -> bool:
        return self.priorValue() <= other.priorValue()
    def __gt__(self, other) -> bool:
        return self.priorValue() > other.priorValue()
    def __ge__(self, other) -> bool:
        return self.priorValue() >= other.priorValue()

@dataclass
class taskQueue:
    lock: asyncio.Lock = field(default_factory=asyncio.Lock)
    queue: PriorityQueue[Task] = field(default_factory=PriorityQueue)
    running: dict[str, Task] = field(default_factory=dict)
    waiting: dict[str, Task] = field(default_factory=dict)
    tobeCancel: dict[str, Task] = field(default_factory=dict)
    canceled: dict[str, Task] = field(default_factory=dict)
    completed: dict[str, Task] = field(default_factory=dict)

    async def add(self, task: Task) -> None:
        async with self.lock:
            self.waiting[task.id] = task
        await self.queue.put(task)

    async def cancel(self, task: Task) -> bool:
        if task.id in self.waiting:
            async with self.lock:
                self.tobeCancel[task.id] = task
            return True
        if task.id in self.running:
            await self.task_done(task)
            async with self.lock:
                self.canceled[task.id] = task
            return True
        else:
            return False
        
    async def modify(self, task: Task) -> bool:
        canceled: bool = await self.cancel(task)
        if canceled:
            await self.add(task)
            return True
        else:
            return False
    
    async def get(self) -> Task | None:
        task: Task = await self.queue.get()
        async with self.lock:
            del self.waiting[task.id]
            if task.id in self.tobeCancel:
                del self.tobeCancel[task.id]
                self.canceled[task.id] = task
                return None
            else:
                self.running[task.id] = task
                return task
        
    async def task_done(self, task: Task) -> bool:
        async with self.lock:
            if task.id in self.running:
                del self.running[task.id]
                self.completed[task.id] = task
                self.queue.task_done()
                return True
            else:
                return False

    async def join(self) -> None:
        await self.queue.join()

funcs: dict[str, Callable] = {'uv': sum}
async def worker(queue: taskQueue):
    while True:
        task: Task | None = await queue.get()
        if (task is None) or (task.func is None):
            pass
        else:
            result: Any = funcs[task.func](task.data)
            task.result = result
            await queue.task_done(task)

async def create_queue(app: FastAPI):
    # 程序启动时执行的钩子函数:
    # 1. 载入uv数据
    with open('./uv_data.pkl','rb') as f:
        app.state.uv_data = pickle.load(f)

    # 2. 启动任务队列 
    queue: taskQueue = taskQueue()
    app.state.queue = queue
    app.state.workers = [asyncio.create_task(worker(queue)) for _i in range(1)]

async def destroy_queue(app: FastAPI):
    # 程序结束时执行的钩子函数
    # 1.持久化uv_data.pkl
    with open('./uv_data.pkl','wb') as f:
        pickle.dump(app.state.uv_data, f)

    # 2.关闭任务队列
    workers:list[asyncio.Task[Any]] = app.state.workers
    order_queue: taskQueue = app.state.queue
    print("等待尚未完成的任务执行完毕, 但只有 60 秒的机会")
    try:
        await asyncio.wait_for(order_queue.join(), timeout=60)
    except asyncio.TimeoutError:
        print("程序结束, 但还有任务尚未完成, 这里直接取消")
    finally:
        [worker.cancel() for worker in workers]