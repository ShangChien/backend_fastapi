import asyncio

async def my_coroutine(lock, id):
    print(f"Coroutine {id} has started.")
    await asyncio.sleep(5)
    async with lock:
        # 在这个区块中，同时只有一个协程可以执行
        print(f"Coroutine {id} has acquired the lock.")
        await asyncio.sleep(1)  # 模拟IO操作
        print(f"Coroutine {id} is releasing the lock.")

async def main():
    # 创建一个异步锁
    lock = asyncio.Lock()

    # 启动三个协程，它们都会尝试获取锁
    await asyncio.gather(
        my_coroutine(lock, "A"),
        my_coroutine(lock, "B"),
        my_coroutine(lock, "C"),
    )

# 运行主函数
asyncio.run(main())
