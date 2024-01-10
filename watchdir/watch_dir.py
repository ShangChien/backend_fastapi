import time
import logging
from pathlib import Path
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler
from collections import OrderedDict
from typing import List,Any,Union
from data_api import single_db_process

class Caches:
    def __init__(self, kv:List[Any]=[], max_size:int=10) -> None:
        self._max_size = max_size
        self._items = OrderedDict()
        if kv:
            for i in kv:
                self.add(i[0], i[1])

    @property
    def items(self) -> OrderedDict:
        return self._items

    @property
    def max_size(self) -> int:
        return self._max_size
        
    @max_size.setter
    def max_size(self,max_size:int) -> None:
        if max_size <= 0:
            print("max_size must be positive")
            return
        self._max_size = max_size
        if max_size < len(self._items):
            self._items=self._items[-max_size:]

    def __getitem__(self, index) ->Any:
        if isinstance(index, slice):
            start, stop, step = index.start, index.stop, index.step
            return list(self._items.items())[start:stop:step]
        else:
            return list(self._items.items())[index]
            
    def add(self, key:Union[int,str], value:Any) -> None:
        if key in self._items:
            logging.warning(f"{key} already exists in the container")
            return
        if len(self._items) >= self.max_size:
            _, _history = self._items.popitem(last=False)
            with open('./.history.log', 'a') as f:
                f.write(f'{_history}\n')
        self._items[key] = value

    def update_kv(self, key:Union[int,str], value:Any) -> None:
        if key not in self._items:
            logging.warning(f"{key} doesn't exist in the container")
            return
        self._items[key] = value

    def rm(self, key) -> None:
        if key not in self._items:
            logging.warning(f"{key} doesn't exist in the container")
            return
        del self._items[key]

class subHandler(FileSystemEventHandler):
    def on_created(self, event) -> None:
        if event.is_directory:
            pass
        else:
            path=Path(event.src_path)
            if path.suffix == '.log':
                print('Sub:',path.as_posix())
                try:
                    single_db_process(path.as_posix())
                    logging.info(f'success to single_db_process {path.as_posix()}')
                except Exception:
                    logging.error(f'fail to single_db_process {path.as_posix()}')

def update_subqueue(sub_queue:Caches) -> None:
    folders = [p for p in dir_path.glob('*') if p.is_dir()]
    last_updated_dir=sorted(folders, key=lambda p: p.stat().st_mtime, reverse=True)[:3]
    for dir in last_updated_dir:
        str_dir = dir.as_posix()
        if str_dir in sub_queue.items:
            continue
        else:
            ##stop & join older observer if queue is full
            if sub_queue.max_size == len(sub_queue.items):
                sub_queue[0][1].stop()
                sub_queue[0][1].join(30)
            ##add observer(auto delete older observer)
            print("add observer", str_dir)
            sub_handler = subHandler()
            ##better way is Observer()
            observer = PollingObserver(timeout=300)
            observer.schedule(sub_handler, str_dir, recursive=True)
            sub_queue.add(str_dir,observer)
            observer.start()

class mainHandler(FileSystemEventHandler):
    def on_created(self, event) -> None:
        if event.is_directory:
            logging.info(f"main hook Directory modied: {event.src_path}")
            update_subqueue(sub_queue)

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename='error.log', filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')
    dir_path =   Path('Z:/Data-month-statistics')
    sub_queue  =  Caches([],3)
    update_subqueue(sub_queue)
    main_handler = mainHandler()
    ##better way is Observer()
    main_observer = PollingObserver(timeout=300)
    main_observer.schedule(main_handler, dir_path.as_posix(), recursive=False)
    main_observer.start()
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        for _key,observer in list(sub_queue.items.items()):
            try:
                observer.stop()
                observer.join()
            except Exception:
                logging.error(f'fail after activate ctrl+c stop & join sub_observer: key-{_key}; value-{observer}')
        main_observer.stop()
        main_observer.join()




