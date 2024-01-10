import os,json,requests,asyncio,aiohttp,logging    # noqa: F811, E401
from typing import List,Union
from datetime import timedelta
from prefect import flow,task,Task
from prefect.client.schemas import TaskRun,FlowRun,State 
from prefect.task_runners import SequentialTaskRunner
from prefect.tasks import task_input_hash
import psycopg as pg
from db import get_login_str,cur_process,get_molData

def task_hook(task: Task, task_run: TaskRun, state: State) -> None:
    logging.error("task error occurred: %s  data: %s  state: %s", task,task_run,state)

def flow_hook(flow, flow_run: FlowRun, state: State) -> None:
    logging.error("task error occurred: %s  data: %s  state: %s", flow,flow_run,state)

@task(retries=3,retry_delay_seconds=5,cache_key_fn=task_input_hash,cache_expiration=timedelta(days=0.05),on_failure=[task_hook])
def get_cookie() -> str:
    with open(os.path.join(os.path.dirname(__file__), '.access.json'), 'r') as f:
        access = json.load(f)
    username = access['cookie']['username']
    password =  access['cookie']['password']
    url='http://192.168.0.238/sunera-dbsys/M0001_01'
    headers = {
        'Accept': '*/*;q=0.8',
        'Accept-Language': 'zh-CN,zh;q=0.9',
        'Cache-Control': 'no-cache',
        'Pragma': 'no-cache',
        'Proxy-Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/112.0.0.0 Safari/537.36'
    }
    response = requests.get(url + '/init', headers=headers, verify=False)
    cookie = response.headers['Set-Cookie'].split(";")[0]
    headers['Content-Type'] = 'application/x-www-form-urlencoded'
    headers['Cookie'] = cookie
    headers['Referer']= 'http://192.168.0.238/sunera-dbsys/A0001_01/init?init=true'
    data = {
      'regionUser': username,
      'password': password
    }
    requests.post(url + '/regist', headers=headers, data=data, verify=False)
    return cookie

async def get_counter_name(calc_name:str,mat_name:str,cookie:str)->Union[List[str],None]:
    url = 'http://192.168.0.238/sunera-dbsys/A0001_01/search'
    data={
      'length': '50',
      'start': '0',
      'draw': '1',
      'jieGouMingCheng': '',
      'gongNengCeng': '',
      'leiXing': '',
      'caiLiaoMingCheng': '',
      'jiSuanZhuangTai': '',
      'smiles': '',
      'zhuanLiJieLun': '',
      'jieGouMingChengList': '',
      'search': 'true',
    }
    if calc_name=='':
        data['jieGouMingCheng']=mat_name
    if mat_name=='':
        data['caiLiaoMingCheng']=calc_name
    if calc_name=='' and mat_name=='':
        print('no name provide!')
        return None
    headers = {
      'Accept': 'application/json, text/javascript, */*; q=0.01',
      'Accept-Language': 'zh-CN,zh;q=0.9',
      'Cache-Control': 'no-cache',
      'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8',
      'Cookie': cookie,
      'Origin': 'http://192.168.0.238',
      'Pragma': 'no-cache',
      'Proxy-Connection': 'keep-alive',
      'Referer': 'http://192.168.0.238/sunera-dbsys/A0001_01/init?init=true',
      'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/112.0.0.0 Safari/537.36',
      'X-Requested-With': 'XMLHttpRequest',
    }
    try:
      async with aiohttp.ClientSession() as session:
          async with session.post(url, data=data, headers=headers,timeout=10) as response:
              res_data = await response.json()
              mat_name = res_data['data'][0]['caiLiaoMingCheng']
              calc_name = res_data['data'][0]['jieGouMingCheng']
              return [mat_name,calc_name] if (mat_name !='')&(calc_name !='') else None
    except Exception as e:
        logging.error("error occurred: %s.  data: %s", e, data)

@flow(task_runner=SequentialTaskRunner(),on_completion=[flow_hook])
async def get_name_pair(names:List[str])->List[str]:#->List[[mat_name,calc_name]]
    cookie = get_cookie()
    tasks = [asyncio.create_task(get_counter_name(calc_name='',mat_name=name,cookie=cookie)) for name in names]
    results = await asyncio.gather(*tasks)
    pair=[x for x in results if x is not None]
    return pair

@task(retries=3,retry_delay_seconds=30,on_failure=[task_hook])
def db_process(file_p:str,cur) -> str:
  try:
    molData	=	get_molData(file_p)
    try:
      cur_process(molData,cur)
    except Exception as e:
      logging.info("%s: sql operate error --(%s)", file_p, e)
  except Exception as e:
    logging.info("%s: parse moldata error --(%s)", file_p, e)
  return file_p+'db process done'

@flow(task_runner=SequentialTaskRunner(),on_failure=[flow_hook])
def batch_db_process(path_files:List[str]) -> None:
  loggin_str=get_login_str(file_p='.access.json')
  with pg.connect(loggin_str) as conn:
    with conn.cursor() as cur:
      for file_p in path_files:
        db_process(file_p,cur)

@flow(task_runner=SequentialTaskRunner(),on_failure=[flow_hook])
def single_db_process(file_p:str) -> None:
  loggin_str=get_login_str(file_p='.access.json')
  with pg.connect(loggin_str) as conn:
    with conn.cursor() as cur:
        db_process(file_p,cur)



