from db.db import get_login_str,pg,MolData
from pydantic import BaseModel
from psycopg import sql
from psycopg.rows import dict_row
from psycopg.sql import Literal
from typing import TypeGuard, Union,List,Dict,Any
from enum import Enum

class Logic(str, Enum):
    AND = "and"
    OR = "or"
    NOT = "not"

class searchItem(BaseModel):
    key: str
    logic: str
    value: str| List[Union[float,str]]

class molReturn(MolData):
    id:int|None=None
    bfp:str|None=None
    sfp:str|None=None
    timestamp:str|None=None
    mole_weight:float|None=None

def is_str(val: Any) -> TypeGuard[str]:
    return isinstance(val, str)

def is_str_list(val: Any) -> TypeGuard[list[str]]:
    return all(isinstance(x, str) for x in val)

def convert2SQLstr(item:searchItem) -> str:
  lower:float = 0.0
  upper:float = 0.0
  if isinstance(item.value, List) & isinstance(item.value[0], float) & isinstance(item.value[1], float) :
    item_value: list[float] = [float(i) for i in item.value[:2]]
    lower = item_value[0] if item_value[0] <= item_value[1] else item_value[1]
    upper = item_value[1] if item_value[0] <= item_value[1] else item_value[0]
  itemStr=''
  if (item.key=='substructure'):
    operate = ''
    if item.value[0]=='>':
      operate = '@>'
    elif item.value[0]=='>':
      operate = '<@'
    elif item.value[0]=='=':
      operate = '@='
    itemStr = f"(m{operate}'{item.value[1]}')"
  elif(item.key=='similarity'):
    itemStr = f'(tanimoto_sml(morganbv_fp(mol_from_smiles({item.value[2]}),2), bfp) BETWEEN {lower} AND {upper})'
  elif (item.key=='mw'):
    itemStr = f'(mol_amw(m) BETWEEN {lower} AND {upper})'
  elif (item.key=='date'):
    itemStr = f'(timestamp BETWEEN to_timestamp({lower/1000}) AND to_timestamp({upper/1000}))'
  elif (item.key in ['name_mat', 'name_calc']) and is_str(item.value):
    valueList = [i[1:] if i.startswith('/') else i for i in item.value.split(' ')]
    valueStr = "','".join(valueList)
    itemStr = f"({item.key} IN ('{valueStr}'))"
  elif (item.key in ['types', 'labels']) and is_str_list(item.value):
    valueStr = "','".join(item.value)
    itemStr = f"(ANY({item.key}) IN ('{valueStr}'))"
  else:
    pass
  return itemStr

def query(data:List[searchItem]) -> str:
  and_arr = []
  or_arr  = []
  not_arr = []
  for i in data:
    i_str = convert2SQLstr(i)
    if i.logic == 'and':
      and_arr.append(i_str)
    elif i.logic == 'not':
      not_arr.append(i_str)
    elif i.logic == 'or':
      or_arr.append(i_str)
  and_str ='(' + (' AND '.join(and_arr))+ ')' if len(and_arr) > 0 else '' 
  or_str  ='OR (' + (' AND '.join(or_arr)) + ')' if len(or_arr) > 0 else '' 
  not_str ='AND NOT (' + (' OR '.join(not_arr))+ ')' if len(not_arr) > 0 else '' 

  sql_str = f'({and_str} {or_str}) {not_str}'
  return sql_str

def db_operate(sql_str:str) -> List[Union[Dict[str,Union[int,str]],None]]:
  loggin_str=get_login_str(file_p='.access.json')
  with pg.connect(loggin_str) as conn:
    with conn.cursor() as cur:
      query_all= f"SELECT id,smiles,name_calc,name_mat FROM moldata.main WHERE {sql_str}"
      cur.execute(query_all) # type: ignore
      all=[{
        'id':i[0],
        'smiles':i[1],
        'name_calc':i[2] if i[2] else 'null',
        'name_mat':i[3] if i[3] else 'null',
      } for i in cur.fetchall()]
      ## print(all[:2])
      return all  # type: ignore
    
def get_mol_by_id(id:int) -> Union[molReturn,None]:
  loggin_str=get_login_str(file_p='.access.json')
  with pg.connect(loggin_str) as conn:
    with conn.cursor(row_factory=dict_row) as cur:
      query_all= f"SELECT * FROM moldata.main WHERE id={id}"
      cur.execute(query_all) # type: ignore
      item=cur.fetchone()
      if item is not None and item["timestamp"] is not None:
        item['timestamp']=item['timestamp'].strftime("%Y-%m-%d %H:%M:%S")
        del item['m']
      return item  # type: ignore

