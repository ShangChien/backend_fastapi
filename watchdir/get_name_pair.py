from typing import List,Union
import asyncio
import aiohttp
import psycopg as pg
import logging

def get_cookie(username, password,url='http://192.168.0.238/sunera-dbsys/M0001_01'):
	import requests
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

async def get_material_name(
	calculate_name:str,
	cookie:str
)->Union[List[str],None]:
	url = 'http://192.168.0.238/sunera-dbsys/A0001_01/search'
	data = {
    'length': '50',
    'start': '0',
    'draw': '1',
    'jieGouMingCheng': calculate_name,
    'gongNengCeng': '',
    'leiXing': '',
    'caiLiaoMingCheng': '',
    'jiSuanZhuangTai': '',
    'smiles': '',
    'zhuanLiJieLun': '',
    'jieGouMingChengList': '',
    'search': 'true',
  }
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
				material_name = res_data['data'][0]['caiLiaoMingCheng']
				return [calculate_name, material_name] if material_name !='' else None
	except Exception as e:
		logging.error("error occurred: %s.  data: %s", e, data['jieGouMingCheng'])

async def main(cal_names:List[str],cookie:str)->List[str]:
	tasks = [asyncio.create_task(get_material_name(cal_name,cookie)) for cal_name in cal_names]
	results = await asyncio.gather(*tasks)
	pair=[x for x in results if x is not None]
	return pair

if __name__ == '__main__':
	logging.basicConfig(level=logging.DEBUG, filename='error.log', filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')
	cookie=get_cookie('chensq','Sunera@2013')
	raws=[]
	with pg.connect("dbname=emolecules user=postgres host=192.168.2.233 port=5432 connect_timeout=10 	password=sunera") as conn:
		with conn.cursor() as cur:
			cur.execute("SELECT name_calc FROM moldata.main WHERE smiles is not Null")
			raws=cur.fetchall()
	cal_names=[s[0] for s in raws]
	pairs=[]
	slice_cal_names=[]
	for i in range(0,len(cal_names),100):
		try:
			slice_cal_names = cal_names[i:i+100]
		except IndexError:
			slice_cal_names = cal_names[i:]
		name_pair=asyncio.run(main(slice_cal_names,cookie))
		pairs.extend(name_pair)
	with open('name_pair.txt','w') as f:
		for pair in pairs:
			f.write(':'.join(pair)+'\n')