import re
import json
from datetime import datetime,timedelta
from pathlib import Path
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.electronic_structure.core import Spin
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds 	# type: ignore
import psycopg as pg
from psycopg.sql import SQL, Identifier, Placeholder
from typing import Union,List,Any
from pydantic import BaseModel
import logger_config
from logging import Logger

log: Logger = logger_config.get_logger(__name__)

class DirPath(BaseModel):
	ground:str|None = None
	udft:str|None = None
	tdft:str|None = None
	polar:str|None = None
	uv:str|None = None
	
class Corrections(BaseModel):
	zeroPoint:float|None=None
	energy: float|None=None
	enthalpy: float|None=None
	gibbsFreeEnergy: float|None=None

class ElectronCrossing(BaseModel):
	energy:list[float]|None=None
	oscillator:list[float]|None=None

class Structures(BaseModel):
	ground:str|None = None
	udft:str|None = None
	tdft:str|None = None

class Polar(BaseModel):
	iso:float|None=None
	aniso:float|None=None
	polarizability:list[float]|None=None

class TDFT(BaseModel):
	s0_s1:ElectronCrossing|None=None
	s0_t1:ElectronCrossing|None=None

class UDFT(BaseModel):
	energy			: float|None=None
	corrections	:	Corrections|None=None

class Spectrum(BaseModel):
	uv:ElectronCrossing|None=None
	uv_low:ElectronCrossing|None=None

class MolData(BaseModel):
	timeStamp					:	datetime|None=None
	name_calc					:	str|None = None
	smiles						:	str|None = None
	functional					:	str|None = None
	basisSet					:	str|None = None
	charge						:	int|None = None
	spinMultiplicity			:	int|None = None
	numBasisFunc				:	int|None = None
	electrons					:	list[int]|None = None
	isPcm						:	bool|None = None
	isSpin						:	bool|None = None
	stationaryType				:	str|None = None
	extraInfo					:	str|None = None
	energy						:	float|None=None
	homo						:	float|None=None
	lumo						:	float|None=None
	eg							:	float|None=None
	polar						:	Polar|None=None
	dirPath						:	DirPath|None=None
	corrections					: 	Corrections|None=None
	structures					:	Structures|None=None
	tdft						:	TDFT|None=None
	udft						:	UDFT|None=None
	spectrum					:	Spectrum|None=None

def get_timestrp(filePath: str) -> Union[datetime,None]:
	with open(filePath,'r') as file:
		for line in file:
			if "Leave Link" in line :
				s = line.split()
				s[8] = s[8][:-1]
				time_str = ' '.join(s[4:9])
				strptime = datetime.strptime(time_str, "%a %b %d %H:%M:%S %Y")
				return strptime
			elif 'Normal termination' in line:
				s = line.split()
				s[10] = s[10][:-1]
				time_str = ' '.join(s[6:11])
				strptime = datetime.strptime(time_str, "%a %b %d %H:%M:%S %Y")
				return strptime
	log.info("%s: Timestrp data not found", filePath)
	raise Exception("Timestrp data no found")
	
def xyz2mol(xyz_str:str) -> Any:
	raw_mol = Chem.MolFromXYZBlock(xyz_str)  # type: ignore
	conn_mol = Chem.Mol(raw_mol)  # type: ignore
	try:
		rdDetermineBonds.DetermineBonds(conn_mol,charge=0)
	except Exception:
		rdDetermineBonds.DetermineConnectivity(conn_mol)
	return conn_mol

def get_corrections(filePath: str) -> Union[Corrections,None]:
	with open(filePath,'r') as file:
		lines=file.readlines()
		for n in range(len(lines)-1,0,-1):
			if 'Zero-point correction' in lines[n]:
				corrections	= Corrections(
					zeroPoint				= float(lines[n].split()[-2]),
					energy 					= float(lines[n+1].split()[-1]),
					enthalpy 				= float(lines[n+2].split()[-1]),
					gibbsFreeEnergy	= float(lines[n+3].split()[-1])
				)
				return corrections
	log.info("%s: Corrections data not found", filePath)
	raise Exception("Corrections data no found")

def get_crossing(filePath: str, target:str='Triplet') -> Union[ElectronCrossing,None]:
	with open(filePath,'r') as file:
		for line in file:
			if target in line:
				strs=line.split()
				crossing = ElectronCrossing(
					energy			=	[float(strs[4])],
					oscillator	=	[float(strs[-2][2:])]
				)
				return crossing
	log.info("%s: ElectronCrossing data not found", filePath)
	raise Exception("ElectronCrossing data no found,target:"+target)

def get_spectrum(filePath: str, target:str='Singlet') -> Union[ElectronCrossing,None]:
	with open(filePath,'r') as file:
		energy			=	[]
		oscillator	=	[]
		for line in file:
			if target in line:
				strs=line.split()
				e=float(strs[4])
				o=float(strs[-2][2:])
				energy.append(e)
				oscillator.append(o)
	if energy==[] or oscillator == []:
		log.info("%s: Spectrum data not found", filePath)
		raise Exception("Spectrum data not found,target:"+target)
	EC = ElectronCrossing(
		energy			=	energy,
		oscillator	=	oscillator
	)
	return EC

def get_polar(filePath:str) -> Union[Polar,None]:
	with open(filePath,'r') as f:
		for line in f:
			if "Exact polarizability" in line:
				float_str=line.split(':')[1]
				xx=float(float_str[:8])				
				yx=float(float_str[8:16])
				yy=float(float_str[16:24])
				zx=float(float_str[24:32])
				zy=float(float_str[32:40])
				zz=float(float_str[40:48])
				xyz=[xx,yx,yy,zx,zy,zz]
				iso	=	(xx+yy+zz)/3
				aniso	=	(((xx-yy)**2+(xx-zz)**2+(zz-yy)**2+6*(yx**2+zy**2+zx**2))/2)**(0.5)
				polar=Polar(
					iso 	= round(iso,3),
					aniso = round(aniso,3),
					polarizability	=	xyz
				)
				return polar
	log.info("%s: Polar data not found", filePath)
	raise Exception("Polar data not found")

def mol2smiles(mol:Any,filePath:str) -> Union[str,None]:
	try:
		smiles:str	=	Chem.CanonSmiles(Chem.MolToSmiles(mol))	# type: ignore
		return smiles
	except Exception:
		log.info("%s: parse smiles failed",filePath)
		raise Exception("parse smiles failed")

def get_name(name:str,target:str)	->	str:
	out=re.sub(target,'',name)
	return out if out[-1]!='-' else out[:-1]

def get_molData(filePath:str)	->	MolData:
	moldata			=	MolData(**{})
	path				=	Path(filePath)
	filename		=	path.stem
	upper_name	= filename.upper()
	gaussian  	= GaussianOutput(path.as_posix())
	if gaussian.properly_terminated:
		xyz_str					=	gaussian.structures[-1].to(fmt="xyz", filename=None)
		mol_rdkit				=	xyz2mol(xyz_str)
		m_Block					=	Chem.MolToMolBlock(mol_rdkit)	# type: ignore
		moldata.smiles	= mol2smiles(mol=mol_rdkit,filePath=filePath)
		if 'UDFT' in upper_name:
			moldata.name_calc = get_name(upper_name,'UDFT')
			moldata.dirPath				=	DirPath(udft=path.as_posix()[3:], **{})
			moldata.structures		=	Structures(udft=m_Block, **{})
			moldata.udft					=	UDFT(
				corrections	=	get_corrections(filePath),
				energy			=	gaussian.final_energy
			)
		elif 'S1' in upper_name or 'UV' in upper_name:
			if 'T1' in upper_name:#标准tdft
				name	= get_name(upper_name,'T1')
				name	= get_name(name,'S1')
				moldata.name_calc =	name
				moldata.dirPath				=	DirPath(tdft=path.as_posix()[3:], **{})
				moldata.spectrum			=	Spectrum(uv_low=get_spectrum(filePath, target='Singlet'),**{})
				moldata.tdft 					= TDFT(
					s0_s1=get_crossing(filePath,target='Singlet'),
					s0_t1=get_crossing(filePath,target='Triplet'),
				)
			else:#仅包含S1
				moldata.name_calc = get_name(upper_name,'S1') if 'S1' in upper_name else get_name(upper_name,'UV')
				moldata.dirpath				=	DirPath(uv=path.as_posix()[3:], **{})
				moldata.tdft 					= TDFT(s0_s1=get_crossing(filePath,target='Singlet'), **{})
				moldata.spectrum			=	Spectrum(uv=get_spectrum(filePath, target='Singlet'), **{})
		elif 'POLAR' in upper_name:
			moldata.name_calc = get_name(upper_name,'POLAR')
			moldata.dirPath				=	DirPath(polar=path.as_posix()[3:],**{})
			moldata.polar					= get_polar(filePath)
		else: ## get ground calculation info
			moldata.name_calc = upper_name	
			moldata.timestamp			=	get_timestrp(filePath)
			moldata.basisset			=	gaussian.basis_set
			moldata.charge				=	gaussian.charge
			moldata.functional		=	gaussian.functional
			moldata.electrons			=	list(gaussian.electrons)
			moldata.ispcm					=	gaussian.is_pcm
			moldata.isspin				=	gaussian.is_spin
			moldata.energy				= gaussian.final_energy
			moldata.homo 					= gaussian.eigenvalues[Spin.up][gaussian.electrons[0]-1] # type: ignore
			moldata.lumo 					= gaussian.eigenvalues[Spin.up][gaussian.electrons[0]] # type: ignore
			moldata.eg						= moldata.lumo - moldata.homo
			moldata.dirpath				=	DirPath(ground=path.as_posix()[3:], **{})
			moldata.structures		= Structures(ground=m_Block, **{})
			moldata.corrections		=	get_corrections(filePath)
			moldata.numbasisfunc	=	gaussian.num_basis_func
			moldata.spinmultiplicity	=	gaussian.spin_multiplicity
			moldata.stationarytype		=	gaussian.stationary_type
	else:
		moldata.extrainfo	=	"error: "+str(gaussian.errors)
		raise Exception("job error")
	return moldata

def cur_process(moldata:MolData,cur:pg.Cursor[Any])->Any:
	name			= moldata.name_calc
	dirPath		=	moldata.dirpath.model_dump() # type: ignore
	fileType 	= [k for k in dirPath if dirPath[k] is not None][0]
	mol_without_none	=	{k: v for k, v in moldata.model_dump().items() if v is not None}
	mol_dict 					= {k: json.dumps(v) if isinstance(v, dict) else v for k, v in mol_without_none.items()}
	cur.execute("SELECT COUNT(*) FROM moldata.main WHERE name_calc = %s", (name,))
	#如果库中没有直接插入
	if  cur.fetchone()[0] == 0:# type:ignore
		query = SQL("INSERT INTO moldata.main ({}) VALUES ({})").format(
			SQL(', ').join(map(Identifier, [s.lower() for s in mol_dict.keys()])),
  	  SQL(', ').join([Placeholder() for i in range(len(mol_dict))])
		)
		#print('new data',query.as_string(cur))
		cur.execute(query,list(mol_dict.values()))
	##已经存在则 updata data
	else: 
		common_field = {k: v for k, v in mol_without_none.items() if not isinstance(v, dict)}
		query = SQL('UPDATE moldata.main SET {} WHERE name_calc = %s').format(
			SQL(', ').join([SQL('{} = %s').format(Identifier(k.lower())) for k in list(common_field.keys())])
		)	
		#print('update',query.as_string(cur))
		#print(query.as_string(cur))
		cur.execute(query, list(common_field.values()) + [name])  # type: ignore
		if fileType=='ground':
			cur.execute('UPDATE moldata.main SET structures = structures || %s WHERE name_calc = %s',
				[json.dumps({'ground':moldata.structures.ground}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET dirPath = dirPath || %s WHERE name_calc = %s',
				[json.dumps({'ground':moldata.dirpath.ground}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET corrections = corrections || %s WHERE name_calc = %s',
				[json.dumps(moldata.corrections.model_dump()),name]) #type:ignore
		elif fileType=='polar':
			cur.execute('UPDATE moldata.main SET dirPath = dirPath || %s WHERE name_calc = %s',
				[json.dumps({'polar':moldata.dirpath.polar}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET polar = polar || %s WHERE name_calc = %s',
				[json.dumps(moldata.polar.model_dump()),name]) #type:ignore
		elif fileType=='tdft':
			cur.execute('UPDATE moldata.main SET dirPath = dirPath || %s WHERE name_calc = %s',
				[json.dumps({'tdft':moldata.dirpath.tdft}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET tdft = tdft || %s WHERE name_calc = %s',
				[json.dumps(moldata.tdft.model_dump()),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET spectrum = spectrum || %s WHERE name_calc = %s',
				[json.dumps({'uv_low':moldata.spectrum.uv_low.model_dump()}),name]) #type:ignore
		elif fileType=='uv':
			cur.execute('UPDATE moldata.main SET dirPath = dirPath || %s WHERE name_calc = %s',
				[json.dumps({'tdft':moldata.dirpath.uv}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET spectrum = spectrum || %s WHERE name_calc = %s',
				[json.dumps({'uv':moldata.spectrum.uv.model_dump()}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET tdft = tdft || %s WHERE name_calc = %s',
				[json.dumps({'s0_s1':moldata.tdft.s0_s1.model_dump()}),name]) #type:ignore
		elif fileType=='udft':
			cur.execute('UPDATE moldata.main SET dirPath = dirPath || %s WHERE name_calc = %s',
				[json.dumps({'udft':moldata.dirpath.udft}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET structures = structures || %s WHERE name_calc = %s',
				[json.dumps({'udft':moldata.structures.udft}),name]) #type:ignore
			cur.execute('UPDATE moldata.main SET udft = udft || %s WHERE name_calc = %s',
				[json.dumps(moldata.udft.model_dump()),name]) #type:ignore
		else:
			log.info('unknown,filetype:',fileType,dirPath)

def get_login_str(file_p:str) -> str:
	with open(file_p) as f:
		data = json.load(f)
	login=data['pg']
	login_str = 'dbname={dbname} user={user} host={host} port={port} connect_timeout=10 password={password}'.format(
    dbname=login['dbname'],user=login['user'],host=login['host'],port=login['port'],password=login['password']
	)
	return login_str

def batch_db_process(path_files:List[str]) -> None:
	loggin_str=get_login_str(file_p='.access.json')
	with pg.connect(loggin_str) as conn:
		with conn.cursor() as cur:
			for file_p in path_files:
				try:
					molData	=	get_molData(file_p)
					print(molData)
					try:
						cur_process(molData,cur)
					except Exception:
						log.info("%s: sql operate error", file_p)
				except Exception:
					log.info("%s: parse moldata error", file_p)

if __name__ == "__main__":
	start_time = datetime.now()
	logname=start_time.strftime('%Y-%m-%d-%H-%M-%S')+"_runINFO.log"
	log.info("程序开始时间: %s", start_time.strftime('%Y-%m-%d %H:%M:%S'))
	dir_path = Path('Z:/Data-month-statistics/202305/')
	i=0#统计次数
	loggin_str=get_login_str(file_p='./.access.json')
	#path_files=[p.as_posix() for p in dir_path.glob('**/*.log')]
	#batch_db_process(loggin_str=loggin_str,path_files=path_files)
	with pg.connect(loggin_str) as conn:
		with conn.cursor() as cur:
			for file_p in dir_path.glob('**/*.log'):
				try:
					molData	=	get_molData(file_p.as_posix())
					try:
						log.info("processing: %s", file_p)
						cur_process(molData,cur)
						log.info("finished: %s", file_p)
						i=i+1
					except Exception:
						log.info("%s: sql operate error", file_p.as_posix())
				except Exception:
					log.info("%s: parse moldata error", file_p.as_posix())		
	end_time = datetime.now()
	elapsed_time = end_time - start_time
	elapsed_seconds = int(elapsed_time.total_seconds())
	elapsed_time_str = str(timedelta(seconds=elapsed_seconds))
	log.info("程序结束时间: %s", end_time.strftime('%Y-%m-%d %H:%M:%S'))
	log.info("程序耗时: %s ,插入次数: %s",elapsed_time_str,str(i))