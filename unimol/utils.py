from rdkit.Chem import rdDepictor,rdDistGeom,MolFromMolBlock,MolFromSmiles,AddHs,rdMolDescriptors # type: ignore
from ase.optimize import FIRE
from pathlib import Path
from ase import Atoms
from tblite.ase import TBLite
from typing import Any, Generic, Optional, TypeVar
from pydantic import BaseModel
import logging

logger = logging.getLogger(__name__)

rdDepictor.SetPreferCoordGen(True)
etkdg:Any = rdDistGeom.ETKDGv3()
etkdg.randomSeed = 0xa700f
etkdg.verbose = False
etkdg.numThreads = 0
etkdg.optimizerForceTol = 0.0135
etkdg.useRandomCoords = True

class dataUnimol(BaseModel):
    models: list[str] = []
    names : list[str] = []
    smiles: list[str] = []
    molBlocks: list[str] = []
    atoms : list[list[str]] = [] 
    coordinates : list[list[list[float]]] = []
    results : dict[str, list[int|float|str]] = {} ## { 'model1': [1,2,3] }

T = TypeVar('T')
class RES(BaseModel, Generic[T]):
    success: bool = True
    data: Optional[T] = None
    error: str= ''

def calculator(aseMol:Atoms, method="GFN2-xTB", criteria=0.05, steps=100) -> RES[Atoms]:
    try:
        #如果不存在opt.log文件，则使用pathlib创建该文件
        if not Path('opt.log').exists():
            Path('opt.log').touch()
        aseMol.calc = TBLite(method=method, verbosity=0)
        FIRE(aseMol, logfile='opt.log').run(fmax = criteria,steps = steps)
        return RES[Atoms](data = aseMol)
    except Exception as e:
        logger.error(e)
        return RES[Atoms](success=False, data = aseMol, error = str(e))
    
def xTBconvertor(molBlock:str|None=None, smi:str|None=None) -> RES[dict]: ## 2D(molBlock) to 3D(dict)
    aseList:list[Atoms] = []
    str_atoms:list[str] = []
    coords:list[list[float]] = []
    RES_result = RES[dict]()
    try:
        mol: Any = MolFromMolBlock(molBlock) if molBlock else MolFromSmiles(smi)
        mol: Any = AddHs(mol)
        nRBonds:int = rdMolDescriptors.CalcNumRotatableBonds(mol)
        nRBonds = nRBonds if nRBonds>1 else 1 
        cids: list[int] = rdDistGeom.EmbedMultipleConfs(mol, numConfs = 3*nRBonds, params=etkdg)
        for cid in cids:
            conf: Any = mol.GetConformer(cid)
            str_atoms = [i.GetSymbol() for i in mol.GetAtoms()]
            coords = conf.GetPositions()
            aseMol = Atoms(str_atoms,coords)
            res_calc: RES[Atoms] = calculator(aseMol,method='GFN1-xTB',criteria=0.5,steps=60)
            aseMol: Atoms = res_calc.data if res_calc.data else aseMol
            if not res_calc.success:
                error_str:str = f'Preliminary screening error: fail to get 3D-info from molBlock {molBlock}\n {res_calc.error}'
                logger.error(error_str)
                RES_result.error = RES_result.error + error_str
            aseList.append(aseMol)
        energy:list[float] = [i.get_potential_energy() for i in aseList]
        index: int = energy.index(min(energy))
        aseMol = aseList[index]
        res_calc: RES[Atoms] = calculator(aseMol,method='GFN2-xTB',criteria=0.05,steps=200)
        if not res_calc.success:
            error_str:str = f'final 3D mol_conf optimization error: fail to get 3D-info from molBlock {molBlock}\n {res_calc.error}'
            logger.error(error_str)
            RES_result.error = RES_result.error + error_str
        str_atoms = aseMol.get_chemical_symbols()
        coords= aseMol.get_positions()
        RES_result.data = {'atoms': str_atoms, 'coords': coords}
        return RES_result
    except Exception as e:
        error_str: str = f'error: fail to get 3D-info from molBlock {molBlock}\n {e}'
        logger.error(error_str)
        RES_result.error = RES_result.error + error_str
        RES_result.success = False
        return RES_result

def gen3D(molBlocks: list[str] | None, smiles: list[str] | None) -> RES[dataUnimol]:
    structures = dataUnimol()
    structures.atoms = []
    structures.coordinates = []
    RES_results = RES[dataUnimol]()
    if (molBlocks is not None):
        logger.info('before generate 3D from molBlock')
        for molBlock in molBlocks:
            lines = molBlock.split('\n')

            ### 特殊处理垃圾GaussView的格式不规范
            if 'GaussView' in lines[2]:
                # 规范化文件格式版本
                lines[3] = lines[3].rstrip()[:-9] + '0999 V2000'
                molBlock = '\n'.join(lines)
                # 替换不能识别键的类型为单键
                molBlock = molBlock.replace(' 4  0  0  0', ' 1  0  0  0')

            isPlaneStruct = []
            singleMolAtoms = []
            singleMolCoordinates = []
            for line in lines[4:-2]:
                word = line.split()
                if len(word)>3 and word[3].isalpha():
                    isPlaneStruct.append(float(word[2]) == 0)
                    singleMolAtoms.append(word[3])
                    singleMolCoordinates.append([float(word[0]),float(word[1]),float(word[2]),])
            ## if mol is plane structure, need to generate 3D conformer
            if isPlaneStruct.count(True) > isPlaneStruct.count(False): 
                res_xtb:RES[dict] = xTBconvertor(molBlock=molBlock)
                if res_xtb.success and res_xtb.data is not None:
                    struct3d:dict = res_xtb.data
                    singleMolAtoms = struct3d['atoms']
                    singleMolCoordinates = struct3d['coords']
                else:
                    error_str:str = f'error: fail to get 3D-info from molBlock {molBlock}\n {res_xtb.error}'
                    RES_results.error = RES_results.error + error_str
                    logger.error(error_str)
            structures.atoms.append(singleMolAtoms)
            structures.coordinates.append(singleMolCoordinates)
    elif (smiles is not None):
        logger.info('before generate 3D from smi')
        for smi in smiles:
            res_xtb:RES[dict] = xTBconvertor(smi=smi)
            if res_xtb.success and res_xtb.data is not None:
                struct3d:dict = res_xtb.data
                structures.atoms.append(struct3d['atoms'])
                structures.coordinates.append(struct3d['coords'])
            else:
                error_str:str = f'error: fail to get 3D-info from molBlock {smi}\n {res_xtb.error}'
                RES_results.error =  RES_results.error + error_str
                logger.error(error_str)
    else:
        logger.error("Invalide input: mol_block and smi are all None")
    RES_results.data = structures
    return RES_results

