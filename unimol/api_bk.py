import logging
logger = logging.getLogger("__main__")
from pydantic import BaseModel
from typing import Any
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdDistGeom
rdDepictor.SetPreferCoordGen(True)
etkdg = rdDistGeom.ETKDGv3()
etkdg.randomSeed = 0xa700f
etkdg.verbose = False
etkdg.numThreads = 0
etkdg.optimizerForceTol = 0.0135
etkdg.useRandomCoords = False

# clf = MolTrain(task='classification', 
#                 data_type='molecule', 
#                 epochs=10, 
#                 batch_size=16, 
#                 metrics='auc',
#                 )
# pred = clf.fit(data = data)
# currently support data with smiles based csv/txt file, and
# custom dict of {'atoms':[['C','C],['C','H','O']], 'coordinates':[coordinates_1,coordinates_2]}

# clf = UniMolRepr(data_type='molecule')
# smiles = ['CCO', 'CCC', 'CCCC']
# reprs = clf.get_repr(smiles)
# XfeatureMol = reprs["cls_repr"]

from unimol_tools import MolPredict, UniMolRepr
from pathlib import Path as P
models={
    'homo': MolPredict(load_model=P(P(__file__).parent, f'./jobs/homo/exp')),
    'lumo': MolPredict(load_model=P(P(__file__).parent, f'./jobs/lumo/exp')),
    'eg'  : MolPredict(load_model=P(P(__file__).parent, f'./jobs/eg/exp')),
    't1'  : MolPredict(load_model=P(P(__file__).parent, f'./jobs/t1/exp')),
    's1'  : MolPredict(load_model=P(P(__file__).parent, f'./jobs/s1/exp')),
    'repr': UniMolRepr(data_type='molecule')
}


class dataUnimol(BaseModel):
    models: list[str] | None = []
    names : list[str] | None = []
    smiles: list[str] | None = []
    molBlocks: list[str] | None = []
    atoms : list[list[str]] | None = [] 
    coordinates : list[list[list[float]]] | None = []
    results : dict[str, list[int|float]] | None = {} ## { 'model1': [1,2,3] }
    
def gen3DbySmiles(smiles: list[str]) -> dataUnimol:
    structures = dataUnimol()
    for v in smiles:
        try:
            mol = Chem.MolFromSmiles(v)
            mol = Chem.AddHs(mol)
            _cid = rdDistGeom.EmbedMultipleConfs(mol, numConfs=1, params=etkdg)
            symbols= [i.GetSymbol() for i in mol.GetAtoms()]
            conf=mol.GetConformer()
            positions=conf.GetPositions()
            structures.atoms.append(symbols)
            structures.coordinates.append(positions)
        except Exception:
            logger.error(f'error: fail to prase smile to 3D  {v}')
    return structures

def get3DbyRdkit(molBlock) -> dict[str, list]:
    singleMolAtoms = []
    singleMolCoordinates = []
    try:
        mol = Chem.MolFromMolBlock(molBlock)
        mol = Chem.AddHs(mol)
        _cid = rdDistGeom.EmbedMultipleConfs(mol, numConfs=1, params=etkdg)
        singleMolAtoms = [i.GetSymbol() for i in mol.GetAtoms()]
        conf = mol.GetConformer()
        singleMolCoordinates = conf.GetPositions().tolist()   
    except Exception:
        logger.error(f'error: fail to get 3D-info from molBlock {molBlock}')
    return {'singleMolAtoms':singleMolAtoms, 'singleMolCoordinates':singleMolCoordinates}

def gen3DbyMolBlock(MolBlocks: list[str]) -> dataUnimol:
    structures = dataUnimol()
    for v in MolBlocks:
        lines = v.split('\n')
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
            struct3d = get3DbyRdkit(molBlock = v)
            singleMolAtoms = struct3d['singleMolAtoms']
            singleMolCoordinates = struct3d['singleMolCoordinates']
        structures.atoms.append(singleMolAtoms)
        structures.coordinates.append(singleMolCoordinates)
    return structures

def inference(data: dataUnimol) -> dataUnimol| Any:
    global models
    inputData = {}
    data.results = {}
    if (data.atoms is not None) and (data.coordinates is not None):
        logger.info(f'process atoms and coord')
        inputData['atoms'] = data.atoms
        inputData['coordinates'] = data.coordinates
    elif (data.molBlocks is not None):
        try:
            logger.info(f'before process molblocks')
            data3D = gen3DbyMolBlock(data.molBlocks)
            inputData['atoms'] = data3D.atoms
            inputData['coordinates'] = data3D.coordinates
        except Exception:
            data.results = {'molblocks':['process error']}
            return data
    elif (data.smiles is not None):
        logger.info(f'before prase smiles to 3D')
        try:
            data3D = gen3DbySmiles(data.smiles)
            inputData['atoms'] = data3D.atoms,
            inputData['coordinates'] = data3D.coordinates
        except Exception:
            data.results = {'smiles':['process error']}
            return data
    else:
        logger.info(f'input data wrong')
        data.results= {'inputData': ['data is invalid!']} 
        return data
    for model in data.models:
        model = model.lower()
        results = None
        if model == 'repr':
            reprs = models['repr'].get_repr(inputData)
            results = reprs["cls_repr"]
        elif model in models.keys():
            results = models[model].predict(data=inputData)
            results = results.flatten().tolist()
        else:
            logger.error(f'unimol model type error: {model}')
        data.results[model] = results
        logger.info(f'unimol res: {data.results[model]}')
    return data