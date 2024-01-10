#!/home/sanyue/.conda/envs/rdkit/bin/python
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdDistGeom
from ase.optimize import FIRE
from ase import Atoms
from xtb.ase.calculator import XTB
from pathlib import Path as P
from typing import Any

rdDepictor.SetPreferCoordGen(True)
etkdg:Any = rdDistGeom.ETKDGv3()
etkdg.randomSeed = 0xa700f
etkdg.verbose = False
etkdg.numThreads = 0
etkdg.optimizerForceTol = 0.0135
etkdg.useRandomCoords = False

def initAseMol(molBlock:str) -> Atoms:
    mol: Any = Chem.MolFromMolBlock(molBlock)
    mol: Any = Chem.AddHs(mol)
    _cid: int = rdDistGeom.EmbedMultipleConfs(mol, numConfs=1, params=etkdg)
    symbols: list[str] = [i.GetSymbol() for i in mol.GetAtoms()]
    conf: int = mol.GetConformer()
    positions: list[list[float]] = conf.GetPositions()
    aseMol = Atoms(symbols,positions)
    return aseMol
        
def calculator(aseMol:Atoms, method="GFN2XTB", criteria=0.05, steps=100) -> Atoms:
    aseMol.calc = XTB(method=method)
    FIRE(aseMol).run(fmax = criteria,steps = steps)
    # remove tmp files
    P(P.cwd(),'gfnff_adjacency').unlink(missing_ok=True)
    P(P.cwd(),'gfnff_topo').unlink(missing_ok=True)
    return aseMol
    
def optimized_Mol(mol)-> dict[str, Any]:
    aseMol: Atoms = calculator(mol,method='GFN2XTB',criteria=0.1,steps=100)
    symbols: list = aseMol.get_chemical_symbols()
    coords: Any = aseMol.get_positions()
    return {'atoms': symbols, 'positions': coords}
