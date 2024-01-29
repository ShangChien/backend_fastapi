from pydantic import BaseModel
from typing import Iterable,Any
from pycdxml import cdxml_slide_generator, cdxml_converter # type: ignore
from rdkit import Chem
from rdkit.Chem import rdMolEnumerator,rdDepictor # type: ignore
from itertools import combinations,product
from joblib import Parallel, delayed
rdDepictor.SetPreferCoordGen(True)

## enumData,enum_atoms_smiles,layout_cdxml,smi2cdxml
class enumSetting(BaseModel):
    array:list[int]=[]
    range:list[int]=[]
    connect2index:list[int]=[]
    keepSame2Index:list[int]=[]

class core(BaseModel):
    id: int
    smiles: str
    enumAtoms:dict[int, enumSetting]={}
    enumBonds:dict[int, enumSetting]={}

class molecule(BaseModel):
  id: int
  smiles: str
  atoms: dict[int,list[int]]={}
  bonds: dict[int,list[int]]={}
    
class enumData(BaseModel):
  core: core
  ligands: list[molecule]

class layout(BaseModel):
  cols:int=10
  rows:int=20
  scale:float=1.0
  singlePage:bool=False

class smi2cdx(BaseModel):
  smiles:list[str]=[]
  layout:layout

def add_dummy2smi(mol:Any,i:int)->list[str|int]:
    dummy = Chem.Atom(0)  # type: ignore # 0 is the atomic number 
    rwmol = Chem.RWMol(mol)  # type: ignore
    rwmol.AddBond(rwmol.AddAtom(dummy),i,Chem.rdchem.BondType.SINGLE)
    nums=rwmol.GetNumAtoms()
    smi=Chem.MolToSmiles(rwmol)  # type: ignore
    return [smi,nums]

def split_bond2smi(mol:Any,i:int)->list[str|int]:
    ibond=mol.GetBondWithIdx(i)
    begin=ibond.GetBeginAtomIdx()
    end=ibond.GetEndAtomIdx()
    rwmol=Chem.RWMol(mol)  # type: ignore
    rwmol.RemoveBond(begin,end)
    rwmol.ReplaceAtom(begin,Chem.Atom(0))  # type: ignore
    rwmol.ReplaceAtom(end,Chem.Atom(0))  # type: ignore
    smi=Chem.MolToSmiles(rwmol)  # type: ignore
    nums=rwmol.GetNumAtoms()
    dummy_index=[
        atom.GetIdx()+1 
        for atom in Chem.MolFromSmiles(smi,sanitize=False).GetAtoms()   # type: ignore
            if atom.GetSymbol() == "*"
    ]
    return [smi,nums,*dummy_index]

def enum_core_site(core:core)->Iterable[dict[int,list[int]]]:
    rateCombo={}
    for k in core.enumAtoms:
        arr=core.enumAtoms[k].array
        rate=sorted(core.enumAtoms[k].range)
        rate[-1]=rate[-1]+1
        rateCombo[k]=[]
        for n in range(*rate):
            rateCombo[k].extend([list(i) for i in combinations(arr, n)])
    combo=[[*i] for i in product(*rateCombo.values())]
    colorIndex=rateCombo.keys()
    combo=map((lambda x : {v:x[i] for i,v in enumerate(colorIndex)}), combo)
    return combo

def getAtomCombos(core:core,ligands:dict[str,dict[int,list]])->list[str]:
    #组合不同颜色的位点库
    mol=Chem.MolFromSmiles(core.smiles)  # type: ignore
    nums=mol.GetNumAtoms()
    colorbase={}
    for i in core.enumAtoms.keys():
        colorbase[i]=[]
        for i1 in core.enumAtoms[i].connect2index:
            colorbase[i].extend(ligands['atoms'][i1])
    #组合带有虚原子的smi
    t_combos=[]
    combos=enum_core_site(core)
    for combo in combos:
        keys=combo.keys()
        keys_combo=[]
        for k in keys:
            keys_combo.append([])
            length=len(combo[k])
            for combo_i in product(colorbase[k],repeat=length):
                _item=[]#[[smi,totalAtoms,site]]
                for i,_v in enumerate(combo_i): 
                    _convert=_v.copy()
                    _convert.append(combo[k][i])
                    _item.append(_convert)
                keys_combo[-1].append(_item)
        for combo_color in product(*keys_combo):
            i=[c for com in combo_color for c in com]
            smis=core.smiles
            site_n=nums
            link='m:'
            for item in i:
                smis=smis+'.'+item[0]
                link=link+str(site_n)+':'+str(item[-1])+','
                site_n=site_n+item[1]
            cxsmi=smis+' '+'|'+link+'|'
            if not (link=='m:'):
                t_combos.append(cxsmi) 
    return t_combos

def convert_smi(i:str)->str:
    mol=Chem.MolFromSmiles(i)  # type: ignore
    smi=Chem.MolToSmiles(rdMolEnumerator.Enumerate(mol)[0])  # type: ignore
    Chem.CanonSmiles(smi)
    return Chem.CanonSmiles(smi)

def pre_ligand(arr:list[molecule])->dict[str,dict[int,list]]:
    from rdkit import Chem
    colorDummyAtomsIndex:dict[int,list]={}
    colorDummyBondsIndex:dict[int,list]={}
    for ligand in arr:
        mol=Chem.MolFromSmiles(ligand.smiles)  # type: ignore
        for color in ligand.atoms:
            if color not in colorDummyAtomsIndex.keys():
                colorDummyAtomsIndex[color]=[]
            colorDummyAtomsIndex[color].extend(
                list(map(lambda i:add_dummy2smi(mol,i),ligand.atoms[color]))
            )
        for color in ligand.bonds:
            if color not in colorDummyBondsIndex.keys():
                colorDummyBondsIndex[color]=[]
            colorDummyBondsIndex[color].extend(list(map(lambda i:split_bond2smi(mol,i),ligand.bonds[color])))
    return {'atoms':colorDummyAtomsIndex,'bonds':colorDummyBondsIndex}

def parallel_convert_smis(smis_link:list[str])->list[str]:
    smis=Parallel(n_jobs=8)(delayed(convert_smi)(i) for i in smis_link)
    return smis # type: ignore

async def enum_atoms_smiles(enumData:enumData)->list[str]:
    ligands=pre_ligand(enumData.ligands)
    smi_link=getAtomCombos(core=enumData.core,ligands=ligands)
    smis_p= parallel_convert_smis(smi_link)
    smis=list(set(smis_p))
    return smis

def smi2cdxml(smi):
    mol = Chem.MolFromSmiles(smi)  # type: ignore
    rdDepictor.Compute2DCoords(mol)
    return cdxml_converter.mol_to_document(mol).to_cdxml()

def layout_cdxml(cdxmls:list[Any],cols:int,rows:int=-1,single_page:bool=True,scale:float=1.0)->list[Any]:
    ##when single_page=false, must explicitly pass rows with value>0
    import math
    length = len(cdxmls)
    props = [[cdxml_slide_generator.TextProperty('index', i+1)] for i in range(length)]
    real_rows = math.ceil(length/cols) if single_page else rows
    #logic
    slides=[]
    if single_page:
        sg = cdxml_slide_generator.CDXMLSlideGenerator(
            style="ACS 1996", 
            number_of_properties=1, 
            columns=cols, rows=real_rows, 
            slide_width=cols*10*scale, 
            slide_height=real_rows*10*scale
        )
        slides.append(sg.generate_slide(cdxmls, props))
    else:
        if rows<0:
            print("Value Error: rows=", rows)
        else:
            per_page_nums=cols*rows
            total_pages=math.ceil(length/per_page_nums)
            def getSlide(pageIndex:int):
                ##必须将初始化sg放到并行函数内部才能运行
                sg = cdxml_slide_generator.CDXMLSlideGenerator(
                    style="ACS 1996", 
                    number_of_properties=1, 
                    columns=cols, rows=real_rows, 
                    slide_width=cols*10*scale, 
                    slide_height=real_rows*10*scale
                )
                start = pageIndex * per_page_nums
                end   = start + per_page_nums
                cdxml = cdxmls[start:end]
                prop = props[start:end]
                slide = sg.generate_slide(cdxml, prop)
                return slide
            #sub_cdxmls=[cdxmls[i*per_page_nums:(i+1)*per_page_nums] for i in range(total_pages)]
            #sub_props=[props[i*per_page_nums:(i+1)*per_page_nums] for i in range(total_pages)]
            #slides=[sg.generate_slide(sub_cdxmls[i], sub_props[i]) for i,_v in enumerate(sub_props)]
            slides=Parallel(n_jobs=10)(delayed(getSlide)(i) for i in range(total_pages))
    return slides # type: ignore