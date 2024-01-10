from rdkit import Chem
from pycdxml import cdxml_slide_generator, cdxml_converter
from rdkit.Chem import rdDepictor
from typing import Union,Iterable,List,Set,Dict,Any
rdDepictor.SetPreferCoordGen(True)
smis1=["C1=CC1c1ccc2c(oc3c(C4CCC4)c(C4CCC4)ccc32)c1C1C=C1",
    "C1=CC1c1cc(C2C=C2)c2oc3c(C4CCC4)c(C4C=C4)c(C4C=C4)cc3c2c1",
    "C1=CC1c1cc(C2CCC2)cc2c1oc1c(C3CCC3)c(C3CCC3)cc(C3CCC3)c12",
    "C1=CC1c1c(C2CCC2)cc2oc3cc(C4CCC4)c(C4CCC4)c(C4CCC4)c3c2c1C1CCC1",
    "C1=CC1c1cc(C2C=C2)c2oc3c(C4CCC4)c(C4CCC4)c(C4C=C4)c(C4CCC4)c3c2c1",
    "C1=CC1c1c(C2CCC2)c(C2C=C2)c2oc3c(C4CCC4)c(C4CCC4)c(C4C=C4)c(C4CCC4)c3c2c1C1CCC1",
    "C1=CC1c1ccc2c(oc3c(C4CCC4)ccc(C4C=C4)c32)c1C1CCC1",
    "C1=CC1c1ccc(C2CCC2)c2oc3c(C4CCC4)c(C4CCC4)c(C4CCC4)cc3c12",
    "C1=CC1c1ccc(C2C=C2)c2c1oc1c(C3CCC3)c(C3CCC3)c(C3CCC3)c(C3CCC3)c12",
    "C1=CC1c1ccc(C2C=C2)c2c1oc1ccc(C3CCC3)c(C3CCC3)c12",
    "C1=CC1c1cc2c(oc3c(C4CCC4)c(C4CCC4)cc(C4C=C4)c32)c(C2CCC2)c1C1C=C1",
    "C1=CC1c1c(C2C=C2)c(C2CCC2)c2c(oc3c(C4CCC4)c(C4CCC4)cc(C4C=C4)c32)c1C1C=C1",
    "c1c(C2CCC2)cc2c(oc3c(C4CCC4)ccc(C4CCC4)c32)c1C1CCC1",
    "C1=CC1c1cc2c(oc3cc(C4CCC4)cc(C4C=C4)c32)c(C2CCC2)c1C1C=C1",
    "C1=CC1c1c(C2CCC2)c(C2CCC2)c2c(oc3c(C4CCC4)ccc(C4CCC4)c32)c1C1CCC1",
    "C1=CC1c1cc2c(oc3c(C4C=C4)c(C4C=C4)c(C4C=C4)c(C4CCC4)c32)c(C2CCC2)c1C1C=C1",
    "C1=CC1c1c(C2CCC2)c(C2CCC2)c2oc3c(C4C=C4)c(C4CCC4)cc(C4CCC4)c3c2c1C1CCC1",
    "C1=CC1c1cc2oc3cc(C4CCC4)c(C4C=C4)c(C4CCC4)c3c2cc1C1C=C1",
    "C1=CC1c1cc2oc3ccc(C4CCC4)c(C4C=C4)c3c2cc1C1CCC1",]

#joblib并行
from joblib import Parallel,delayed
def smi2cdxml(smi):
	mol = Chem.MolFromSmiles(smi)
	return cdxml_converter.mol_to_document(mol).to_cdxml()
cdxmls=Parallel(n_jobs=6)(delayed(smi2cdxml)(smi) for smi in smis1)

def layout_cdxml(cdxmls:List[Any],cols:int,rows:int=-1,single_page:bool=True)->List[Any]:
	##when single_page=false, must explicitly pass rows with value>0
	import math
	length = len(cdxmls)
	props = [[cdxml_slide_generator.TextProperty('index', i+1)] for i in range(length)]
	real_rows = math.ceil(length/cols) if single_page else rows
	sg = cdxml_slide_generator.CDXMLSlideGenerator(
		style="ACS 1996", 
		number_of_properties=1, 
		columns=cols, rows=real_rows, 
		slide_width=cols*10, 
		slide_height=real_rows*10
	)
	#logic
	slides=[]
	if single_page:
		slides.append(sg.generate_slide(cdxmls, props))
	else:
		if rows<0:
			print("Value Error: rows=", rows)
		else:
			per_page_nums=cols*rows
			total_pages=math.ceil(length/per_page_nums)
			sub_cdxmls=[cdxmls[i*per_page_nums:(i+1)*per_page_nums] for i in range(total_pages)]
			sub_props=[props[i*per_page_nums:(i+1)*per_page_nums] for i in range(total_pages)]
			slides=[sg.generate_slide(sub_cdxmls[i], sub_props[i]) for i,_v in enumerate(sub_props)]
	return slides

slides=layout_cdxml(cdxmls=cdxmls,cols=6,rows=30)
len(slides)
for index,slide in enumerate(slides):
	with open("slide-"+str(index)+".cdxml", "w", encoding='UTF-8') as xf:
		xf.write(slide)
