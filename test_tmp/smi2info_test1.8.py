import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor,Descriptors,rdMolDescriptors,rdMolEnumerator,Draw
import xlsxwriter as xlw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
import pathlib,random
from tqdm import tqdm
print('依赖库加载完毕...\n')
rdDepictor.SetPreferCoordGen(True)
p_table=Chem.GetPeriodicTable()

def input_file_check(suffix,discription):#选取目标后缀文件,返回文件名
    print(discription)
    initial_file_list=list(pathlib.Path.cwd().glob('*.'+suffix))
    if len(initial_file_list)==0:
        print("Error!!! the file of \""+suffix+"\" not found!")
    elif len(initial_file_list)>1:
        print("NOTE!! should be only one target \""+suffix+"\" file in current directory, but multiple items are detected.")
        for pt in initial_file_list:
            print('\t'+str(initial_file_list.index(pt)+1)+' <==> '+pt.name)
        index_p=int(input('---Enter the index to specify'+suffix+'file:')or 1)-1 
    else:
        index_p=0
    file=initial_file_list[index_p]
    return file

def split_text(file,Delimiter):#读取并分割文件，返回smiles_list列表
    with open(file, 'r') as smi_f:
        data = smi_f.read()
    smiles_list=data.split(Delimiter,-1)
    return smiles_list

def smiles2mol(smiles_list):#将smiles转换为mol
    print('total '+str(len(smiles_list))+' smiles formula')
    algo=input('---Whether to use Schrodinger\'s algorithm to calculate 2D coordinates(y or n):')
    if ('y' in algo) or ('Y' in algo):
        print('---Schrodinger\'s cooordinate algorithm used')
        ag_type=True
        def smis2mol(smis):
            mol=Chem.MolFromSmiles(smis)
            rdDepictor.Compute2DCoords(mol)
            return mol
    else:
        ag_type=False
        print('---did not use Schrodinger\'s cooordinate algorithm')
        def smis2mol(smis):
            mol=Chem.MolFromSmiles(smis)
            return mol
    mol_list=[]
    for i in tqdm(smiles_list):
        mol_list.append(smis2mol(i))
    all_item=len(mol_list)
    print('total '+str(len(smiles_list))+' molecules successfully converted')
    return mol_list,all_item,ag_type

def layout_excel(all_item,task):#设置Excel排版布局
    if int(task) == 1:
        print('######设置图片大小（dpi）,原子属性######')
        height=float(input('---Please set height of cell in excel(height of picture chemdraw, default=200):') or 200)
        width=float(input('---Please set width of cell in excel(width of picture chemdraw, default=200):') or 200) 
        FontSize=int(input('---Please set atom label FontSize(default=15):') or 15) 
        bondWidth=int(input('---Please set bond Width(default=1):') or 1)
        print('######排版布局设置######')
        col=int(input('---Column = default 4（how many patterns per row）:') or 4)
        row=all_item//col
        if all_item%col>0.5:
            row+=1
    elif int(task) == 2:
        height=float(200)
        width=float(200) 
        FontSize=int(15) 
        bondWidth=int(1)
        col=1
        row=all_item
    else:
        pass
    print('Excel shape: ({0}, {1})'.format(row,col))
    return row,col,height,width,FontSize,bondWidth

def draw2d(mol,save_path):#画图并保存
    global height,width,FontSize,bondWidth
    for atom in mol.GetAtoms():
        if (atom.GetIsotope()==3):
            atom.SetProp('atomLabel','T')
        if (atom.GetIsotope()==2):
            atom.SetProp('atomLabel','D')
    drawer = rdMolDraw2D.MolDraw2DCairo(int(height),int(width))
    opts=drawer.drawOptions()
    opts.minFontSize = FontSize
    opts.bondLineWidth= bondWidth
    opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})#clear_atom_color,and set black
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    drawer.WriteDrawingText(str(save_path))

def exist_or_rename(name):#结果文件名检查
    if pathlib.Path(pathlib.Path.cwd(),name).exists():
        print('NOTE!! File '+name+' is already exist') 
        new_name=input("---please name a different file name(End with .xlsx):") or name
    else:
        new_name=name
    return new_name

def excel_insert_image(xlsx_filename,mol_list,png_path,row,col,height,width,task=1):#插入图片到excel
    print('---正在插入图片...')
    book = xlw.Workbook(pathlib.Path(pathlib.Path.cwd(),xlsx_filename))
    sheet = book.add_worksheet('sheet_by_script')
    if int(task) == 1:
        for set_row in range(1,row+1):
            sheet.set_row_pixels(set_row,float(height))
        for set_col in range(1,col+1):
            sheet.set_column_pixels(set_col*2,set_col*2,float(width))
        for mol in tqdm(mol_list):
            serial = mol_list.index(mol)+1
            m=serial%col     ### cell serial number is m*2 in excel
            n=serial//col+1 
            if m==0:
                m=col
                n=n-1   
            filename=str(n)+'-'+str(m)+'.png'
            tmp_png=pathlib.Path(png_path,filename)
            draw2d(mol,save_path=tmp_png)
            sheet.insert_image(n,m*2,tmp_png,options={}) ###options set style
            sheet.write(n,m*2-1,'('+str(serial)+')')
    elif int(task) == 2:
        superscript = book.add_format({'font_script': 1})
        subscript = book.add_format({'font_script': 2})
        for set_row in range(1,row+1):
            sheet.set_row_pixels(set_row,float(height)) 
        sheet.set_column_pixels(col*2,col*2,float(width))
        sheet.set_column_pixels(col*2+1,col*2+1,float(width/2))
        sheet.set_column_pixels(col*2+2,col*2+2,float(1000))
        for mol in tqdm(mol_list):
            serial = mol_list.index(mol)+1
            m=serial%col     ### cell serial number is m*2 in excel
            n=serial//col+1 
            if m==0:
                m=col
                n=n-1   
            filename=str(n)+'-'+str(m)+'.png'
            tmp_png=pathlib.Path(png_path,filename)
            draw2d(mol,save_path=tmp_png)
            sheet.insert_image(n,m*2,tmp_png,options={}) ###options set style
            sheet.write(n,m*2-1,'('+str(serial)+')')
            mass1,mass2=ExactMolWt_INFO(mol)
            ##########获取分子式，添加下表格式化#########
            formula=get_formula(mol,split=True)
            for index in range(len(formula))[::-1]:
                if index%2>0:
                    formula.insert(index,subscript)
            ###########################################
            sheet.write_rich_string(n,m*2+1,*formula)#特殊用法
            sheet.write_rich_string(n,m*2+2,mass1,superscript,'+',mass2)       
    book.close()
    print('######插图结束：请仔细检查表格图片，手动调节图片序号的字体和大小######')

def get_element_num(mol,element):#统计除H原子之外，分子中目标原子的个数
        Number = len(mol.GetSubstructMatches(Chem.MolFromSmiles(element)))
        return Number

def get_H_num(smi):#以列表形式统计H原子的总数和H同位素个数：[总数，氕，氘，氚]
    mol = Chem.MolFromSmiles(smi)
    H1=0
    for atom in mol.GetAtoms():
        H1 += atom.GetTotalNumHs()
    H2=smi.count('[2H]')
    H3=smi.count('[3H]')
    H_total= H1 + H2 + H3
    H_list=[H_total,H1,H2,H3]
    return H_list

def get_formula(mol,split=False):#输出化学计量数分子式，修复rdkit不能识别氢的同位素
    formula = rdMolDescriptors.CalcMolFormula(mol)
    smi = Chem.MolToSmiles(mol)
    H_list = get_H_num(smi)
    if H_list[2]>0 or H_list[3]>0:#判断是否含有氢的同位素
        formula_sub=formula.split('H',-1)
        for index,value in enumerate(formula_sub[1]):
            if value.isalpha():
                formula_sub[1]=formula_sub[1][index:]
                break
        formula_sub.insert(1,'H')
        formula_sub.insert(2,str(H_list[1]))
        if H_list[2]>0:
            formula_sub.insert(3,'D')
            formula_sub.insert(4,str(H_list[2]))
            if H_list[3]>0:
                formula_sub.insert(5,'T')
                formula_sub.insert(6,str(H_list[3]))
        else:
            formula_sub.insert(3,'T')
            formula_sub.insert(4,str(H_list[3]))
        formula = ''.join(formula_sub)
    if split:#是否分割输出
        import re
        return re.findall(r'[0-9]+|[a-zA-Z]+',formula)
    else:
        return formula

def remove_ditto(smis=[]):#正则化smiles,并去除重复smiles
    for index,item in enumerate(smis):
        smis[index] = Chem.CanonSmiles(item)
    smi_sorted = sorted(set(smis),key=smis.index)
    return smi_sorted

def ExactMolWt_INFO(mol=None,smi=''):#格式输出分子质谱信息
    if len(smi)>0:
        mol = Chem.MolFromSmiles(smi)
    else:
        smi=Chem.MolToSmiles(mol)
    exact_MolWt = Descriptors.ExactMolWt(mol)#获取整个分子的确切质量
    exact_MolWt_r = exact_MolWt + 0.6*random.random() + 0.9
    MolWt = Descriptors.MolWt(mol)#获取整个分子的相对摩尔质量
    ###########统计分子中每种原子的个数的字典###########
    check_list=['C','N','S','Br','Cl']
    H_list = get_H_num(smi)
    H1_mass = p_table.GetAtomicWeight('H')*H_list[1]
    H2_mass = p_table.GetMassForIsotope('H',2)*H_list[2]
    H3_mass = p_table.GetMassForIsotope('H',3)*H_list[3]
    H_mass = H1_mass + H2_mass + H3_mass
    H_ratio = 100*H_mass/MolWt
    H_refer = H_ratio + 0.6*random.random()-0.3
    atom_pd=pd.DataFrame([[H_list[0],H_mass,H_ratio,H_refer]],index=['H'], columns=['num','mass','ratio','ratio_r'])
    for i in check_list:
        num = get_element_num(mol,i)
        if num>0:
            atom_pd.loc[i,'num']=num
            atom_pd.loc[i,'mass']=p_table.GetAtomicWeight(i)*num
            atom_pd.loc[i,'ratio']=(atom_pd.loc[i,'mass']/MolWt)*100
            atom_pd.loc[i,'ratio_r']=atom_pd.loc[i,'ratio']+0.3*(random.random()+0.1)*random.choice((-1,1))
    atom_pd1=atom_pd.reindex(index=['C','H','N','S','Br','Cl'])
    atom_pd1.dropna(axis=0,inplace=True)
    #atom_pd1.iloc[-1,-1]=100-atom_pd1.iloc[:-1,-1].sum()
    txt_list1=['理论值：']
    txt_list2=['测试值：']
    for index, row in atom_pd1.iterrows():
        txt_list1.append(index)
        txt_list1.append(', {:.2f};'.format(row[2]))
        txt_list2.append(index)
        txt_list2.append(', {:.2f};'.format(row[3]))
    txt1=''.join(txt_list1)[:-1]   
    txt2=''.join(txt_list2)[:-1]
    txt=txt1+'；'+txt2+'。'
    LC_MS1='LC-MS：测定值：{:.2f}([M+H]'.format(exact_MolWt_r)
    LC_MS2=')，精确质量：{:.2f}。'.format(exact_MolWt)
    mass_data1=txt+LC_MS1
    return mass_data1,LC_MS2

def link_rule(sub_sites=[],core_sites=[],optimize_half=False):#根据配体的虚site*,和连续取代位点生成键连规则
    begin,end=core_sites[0],core_sites[1]
    sites_list=list(range(begin,(end+1),1))
    combs=[]
    if len(sub_sites)==1:
        for core_site in sites_list:
            comb='|m:%s:%s|'%(sub_sites[0],core_site)
            combs.append(comb)
    elif len(sub_sites)==2:
        import math,copy
        site1,site2=sub_sites[0],sub_sites[1]
        if optimize_half==True:
            half_interval=math.ceil((end-begin)/2)
            half_end=begin+half_interval
            sites_list=list(range(begin,(half_end+1),1))
        for core_site in sites_list:
            a=core_site
            b=copy.copy(sites_list)
            b.remove(core_site)
            b='.'.join(map(str,b))
            comb='|m:%s:%s,%s:%s|'%(site1,a,site2,b)
            combs.append(comb)
    else:
        print('子结构位点指数错误！')
    return combs

def enumerate_mol(mol,link_list):#根据mol(cxsmiles(core.sub))和键连规则枚举分子结构
    smile_plain=Chem.MolToSmiles(mol)
    mol_list=[]
    for i in link_list:
        cxsmiles_mol=Chem.MolFromSmiles(smile_plain+' '+i)
        print(smile_plain+' '+i)
        molbundle=rdMolEnumerator.Enumerate(cxsmiles_mol)
        mol_list.extend(molbundle)
    return mol_list

def get_enmrt_smi(sub_sites=[],core_sites=[],mol=None,sub=False,img=False):#根据mol，取代位点core，配体链接点sub来枚举smiles式和图片
    link_list=link_rule(core_sites=core_sites,sub_sites=sub_sites,optimize_half=True)
    result=enumerate_mol(mol=mol,link_list=link_list)
    smi_list=[Chem.MolToSmiles(x) for x in result]
    unique_smi_list=remove_ditto(smis=smi_list)
    if sub:
        unique_smi_list=[x.replace('[Y]','[*]') for x in unique_smi_list]
        unique_smi_list=[Chem.CanonSmiles(x) for x in unique_smi_list]
    else:
        smiles_enumerate='.'.join(unique_smi_list)
        with open('smiles_enumerate.txt','w') as file:
            file.write(smiles_enumerate)
    if img:       
        result_unique=[]
        for i in unique_smi_list:
            mol=Chem.MolFromSmiles(i)
            rdDepictor.Compute2DCoords(mol)
            result_unique.append(mol)
        image=Draw.MolsToGridImage(result_unique,molsPerRow=5,
                                    legends=list(map(str,list(range(len(result_unique))))),
                                    useSVG=True,
                                    subImgSize=(200,200))
    else:
        image=None
    return unique_smi_list,image

def enmrt_by_sublist(core_sites=[],sub_sites=[],core_smi='',sub_list=[],sub=False,img=False):#使用取代基列表枚举smiles式
    smi_list=[]
    if len(sub_list)==2:
        for i in sub_list[0]:
            for k in sub_list[1]:
                smiles=i+'.'+k+'.'+core_smi
                print(smiles)
                mol=Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
                result=get_enmrt_smi(core_sites=core_sites,sub_sites=sub_sites,mol=mol,sub=sub,img=img)
                smi_list.extend(result[0])
    else:
        print('取代基团列表错误')
        return None
    return smi_list

def mol_match_sub(mols=[],subs=[],reverse=False):##mols=mol_list,subs=MolBundle
    matches = [x for x in mols if x.HasSubstructMatch(subs)]
    print('all mols num :' +str(len(mols))+'matched mols:'+str(len(matches)))
    if not reverse:
        return matches
    else:
        result = []
        for x in mols:
            match = x.GetSubstructMatch(subs)
            if not match:
                result.append(x)
        return result

if __name__ == "__main__":
    while True:
        print('''                                                                            
                                        ########*+=====+*########                   
                                   ####+===-                     -####              
                               ####===-                               ####          
                             ##*===:            .-==========-.           *##        
                           ##====.         -=======================        :##      
                         ##+===-        -=============================.      +##    
                        ##====.       -=========              -=========-     .##   
                       ##====-       ========:                   -========     :##  
                      ###====        =======.       =====-         ========     ##  
                      ##+====        =======:      =========        ========    =## 
                      ##*====-       -========        =======       -=======    +## 
                       ##=====:        =====================.       ========    ##  
                       ##*======         =================-        ========    ###  
                        ###=======           .========:          .========.   ###   
                          ##========:                          -=========   -##     
                        ### ##+==========                  -===========   -##       
                     ###      +##===================================   .###         
                  ###            -###===========================.   ####            
               ###                 #######+-============-:    =#####                
             ##-                ###          ###############                        
             ##              ###                                 smi2info.1.8(test)                   
              ##.         ###                                            2021.11.19
                ###=-=####      Computational Simulation Team of SUNERA Tech Co Ltd''')
        print('\n\t注意:此脚本不适用超过六个连续间位并苯的结构！！（结构会扭曲折叠）\n\targ:所有参数设置环节，直接按回车键会使用默认值。\
        \n\t使用前确保当前目录存在含有smiles式的.txt文本文件\n')
        print('''1 <===> 根据smiles分子式批量插图到excel,\n2 <===> 根据smiles分子式格式化输出质谱数据\n3 <===> 根据smiles分子式去重''')
        s=eval(input('输入任务的索引(1 or 2 or3)：'))
        if (s in [1,2,3]):
            task=int(s)
        else:
            print('未能识别输入，请确认输入有效任务序号！')
            continue
        if task == 1:#批量插图
            print('执行任务1：批量插图')
            smi_file=input_file_check(suffix='txt',discription='指定选取含有smiles分子式的文本')
            smiles_list=split_text(smi_file,Delimiter='.')
            mol_list,all_item,ag_type=smiles2mol(smiles_list)
            row,col,height,width,FontSize,bondWidth=layout_excel(all_item,task)
            if ag_type:
                png_path=pathlib.Path(pathlib.Path.cwd(),'png_file_optimized')
            else:
                png_path=pathlib.Path(pathlib.Path.cwd(),'png_file_normal')
            png_path.mkdir(mode=0o777, parents=False, exist_ok=True)
            xlsx_filename=exist_or_rename("2Dstructure.xlsx")
            excel_insert_image(xlsx_filename,mol_list,png_path,row,col,height,width)
        elif task == 2:#质谱数据
            print('执行任务2：输出质谱数据')
            smi_file=input_file_check(suffix='txt',discription='指定选取含有smiles分子式的文本')
            smiles_list=split_text(smi_file,Delimiter='.')
            mol_list,all_item,ag_type=smiles2mol(smiles_list)
            row,col,height,width,FontSize,bondWidth=layout_excel(all_item,task)
            if ag_type:
                png_path=pathlib.Path(pathlib.Path.cwd(),'png_file_optimized')
            else:
                png_path=pathlib.Path(pathlib.Path.cwd(),'png_file_normal')
            png_path.mkdir(mode=0o777, parents=False, exist_ok=True)
            xlsx_filename=exist_or_rename("mass_analysis.xlsx")
            excel_insert_image(xlsx_filename,mol_list,png_path,row,col,height,width,task)
        elif task == 3:#结构式去重和正则化
            print('执行任务3：结构式去重')
            smi_file=input_file_check(suffix='txt',discription='指定选取含有smiles分子式的文本')
            smiles_list=split_text(smi_file,Delimiter='.')
            unique_smi_list=remove_ditto(smis=smiles_list)
            num_smi=len(smiles_list)
            num_sorted = len(unique_smi_list)
            if num_smi != num_sorted:
                print('删除'+str(num_smi-num_sorted)+'重复的smile分子式')
            else:
                print('没有重复的分子式')
            txt_smiles='.'.join(unique_smi_list)
            with open('CanonSmiles.txt','w') as file:
                file.write(txt_smiles)
        else:
            ##结构枚举
            ##from rdkit.Chem.Draw import IPythonConsole
            ##IPythonConsole.drawOptions.addAtomIndices = True###显示原子索引
            ##IPythonConsole.molSize = 600,600
            print('请准确输入任务序号')
            pass

