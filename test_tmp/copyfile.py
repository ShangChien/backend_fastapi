from  pathlib import Path
import re,shutil

#定义遍历和复制的函数，需要传入的参数：file=不带后缀的文件名;date=四个字符年月（一般从file中提取）
def find_copy(file='',date=''):
    #使用pathlib模块中glob函数遍历所有匹配的文件，生成可以用for循环切片的迭代器（_generator）
    _generator=Path(Path('Z:/Data-month-statistics'),'20'+date).glob('./*/'+file[:15]+'/'+file+'.gjf/'+file+'.log')
    n=0 #定义一个n来判断for循环内部是否执行
    for item in _generator:
        n+=1
        with open(item) as data:
            return data
        ##shutil.copyfile(item,Path(Path.cwd(),'in',file+'.log'))
    ##if n==0:#如果n为0则上面的for循环内部没有执行，即迭代器的元素为0没有找到文件
        #改变glob函数的正则表达式再次匹配，生成迭代器（_generator）
    _generator=Path(Path('Z:/Data-month-statistics'),'20'+date).glob('./*/'+file[:15]+'.gjf/'+file+'.log')
    for item in _generator:
        n+=1
        with open(item) as data:
            return data
        ##shutil.copyfile(item,Path(Path.cwd(),'in',file+'.log'))
    ##return True if n > 0 else False #根据n是否为0判断是否成功，成功返回True,失败返回False

if __name__== "__main__":                                   #如果是主程序，则以下命令运行
    Path('./in').mkdir(exist_ok=True)                       #如果当前文件夹下不存在in目录则新建
    source_root = Path('Z:/Data-month-statistics')          #源文件的根目录

    with open('input.txt') as a:                            #打开输入文件
        lines = a.readlines()                               #读取每一行文本，以字符串格式保存在列表中
    files = [x.strip() for x in lines if re.match(r'.*\d{6}-\w{3}-\d{4}.*',x)]
    #files = [re.findall(r'\d{6}-\w{3}-\d{4}',x)[0] for x in lines if re.match(r'.*\d{6}-\w{3}-\d{4}.*',x)]
    #匹配每一行是否含有文件名。如果是则截断提取文件名，将所有文件名以字符串元素形成列表。
    #python内置re模块，具体正则表达式匹配规则可以参考菜鸟教程

    for file in files:
        if find_copy(file=file,date=file[:4]):#文件在计算编号的月份文件夹中
            print(file,'copied successfully')
        elif find_copy(file=file,date=str(int(file[:4])+1)):#文件在计算编号的月份+1文件夹中
            print(file,'copied successfully')
        elif file[2:4] == '12':
            if find_copy(file=file,date=str(int(file[:2])+1)+'01'):#文件计算编号跨年
                print(file,'copied successfully')
            else:
                print(file,'copied failed')
        else:
            print(file,'copied failed')
    input('复制完成，输入任何字符，回车关闭窗口')



