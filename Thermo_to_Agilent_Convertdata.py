import os,csv,xlwt
# import shutil
from tkinter import *
from tkinter.filedialog import askdirectory
import tkinter.messagebox

import time


def makenew_dir():
        #为每个oigin_dir目录下的csv文件建立对应目录
    sample_Name=[]
    machine_Name=[]
    FileName=[]
    seq=[]         
    for root,dirs,files in os.walk(origin_dir):
        for file in files:
            info=open(origin_dir+'\\'+file,encoding='utf-8')
            headreader=csv.reader(info)
            for i,row in enumerate(headreader):
                 if i==0:
                    ss=row[0]
                    ss=ss.split(':')                
                    sample_Name.append(ss[0])
                 elif i==2:
                    mm=row[0].split('=')
                    mm=mm[1].split(';')
                    machine_Name.append(mm[0])                          
        for name in files:
            fname,fext=os.path.splitext(name)
            seq.append(fname.split('_')[-1])                       
            fname=fname.replace('_','')
            FileName.append(fname)            
            os.mkdir(os.path.join(output_dir,fname+'.D'))
            
#对原csv文件前12行删除，并在第一行增加转换后目录。
            newpath=os.path.join(output_dir+'\\'+fname+'.D',fname+'.csv')
            originpath=os.path.join(origin_dir+'\\'+name)
            o=open(originpath,'r')
            h=open(newpath,'w',newline='')
            s=os.path.dirname(output_dir)
            s=[''.join(s)]
            writer=csv.writer(h)
            writer.writerow(s)
            reader=csv.reader(o)
            for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]:
                next(reader)
            lineNo=0
            for line in reader:
                lineNo=lineNo+1
                if  lineNo==3:
                    pass
                else :
                    writer.writerow(line)
            h.close()                   
    return sample_Name,machine_Name,FileName,seq
    
def built_LIST(sample_Name,machine_Name,FileName,seq):
    
    workbook = xlwt.Workbook()
    worksheet=workbook.add_sheet('Sheet1')
    style = xlwt.XFStyle() # 初始化样式
    font = xlwt.Font() # 为样式创建字体
    font.name = 'Times New Roman' 
    font.bold = True # 黑体
    font.underline = True # 下划线
    font.italic = True # 斜体字
    style.font = font # 设定样式          
    worksheet.write(0, 2, 'SampleName')
    worksheet.write(0, 1, 'FileName')
    worksheet.write(0, 3, 'MachineName')
    worksheet.write(0, 0, 'sequence Num')
    FILE=FileName[0][:-1]
    for x in range(1,len(sample_Name)+1):
        
        worksheet.write(x,2,sample_Name.pop())
        worksheet.write(x,3,machine_Name.pop())
        worksheet.write(x,1,FileName.pop())
        worksheet.write(x,0,eval(seq.pop()))
    os.chdir(output_dir)
    
    
    workbook.save(FILE+'_LIST.xls')  


        

def convertdata():
    #global output_dir
    #global origin_dir 
    
    sample_Name,machine_Name,FileName,seq=makenew_dir()
    
    built_LIST(sample_Name,machine_Name,FileName,seq)
    tkinter.messagebox.showinfo(title='提示信息！', message='转换已完成！')


def main():    
    origin_dir = ''
    
    output_dir=''
    def selectPath_in():
        global origin_dir        
        origin_dir=askdirectory(title=u'选择输入文件夹',initialdir=(os.path.expanduser('H:/')))
        print('打开文件：', origin_dir)
        if origin_dir is not None:
            path1.set(origin_dir)
             
        
    
    def selectPath_out():
        global output_dir
        output_dir=askdirectory(title=u'选择目标文件夹',initialdir=(os.path.expanduser('H:/')))
        print('保存文件夹：',output_dir)
        if output_dir is not None:
            
            path2.set(output_dir)
        
            
    
        
        


    t=Tk()
    t.title('数据转换程序2020')  
    
    path1=StringVar()
    Label(t,text = "输入路径:").grid(row = 0, column = 0)
    Entry(t, textvariable = path1).grid(row = 0, column = 1)
    Button(t, text = "浏览", command = selectPath_in).grid(row = 0, column = 2)
    

    path2=StringVar()
    Label(t,text = "输出路径:").grid(row = 1, column = 0)
    Entry(t, textvariable =path2).grid(row = 1, column = 1)
    Button(t, text = "浏览", command = selectPath_out).grid(row = 1, column = 2)
    
    
    

    Button(t, text = "转换", command = convertdata).grid(row = 0, column = 3,rowspan=2)

    
    

    t.mainloop()
  
   
    

if __name__=='__main__':
    main()
    
            
           
           



  




        
