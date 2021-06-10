###############################
import urllib.request
import urllib.parse
from urllib.request import urlopen
import copy
import time
import csv
import gzip
import ftplib
import glob
import os
import shutil
import xarray
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal

###################################

def NOAA_gage_data_request(begin_date, 
                           end_date, 
                           station,
                           product="water_level",
                           interval="60", 
                           units="english", 
                           time_zone="gmt", 
                           datum="NAVD"):
    option={}
    option['begin_date']=begin_date
    option['end_date']=end_date
    option['station']=station
    option['product']=product
    option['interval']=interval
    option['units']=units
    option['time_zone']=time_zone
    option['datum']=datum
    option['application']='web_services'
    option['format']='csv'
    url= 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter'
    url_values= urllib.parse.urlencode(option)
    full_url= url+'?'+url_values
    data=urllib.request.urlretrieve(full_url)
    urllib.request.urlretrieve(full_url, "gage_data/"+str(station)+'.csv')

def USGS_gage_data_request(begin_date,end_date,site_no):
    from urllib.request import urlopen
#    data=urlopen("https://waterdata.usgs.gov/nwis/uv?cb_00060=on&&format=rdb&site_no="+site_no+"&period=&begin_date="+begin_date+"&end_date="+end_date)
    data=urlopen("https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00060=on&&format=rdb&site_no="+site_no+"&period=&begin_date="+begin_date+"&end_date="+end_date) 
    lines=data.read().decode("utf8").split("\n")
    dlist_usgs=[]
    vlist_usgs=[]
    for line in lines:
        if (line ==''):
            pass
        elif (line[0]=='#'):
            pass
        elif (line[0]=='a'):
            pass
        else:
            try:
                float(line.split('\t')[4])
                dlist_usgs.append(line.split('\t')[2])
                vlist_usgs.append(line.split('\t')[4])
            except:
                pass
    np.savetxt("gage_data/"+site_no+".csv", np.c_[dlist_usgs,vlist_usgs], fmt="%s,%s")
    return dlist_usgs, vlist_usgs
    

def csvParseToLists(csvFile):
    dList = []
    vList = []
    with open(csvFile) as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            #print (row)
            dList.append(row["Date Time"])
            vList.append(row[" Water Level"])    
    return dList, vList#map(float, vList)
                       
def hecTimeParser(timesList):
    hecTimes = []
    for time in timesList:
        try:
            hecTime=HecTime()
            hecTime.set(time)
            hecTimes.append(hecTime.value())
        except Exception as e :
            print(e)    
    return hecTimes



def getQPF(interval=None,destination='precip/qpf/'):
    try:
        os.mkdir(destination)
    except:
        print(destination+" already exists.")
    ftp = ftplib.FTP('ftp.wpc.ncep.noaa.gov')
    ftp.login()
    ftp.cwd('/5km_qpf/')
    filenames=ftp.nlst()
    filenames[:] = [x for x in filenames if "p06m" in x]
    local_filenames=[os.path.basename(x) for x in glob.glob(destination+'*')]
    local_filenames[:]=[filename for filename in local_filenames if "idx" not in filename]
    filenames = list(set(filenames)^set(local_filenames))
    filenames = list(set(filenames)-set(local_filenames))
    for i in filenames:
        with open(destination+i, 'wb') as localfile:
            try:
                ftp.retrbinary('RETR '+i, localfile.write)
            except:
                try:
                    localfile.close()
                    os.remove(destination+'gz/'+i)
                    ftp.retrbinary('RETR '+i, localfile.write)
                except:
                    localfile.close()
                    print(i," had problem")
        localfile.close()
    ftp.close()
    return filenames

def getQPE(interval=None,destination='precip/qpe/'):
    try:
        os.mkdir(destination)
    except:
        print(destination+" already exists.")
    try:
        os.mkdir(destination+"gz")
    except:
        print(destination+"gz already exists.")
        
    ftp = ftplib.FTP('tgftp.nws.noaa.gov')
    ftp.login()
    ftp.cwd('/data/rfc/lmrfc/xmrg_qpe/')
    filenames=ftp.nlst()
    filenames[:] = [x for x in filenames if "xmrg" not in x]
    local_filenames=[os.path.basename(x) for x in glob.glob(destination+'gz/*')]
    filenames = list(set(filenames)^set(local_filenames))
    filenames = list(set(filenames)-set(local_filenames))
    for i in filenames:
        with open(destination+'gz/'+i, 'wb') as localfile:
            try:
                ftp.retrbinary('RETR '+i, localfile.write)
            except:
                try:
                    localfile.close()
                    os.remove(destination+'gz/'+i)
                    ftp.retrbinary('RETR '+i, localfile.write)
                except:
                    localfile.close()
                    print(i," had problem")
        localfile.close()
    ftp.close()
    
    return filenames

def QPE_unzip(interval=None,destination='precip/qpe/'):
    try:
        os.mkdir(destination+"grib")
    except:
        print(destination+"grib already exists.")
        
    try:
        os.mkdir(destination+"../temp")
    except:
        print(destination+"../temp already exists.")
        
    zipped_filenames=[os.path.basename(x) for x in glob.glob(destination+'gz/*')]
    filenames=[os.path.basename(x) for x in glob.glob(destination+'grib/*')]
    filenames = [x+".gz" for x in filenames if "gz" not in x]
    filenames = list(set(zipped_filenames)^set(filenames))
    for i in filenames:
        shutil.copyfile(destination+'gz/'+i, destination+'../temp/'+i)
        with gzip.open(destination+'../temp/'+i, 'rb') as f_in:
            with open(destination+'grib/'+i[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            f_out.close()
        f_in.close()
        os.remove(destination+'../temp/'+i)
    return filenames

def stack_gribs_qpe(grib_dir='precip/qpe/grib/'):
    files=[os.path.basename(x) for x in glob.glob(grib_dir+'*')]
    files[:]=[filename for filename in files if "idx" not in filename]
    latitude=xarray.load_dataset(grib_dir+files[0],engine='cfgrib').latitude   
    longitude=xarray.load_dataset(grib_dir+files[0],engine='cfgrib').longitude-360   
    precip=[]
    for filename in files:
        with xarray.open_dataset(grib_dir+filename,engine='cfgrib') as ds:
            precip.append(ds.load().tp)
        ds.close()
    clean=files=[os.path.basename(x) for x in glob.glob(grib_dir+'*idx')]
    for junk in clean:
        os.remove(grib_dir+junk)
    return [longitude,latitude,precip]

def stack_gribs_qpf(grib_dir='precip/qpf/',start_time=None):        
    files=[os.path.basename(x) for x in glob.glob(grib_dir+'*')]
    files[:]=[filename for filename in files if "idx" not in filename]
    if start_time==None:
        start_time=max([int(filename[5:-8]) for filename in files ])
    files[:]=[filename for filename in files if str(start_time) in filename]
    print(files)
    latitude=xarray.load_dataset(grib_dir+files[0],engine='cfgrib').latitude   
    longitude=xarray.load_dataset(grib_dir+files[0],engine='cfgrib').longitude-360   
    precip=[]
    for filename in files:
        with xarray.open_dataset(grib_dir+filename,engine='cfgrib') as ds:
            precip.append(ds.load().tp)
        ds.close()
    clean=files=[os.path.basename(x) for x in glob.glob(grib_dir+'*idx')]
    for junk in clean:
        os.remove(grib_dir+junk)
    return [longitude,latitude,precip]

def plot_precip(data,workdir):
    import matplotlib.pylab as pl
    from matplotlib.colors import ListedColormap
    cmap = pl.cm.gist_ncar
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)
    files = glob.glob('//precip//png//*')
    for f in files:
        os.remove(f)

    for i in range(len(data[2])):
        plt.pcolormesh(data[0]-360, data[1],data[2][i],cmap=my_cmap)
        plt.axis('off')
        plt.savefig('precip//png//test'+f"{i:04d}.png"+'.png',transparent=True,bbox_inches='tight', pad_inches=0,dpi=100)
        plt.clf()
    try:
        stream = os.remove(workdir+'\output.mov')
    except:
        pass
    try:
        stream = os.remove(workdir+'\output.gif')
    except:
        pass
    try:
        stream = os.popen('ffmpeg -i '+workdir+'//precip//png//test%4d.png.png -vcodec prores_ks -pix_fmt yuva444p10le -profile:v 4444 -q:v 50 -filter:v "setpts=2.5*PTS" '+workdir+'//output.mov')
    except:
        pass
    try:
        stream = os.popen('ffmpeg -i '+workdir+'//output.mov'+' -lavfi split[v],palettegen,[v]paletteuse '+workdir+'//output.gif')
    except:
        pass

def run_ras(path,plan=None):
    import rascontrol
    rc = rascontrol.RasController(version='5X')
    rc.open_project(path)
    if plan == None:
        rc.run_current_plan()
    else:
        rc.set_plan(plan)
        rc.run_current_plan()
    rc.close()

def nan_filter(data):
    for i in range(len(data)):
        if data[i] == '':
            data[i]=data[i-1]
    return data

def convert_data_to_RAS_string(data):
    data=round(data,3)
    if len(str(data))>8:
        if str(data)[0]=='-':
            return '%.1E' % Decimal(str(data))
        else:
            return '%.2E' % Decimal(str(data))
    else:
        return (8-len(str(data)))*' '+str(data)

def write_BC(data):
    line=''
    count=0
    for i in range(len(data)):
        line+=convert_data_to_RAS_string(data[i])
        if (count+1)%10==0:
            line+='\n'
        count+=1
    line+='\n'
    return line

def flow_gage_2_unsteady_flow_file(gage,unsteady_flow_file,bc_name):   
    rawdata=np.loadtxt("gage_data\\"+gage+".csv", dtype=str, delimiter=',')
    data=[float(item) for item in rawdata[:,1].tolist()]
    
    f = open(unsteady_flow_file, "r")
    lines=f.readlines()
    switch=False
    for number,line in enumerate(lines):
        if ("Boundary Location" in line) and (bc_name in line):
            startline=number+2
            switch=True
        elif ("Stage Hydrograph TW Check" in line) and (switch==True):
            endline=number
            break
    f.close()
    new_lines=[]
    for i,line in enumerate(lines):
        if (i> startline) and (i<endline):
            pass
        else:
            new_lines.append(line)

    new_lines[startline]="Flow Hydrograph="+str(len(data))+'\n'
    new_lines.insert(startline+1,write_BC(data))
    with open(unsteady_flow_file, "w") as f:
        for i,line in enumerate(new_lines):  
            f.writelines(line)
    f.close()

def stage_gage_2_unsteady_flow_file(gage,unsteady_flow_file,bc_name):   
    rawdata=np.loadtxt("gage_data\\"+gage+".csv", dtype=str, delimiter=',')
    rawdata=nan_filter(rawdata[1:,1].tolist())
    data=[float(item) for item in rawdata]
    
    f = open(unsteady_flow_file, "r")
    lines=f.readlines()
    switch=False
    for number,line in enumerate(lines):
        if ("Boundary Location" in line) and (bc_name in line):
            startline=number+2
            switch=True
        elif ("Stage Hydrograph TW Check" in line) and (switch==True):
            endline=number-1
            break
    f.close()
    new_lines=[]
    for i,line in enumerate(lines):
        if (i> startline) and (i<endline):
            pass
        else:
            new_lines.append(line)

    new_lines[startline]="Stage Hydrograph="+str(len(data))+'\n'
    new_lines.insert(startline+1,write_BC(data))
    
    with open(unsteady_flow_file, "w") as f:
        for i,line in enumerate(new_lines):  
            f.writelines(line)
    f.close()

def date_2_plan_file(start_date,end_date,plan_file):   
    f = open(plan_file, "r")
    lines=f.readlines()
    for number,line in enumerate(lines):
        if "Simulation Date" in line:
            line_number=number
            old_date=line.split(",")
            new_line=old_date[0][:-9]+start_date+","+old_date[1]+","+end_date+","+old_date[3]

    f.close()
    new_lines=[]
    for i,line in enumerate(lines):
        if (i==line_number):
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    
    with open(plan_file, "w") as f:
        for i,line in enumerate(new_lines):  
            f.writelines(line)
    f.close()

def find_u_file(plan_file):
    with open(plan_file, "r") as f:
        lines=f.readlines()
    f.close()
    for line in lines:
        if ("Flow File" in line):
            u=line[-4:-1]
    flow_path=glob.glob('\\'.join(plan_file.split("\\")[:-1])+'\\*'+u)
    return flow_path[0]

def BC_list(unsteady_flow_file):
    f = open(unsteady_flow_file, "r")
    lines=f.readlines()
    f.close()
    BC_list=[]
    for line in lines:
        if ("Boundary Location" in line):
            BC_list.append(line.split()[-1][1:])
    return BC_list
