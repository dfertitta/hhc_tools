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
from osgeo import gdal
import html.parser
from os import makedirs
import h5py
import math
import progressbar
import xml.etree.ElementTree as ET
import netCDF4 as nc
from scipy import spatial
import datetime
import pandas as pd

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

def USACE_gage_data_request(begin_date, end_date, site_no,rating_curve, variable="HG"):
    """downloads and processes xml data into time series"""
    import ssl
    import urllib.request

    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        USACE_xml=urlopen("https://rivergages.mvr.usace.army.mil/watercontrol/webservices/rest/webserviceWaterML.cfc?meth=getValues&site="+site_no+"&location="+site_no+"&&variable="+variable+"&beginDate="+begin_date+"0:00&endDate="+end_date+"0:00&authToken=RiverGages&method=RGWML")
    except:
        try:
            import ssl
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            USACE_xml=urlopen("https://rivergages.mvr.usace.army.mil/watercontrol/webservices/rest/webserviceWaterML.cfc?meth=getValues&site="+site_no+"&location="+site_no+"&&variable="+variable+"&beginDate="+begin_date+"0:00&endDate="+end_date+"0:00&authToken=RiverGages&method=RGWML",context=ctx)
        except:
            raise
    root = ET.parse(USACE_xml).getroot()
    dlist_usace=[]
    vlist_usace=[]
    for value in root[1][2][:]:
        dlist_usace.append(value.attrib["dateTime"][:-6])
        vlist_usace.append(stage_2_flow_rating(float(value.text),rating_curve))

    np.savetxt("gage_data/"+site_no+".csv", np.c_[dlist_usace,vlist_usace], fmt="%s,%s")

    return dlist_usace, vlist_usace

def NWS_gage_data_request_forecast(site_name,site_no,rating_curve):
    """only forecasted stages are provided, we are converting them to flows here with our own rating curves"""
    nws_xml=urlopen("https://water.weather.gov/ahps2/hydrograph_to_xml.php?gage="+site_name+"&output=xml&time_zone=cdt") 
    root = ET.parse(nws_xml).getroot()
    dlist_nws_forecast=[]
    vlist_nws_forecast=[]
    for type_tag in root.findall('forecast/datum'):
        dlist_nws_forecast.append(type_tag.find('valid').text[:-15]+" "+type_tag.find('valid').text[-14:-12]+":00")
        vlist_nws_forecast.append(stage_2_flow_rating(float(type_tag.find('primary').text),rating_curve))

    np.savetxt("gage_data/"+site_no+"_forecast.csv", np.c_[dlist_nws_forecast,vlist_nws_forecast], fmt="%s,%s")

    return dlist_nws_forecast, vlist_nws_forecast

        
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
        try:
            shutil.copyfile(destination+'gz/'+i, destination+'../temp/'+i)
            with gzip.open(destination+'../temp/'+i, 'rb') as f_in:
                with open(destination+'grib/'+i[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                f_out.close()
            f_in.close()
            os.remove(destination+'../temp/'+i)
        except:
            print("could not process "+i)
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
    try:
        rc = rascontrol.RasController(version='60')
    except:
        try:
            rc = rascontrol.RasController(version='5X')
        except:
            raise
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
            
    new_lines[startline-1]="Interval="+ras_timestep(gage)+'\n'
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
            
    new_lines[startline-1]="Interval="+ras_timestep(gage)+'\n'
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

def extract_2D_RAS(filename):

    f1 = h5py.File(filename,'r+')
    points=f1['Geometry']['2D Flow Areas']['test']['Cells Center Coordinate'][()]
    maxSWL=f1['Results']['Unsteady']['Output']['Output Blocks']['Base Output']['Summary Output']['2D Flow Areas']['test']['Maximum Water Surface'][()]
    cellmin=f1['Geometry']['2D Flow Areas']['test']['Cells Minimum Elevation'][()]
    f1.close()
    maxDEP=maxSWL[0]-cellmin
    xyz=np.concatenate((points,np.array([maxDEP.tolist()]).T),axis=1)
    xyz_nan=[]
    for i in xyz:
        if (math.isnan(i[2])) or (i[2]==0.0):
            pass
        else:
            xyz_nan.append(i.tolist())

    xyz_nan=np.array(xyz_nan)
    np.savetxt(filename+".xyz",xyz_nan,delimiter=',')

    
def xyz2tif(filename):

    kwargs_1 = {
        'format': 'GTiff',
        'outputType': gdal.GDT_Int16
    }
    kwargs_2 = {    
        'srcSRS': 'ras\Eq_Albs_USGS.prj',
        'dstSRS': 'EPSG:3857'
    }

    xyz_file = filename+'.hdf.xyz'
    tif_file1 = filename+'1.hdf.tif'
    tif_file2 = filename+'2.hdf.tif'


    step_1 = gdal.Translate(tif_file1, 
                            xyz_file,
                            **kwargs_1)

    step_1 = None

    step_2 = gdal.Warp(tif_file2,
                       tif_file1,
                       **kwargs_2)
    step_2 = None

def tif2png(filename):
    kwargs = {
        'colorFilename': 'cmap.txt',
        'addAlpha':'True',
    }

    
    tif_file2 = filename+'2.hdf.tif'
    png_file = filename+'.hdf.png'
    step_3 = gdal.DEMProcessing(png_file,
                                tif_file2,
                                "color-relief",
                                **kwargs)

    step_3 = None

def href_list(url):
    link_list=[]
    with urlopen(url) as website:
        for i in str(website.read()).split():
            if 'href' in i:
                if i.split('"')[1] == '../':
                    pass
                else:
                    link_list.append(i.split('"')[1][:-1])
    website.close()
    return link_list

def show_progress(block_num, block_size, total_size):
    global pbar
    if pbar is None:
        pbar = progressbar.ProgressBar(maxval=total_size)
        pbar.start()

    downloaded = block_num * block_size
    if downloaded < total_size:
        pbar.update(downloaded)
    else:
        pbar.finish()
        pbar = None

def retrieve_recent_advisory_LSU(year,event,mesh,files=["maxele.63.nc"]):
    global pbar
    loni=False
    for adv in reversed(href_list("https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/")):
        if mesh in href_list("https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv):
            try:
                sim=href_list("https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv+"/"+mesh+"/supermic.hpc.lsu.edu/")
                tracks=href_list("https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv+"/"+mesh+"/supermic.hpc.lsu.edu/"+sim[0]+"/")
                break
            except:
                print ("no supermic.hpc.lsu.edu run, try loni")
                try:
                    sim=href_list("https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv+"/"+mesh+"/qbc.loni.org/")
                    tracks=href_list("https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv+"/"+mesh+"/qbc.loni.org/"+sim[0]+"/")
                    loni=True
                    break
                except:
                    pass
        else:
            pass
    for track in tracks:
        if loni == False:
            files_path="https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv+"/"+mesh+"/supermic.hpc.lsu.edu/"+sim[0]+"/"+track
        else:
            files_path="https://fortytwo.cct.lsu.edu/thredds/fileServer/"+year+"/"+event+"/"+adv+"/"+mesh+"/qbc.loni.org/"+sim[0]+"/"+track            
        print(files_path)
        makedirs("./surge/"+year+"/"+event+"/"+adv+"/"+track+"/", exist_ok=True)    
        for file in files:
            try:
                print("downloading "+file)
                pbar = None
                urllib.request.urlretrieve(files_path+"/"+file,"surge/"+year+"/"+event+"/"+adv+"/"+track+"/"+file,show_progress)
            except:
                print("could not retrieve data for "+file)
    return "surge/"+year+"/"+event+"/"+adv+"/"+track+"/"

def stage_2_flow_rating(stage,polynomial_coefficents):
    return np.round(polynomial_coefficents[0]*stage**2+polynomial_coefficents[1]*stage+polynomial_coefficents[2],2)


def construct_adcirc_date(base_date,time_stamps_from_nc):
    converted_ts=[]
    for time_stamp in time_stamps_from_nc:
        converted_ts.append((base_date+datetime.timedelta(seconds=time_stamp)).strftime('%Y-%m-%d %H:%M:%S'))
        
    return converted_ts

def get_adcirc_time_series(latitude,longitude,fort_63_path,units='ft'):
    fn = fort_63_path
    ds = nc.Dataset(fn)
    nodes = [list(a) for a in zip(ds['y'][:],ds['x'][:])]
    pt=[latitude,longitude]
    node_index=spatial.KDTree(nodes).query(pt)[1]
    time=ds['time'][:].tolist()
    zeta_time=ds['zeta'][:].data[:,node_index].tolist()
    if units == 'ft':
        zeta_time=ds['zeta'][:].data[:,node_index]*3.2808
        zeta_time=zeta_time.tolist()
    else:
        zeta_time=ds['zeta'][:].data[:,node_index].tolist()
    new_time=construct_adcirc_date(datetime.datetime.strptime(ds['time'].base_date, '%Y-%m-%d %H:%M:%S'),time)
    return new_time, zeta_time

def extract_grib_precip_list(grib_path):
    with xarray.open_dataset(grib_path,engine='cfgrib') as ds:
        parsed_data=np.flip(ds.tp.data,axis=0).flatten()
    ds.close()
    return parsed_data

def extract_grib_precip_list_qpf(grib_path):
    return np.flip(transcribe_qpf_to_qpe_grid(grib_path),axis=0).flatten()

def extract_grib_timestamp(grib_path):
    with xarray.open_dataset(grib_path,engine='cfgrib') as ds:
        grib_timestamp=str(ds.time.data).split("T")
    ds.close()
    grib_timestamp[1]=grib_timestamp[1].split(".")[0]
    return grib_timestamp[0]+" "+grib_timestamp[1]    
    
def build_timestamp_list(grib_path_list):
    timestamp_list=np.array([],dtype='|S19')
    for path in grib_path_list:
        timestamp_list=np.append(timestamp_list,extract_grib_timestamp(path).encode())
    return timestamp_list

def build_timestamp_list_qpf(grib_path_list):
    timestamp_list=np.array([],dtype='|S19')
    
    file_times=[]
    forecast_times=[]
    for i in grib_path_list:
        file_times.append(int(i[:-4].split("p06m_")[1].split("f")[0]))
        forecast_times.append(i[:-4].split("p06m_")[1].split("f")[1])
    file_times=list(set(file_times))
    file_times.sort()
    
    base_date=datetime.datetime.strptime(str(file_times[-1]),'%Y%m%d%H')

    for times in forecast_times:
        timestamp_list=np.append(timestamp_list,(base_date+datetime.timedelta(hours=int(times))).strftime('%Y-%m-%d %H:%M:%S').encode())
    return timestamp_list

def build_precip_data_array(grib_path_list,qpf_path_list=None):
    """this builds a cummulative rainfall time series out of the grib datasets"""
    precip_array=np.array([],dtype="float32")
    for path in grib_path_list:
        if len(precip_array) == 0:
            precip_array=np.array([np.zeros(len(extract_grib_precip_list(path)))],dtype="float32")
        else:
            precip_array=np.concatenate((precip_array,[extract_grib_precip_list(path)+precip_array[-1]]))
            
    if qpf_path_list!=None:
        for path in qpf_path_list:
            precip_array=np.concatenate((precip_array,[extract_grib_precip_list_qpf(path)+precip_array[-1]]))
    
    return precip_array

def write_unsteady_hdf(timelist,precip_array,file_name):
    shutil.copyfile('templates/event_conditions_template.hdf', file_name)
    f1=h5py.File(file_name,'r+')
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values'].resize(0,axis=0)
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values'].resize(len(timelist),axis=0)
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values'][:]=precip_array
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values (Vertical)'].resize(0,axis=0)
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values (Vertical)'].resize(len(timelist),axis=0)
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values (Vertical)'][:]=precip_array
    f1['Event Conditions']['Meteorology']['Precipitation']['Imported Raster Data']['Values'].attrs['Times']=timelist
    f1.close()
    
def qpe_to_hdf(folder_path,hdf_path):
    times=build_timestamp_list(glob.glob(folder_path+'\*grib'))
    pdata=build_precip_data_array(glob.glob(folder_path+'\*grib'))
    write_unsteady_hdf(times,pdata,hdf_path)
    
def find_recent_qpf(qpf_dir):
    qpf_files=glob.glob(qpf_dir+'\*.grb')
    file_times=[]
    forecast_times=[]
    for i in qpf_files:
        file_times.append(int(i[:-4].split("p06m_")[1].split("f")[0]))
        forecast_times.append(i[:-4].split("p06m_")[1].split("f")[1])
    file_times=list(set(file_times))
    file_times.sort()
    forecast_times=list(set(forecast_times))
    forecast_times.sort()
    return file_times[-1],forecast_times,glob.glob(qpf_dir+'\\*'+str(file_times[-1])+'*.grb')

def build_timestamp_list_qpf(grib_path_list):
    timestamp_list=np.array([],dtype='|S19')
    
    file_times=[]
    forecast_times=[]
    for i in grib_path_list:
        file_times.append(int(i[:-4].split("p06m_")[1].split("f")[0]))
        forecast_times.append(i[:-4].split("p06m_")[1].split("f")[1])
    file_times=list(set(file_times))
    file_times.sort()
    
    base_date=datetime.datetime.strptime(str(file_times[-1]),'%Y%m%d%H')

    for times in forecast_times:
        timestamp_list=np.append(timestamp_list,
                                 (base_date+datetime.timedelta(hours=int(times))).strftime('%Y-%m-%d %H:%M:%S').encode())
    return timestamp_list

def transcribe_qpf_to_qpe_grid(qpf_file):
    qpe_to_qpf_map=np.loadtxt('./precip/qpe_to_qpf_map.csv')    
    with xarray.open_dataset(qpf_file , engine='cfgrib') as qpf:
        qpf
    qpf_precip=qpf.tp.data.flatten()
    forecast_on_qpe_grid=np.empty(shape=(419*419))
    for i in range(len(forecast_on_qpe_grid)):
        forecast_on_qpe_grid[i]=qpf_precip[int(qpe_to_qpf_map[i])]
    qpf.close()
    return forecast_on_qpe_grid.reshape(419, 419)
    
def trim_qpf_files(latest_observed_timestamp,qpf_file_list):
    latest_qpe_time=datetime.datetime.strptime(latest_observed_timestamp,'%Y-%m-%d %H:%M:%S')
    qpf_timestamp_list=build_timestamp_list_qpf(qpf_file_list)
    trimmed_file_list=[]
    for i in range(len(qpf_timestamp_list)):
        if datetime.datetime.strptime(qpf_timestamp_list[i].decode(),'%Y-%m-%d %H:%M:%S')>latest_qpe_time:
            trimmed_file_list.append(qpf_file_list[i])
    return trimmed_file_list
    
def full_precip_to_hdf(qpe_path,qpf_path,hdf_path):
    qpe_times=build_timestamp_list(glob.glob(qpe_path+'\*grib'))
    qpf_files=find_recent_qpf(qpf_path)[2]  
    qpf_files=trim_qpf_files(qpe_times[-1].decode(),qpf_files)
    qpf_times=build_timestamp_list_qpf(qpf_files)
    
    times=np.append(qpe_times,qpf_times)
    pdata=build_precip_data_array(glob.glob(qpe_path+'\*grib'),qpf_files)
    write_unsteady_hdf(times,pdata,hdf_path)
    
def ras_timestep(gage_id):
    ras_timestep_dict={str(360):'6MIN',
                  str(600):'10MIN',
                  str(900):'15MIN',
                  str(1800):'30MIN',
                  str(3600):'1HOUR',
                  str(21600):'6HOUR',
    }
    
    rawdata=np.loadtxt("gage_data\\"+gage_id+".csv", dtype=str, delimiter=',')
    tstep=(datetime.datetime.strptime(rawdata[:,0][-1],'%Y-%m-%d %H:%M')
           -datetime.datetime.strptime(rawdata[:,0][-2],'%Y-%m-%d %H:%M')).seconds
    return ras_timestep_dict[str(tstep)]

def maxele2csv(ncfile, bbox=None):
    """Inputs:
    ncfile = path and filename of netCDF for maxele
    bbox = bounding box for area, [lon1, lon2, lat1, lat2], lat/lon listed respectively in ascending order. 
    Input for bbox is optional is no box (default)."""
    project_dir = os.path.dirname(ncfile)

    sim0 = []        

    surge = nc.Dataset( ncfile)

    sim0.append(surge.variables['x'][:])
    sim0.append(surge.variables['y'][:])
    sim0.append(surge.variables['zeta_max'][:])

    dum0 = np.asarray(sim0) # put into array

    # Put list into dataframe
    lonlat = pd.DataFrame(np.transpose(dum0[0:2,:]))

    dum0 = np.nanmax(dum0[2:,:],axis=0) # find max of zetamax at each location for all the simulations
    dum0[dum0 == -99999] = np.nan # set -99999 values to be NaN
    zetamax = dum0*3.2808399 # convert to ft
    zetamax[np.isnan(zetamax)] = -999999.
    zetamax = pd.DataFrame(zetamax)


    # combine the lon/lat columns with the max elevation to one dataframe
    surgemax = pd.concat([lonlat,zetamax],axis=1) 
    # name the columns of the dataframe for csv output
    surgemax.columns = ['Lon','Lat','SWL'] 

    if bbox == None:
        pass
    else:    
        # drops the values so the comparison signs are opposite of the matlab source code
        surgemax.drop(surgemax[surgemax['Lon'] < bbox[0]].index, inplace = True) 
        surgemax.drop(surgemax[surgemax['Lon'] > bbox[1]].index, inplace = True)
        surgemax.drop(surgemax[surgemax['Lat'] < bbox[2]].index, inplace = True)
        surgemax.drop(surgemax[surgemax['Lat'] > bbox[3]].index, inplace = True)
    
    pd.DataFrame.to_csv(surgemax, project_dir + "\\Surgeout" + ".csv", index=False)