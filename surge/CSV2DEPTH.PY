print("Importing arcpy...")
import arcpy
import os
arcpy.env.overwriteOutput = True
print("Done!")

wrkdir=os.getcwd()

def csv_files():
    global csv_files
    csv_files=[]
    for file in os.listdir("."):
        if file.endswith(".csv"):
        	csv_files.append(file)
    return csv_files

def find_prj():
    global prjctns
    prjctns=[]
    num=0
    for file in os.listdir("."):
        if file.endswith(".prj"):
		num+=1
        	prjctns.append([num,file])
    return prjctns

def get_prj():
    global prj_path
    prj_ = open('prj_path.txt',"r+")
    prj_path=prj_.readline()
    prj_.close()
    return prj_path

def find_shape():
	shps=[]
	try:
        	for file in os.listdir(".\\"):
			if file.endswith(".shp"):
        			shps.append(file)
	except:
		print("Copy over point shape")
    	return shps



try:
	prjpath=get_prj()
except:
	print("Projection not set!")
	find_prj()
	print prjctns
	prj_ind = int(raw_input("Which projection to use for ADCIRC data? (type = 1,2,3..)"))-1
	prj=prjctns[prj_ind][1]
	try:
		with open('prj_path.txt', 'a') as the_file:
    			the_file.write(str(wrkdir+"\\"+prj))
		prjpath=get_prj()
	except:
		print("failed")

def make_TIN(in_,prj):
    arcpy.MakeXYEventLayer_management(in_,"Lon","Lat","temp",prj,"SWL")
    arcpy.CreateTin_3d(out_tin=wrkdir+"\\"+in_[:-4],
                       spatial_reference=wrkdir+"\\Eq_Albs_USGS.prj",
                       
                       in_features="temp SWL Mass_Points <None>",
                       constrained_delaunay="CONSTRAINED_DELAUNAY")

def TIN2TIF(tin):
    arcpy.env.extent = arcpy.Extent(624127.347083, 2095613.485647, 2684127.347083, 3259613.485647)
    arcpy.TinRaster_3d(in_tin=wrkdir+"\\"+tin,
                       out_raster=wrkdir+"\\"+tin+".tif",
                       data_type="FLOAT",
                       method="LINEAR",
                       sample_distance="CELLSIZE 200",
                       z_factor="1")

def Calculate_depth(surge,terrain):
    arcpy.Minus_3d(in_raster_or_constant1=surge,
                in_raster_or_constant2=terrain,
                out_raster=wrkdir+"\\"+name+"_depth.tif")    

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

csv_files()
name=csv_files[0][:-4]

print (str(len(csv_files))+" files detected")

print("converting "+name+" from shape to TIN...")
make_TIN(name+".csv",prjpath)

print("converting "+name+" from TIN to TIf...")
TIN2TIF(name)

print("Calculating depth...")
Calculate_depth(wrkdir+"\\"+name+".tif",wrkdir+"\\SlaMM200ft.tif")

print("Removing Nulls...")
remove_null=arcpy.sa.SetNull(wrkdir+"\\"+name+"_depth.tif",wrkdir+"\\"+name+"_depth.tif","VALUE < 0.05")

print("Exporting final depth as: "+name+"_depth_final.tif")
remove_null.save(wrkdir+"\\"+name+"_depth_final.tif")

arcpy.CheckInExtension("3D")
arcpy.CheckInExtension("Spatial")
