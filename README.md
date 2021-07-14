# hhc_tools
Engineering tools for acquisition and processing of hydraulic data.



Setup instructions unsing anaconda:

conda create -n "hhc_tools" python=3.9 ipython
conda activate hhc_tools

#The lines below can be placed into a batch (.bat) file and run from the hhc_tools environment.

echo y|conda install jupyter
echo y|conda install numpy
echo y|conda install scipy
echo y|conda install ipyleaflet
echo y|conda install -c conda-forge ipyleaflet
echo y|conda install -c conda-forge voila
echo y|conda install matplotlib
echo y|conda install -c conda-forge bqplot
echo y|conda install -c conda-forge pygrib
echo y|conda install -c conda-forge xarray
echo y|conda install -c conda-forge cfgrib
echo y|conda install -c conda-forge h5netcdf
echo y|conda install -c conda-forge progressbar
echo y|conda install -c conda-forge netcdf4
echo y|conda install -c menpo wget
echo y|conda install -c conda-forge gdal
echo y|conda install -c conda-forge ffmpeg
echo y|conda install -c conda-forge certifi

git clone https://github.com/mikebannis/rascontrol.git
cd rascontrol
pip install .
cd ..
