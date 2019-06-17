# Exercise 5
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import wget
import netCDF4 

# import functions from utils here

src = 'C:/Users/Ruth/Desktop/Eigene_Dateien/01_Studium/PhyGeo_Msc/Profilmodule/Python/exercise-5-RuthNvll/'
input_dir = src + "data/"
output_dir = src + "solution/"


# 1. Go to http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php 
#    datafiles and download the 0.25 deg. file for daily mean temperature. Save the 
#    file into the data directory but don't commit it to github. [2P]
url= "https://www.ecad.eu/download/ensembles/data/Grid_0.25deg_reg_ensemble/tg_ens_mean_0.25deg_reg_v19.0e.nc"

#----------------------------------------------------------------------
# Everything below this line could not be executed, due to a problem with loading the data.
# So everything below is just theory. 

wget.download(url, out = input_dir/ "0.25deg_reg_v19.0e.nc"  )

# 2. Read the file using xarray. Get to know your data. What's in the file?
#    Calculate monthly means for the reference periode 1981-2010 for Europe 
#    (Extent: Lon_min:-13, Lon_max: 25, Lat_min: 30, Lat_max: 72). [2P]

tab = xr.open_dataset(input_dir / "tg_ens_mean_0.25deg_reg_v19.0e.nc")
print tab.dimensions.keys()

tab_sel1 = tab.sel(time = slice ("1981-01-01", "2010-12-31"), 
                   latitude = slice (30,72), 
                   longitude = slice (-13,25))
tab_sel1_mean = tab_sel1.groupby("time.month").mean("time")


# 3. Calculate monthly anomalies for 2018 for the reference period and extent in #2.
#    Make a quick plot of the anomalies for the region. [2P]

tab_sel2 = tab.sel(time = slice ("2018", "2018"), 
                   latitude=slice(30,72), 
                   longitude=slice(-13,25))
tab_sel2_mean = tab_sel2.groupby("time.month").mean("time")
anom2018 = tab_sel2_mean - tab_sel1_mean

anom2018["tg"].plot(x = "longitude", col = "month")


# 4. Calculate the mean anomaly for the year 2018 for Europe and compare it to the 
#    anomaly of the element which contains Marburg. Is the anomaly of Marburg lower 
#    or higher than the one for Europe? [2P] 

meanAnamalie = anom2018.mean()["tg"]
marburg = [8.77, 50.8]
anomalieMarburg = anom2018.sel(latitude = marburg[1],
                                    longitude = marburg[0],
                                    method = "nearest").mean()

if meanAnamalie > anomalieMarburg:
    print("The mean anomalie for Europe is higher than in Marburg")
else:
    print("The mean anomalie for Marburg is higher than for Europe")


# 5. Write the monthly anomalies from task 3 to a netcdf file with name 
#    "europe_anom_2018.nc" to the solution directory.
#    Write the monthly anomalies for Marburg to a csv file with name 
#    "marburg_anom_2018.csv" to the solution directory. [2P]

anom2018.to_netcdf(output_dir / "europe_anomalie_2018.nc")
Marburg2018 = anom2018.sel(latitude = 50.81, 
                           longitude = 8.77, 
                           method = "nearest").to_dataframe()["tg"]
Marburg2018.to_csv(output_dir / "marburg_anom_2018.csv")
