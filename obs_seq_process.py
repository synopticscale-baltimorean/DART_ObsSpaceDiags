#~#~#~#~~#~#~#~#~#~#~#~~#~#~#~#
#obs_seq.py 

#module for loading and filtering obs_seq files 
#and outputing data to pandas df or np arrays

#dependencies: 

#shell files
#input obs_seq files: process_obs_seq.sh
#processed obs_seq.final files w/ mean output only: process_obs_seq_final_outputmem.sh, process_obs_seq_final_outputmem20.sh
#processed obs_seq.final files / member output: process_obs_seq_final_member.sh, process_obs_seq_final_allmems.sh

#~#~#~#~~#~#~#~#~#~#~#~~#~#~#~#

#INSTRUCTIONS: --set system paths in filepaths.txt--

#directory with stored ensemble runs:
#top_dir=/glade/scratch/jmccurry/WOF/realtime
#directory with stored forecast output
#forecast_dir=/glade/campaign/univ/umcp0011/WOF_FORECAST_ARCHIVE
#directory where you want to do data analysis
#work_dir=/glade/u/home/jmccurry/colab_scripts
####################################
#import required modules -if you run into any 'module not found' errors 
#on cheyenne or casper make sure you did 'module load ncarenv' and then 'ncar_pylib'
from pylab import *
import numpy as np 
import sys
from os import path 
from os import system
import os
import subprocess
from scipy.interpolate import interp1d
import math
from netCDF4 import Dataset as netcdf_dataset
import wrf
from wrf import to_np, vertcross, CoordPair
#######################################
class obs_seq:
   #class for input,unassimilated obs sequence files
   def __init__(self,timestamp,obs_seq):
    with open('user_filepaths.txt') as f:
        lines = f.read().split('\n')
    self.filter_flag=0 #has this obs_seq object been filtered

    self.obs_dir=lines[3]
    self.obs_types = {'ACARS_U_WIND_COMPONENT':16, 'ACARS_V_WIND_COMPONENT':17,'ACARS_TEMPERATURE':18, \
       'DOPPLER_RADIAL_VELOCITY':36,'COMBINED_2M_TEMPERATURE':[42,27], \
       'RADAR_REFLECTIVITY':37, 'RADAR_CLEARAIR_REFLECTIVITY':38,'COMBINED_SFC_U_WIND_COMPONENT':[40,25]  \
       ,'LAND_SFC_ALTIMETER':74, 'METAR_ALTIMETER':75,'COMBINED_SFC_V_WIND_COMPONENT':[26,41]}
    self.orig_data = load_obs_seq(timestamp,obs_seq,obsdir=self.obs_dir)
    self.varnames = dict([('type',0),('obs_type',0),('value',1),('vals',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('radar_xloc',7),('radar_yloc',8)])
    #dictionary containing radar locations in dart radian coords 
    #if you want to do selection of radar obs for only one radar then 
    #make sure it's added onto this list
    self.radar_locs = dict([('KVNX',[4.57053374392750,0.6412447157008800]),('KOAX',[4.601266986004480,0.7211764997633801]),('KEAX',[4.637959886345540,0.6773666343042200]),('KDVN',[4.702251524537250,0.7262606074423900])])
   def data_check(self):
      if(self.filter_flag==0):
        self.data = np.copy(self.orig_data)
      else:
        pass
   def filter_data_type(self,obs_type):
   #select data of a single observational type
   #usage: df = obs_seq.filter_data_type(36)
      self.data_check()
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      self.filter_flag = 1
      return self 
   def filter_data_type_list(self,obs_type_list):
   #select data of multiple observational types
   #usage: df = obs_seq.filter_data_type([36,45])
      self.data_check()
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      self.filter_flag = 1
      return self 
   def filter_outliers(self,var,low_thresh,high_thresh):
   #select data within a specified range for a particular obs variable
   #make sure to filter data type first 
   #usage: df = obs_seq.filter_data_type(37).filter_outliers('value',40,50) -- selects radar obs between 40 and 50 dbz
      self.data_check()
      self.data[:,self.varnames[str(var)]][self.data[:,self.varnames[str(var)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(var)]][self.data[:,self.varnames[str(var)]]> high_thresh] = 'NaN' 
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(var)]])]
      self.filter_flag = 1
      return self 
   def filter_radar_location(self,radar):
   #select data from a single radar location
      self.data_check()
      self.data[:,self.varnames['radar_xloc']][self.data[:,self.varnames['radar_xloc']]!=self.radar_locs[radar][0]] = 'NaN'  
      self.data[:,self.varnames['radar_yloc']][self.data[:,self.varnames['radar_yloc']]!=self.radar_locs[radar][1]] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_xloc']])]
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_yloc']])]
      self.filter_flag = 1

      return self 
   def return_np(self,var):
   #return 1D np array containing desired variable
   #usage: radar_vals_arr = obs_seq.filter_data_type(37).filter_outliers('value',40,50).return_np('value')
      self.data_check()

      data = self.data[:,self.varnames[str(var)]]
      return data 
      self.filter_flag = 0

class obs_seq_final:
    #class for obs_seq.final files 
    #For cycle means: obs_seq_final(timestamp,obs sequence file specifier eg. 'combined_full_cropped'\
    #,experiment name eg. '20190703.64_mem_3km_pf_hybrid')
    #For forecast means: obs_seq_final(forecast timestamp,obs seq file specifier,forecast=forecast init timestamp,\
    #forecast_subdir=directory under main forecast archive eg. 'hybrid/25dbz_64_mem_fullloc') NOTE: only works for \
    #forecast ensemble sizes of 20 members currently
    #For forecast members: 
    #for all forecast members:
    #cannot get cycle members currently


   def __init__(self,timestamp,obs_seq,forecast_subdir='',experiment_name='',forecast='no',outputmem='',member='no'):
    with open('user_filepaths.txt') as f:
        lines = f.read().split('\n')
    self.filter_flag=0 #has this obs_seq object been filtered
    self.top_dir=lines[0]
    self.obs_types = {'ACARS_U_WIND_COMPONENT':16, 'ACARS_V_WIND_COMPONENT':17,'ACARS_TEMPERATURE':18, \
   'DOPPLER_RADIAL_VELOCITY':36,'COMBINED_2M_TEMPERATURE':[42,27], \
   'RADAR_REFLECTIVITY':37, 'RADAR_CLEARAIR_REFLECTIVITY':38,'COMBINED_SFC_U_WIND_COMPONENT':[40,25]  \
   ,'LAND_SFC_ALTIMETER':74, 'METAR_ALTIMETER':75,'COMBINED_SFC_V_WIND_COMPONENT':[26,41]}
    display(self.top_dir)
    self.forecast_dir=lines[1]
    self.work_dir=lines[2]
    self.varnames = dict([('type',0),('obs_type',0),('value',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('prior',7),('posterior',8),('radar_xloc',9),('radar_yloc',10),('obs_err',11)])
    if ( forecast=='no'):
     display('{}/{}/{}/obs_seq.final.{}.obs_seq.{}.{}'.format(self.top_dir,experiment_name,timestamp,timestamp,obs_seq,timestamp))
     self.orig_data = load_obs_final_cycle_mean(timestamp,obs_seq,self.top_dir,experiment_name)
    #need a few lines to handle user not input-ing experiment name when forecast=='no'
    
    else:
        if ( member=='no'):
           self.orig_data = load_obs_final_forecast_mean(timestamp,obs_seq,forecast,forecastdir=self.forecast_dir,forecastsubdir=forecast_subdir)
        elif( member=='all'):
           self.orig_data = load_obs_final_forecast_allmems(timestamp,obs_seq,forecast,forecastdir=self.forecast_dir,forecastsubdir=forecast_subdir) 
           self.varnames = dict([('obs_type',0),('value',1),('vals',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('prior',7),('post',8),('prior_sp',9),('post_sp',10),('obs_err',11),('priorA',7),('priorB',8),('priorC',9),('priorD',10),('priorE',11),('priorF',12),('priorG',13),('priorH',14),('priorI',15),('priorJ',16),('priorK',17),('priorL',18),('priorM',19),('priorN',20),('priorO',21),('priorP',22),('priorQ',23),('priorR',24),('priorS',25),('priorT',26)])
        else:
           self.orig_data = load_obs_final_forecast_member(timestamp,obs_seq,forecast,(4+2*int(member)),outputmem,forecastdir=self.forecast_dir,forecastsubdir=forecast_subdir) 
   def data_check(self):
      if(self.filter_flag==0):
        self.data = np.copy(self.orig_data)
      else:
        pass
   def filter_data_QC(self,QC):
      self.data_check()
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,5]==QC)]
      self.filter_flag=1

      return self
   def filter_data_QC_list(self,QC):
      self.data_check()
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,5],QC))]
      self.filter_flag=1

      return self 
   def filter_data_type(self,obs_type):
      self.data_check()
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      self.filter_flag=1
      return self 
   def filter_outliers(self,var,low_thresh,high_thresh):
      self.data_check()
      self.data[:,self.varnames[str(var)]][self.data[:,self.varnames[str(var)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(var)]][self.data[:,self.varnames[str(var)]]> high_thresh] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(var)]])]
      self.filter_flag=1

      return self 
   def filter_outliers_range(self,z1,z2,low_thresh,high_thresh):
      self.data_check()
      for z in range(z1,z2):
          self.data[:,z][self.data[:,z]< low_thresh] = 'NaN'  
          self.data[:,z][self.data[:,z]> high_thresh] = 'NaN'  
          self.data = self.data[~np.isnan(self.data[:,z])]
      self.filter_flag=1

      return self 
   def return_np(self,var):
      self.data_check()
      data = self.data[:,self.varnames[str(var)]]
      self.filter_flag=0
      return data 
   def crop(self):
      pass
   def thin_z(self):
      pass
   def thin_xy(self):
      pass
   def write_to_obs(self):
      pass 

   
def load_obs_final_cycle_mean(timestamp,obs_seq,topdir,experiment_name):

   procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
   if (1==1):
        path='{}/{}/{}/obs_seq.final.{}.obs_seq.{}.{}'.format(topdir,experiment_name,timestamp,timestamp,obs_seq,timestamp)
        subprocess.call(['{}/process_obs_final_outputmem.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')
   return data
def load_obs_final_forecast_mean(timestamp,obs_seq,forecast_init,forecastdir,forecastsubdir):

   procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
   if (1==1):
        path='{}/{}/WRFOUTS_FCST{}/{}/obs_seq.verify.{}.obs_seq.{}.{}'.format(forecastdir,forecastsubdir,forecast_init,timestamp,timestamp,obs_seq,timestamp)
        
        subprocess.call(['{}/process_obs_final_outputmem20.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')
   return data
   
def load_obs_final_forecast_member(timestamp,obs_seq,forecast_init,line_member,outputmem,forecastdir,forecastsubdir):

   procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
   if (1==1):
        path='{}/{}/WRFOUTS_FCST{}/{}/obs_seq.verify.{}.obs_seq.{}.{}'.format(forecastdir,forecastsubdir,forecast_init,timestamp,timestamp,obs_seq,timestamp)

        subprocess.call(['{}/process_obs_final_member.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),'{}'.format(line_member),'{}'.format(line_member + 1)])
   data = np.loadtxt(procname,delimiter=',')
   return data
def load_obs_final_forecast_allmems(timestamp,obs_seq,forecast_init,forecastdir,forecastsubdir):
   procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
   if (1==1):
        path='{}/{}/WRFOUTS_FCST{}/{}/obs_seq.verify.{}.obs_seq.{}.{}'.format(forecastdir,forecastsubdir,forecast_init,timestamp,timestamp,obs_seq,timestamp)
        subprocess.call(['{}/process_obs_final_allmems.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')

   return data
def load_obs_seq(timestamp,obs_seq,obsdir):
   procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),obs_seq,timestamp)
   if (1==1):
        path='{}/obs_seq.{}.{}'.format(obsdir,obs_seq,timestamp)
        subprocess.call(['{}/process_obs_seq.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')
   return data
def main():
   pass
if __name__ == "__main__":
    main()
