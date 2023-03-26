## This script contains the classes and functions for the execution of the Arusha
## Assignment in PCRaster (www.pcraster.geo.uu.nl)

#In the exercise markdown document, there may contain simple functions necessary.


#All the input files were previously translated to .map (PCraster format) and adjusted to coincide.
from pcraster import *
from pcraster.framework import *
import numpy as np
import pandas as pd
import datetime as datetime

class RechargeModel(DynamicModel):
    def __init__(self, land_use_map, start_date, final_date):
        DynamicModel.__init__(self)
        
        setclone(land_use_map)
        self.land_use_map = land_use_map
        self.start_date = start_date
        self.final_date = final_date
        
    
    def initial(self):
        land_use_folder = "LULC Input"
        other_folder = "Other input"
        self.land_use = self.readmap(os.path.join(land_use_folder, "landuse"))
        self.soil = self.readmap(os.path.join(other_folder, "soil"))
        self.dem = self.readmap(os.path.join(other_folder,"elevation"))
        self.slope = self.readmap(os.path.join(other_folder,"slope"))
        self.threshold = lookupscalar(os.path.join(other_folder,"land_use_threshold.tbl"),self.land_use)
        self.soilcapacity = lookupscalar(os.path.join(other_folder, "soil_cap.tbl"),self.soil)
        self.soilextdepth = lookupscalar(os.path.join(other_folder,"soil_ext_depth.tbl"),self.land_use)
        self.cropcoef = lookupscalar(os.path.join(other_folder,"crop_coef.tbl"),self.land_use)
        self.initial_su = lookupscalar(os.path.join(other_folder,"initial_sto_pct.tbl"),self.land_use)
        
        #Initial Soil Storage computation:
        self.su = self.initial_su * self.soilcapacity * self.soilextdepth
        self.recharge_list = []
        self.suet = self.su
                
        #Initial Runoff (= 0):
        self.runoff_acc = scalar(0)
        
        #Initial accumulated recharge (= 0):
        self.recharge_acc = scalar(0)
        
        #Initial accumulated actual evap:
        self.evap_acc = scalar(0)
        
        #Auxiliar for reporting the annual data (every 12 months):
        self.runoff_acc_ann_bef = scalar(0)
        self.recharge_acc_ann_bef = scalar(0)
        self.evap_acc_ann_bef = scalar(0)
        
        #%%
        # Now we set one counter wich will be useful for saving the results:
        self.counter = 0 # Auxiliar for saving the results every 12 months:
        self.day = 1
        self.month = 1
        
        # Here the dates are set:
        # Daily:
        self.dates = pd.date_range(start = self.start_date, end = self.final_date, freq='D')
        
        # Monthly:
        self.datesmonthly = pd.date_range(start = self.start_date, end = self.final_date, freq='M')
        
    def dynamic(self):
        
        # The year is computed:        
        self.year = self.dates[self.day - 1].year
        
        climate_folder = "Climate data"
        prec = timeinputscalar(os.path.join(climate_folder,"precipitation.tss"), 1)
        refET = timeinputscalar(os.path.join(climate_folder,"et.tss"), 1)
        potET = refET * self.cropcoef
        
        #%%
        # Water-balance computation:
        # Effective Precipitation:
        eff_prec = min(prec, self.threshold)
        
        # Runoff:
        self.runoff = prec - eff_prec
        self.runoff_acc = self.runoff_acc + self.runoff
        
        # Soil storage:
        self.su = self.su + eff_prec
        
        # Actual evapotranspiration:
        # AET is calculated based on model by Allen et al. (1998). It is computed based on PET and limited by soil storage.
        
        self.actET = ifthenelse(self.su >= self.suet, potET, potET * self.su/self.suet)
        self.su = self.su - self.actET
        
        self.evap_acc = self.evap_acc + self.actET
        
        
        # Recharge:
        self.recharge = max(0, (self.su - self.soilcapacity * self.soilextdepth))        
        self.recharge_acc += self.recharge
        
        self.su = self.su - self.recharge
        
        #%%
        # Monthly report (accumulated by month):
        if self.dates[self.day - 1] == self.datesmonthly[self.month - 1]:
            if self.month == 1:
                
                # The current accumulation is saved:
                self.runoff_acc_mon_bef = self.runoff_acc
                self.recharge_acc_mon_bef = self.recharge_acc
                self.evap_acc_mon_bef = self.evap_acc
                
                # Report:
                self.report(self.runoff_acc, os.path.join("output/rum" + str(self.year)))      
                self.report(self.recharge_acc, os.path.join("output/rem" + str(self.year)))
                self.report(self.evap_acc, os.path.join("output", "evm" + str(self.year))) 
                
                self.month = self.month + 1
                self.day = self.day + 1
                
                self.counter = self.counter + 1
            else:
                
                # Current month:
                self.runoff_acc_mon =  self.runoff_acc - self.runoff_acc_mon_bef
                self.recharge_acc_mon = self.recharge_acc - self.recharge_acc_mon_bef
                self.evap_acc_mon = self.evap_acc - self.evap_acc_mon_bef
                
                # The current accumulation is saved:
                self.runoff_acc_mon_bef = self.runoff_acc
                self.recharge_acc_mon_bef = self.recharge_acc
                self.evap_acc_mon_bef = self.evap_acc
                
                # Report:
                self.report(self.runoff_acc_mon, os.path.join("output/rum" + str(self.year)))      
                self.report(self.recharge_acc_mon, os.path.join("output/rem" + str(self.year)))
                self.report(self.evap_acc_mon, os.path.join("output", "evm" + str(self.year))) 
                
                self.month = self.month + 1 
                self.day = self.day + 1
                
                self.counter = self.counter + 1
        else:
            self.day = self.day + 1
        
        #%%
        # Yearly report (accumulated each 12 months):
        if self.counter == 12:
            
            # Current month:
            self.runoff_acc_ann =  self.runoff_acc - self.runoff_acc_ann_bef
            self.recharge_acc_ann = self.recharge_acc - self.recharge_acc_ann_bef
            self.evap_acc_ann = self.evap_acc - self.evap_acc_ann_bef
                
            # The current accumulation is saved:
            self.runoff_acc_ann_bef = self.runoff_acc
            self.recharge_acc_ann_bef = self.recharge_acc
            self.evap_acc_ann_bef = self.evap_acc
                
            # Report:
            self.report(self.runoff_acc_ann, os.path.join("output/ruy" + str(self.year)))      
            self.report(self.recharge_acc_ann, os.path.join("output/rey" + str(self.year)))
            self.report(self.evap_acc_ann, os.path.join("output", "evy" + str(self.year))) 
                
            self.counter = 0 

        #%%
        ## Total daily recharge table report:
        recharge2 = (maptotal(self.recharge)*100**2) / maparea(self.recharge)
        recharge2 = pcr2numpy(recharge2,-1)
        recharge2 = recharge2[50,50]
        self.recharge_list.append(recharge2)
        np.savetxt("output/totalrechargedaily.csv", self.recharge_list, delimiter=",")
    
