#!/usr/bin/python
import pandas as pd
import numpy as np
import scipy.stats
import math

#Read Text file"DecayTimecourse.txt" via pd Dataframe and replace nan values with zero.

Data=pd.read_csv("DecayTimecourse.txt",delimiter=('\t'),index_col=False)
Data=Data.replace(np.nan, 0)
Data.to_csv("DecayTimecourse.csv")

# Extracted the First_time_course data  from Data and find out value of decay constant (slope) through scipy stats and calculated half life.

First_Timecourse_Data=Data.iloc[0:,1:10]
Timecourse_Data_Regr1= pd.DataFrame() # created database to store stats calculated values.
for i,row in First_Timecourse_Data.iloc[1:,:].iterrows(): #interate row one by one
      r= scipy.stats.linregress(First_Timecourse_Data.iloc[0,:],row) #linear regression of each gene expression row with first row(Time course).
      Timecourse_Data_Regr1=pd.concat([Timecourse_Data_Regr1,pd.DataFrame(r._asdict(),index=[i])])# all stats values stored in the database.
Slope_Data1 = Timecourse_Data_Regr1[['slope']].copy()# copied the value of slope to Slope_Data1 from Timecourse_Data_Regr1 database
Slope_Halflife_Data1=Slope_Data1.assign(Halflife_Timecourse1=math.log(.5)/Slope_Data1.slope)#Calculated the value of halflife using decay constant(slope)
Slope_Halflife_Data1.rename(columns={'slope': 'Decay constant(lambda)_1'}, inplace=True)# changed the name of column from slope to decay constant.


# Extracted the Second_time_course data  from Data  and find out value of decay constant (slope) through scipy stats and calculated half life.
Second_Timecourse_Data=Data.iloc[0:,10:19]
Timecourse_Data_Regr2= pd.DataFrame() # created database to store stats calculated values.
for i,row in Second_Timecourse_Data.iloc[1:,:].iterrows():#interate row one by one
      r= scipy.stats.linregress(Second_Timecourse_Data.iloc[0,:],row)#linear regression of each gene expression row with first row(Time course).
      Timecourse_Data_Regr2=pd.concat([Timecourse_Data_Regr2,pd.DataFrame(r._asdict(),index=[i])])# all stats values stored in the database.
Slope_Data2= Timecourse_Data_Regr2[['slope']].copy()# copied the value of slope to Slope_Data1 from Timecourse_Data_Regr1 database
Slope_Halflife_Data2=Slope_Data2.assign(Halflife_Timecourse2=math.log(.5)/Slope_Data2.slope)#Calculated the value of halflife using decay constant(slope)
Slope_Halflife_Data2.rename(columns={'slope': 'Decay constant(lambda)_2'}, inplace=True)#Changed the name of column from slope to decay constant.

# Extracted the Third_time_course data  from Data  and find out value of decay constant (slope) through scipy stats and calculated half life.
Third_Timecourse_Data=Data.iloc[0:,19:29]
Timecourse_Data_Regr3= pd.DataFrame()# created database to store stats calculated values.
for i,row in Third_Timecourse_Data.iloc[1:,:].iterrows():#interate row one by one
      r= scipy.stats.linregress(Third_Timecourse_Data.iloc[0,:],row)#linear regression of each gene expression row with first row(Time course).
      Timecourse_Data_Regr3=pd.concat([Timecourse_Data_Regr3,pd.DataFrame(r._asdict(),index=[i])])# all stats values stored in the database.
Slope_Data3 = Timecourse_Data_Regr3[['slope']].copy()# copied the value of slope to Slope_Data1 from Timecourse_Data_Regr1 database
Slope_Halflife_Data3=Slope_Data3.assign(Halflife_Timecourse3=math.log(.5)/Slope_Data3.slope)#Calculated the value of halflife using decay constant(slope)
Slope_Halflife_Data3.rename(columns={'slope': 'Decay constant(lambda)_3'}, inplace=True)#Changed the name of column from slope to decay constant.

# Combined Decay constant and halflife data along with transcript names from previous steps.
Combine_Halflife_data = pd.concat([Data.iloc[1:,0],Slope_Halflife_Data1,Slope_Halflife_Data2,Slope_Halflife_Data3], axis=1)
Combine_Halflife_data.replace([np.inf, -np.inf], np.nan, inplace=True)# Removed all the infinity values.
Combine_Halflife_data= Combine_Halflife_data.fillna(0)# filled all nan (empty) values with zero.
# Calculated the Average of 3 time course half life data.
Combine_Halflife_data['Average'] = Combine_Halflife_data[['Halflife_Timecourse1', 'Halflife_Timecourse2','Halflife_Timecourse3']].mean(axis=1)
Combine_Halflife_data.to_csv("FinalHalflifeaverages.csv") # Saved the decay constant and half life information in FinalHalflifeaverages.csv file as output.

# Top 618 genes which are having high life compared to other tanscript (10% of 6184= 618).
Top_halflife = Combine_Halflife_data.nlargest(618, 'Average', keep='first')
Top_halflife.to_csv("Top10%_halflife.csv")# The output file is saved as Top10%_halflife.csv

#Botton 618 genes which are having low half life compared to other transcript ((10% of 6184= 618)
Bottom_halflife=Combine_Halflife_data.nsmallest(618, 'Average', keep='first')
Bottom_halflife.to_csv("Bottom10%_halflife.csv")#The output file is saved as Bottom0%_halflife.csv
