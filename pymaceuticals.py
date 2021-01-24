#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 1) The drugs Capomulin and Ramicane appear to be significantly more effective than Infubinol and Ceftamin. The median final tumor volumes in mice treated with Infubinol and Ceftamin were larger than their mean tumor volumes across all timepoints, while the median final tumor volumes in mice treated with Capomulin and Ramicane were lower than their mean tumor volumes across all timepoints.
# 
# 2) For one mouse treated with Capomulin (Mouse ID s185), tumor volume was reduced by close to 50% over the course of its treatment.
# 
# 3) For mice treated with Capomulin, weight and average tumor volume were proven to be significantly correlated. These two measurements had a correlation factor of 0.84, implying a strong relationship.

# 

# In[1]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
total_data = pd.merge(mouse_metadata, study_results, how = "left", on=["Mouse ID", "Mouse ID"])

# Display the data table for preview
total_data


# In[2]:


# Checking the number of mice.
num_mice = len(total_data['Mouse ID'].unique())
num_mice


# In[3]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mice = total_data[total_data.duplicated(['Mouse ID', 'Timepoint'])]

duplicate_mice


# In[4]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
duplicate_indices = total_data[total_data['Mouse ID'] == 'g989'].index
total_data.drop(duplicate_indices, inplace = True)
total_data


# In[5]:


# Checking the number of mice in the clean DataFrame.
clean_mice_num = len(total_data['Mouse ID'].unique())
clean_mice_num


# ## Summary Statistics

# In[6]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 

# create df to hold results of calculations
summary_stats = pd.DataFrame(total_data.groupby("Drug Regimen").count())

summary_stats["Mean"] = pd.DataFrame(total_data.groupby("Drug Regimen")["Tumor Volume (mm3)"].mean())

summary_stats["Median"] = pd.DataFrame(total_data.groupby("Drug Regimen")["Tumor Volume (mm3)"].median())

summary_stats["Standard Deviation"] = pd.DataFrame(total_data.groupby("Drug Regimen")["Tumor Volume (mm3)"].std())

summary_stats["Variance"] = pd.DataFrame(total_data.groupby("Drug Regimen")["Tumor Volume (mm3)"].var())

summary_stats["SEM"] = pd.DataFrame(total_data.groupby("Drug Regimen")["Tumor Volume (mm3)"].sem())

# # Assemble the resulting series into a single summary dataframe.
summary_stats = summary_stats[["Mouse ID", "Mean", "Median", "Standard Deviation", "Variance", "SEM"]]

summary_stats = summary_stats.rename(columns = {"Mouse ID" : "Trials"})

summary_stats


# ## Bar and Pie Charts

# In[7]:


# Generate a bar plot showing the total number of measurements taken on each drug regimen using pandas.
measurements = summary_stats[['Trials']]
measurements

measurements.plot(kind = "bar", figsize = (7,4), color = "r", legend = False)

plt.title("Trials per Drug")
plt.ylabel("Number of Trials")

plt.ylim(0, max(measurements["Trials"])+20)

plt.show()


# In[8]:


# Generate a bar plot showing the total number of measurements taken on each drug regimen using pyplot.
x_axis = np.arange(len(measurements))
tick_locations = [value for value in x_axis]

plt.figure(figsize = (7,4))
plt.bar(x_axis, measurements['Trials'], color = "r", width = .6)
plt.xticks(tick_locations, measurements.index.values, rotation = "vertical")

plt.title("Trials per Drug")
plt.ylabel("Number of Trials")

plt.xlim(-0.75, len(x_axis)+0.25)
plt.ylim(0, max(measurements["Trials"])+20)

plt.show()


# In[9]:


# Generate a pie plot showing the distribution of female versus male mice using pandas
total_data

# group by sex
gender_data = pd.DataFrame(total_data.groupby(['Sex']).count()).reset_index()

gender_data = gender_data[["Sex", "Mouse ID"]]

# figure dimensions, aspects
plt.figure(figsize = (10,5))
ax1 = plt.subplot(121, aspect = "equal")

gender_data.plot(kind = "pie", ax = ax1, y = "Mouse ID", autopct='%1.1f%%', startangle = 150, shadow = False, labels = gender_data["Sex"], legend = False, fontsize = 15)

plt.title("Percentage of Trials by Sex")
plt.xlabel = ""
plt.ylabel = ""

plt.show()


# In[10]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
gender_numbers = (total_data.groupby(["Sex"])["Mouse ID"].count()).tolist()
gender_numbers


# In[11]:


# pie chart labels
labels = ["Females", "Males"]
colors = ["blue", "orange"]
explode = (0, 0)

# dimensions, etc
plt.pie(gender_numbers, explode = explode, labels = labels, colors = colors, autopct="%1.1f%%", shadow = False, startangle = 150)
plt.axis("equal")

plt.title("Percentage of Trials by Sex")

plt.show()


# ## Quartiles, Outliers and Boxplots

# In[12]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
chronological = total_data.sort_values(["Drug Regimen", "Mouse ID", "Timepoint"], ascending=True)
final_timepoint = chronological.loc[chronological["Timepoint"] == 45]
final_timepoint.reset_index()


# In[13]:


# separate data frames for each of the 4 drugs

capomulin = final_timepoint[final_timepoint["Drug Regimen"] == "Capomulin"].reset_index()

ramicane = final_timepoint[final_timepoint["Drug Regimen"] == "Ramicane"].reset_index()

infubinol = final_timepoint[final_timepoint["Drug Regimen"] == "Infubinol"].reset_index()

ceftamin = final_timepoint[final_timepoint["Drug Regimen"] == "Ceftamin"].reset_index()


# In[14]:


# turning tumor volume data into list for plotting
capo_vol = capomulin.sort_values(["Tumor Volume (mm3)"], ascending=True).reset_index()
capo_vol = capo_vol["Tumor Volume (mm3)"]

rami_vol = ramicane.sort_values(["Tumor Volume (mm3)"], ascending=True).reset_index()
rami_vol = rami_vol["Tumor Volume (mm3)"]

infu_vol = infubinol.sort_values(["Tumor Volume (mm3)"], ascending=True).reset_index()
infu_vol = infu_vol["Tumor Volume (mm3)"]
infu_vol

ceft_vol = ceftamin.sort_values(["Tumor Volume (mm3)"], ascending=True).reset_index()
ceft_vol = ceft_vol["Tumor Volume (mm3)"]


# Calculate the IQR and quantitatively determine if there are any potential outliers. 

# drugs = [capo_vol, rami_vol, infu_vol, ceft_vol]

# for drug in drugs:
#     quartiles = drug.quartile([.25, .5, .75])
#     lowerq = quartiles[0.25]
#     print(lowerq)


capo_quartiles = capo_vol.quantile([.25,.5,.75])
capo_lowerq = capo_quartiles[0.25]
capo_upperq = capo_quartiles[0.75]
capo_iqr = capo_upperq - capo_lowerq

capo_lower_bound = capo_lowerq - (1.5*capo_iqr)
capo_upper_bound = capo_upperq + (1.5*capo_iqr)

print(f"The lower quartile of Capomulin is {capo_lowerq} and the upper quartile is {capo_upperq}, with an IQR of {capo_iqr}.")
print(f"Values below {capo_lower_bound} or above {capo_upper_bound} could be outliers.")
print("------------------------------------------------------------------------------------------")

rami_quartiles = rami_vol.quantile([.25,.5,.75])
rami_lowerq = rami_quartiles[0.25]
rami_upperq = rami_quartiles[0.75]
rami_iqr = rami_upperq - rami_lowerq

rami_lower_bound = rami_lowerq - (1.5*rami_iqr)
rami_upper_bound = rami_upperq + (1.5*rami_iqr)

print(f"The lower quartile of Ramicane is {rami_lowerq} and the upper quartile is {rami_upperq}, with an IQR of {rami_iqr}.")
print(f"Values below {rami_lower_bound} or above {rami_upper_bound} could be outliers.")
print("------------------------------------------------------------------------------------------")

infu_quartiles = infu_vol.quantile([.25,.5,.75])
infu_lowerq = infu_quartiles[0.25]
infu_upperq = infu_quartiles[0.75]
infu_iqr = infu_upperq - infu_lowerq

infu_lower_bound = infu_lowerq - (1.5*infu_iqr)
infu_upper_bound = infu_upperq + (1.5*infu_iqr)

print(f"The lower quartile of Infubinol is {infu_lowerq} and the upper quartile is {infu_upperq}, with an IQR of {infu_iqr}.")
print(f"Values below {infu_lower_bound} or above {infu_upper_bound} could be outliers.")
print("------------------------------------------------------------------------------------------")

ceft_quartiles = ceft_vol.quantile([.25,.5,.75])
ceft_lowerq = ceft_quartiles[0.25]
ceft_upperq = ceft_quartiles[0.75]
ceft_iqr = ceft_upperq - ceft_lowerq

ceft_lower_bound = ceft_lowerq - (1.5*ceft_iqr)
ceft_upper_bound = ceft_upperq + (1.5*ceft_iqr)

print(f"The lower quartile of Ceftamin is {ceft_lowerq} and the upper quartile is {ceft_upperq}, with an IQR of {ceft_iqr}.")
print(f"Values below {ceft_lower_bound} or above {ceft_upper_bound} could be outliers.")

    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    
    
    # add subset 
    
    
    # Determine outliers using upper and lower bounds
    


# In[21]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest

fig, ax1 = plt.subplots()

four_drugs_df = pd.DataFrame({
    'Capomulin': capo_vol, 'Ramicane': rami_vol, 'Infubinol':infu_vol, 'Ceftamin':ceft_vol
})
ax1 = four_drugs_df.boxplot(column=['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin'], figsize=(10,50), grid=False)
ax1.set_title('Final Tumor Volumes')
ax1.set_xlabel('Drug Regimen')
ax1.set_ylabel('Final Tumor Volume (mm3)')

plt.show()


# ## Line and Scatter Plots

# In[24]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
# show all capomulin mice to see available choices
import matplotlib.pyplot as plt

all_capomulin_df = total_data.loc[total_data["Drug Regimen"] == "Capomulin"]
all_capomulin_df = all_capomulin_df.reset_index()
all_capomulin_df.head()

# generate line plot (tumor vol vs. time) for Mouse s185
# Mouse s185 data frame
chosen_mouse = all_capomulin_df.loc[all_capomulin_df['Mouse ID'] == 's185']

chosen_mouse = chosen_mouse.loc[:, ["Timepoint", "Tumor Volume (mm3)"]]

chosen_mouse = chosen_mouse.reset_index(drop=True)

# line plot
chosen_mouse.set_index("Timepoint").plot(figsize=(8,6), linewidth=2.5, color="red")

plt.show()


# In[33]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
all_capomulin_df

capomulin_weight_vol = all_capomulin_df.loc[:, ["Mouse ID", "Weight (g)", "Tumor Volume (mm3)"]]
capomulin_weight_vol

avg_capo_vol = pd.DataFrame(capomulin_weight_vol.groupby(["Mouse ID", "Weight (g)"])["Tumor Volume (mm3)"].mean()).reset_index()
avg_capo_vol

avg_capo_vol = avg_capo_vol.rename(columns = {
    'Tumor Volume (mm3)': 'Average Volume (mm3)'
})

avg_capo_vol

avg_capo_vol.plot(kind="scatter", x="Weight (g)", y="Average Volume (mm3)", grid=True, figsize=(5,5), title="Weight vs. Average Tumor Volume with Capomulin")


plt.show()


# ## Correlation and Regression

# In[60]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
from scipy.stats import linregress

x_values = avg_capo_vol["Weight (g)"]
y_values = avg_capo_vol["Average Volume (mm3)"]

(slope, intercept, rvalue, pvalue, stderr) = linregress(x_values, y_values)
regress_values = x_values * slope + intercept

line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))

plt.scatter(x_values, y_values)

plt.plot(x_values,regress_values,"r-")

plt.annotate(line_eq,(6,10),fontsize=15,color="red")

# plt.xlabel("Weight")
# plt.ylabel("Average Tumor Volume (mm3)")

plt.show()

print(line_eq)

# correlation coefficient (Mouse Weight vs Tumor Volume)
avg_capo_vol

weight = avg_capo_vol["Weight (g)"]
volume = avg_capo_vol["Average Volume (mm3)"]
correlation = st.pearsonr(weight,volume)
print(f"The correlation between mouse weight and average tumor volume is {round(correlation[0],2)}")


# In[ ]:




