import csv
import matplotlib.pyplot as plt
import pandas as pd
import parameters as param

"""
Little snippet of code to read from exported Ansys csv table
Really rudimentary at the moment will get better as I learn to bridge Ansys & Python

Inputs:

    - CSV table data containing as a function of #robtos/module:
        - frame thicknesses [mm]
        - max deformation [m]
        - max stresses [Pa]
        - average stress [Pa]

Output:

    - Figure: frame thickness vs mas def
    - Figure: frame thickness vs max stress
    - TBD: Figure: frame thickness vs average stress (max stress is influenced by stress concentration, good to have average stress to check coherence)

"""
# Fileneame = path to log file from Ansys analyses
# filename = r"C:\Users\rombach\Documents\Astrobots\Spec-s5_workshop\FEA\spec-s5_files\user_files/DesignPointLog.csv"
# filename = "./Results/turbo_table_data6.csv"
material = "stainless_steel"
filename = f"C:/Users/rombach/Documents/Astrobots/Spec-s5_workshop/FEA/Turbo_table_{material}.csv"

comment_prefix = "#"  # Set this to the comment prefix used in your file
separator = ";"

save_plots = False
num_lines_to_skip = 0 # log file logs since beginning of analyses, skip first rows that are now irrelevant

saving_df = {"save_plots": save_plots}
saving = param.SavingResults(saving_df)

data = []

with open(filename, "r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=separator)
    for line_number, row in enumerate(csv_reader):
        if not row or row[0].startswith(comment_prefix) or line_number < num_lines_to_skip:
            continue  # Skip empty rows and comment lines
        data.append(row)

# Now the 'data' list contains the non-comment CSV rows

df = pd.DataFrame(data)
print(df) 
df.drop([0],inplace=True)
df.drop([0,1,3], axis=1, inplace=True)
print(df) 
df.columns = ['Frame thickness', 'def102', 'maxstress102', 'avstress102', 'maxstress75', 'avstress75', 'def75', 'def63', 'maxstress63', 'avstress63']
df.sort_values(by=['Frame thickness'], ascending = True, inplace = True)
print(df) 

# Now you can extract columns from the data and create a plot
x = df['Frame thickness'].astype(float).to_numpy()

figtitle = "Max deformation of loaded focal plate VS plate thickness"
filename = "Max_def_VS_frame_thickness"
plt.figure(figsize=(12,8))
plt.plot(x, df['def63'].astype(float).to_numpy()*10**6, '.-', label='63 robots/module')
plt.plot(x, df['def75'].astype(float).to_numpy()*10**6, '.-', label='75 robots/module')
plt.plot(x, df['def102'].astype(float).to_numpy()*10**6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel(r"Max deformation [$\mu$m]")
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)

figtitle = "Max stress on loaded focal plate VS plate thickness"
filename = "Max_stress_VS_frame_thickness"
plt.figure(figsize=(12,8))
plt.plot(x, df['maxstress63'].astype(float).to_numpy()*10**-6, '.-', label='63 robots/module')
plt.plot(x, df['maxstress75'].astype(float).to_numpy()*10**-6, '.-', label='75 robots/module')
plt.plot(x, df['maxstress102'].astype(float).to_numpy()*10**-6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel("Max stress [MPa]")

# plt.tick_params(which='minor', length=5, width=2)
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)

figtitle = "Average stress on loaded focal plate VS plate thickness"
filename = "Av_stress_VS_frame_thickness"
plt.figure(figsize=(12,8))
plt.plot(x, df['avstress63'].astype(float).to_numpy()*10**-6, '.-', label='63 robots/module')
plt.plot(x, df['avstress75'].astype(float).to_numpy()*10**-6, '.-', label='75 robots/module')
plt.plot(x, df['avstress102'].astype(float).to_numpy()*10**-6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel("Average stress [MPa]")
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)
plt.show()