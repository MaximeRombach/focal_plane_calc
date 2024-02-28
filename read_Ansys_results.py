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
filename = r"C:\Users\rombach\Documents\Astrobots\Spec-s5_workshop\FEA\Turbo_table.csv"

comment_prefix = "#"  # Set this to the comment prefix used in your file
separator = ","

save_plots = True
num_lines_to_skip = 0 # log file logs since beginning of analyses, skip first rows that are now irrelevant

saving_df = {"save_plots": save_plots}
project_name = "Spec-s5"
saving = param.SavingResults(saving_df, project_name)

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
df.drop([0,1],inplace=True) # Drop 1st row with only P#

# df.drop([0, df.columns[len(df.columns)-1]], axis=1, inplace=True) # Drop 1st column with only DP# and last one with redundancy of frame thickness
df.drop([0], axis=1, inplace=True) # Drop 1st column with only DP# and last one with redundancy of frame thickness

print(df) 
df.columns = ['Frame thickness', 'def3', 'maxstress3', 'avstress3','def2', 'maxstress2', 'avstress2','def1_5', 'maxstress1_5', 'avstress1_5', 'defff1_5', 'maxstressff1_5', 'avstressff1_5']
df.sort_values(by=['Frame thickness'], ascending = True, inplace = True)


# Now you can extract columns from the data and create a plot
x = df['Frame thickness'].astype(float).to_numpy()

figtitle = "Max deformation of loaded focal plate VS plate thickness for intermediate frameless solution"
filename = "Max_def_VS_frame_thickness"
# labels = ['1.5mm walls', '2mm walls', '3mm walls', 'Full framed - 1.5mm walls']
labels = ['3mm raft gap', '3.5mm raft gap', '4.5 mm gap', 'Full framed - 3mm gap']
plt.figure(figsize=(12,8))
plt.plot(x, df['def1_5'].astype(float).to_numpy()*10**6, '.-', label=labels[0])
plt.plot(x, df['def2'].astype(float).to_numpy()*10**6, '.-', label=labels[1])
plt.plot(x, df['def3'].astype(float).to_numpy()*10**6, '.-', label=labels[2])
plt.plot(x, df['defff1_5'].astype(float).to_numpy()*10**6, '.-', label=labels[3])
plt.xlabel("Frame thickness [mm]")
plt.ylabel(r"Max deformation [$\mu$m]")
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)

figtitle = "Max stress on loaded focal plate VS plate thickness for intermediate frameless solution"
filename = "Max_stress_VS_frame_thickness"
plt.figure(figsize=(12,8))
plt.plot(x, df['maxstress1_5'].astype(float).to_numpy()*10**-6, '.-', label=labels[0])
plt.plot(x, df['maxstress2'].astype(float).to_numpy()*10**-6, '.-', label=labels[1])
plt.plot(x, df['maxstress3'].astype(float).to_numpy()*10**-6, '.-', label=labels[2])
plt.plot(x, df['maxstressff1_5'].astype(float).to_numpy()*10**-6, '.-', label=labels[3])
# plt.plot(x, df['maxstress102'].astype(float).to_numpy()*10**-6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel("Max stress [MPa]")

# plt.tick_params(which='minor', length=5, width=2)
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)

figtitle = "Average stress on loaded focal plate VS plate thickness for intermediate frameless solution"
filename = "Av_stress_VS_frame_thickness"
plt.figure(figsize=(12,8))
plt.plot(x, df['avstress1_5'].astype(float).to_numpy()*10**-6, '.-', label=labels[0])
plt.plot(x, df['avstress2'].astype(float).to_numpy()*10**-6, '.-', label=labels[1])
plt.plot(x, df['avstress3'].astype(float).to_numpy()*10**-6, '.-', label=labels[2])
plt.plot(x, df['avstressff1_5'].astype(float).to_numpy()*10**-6, '.-', label=labels[3])
# plt.plot(x, df['avstress75'].astype(float).to_numpy()*10**-6, '.-', label='75 robots/module')
# plt.plot(x, df['avstress102'].astype(float).to_numpy()*10**-6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel("Average stress [MPa]")
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)
plt.show()
