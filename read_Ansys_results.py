import csv
import matplotlib.pyplot as plt
import pandas as pd
import parameters as param


filename = "./Results/turbo_table_data3.csv"
comment_prefix = "#"  # Set this to the comment prefix used in your file
separator = ","

save_plots = False
save_frame_as_dxf = False
save_csv = False
save_txt = False

saving_df = {"save_plots": save_plots, "save_dxf": save_frame_as_dxf, "save_csv": save_csv, "save_txt": save_txt}
saving = param.SavingResults(saving_df)

data = []

with open(filename, "r") as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        if not row or row[0].startswith(comment_prefix):
            continue  # Skip empty rows and comment lines
        data.append(row)

# Now the 'data' list contains the non-comment CSV rows

df = pd.DataFrame(data)
df.drop([0,1],inplace=True)
df.drop([0,1,2], axis=1, inplace=True)
df.columns = ['Frame thickness', 'def63', 'def102', 'stress102', 'def75', 'stress75']
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
# plt.plot(x, df['def63'].astype(float).to_numpy()*10**-6, '.-', label='63 robots/module')
plt.plot(x, df['stress75'].astype(float).to_numpy()*10**-6, '.-', label='75 robots/module')
plt.plot(x, df['stress102'].astype(float).to_numpy()*10**-6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel("Max stress [MPa]")
plt.title(figtitle)
plt.grid(True)
plt.legend()
saving.save_figures_to_dir(filename)
plt.show()
