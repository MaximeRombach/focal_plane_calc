import csv
import matplotlib.pyplot as plt
import pandas as pd

filename = "./Results/turbo_table_data2.csv"
comment_prefix = "#"  # Set this to the comment prefix used in your file
separator = ","

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
print(type(x[0]))
plt.figure(figsize=(12,8))
plt.plot(df['Frame thickness'].astype(float).to_numpy(), df['def63'].astype(float).to_numpy()*10**6, '.-', label='63 robots/module')
plt.plot(df['Frame thickness'].astype(float).to_numpy(), df['def75'].astype(float).to_numpy()*10**6, '.-', label='75 robots/module')
plt.plot(df['Frame thickness'].astype(float).to_numpy(), df['def102'].astype(float).to_numpy()*10**6, '.-', label='102 robots/module')
plt.xlabel("Frame thickness [mm]")
plt.ylabel(r'Max deformation [$\mu$m]')
plt.title("Max deformation of loaded focal plate VS plate thickness")
plt.grid(True)
plt.legend()
plt.show()
