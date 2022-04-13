
# Read in Data (Cristina Soriano):
# TODO: Read in Data with pandas (as Dataframe): see: https://datatofish.com/import-csv-file-python-using-pandas/
# TODO: Please delete the character (\n) from data set.
# TODO: As you can see in, e.g. data[2] you have two values in one string please do two columns e.g., https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
# TODO: Make the data look like the table in excel that I uploaded in our folder.
file = open('input.csv')
print(file)
data = []
for i in file:
    data.append(i)

cluster_num = data[0][3]

data[1]=data[1].split(';')
rows = data[1][0]
columns = data[1][1]

