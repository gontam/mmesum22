# This is a test.

file = open('input.csv')
print(file)
data = []
for i in file:
    data.append(i)

cluster_num = data[0][3]

data[1]=data[1].split(';')
rows = data[1][0]
columns = data[1][1]

