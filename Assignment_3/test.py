import pandas as pd

data = pd.read_csv("input.csv", error_bad_lines = False)

print(data)