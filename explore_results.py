import pandas as pd
import sys

input_file = sys.argv[1]

results = pd.read_csv(input_file, sep="\t")

results.head()