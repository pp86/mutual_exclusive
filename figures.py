import pandas as pd

es_threshold = 0.0001

results_raw = pd.read_csv("data/results_all_tumors_DDR_down_others_UP_union_ahdmi_all_tumors.txt", sep="\t")

results = results_raw[results_raw['ES'] < es_threshold]
results = results[results['A'] != "UNION_ATM_BRCAx"]
results = results[results['mutB'] >= 200]
results = results[['A','B']]

conteggi = results.groupby("B").size().reset_index()
conteggi = conteggi[conteggi[0] > 1]

results2 = results.merge(conteggi, left_on='B', right_on='B').reset_index()

print(results2.head())

print(conteggi.head())

nodes = list(set(list(results2['A']) + list(results2['B'])))

print(len(nodes))

with open("data/graph.gdf", 'w') as f:
    f.write("nodedef>name VARCHAR,label VARCHAR,color VARCHAR\n")
    for g in list(set(list(results2['A']))):
        f.write(g+", "+g+","+"\'240,90,60\'"+"\n")
    for g in list(set(list(results2['B'])).difference(set(list(results2['A'])))):
        f.write(g + ", " + g + "," + "\'100,120,190\'" + "\n")
    f.write("edgedef>node1 VARCHAR,node2 VARCHAR\n")
    for e in zip(list(results2['A']), list(results2['B'])):
        f.write(e[0]+","+e[1]+"\n")


for g in set(list(results2['B'])):
    print(g)


