from tqdm import tqdm
import sys

patients_rna = []
patients_cna = []
ddr={'APC','ARID1A','RB1','NF1','BRCA1','BRCA2','PTEN','CDH1','ATM','TP53'}
graph = {}
tot = 0

tumor = sys.argv[1]

with open("raw_data/"+tumor+"/data_RNA_clean.txt") as f:
    patients_rna = next(f).strip().split()[2:]

    for line in tqdm(f):
        values = line.strip().split()
        gname = values[0]
        if gname in ddr:
            scores = [float(x) for x in values[2:]]
            for z in zip(patients_rna, scores):
                if z[1] < -2:
                    tot += 1
                    if gname in graph:
                        graph[gname].add(z[0])
                    else:
                        graph[gname] = set([z[0]])
        else:
            scores = [float(x) for x in values[2:]]
            for z in zip(patients_rna, scores):
                if z[1] > 2:
                    tot += 1
                    if gname in graph:
                        graph[gname].add(z[0])
                    else:
                        graph[gname] = set([z[0]])

print(len(graph.keys()))

with open("raw_data/"+tumor+"/data_CNA.txt") as f:
    patients_cna = next(f).strip().split()[2:]

    for line in tqdm(f):
        values = line.strip().split()
        gname = values[0]
        if gname in ddr:
            scores = [int(x) for x in values[2:]]
            for z in zip(patients_rna, scores):
                if z[1] == -2:
                    if gname in graph:
                        graph[gname].add(z[0])
                    else:
                        graph[gname] = set([z[0]])
        else:
            scores = [float(x) for x in values[2:]]
            for z in zip(patients_rna, scores):
                if z[1] == 2:
                    tot += 1
                    if gname in graph:
                        graph[gname].add(z[0])
                    else:
                        graph[gname] = set([z[0]])

print(len(graph.keys()))
print("minori 2: " + str(tot))


with open("data/brca_ddr_down_other_up.txt", "w") as f:
    for g in tqdm(graph.keys()):
        row = [g]
        for p in list(set(patients_cna).union(patients_rna)):
            if p in graph[g]:
                row.append("1")
            else:
                row.append("0")
        f.write(",".join(row) + "\n")