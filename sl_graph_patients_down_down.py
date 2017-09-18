import random
import sys
import scipy.stats as ss
from tqdm import tqdm
from operator import add

tumor_name = sys.argv[1]

def parser_name(line):
    s = line.strip().split(",")
    return s[0]

S = 1000
ddr_genes = set(["UNION_ATM_BRCAx","PMS1", "PMS2","APC","ARID1A","RB1", "NF1","CTNNB1","DUSP4","DUSP6","MSH2", "MSH3","MSH6","MLH1", "PTEN", "CDH1", "TP53", "ATM", "BRCA1", "BRCA2"])
#ddr_genes = set(["TP53"])

list_name_CNV = [parser_name(x) for x in open("data/"+tumor_name+"_down.txt")] + ["UNION_ATM_BRCAx"]
#list_name_ME = [parser_name(x) for x in open("data/BLCA_ME_filtered")]

def parser_CNV(line):
    s = list([1  if x == 2 or x == -2 else 0 for x in map(int, line.strip().split("\t")[3:])])
    return s

def parser_nodup(line):
    s = list([1 if x == 1 else 0 for x in map(int, line.strip().split(",")[1:])])
    return s

matrix_cnv = [parser_nodup(x) for x in open("data/"+tumor_name+"_down.txt")]
#matrix_me = [parser_ME(x) for x in open("data/BLCA_ME_filtered")]

## ADD fake ATM + BRCA1 + BRCA2
_rows = [x[0] for x in zip(range(len(list_name_CNV)), list_name_CNV) if x[1] in set(["ATM","BRCA1", "BRCA2"])]
union_gene = [1 if x > 0 else 0 for x in list(map(add, matrix_cnv[_rows[0]], map(add, matrix_cnv[_rows[1]], matrix_cnv[_rows[2]]))) ]
print(_rows)
print(union_gene)
print(matrix_cnv[_rows[0]])
print(matrix_cnv[_rows[1]])
print(matrix_cnv[_rows[2]])
matrix_cnv.append(union_gene)

print("len(list_name_CNV)",len(list_name_CNV),
      "\nlen(matrix_cnv)", len(matrix_cnv))


#per il hypergeometric
genes_mutation_count = {}
for i in range(len(list_name_CNV)):
    genes_mutation_count[list_name_CNV[i]] = (sum(matrix_cnv[i]), matrix_cnv[i])
num_patients = len(matrix_cnv[0])


original_bpg = []
for col in range(len(matrix_cnv[0])):
    #genes = []
    genes = set()
    for i in range(len(list_name_CNV)):
        if matrix_cnv[i][col] == 1:
            #genes.append(list_name_CNV[i])
            genes.add(list_name_CNV[i])
    original_bpg.append((col, len(genes), genes))

original_bpg=list(filter(lambda x : x[1] > 0, original_bpg))


print("orignal_bpg, taaac\n")

cooccurences_data = {}
occurences = {}
for s in original_bpg:
    for i in list(set(list(s[2])).intersection(ddr_genes)):
        if i in occurences:
            occurences[i] += 1
        else:
            occurences[i] = 1

        for j in list(s[2]):
            if i != j:
                if (i,j) in cooccurences_data:
                    cooccurences_data[(i,j)] += 1
                else:
                    cooccurences_data[(i,j)] = 1

print("occurrences and cooccurrences, taaac\n")

for o in occurences:
    print(o, occurences[o])

sample_nums = range(len(original_bpg))

cooccurences_null = {}

actual_swaps = 0
for it in tqdm(range(S)):
    actual_swaps = 0
    for _ in range(100000):
        rs = random.sample(sample_nums, 2)
        g1name = random.sample(original_bpg[rs[0]][2],1)[0]
        g2name = random.sample(original_bpg[rs[1]][2],1)[0]
        if (g1name not in original_bpg[rs[1]][2]) and (g2name not in original_bpg[rs[0]][2]):
            original_bpg[rs[0]][2].remove(g1name)
            original_bpg[rs[0]][2].add(g2name)
            original_bpg[rs[1]][2].remove(g2name)
            original_bpg[rs[1]][2].add(g1name)
            actual_swaps += 1

	
    #print(" ".join(["1" if "PTEN" in original_bpg[n][2] else "0" for n in sample_nums]))
    #print("actual_swaps: %d" % actual_swaps )

    i_cooccurences_null = {}
    for s in original_bpg:
        for i in list(set(list(s[2])).intersection(ddr_genes)):
            for j in list(s[2]):
                if i != j:
                    if (i,j) in i_cooccurences_null:
                        i_cooccurences_null[(i,j)] += 1
                    else:
                        i_cooccurences_null[(i,j)] = 1

    for p in i_cooccurences_null:
        if p in cooccurences_null:
            cooccurences_null[p].append(i_cooccurences_null[p])
        else:
            cooccurences_null[p] = [ i_cooccurences_null[p] ]


results = []
for k in cooccurences_null.keys():
    ranv = cooccurences_null[k]
    reav = 0
    if k in cooccurences_data:
        reav = cooccurences_data[k]

    pval = 1.0 - float(len([x for x in ranv if x > reav]))/float(S)

    #calcolo hyper
    g1 = k[0]
    g2 = k[1]
    inters = sum([x[0] * x[1] for x in list(zip(genes_mutation_count[g1][1],genes_mutation_count[g2][1]))])
    hg = ss.hypergeom.cdf(inters, num_patients, genes_mutation_count[g1][0],genes_mutation_count[g2][0])
    #hamming
    hd = float(genes_mutation_count[g1][0] + genes_mutation_count[g2][0] - 3*inters)/float(num_patients)

    if pval < 2:
        results.append((k,reav,ranv,pval,hg,hd))

results = sorted(results, key = lambda x:x[2])

with open("results/results_"+tumor_name+"_down_down_DDR_union.txt", "w") as f:
    f.write("\t".join(["A", "B", "mutA", "mutB", "int", "ES", "hyperg", "hdmi"]) + "\n")
    for r in results:
        f.write("\t".join(map(str, [r[0][0],
                                    r[0][1],
                                    genes_mutation_count[r[0][0]][0],
                                    genes_mutation_count[r[0][1]][0],
                                    r[1],
                                    r[3],
                                    r[4],
                                    r[5]])) + "\n")
