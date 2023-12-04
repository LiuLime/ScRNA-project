datas = {'normal:t1:a': ['g1', 'g2', 'g3'],
         'normal:t1:b': ['g1', 'g2'],
         'dis:t1:c': ['g2', 'g4', 'g5', 'g6'],
         'dis:t1:d': ['g2', 'g5', 'g7'],
         'normal:t2:c': ['g2', 'g4', 'g5', 'g6'],
         'dis:t2:c': ['g2', 'g4', 'g5', 'g6'], }

genes = {}

for key in datas.keys():
    data_list = datas[key]
    key_parts = key.rsplit(':', 1)
    for gene in data_list:
        if gene not in genes.keys():
            genes[gene] = {}
        if key_parts[0] not in genes[gene].keys():
            genes[gene][key_parts[0]] = []
        genes[gene][key_parts[0]].append(key_parts[1])

print(genes)

# 1
for gene_key in genes:
    tissues = genes[gene_key]
    if len(tissues) == 1:
        for tissues_key in tissues.keys():
            if len(tissues[tissues_key]) == 1:
                print("Method 1", gene_key)
print("******************************")
# 2
for gene_key in genes:
    tissues = genes[gene_key]
    if len(tissues) == 1:
        for tissues_key in tissues.keys():
            if len(tissues[tissues_key]) >= 2:
                print("Method 2", gene_key)
print("******************************")
# 3
for gene_key in genes:
    tissues = genes[gene_key]
    if len(tissues) >= 2:
        sets = []
        for tissues_key in tissues.keys():
            for tissue in tissues[tissues_key]:
                sets.append(set(tissue))
        all_equal = True
        base_set = sets[0]
        for s in sets[1:]:
            if s != base_set:
                all_equal = False
                break
        if all_equal:
            print("Method 3", gene_key)
print("******************************")
# 4
for gene_key in genes:
    tissues = genes[gene_key]
    if len(tissues) >= 2:
        sets = []
        for tissues_key in tissues.keys():
            for tissue in tissues[tissues_key]:
                sets.append(set(tissue))
        all_equal = True
        index = 0
        for s in sets:
            base_set = s
            for inner_s in sets[index+1:]:
                if inner_s != base_set:
                    all_equal = False
                    break
            if not all_equal:
                print("Method 4", gene_key)
                break
            index += 1
