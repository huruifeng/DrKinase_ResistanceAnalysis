# %% ================================================================
exists_mutated_hotspot = {}
with open("results/GDSC_drug_if_exists_mutated_hotspot.txt", "r") as f_x:
    for line_x in f_x:
        line_ls_x = line_x[:-1].strip().split("\t")
        kinase_x = line_ls_x[0]
        cell_id_x = line_ls_x[2]
        cell_name_x = line_ls_x[3]
        cancer_x = line_ls_x[4]
        drug_x = line_ls_x[5]

        tissue_x = line_ls_x[9]
        subtissue_x = line_ls_x[10]

        id = ":".join([kinase_x, cancer_x, drug_x, subtissue_x])
        if id in exists_mutated_hotspot:
            exists_mutated_hotspot[id].append(line_x)
        else:
            exists_mutated_hotspot[id] = [line_x]

# %% ================================================================
fp = open("results/GDSC_drug_resistance_candidates.txt","w")
with open("results/GDSC_drug_has_mutated_hotspot.txt","r") as f:
    next(f)
    for line in f:
        print(line)
        line_ls = line.strip().split("\t")
        kinase = line_ls[0]
        cell_id = line_ls[2]
        cell_name = line_ls[3]
        cancer = line_ls[4]
        drug = line_ls[5]
        tissue = line_ls[-4]
        subtissue = line_ls[-3]

        id = ":".join([kinase, cancer, drug, subtissue])
        if id in exists_mutated_hotspot:
            for line_x in exists_mutated_hotspot[id]:
                fp.write(line_x)
fp.close()

# %% ================================================================
fp = open("results/GDSC_drug_resistance_candidates_dedup.txt","w")
line_ls = {}
with open("results/GDSC_drug_resistance_candidates.txt","r") as f:
    for line in f:
        if line in line_ls:
            continue
        else:
            fp.write(line)
            line_ls[line] = 1
fp.close()














