# %% ================================================================
kinase2seq = {}
uniprot2seq = {}
with open("data/kinase_basic_info.txt","r") as f:
    next(f)
    for line in f:
        line_ls = line.strip().split("\t")
        kinase = line_ls[0]
        uniprot_id = line_ls[1]
        seq = line_ls[-1]

        kinase2seq[kinase] = seq
        uniprot2seq[uniprot_id] = seq

# %% ================================================================
cellid2tissue = {}
cellid2subtissue = {}
with open("data/Cell_listTue Jun  9 16_15_07 2020.tsv","r") as f:
    next(f)
    for line in f:
        line_ls = line.strip().split("\t")
        cell_id = line_ls[2].strip('"')
        tissue = line_ls[5].strip('"')
        sub_tissue = line_ls[6].strip('"')

        cellid2tissue[cell_id] = tissue
        cellid2subtissue[cell_id] = sub_tissue

# %% ================================================================
kinase_dict = {}
with open("results/predicted_results_mutation_map.csv","r") as f:
    next(f)
    for line in f:
        line_ls = line.strip().split(",")
        hotspot = line_ls[6]
        dataset = line_ls[7]
        if dataset != "GDSC":
            continue
        uniprot_id = line_ls[0]

        cell = line_ls[9]

        pos = int(line_ls[-3])
        aa = line_ls[-2]

        seq = uniprot2seq[uniprot_id]
        if pos > len(seq):
            continue
        if aa != seq[pos-1]:
            continue

        if uniprot_id in kinase_dict:
            if cell in kinase_dict[uniprot_id]:
                kinase_dict[uniprot_id][cell].append([pos,hotspot])
            else:
                kinase_dict[uniprot_id][cell] = [[pos,hotspot]]
        else:
            kinase_dict[uniprot_id] = {}
            kinase_dict[uniprot_id][cell] = [[pos,hotspot]]

# %% ================================================================
fp = open("results/GDSC_drug_if_exists_mutated_hotspot.txt","w")
fp.write("kinase\tdataset\tcell_id\tcell_name\tcancer\tdrug\tic50\tuniprot_id\texist\ttissue\tsubtissue\thotspots\tpos\n")

fp_dr = open("results/GDSC_drug_has_mutated_hotspot.txt","w")
fp_dr.write("kinase\tdataset\tcell_id\tcell_name\tcancer\tdrug\tic50\tuniprot_id\texist\ttissue\tsubtissue\thotspots\tpos\n")
with open("data/GDSC_drug.txt","r") as f:
    next(f)
    for line in f:
        line_ls = line.strip().split("\t")
        kinase = line_ls[-1]
        cell_id = line_ls[2]
        # print(line)
        tissue = cellid2tissue[cell_id]
        sub_tissue = cellid2subtissue[cell_id]

        pos_ls = []
        hs_ls = []
        if kinase in kinase_dict:
            if cell_id in kinase_dict[kinase]:
                for x_i in kinase_dict[kinase][cell_id]:
                    pos_ls.append(str(x_i[0]))
                    hs_ls.append(str(x_i[1]))

                fp.write(line[:-1] + "\tYes\t" +tissue+"\t"+sub_tissue+"\t" + ",".join(hs_ls) + "\t" + ",".join(pos_ls)+ "\n")
                fp_dr.write(line[:-1] + "\tYes\t" +tissue+"\t"+sub_tissue+"\t" +",".join(hs_ls) + "\t" + ",".join(pos_ls)+ "\n")
            else:

                fp.write(line[:-1] + "\tNo\t"+tissue+"\t"+sub_tissue+"\t-\t-\n")
        else:
            fp.write(line[:-1] + "\tNo\t"+tissue+"\t"+sub_tissue+"\t-\t-\n")

fp.close()
fp_dr.close()












