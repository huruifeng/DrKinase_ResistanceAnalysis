#05102020
#Ruifeng Hu

# %% ================================================================
kinase_mutations = {}
with open("data/all_mutation.txt","r") as f:
    next(f)
    for line in f:
        line_ls = line[:-1].split("\t")

        dataset = line_ls[0].strip()
        kinase = line_ls[1].strip()
        gene_id = line_ls[2].strip()
        hgnc_id = line_ls[3].strip()
        cancer = line_ls[4].strip()
        sample_id = line_ls[5].strip()
        AAchange = line_ls[6].strip()
        mut_type = line_ls[7].strip()
        mut_pos = line_ls[8].strip()
        mut_aa = line_ls[9].strip()
        uniprot_id = line_ls[10].strip()
        mut_x = {
            "dataset": dataset,
            "kinase": kinase,
            "gene_id": gene_id,
            "hgnc_id": hgnc_id,
            "cancer": cancer,
            "sample_id": sample_id,
            "AAchange": AAchange,
            "mut_type": mut_type,
            "mut_pos": int(mut_pos),
            "mut_aa": mut_aa,
        }
        if uniprot_id not in kinase_mutations:
            kinase_mutations[uniprot_id] = [mut_x]
        else:
            kinase_mutations[uniprot_id].append(mut_x)

# %% ================================================================
fp_out = open("results/predicted_results_mutation_map.csv","w")
fp_out.write("Entry,Hit,Seq,Start,End,Score,DR_type,mut_dataset,cancer,sample_id,AAchange,mut_pos,mutAA,mut_type\n")
with open("data/predicted_results.csv","r") as f:
    next(f)
    for line in f:
        print(line)
        line_ls = line[:-1].split(",")

        uniprot_id = line_ls[0].strip()
        DR_type = line_ls[-1].strip()

        if DR_type == "Gatekeeper":
            loci = [int(line_ls[3])]
        else:
            loci = list(range(int(line_ls[3]),int(line_ls[4])+1))

        if uniprot_id not in kinase_mutations:
            continue

        for mut in kinase_mutations[uniprot_id]:
            if mut["mut_pos"] in loci:
                print(mut["mut_pos"], end=", ")
                fp_out.write(",".join(line_ls)
                             + ","
                             + ",".join([mut["dataset"],mut["cancer"],mut["sample_id"],mut["AAchange"],str(mut["mut_pos"]),mut["mut_aa"],mut["mut_type"]])
                             + "\n")
fp_out.close()



