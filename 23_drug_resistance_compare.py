drug_resistance = {}
fp_db = open("GDSC_drug_resistance_compare_dbtable.txt","w")
with open("GDSC_drug_resistance_ls.txt","r") as f_x:
    for line_x in f_x:
        line_ls_x = line_x.strip().split("\t")
        kinase_x = line_ls_x[0]
        cell_id_x = line_ls_x[2]
        cell_name_x = line_ls_x[3]
        cancer_x = line_ls_x[4]
        drug_x = line_ls_x[5]
        ic50 = float(line_ls_x[6])
        uniprot_id = line_ls_x[7]
        Y_N_x = line_ls_x[8]

        tissue_x = line_ls_x[9]
        subtissue_x = line_ls_x[10]

        id = ":".join([kinase_x, drug_x, cancer_x, subtissue_x])

        if id in drug_resistance:
            if Y_N_x == "Yes":
                drug_resistance[id]["Yes"]["cell"].append(cell_id_x + ":" + cell_name_x)
                drug_resistance[id]["Yes"]["ic50"].append(ic50)

                fp_db.write("\t".join([id, "Yes", cell_id_x, cell_name_x, str(ic50)])+"\n")

            if Y_N_x == "No":
                drug_resistance[id]["No"]["cell"].append(cell_id_x + ":" + cell_name_x)
                drug_resistance[id]["No"]["ic50"].append(ic50)

                fp_db.write("\t".join([id, "No", cell_id_x, cell_name_x, str(ic50)])+"\n")
        else:
            drug_resistance[id] = {}
            drug_resistance[id]["Yes"] = {}
            drug_resistance[id]["No"] = {}
            drug_resistance[id]["Yes"]["cell"] = []
            drug_resistance[id]["Yes"]["ic50"] = []
            drug_resistance[id]["No"]["cell"] = []
            drug_resistance[id]["No"]["ic50"] = []
            if Y_N_x == "Yes":
                drug_resistance[id]["Yes"]["cell"].append(cell_id_x + ":" + cell_name_x)
                drug_resistance[id]["Yes"]["ic50"].append(ic50)
                fp_db.write("\t".join([id, "Yes", cell_id_x, cell_name_x, str(ic50)])+"\n")
            if Y_N_x == "No":
                drug_resistance[id]["No"]["cell"].append(cell_id_x + ":" + cell_name_x)
                drug_resistance[id]["No"]["ic50"].append(ic50)
                fp_db.write("\t".join([id, "No", cell_id_x, cell_name_x, str(ic50)])+"\n")
fp_db.close()

n = 0
fp = open("GDSC_drug_resistance_compare.txt","w")
fp_err = open("GDSC_drug_resistance_compare_only_hotspot.txt","w")
for id_i in drug_resistance:
    print(id_i)
    yes_cell = ",".join(drug_resistance[id_i]["Yes"]["cell"])
    yes_ic50 = ",".join([str(round(i, 6)) for i in drug_resistance[id_i]["Yes"]["ic50"]])
    yes_avg_ic50 = sum(drug_resistance[id_i]["Yes"]["ic50"]) / len(drug_resistance[id_i]["Yes"]["ic50"])
    yes_avg_ic50 = str(round(yes_avg_ic50, 6))

    if len(drug_resistance[id_i]["No"]["ic50"]) == 0:
        no_cell = "NULL"
        no_avg_ic50 = "NULL"
        no_ic50 = "NULL"
        n+=1
        fp_err.write(id_i + "\n")
    else:
        no_cell = ",".join(drug_resistance[id_i]["No"]["cell"])
        no_avg_ic50 = sum(drug_resistance[id_i]["No"]["ic50"]) / len(drug_resistance[id_i]["No"]["ic50"])
        no_avg_ic50 = str(round(no_avg_ic50, 6))
        no_ic50 = ",".join([str(round(i, 6)) for i in drug_resistance[id_i]["No"]["ic50"]])

    fp.write("\t".join([id_i,yes_cell,yes_ic50,no_cell,no_ic50,yes_avg_ic50,no_avg_ic50])+"\n")

fp.close()
print(n)













