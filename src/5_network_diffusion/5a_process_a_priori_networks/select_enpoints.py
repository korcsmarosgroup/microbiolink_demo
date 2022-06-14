expressed_genes_H = {}
expressed_genes_P = {}

with open("microbiolink_demo/results/TieDie/1_process_a_priori_networks/ec_ligands_epithelial_healthy.csv", "w") as output_file:
    with open(
            "microbiolink_demo/results/TieDie/1_process_a_priori_networks/ec_ligands_epithelial_periodontitis.csv",
            "w") as output_file2:
        with open("microbiolink_demo/resources/OmniPath/Omnipath_ligands.tsv") as protein_list:
            with open('microbiolink_demo/results/human_proteins/epi_DC_oral_cavity_zscore_10p_avg_exp.csv') as expression:
                output_file.write("Gene" + "," + "Expression" + "\n")
                output_file2.write("Gene" + "," + "Expression" + "\n")

                for line in expression:
                    line = line.strip().split(",")

                    try:
                        if float(line[3]) > 0:
                            expressed_genes_H[line[0]] = line[3]

                    except ValueError:
                        continue

                    try:
                        if float(line[4]) > 0:
                            expressed_genes_P[line[0]] = line[4]

                    except ValueError:
                        continue


                protein_list.readline()
                for line_2 in protein_list:
                    line_2 = line_2.strip().split("\t")
                    if line_2[2] in expressed_genes_H:
                        output_file.write(line_2[1] + "," + expressed_genes_H[line_2[2]] + "\n")

                    if line_2[2] in expressed_genes_P:
                        output_file2.write(line_2[1] + "," + expressed_genes_P[line_2[2]]+ "\n")
