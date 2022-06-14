import requests


files = [(4,'microbiolink_demo/results/human_proteins/pm_periodontitis.fasta'),
         (3,'microbiolink_demo/results/human_proteins/pm_healthy.fasta')]


for i in files:
    with open(
            'microbiolink_demo/results/human_proteins/epi_DC_oral_cavity_zscore_10p_avg_exp.csv') as scdata:
        scdata.readline()

        output_file2 = i[1].split('/pm_')[1]
        output_file2 = output_file2.split('.')[0] + '_expressed_genes.csv'
        output_file2 = 'microbiolink_demo/results/human_proteins/' + output_file2

        with open(output_file2, 'w') as output:
            output.write("Gene" + "," + "Expression" + "\n")

            expressed_genes = {}
            for genes in scdata:
                genes = genes.strip().split(",")

                #Checking whether the gene is expressed or not (NaN)
                try:
                    if float(genes[i[0]]) > 0:
                        expressed_genes[genes[0]] = genes[i[0]]
                        output.write(genes[0] + "," + str(genes[i[0]]) + "\n")

                except ValueError:
                    continue


            with open(i[1], 'w') as output_file1:
                with open(
                        "microbiolink_demo/resources/OmniPath/Omnipath_membrane_protein.tsv") as membraneprotein_list:

                    proteins = []

                    for line in membraneprotein_list:
                        line = line.strip().split(" ")
                        uniprot = line[0]
                        gs = line[1]

                        if uniprot not in proteins and gs in expressed_genes:

                            url = "https://www.uniprot.org/uniprot/?query=" + uniprot + "&format=fasta&limit=10&sort=score"
                            r = requests.get(url)
                            output_file1.write(r.text)
                            proteins.append(line)
