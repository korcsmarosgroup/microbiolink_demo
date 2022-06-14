from pyfasta import Fasta
import re

def rename(fasta_key):
    fasta_key = fasta_key.split("|")
    fasta_key = fasta_key[1]
    return fasta_key

files = [('microbiolink_demo/results/human_proteins/pm_healthy.fasta',
          'microbiolink_demo/resources/microbiome/healthy_strains_domain.tsv',
          'microbiolink_demo/results/HMI/healthy_PPIs_all.csv'),
         ('microbiolink_demo/results/human_proteins/pm_periodontitis.fasta',
          'microbiolink_demo/resources/microbiome/periodontitis_strains_domain.tsv',
          'microbiolink_demo/results/HMI/periodontitis_PPIs_all.csv')]

for i in files:
    # fasta processing -> human.keys() print the keys, human[key_name] print the sequence
    human = Fasta(i[0])

    # elm identifier(key) - regex(value) dictionary
    with open ("microbiolink_demo/resources/ELM/elm_classes_2020.tsv") as motif_table:
        motif_table.readline()
        elm_regex = {}
        for line in motif_table:
            line = line.strip().split("\t")
            elm_regex[line[1]] = line[4]

    # motif(key) - domain(value list) dictionary
    with open ("microbiolink_demo/resources/ELM/elm_interaction_domains_2020.tsv") as motif_domain_table:
        motif_domain_table.readline()
        motif_domain = {}
        for line in motif_domain_table:
            line = line.strip("\n").split("\t")
            if line[0] not in motif_domain:
                motif_domain[line[0]] = []
            motif_domain[line[0]].append(line[1])

    # bacteria: pfam(key) - uniprot(value list) dictionary
    with open (i[1], "r") as domain_table:
        domain_table.readline()
        pfam_uniprot = {}
        for line in domain_table:
            line = line.strip().split("\t")
            if len(line) > 1:
                if ";" in line[1]:
                    line[1] = line[1].split(";")
                    for pfam in line[1]:
                        if pfam not in pfam_uniprot:
                            pfam_uniprot[pfam] = []
                        pfam_uniprot[pfam].append(line[0])

    # uniprot(key) - motif(value list) dictionary
    uniprot_motif = {}
    for key in human.keys():
        for motif in elm_regex:
            match = re.search(str(elm_regex[motif]), str(human[key]))
            if match:
                if rename(key) not in uniprot_motif:
                    uniprot_motif[rename(key)] = []
                uniprot_motif[rename(key)].append((motif,str(match.start()),str(match.end())))


    with open (i[2], "w") as output:
        output.write("Human protein" + "\t" + "Motif" + "\t" + "START" + "\t" + "END" + "\t" + "Domain" + "\t" + "Bacterial protein" + "\n")
        for pfam, uniprot_list in pfam_uniprot.items():
            for uniprot in uniprot_list:
                for motif in motif_domain:
                    if pfam in motif_domain[motif]:
                        for uni, motif_list in uniprot_motif.items():
                            for motif_2 in motif_list:
                                if motif_2[0] == motif:
                                    output.write(uni + "\t" + "\t".join(motif_2) + "\t" + "\t" + pfam + "\t" + uniprot + "\n")
