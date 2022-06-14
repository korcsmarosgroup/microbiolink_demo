import pandas as pd

# importing expression information from healthy and diseased conditions in different cells as pandas dataframes (insert the path for your file)
df_scdata = pd.read_csv('microbiolink_demo/resources/scdata/gingiva_epithelium_DC.tsv', index_col = 0, header=0, sep='\t')
print(df_scdata)

ligands = []
with open ('microbiolink_demo/results/TieDie/ec_ligands_epithelial_healthy.csv') as L_H:
    L_H.readline()
    for line in L_H:
        line = line.strip().split(',')
        ligands.append(line[0])

interactions = []
Ls = []
Rs = []
LRIs = {}
with open ('microbiolink_demo/resources/OmniPath/Omnipath_LR_interaction.tsv') as interaction:
    interaction.readline()
    for line in interaction:
        line = line.strip().split('\t')
        if line[3] in ligands:
            interactions.append([line[8],line[7]])
            if line[8] not in LRIs:
                LRIs[line[8]] = []
            LRIs[line[8]].append(line[7])

cell_types = ['Epithelial_2_Healthy', 'DC_Healthy']

output_name = 'microbiolink_demo/results/intercellular_interactions' +  cell_types[0] + "_" + cell_types[1] +"_LRIs.csv"
source_cell = cell_types[0].split("_")[0]
target_cell = cell_types[1].split("_")[0]

df_target_cell = df_scdata[cell_types[1]]

for target in df_target_cell.index:

    # selecting those ones which play role in intercellular communication as a receiver
    if target in LRIs:
        Rs.append(target)
    else:
        continue


Rs = set(Rs)
Ls = set(LRIs.keys())

all_possible_interactions = []
for l in Ls:
    for r in Rs:
        all_possible_interactions.append([l, r])

ligand_dict = {}
for i in all_possible_interactions:

    if i in interactions:

        if i[0] not in ligand_dict:
            ligand_dict[i[0]] = []
        ligand_dict[i[0]].append([i[1], '1'])
    else:
        if i[0] not in ligand_dict:
            ligand_dict[i[0]] = []
        ligand_dict[i[0]].append([i[1], '0'])

receptor_names = [ln[0] for ln in ligand_dict[list(ligand_dict.keys())[0]]]

ligand_receptor_df = pd.DataFrame.from_dict(ligand_dict, orient='columns')

for col in ligand_receptor_df.columns:
    ligand_receptor_df[col] = ligand_receptor_df[col].apply(lambda x: ','.join(map(str, x))[-1])
ligand_receptor_df.index = receptor_names

ligand_receptor_df.to_csv(output_name)