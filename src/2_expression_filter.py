# Author: Lejla Gul
# Gene expression filtration based on individual cell count and z-score
# z-score method: Hart et al, 2013 and a discussion (https://www.biostars.org/p/94680/)
# Aim: discarding gene expressions with low value
# Inputs
#   1. Metadata to describe barcodes and cell type-condition
#   2. Normalised counts from individual cells
#   3. Average expression table
# Output: CSV file which describes genes expressed more than 10% of cells
#         (in each cell type - condition cluster separately)
#         and with a significant average expression ( z-score > -3)

import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

# Read the metatable with information about cells
with open('metadata_file') as metadata:
    metadata.readline()
    barcode_type = {}
    for line in metadata:
        line = line.strip().split("\t")
        barcode_type[line[0]] = (line[1], line[5])


# Read average expression table
# If it is necessary, please format the header to have the same names like in metadata table (e.g. Goblet_Healthy)
with open("average_expression_file") as average_exp:
    avg_exp = {}
    header = average_exp.readline()
    header = header.strip().split("\t")

    for line in average_exp:
        line = line.strip().split("\t")
        count1 = line[1:]
        gene = line[0].replace("\"", "")
        for index, cell_condition in enumerate(header):
            if cell_condition not in avg_exp:
                avg_exp[cell_condition] = []
            avg_exp[cell_condition].append((gene, count1[index]))


# Going through the normalised count table: calculate the expression rate of genes among individual cells,
# if it is more than 10%, avg exp. value is kept, if it is equal/less, avg.exp = 0.0
with open("normalised_count_file") as counts:
    type_index = {}
    avg_exp_filtered = {}

    cell_barcode = counts.readline()
    cell_barcode = cell_barcode.strip().split(',')[1:]
    for index, barcode in enumerate(cell_barcode):
        if barcode in barcode_type:
            cell_barcode[index] = barcode_type[barcode]
    for index, type_ in enumerate(cell_barcode):
        if type_ not in type_index:
            type_index[type_] = []
        type_index[type_].append(index)


    for line in counts:
        type_count = {}
        line = line.strip().split(",")
        gene = line[0]
        count_info = line[1:]
        for t, i_list in type_index.items():
            if t not in type_count:
                type_count[t] = []
            for i in i_list:
                type_count[t].append(float(count_info[i]))

        for type_ in type_count:
            type_count[type_] = np.array(type_count[type_])
            nonzero_count = np.count_nonzero(type_count[type_] != 0)
            rate = nonzero_count / len(type_count[type_]) * 100

            type_ = str(type_[0]) + "-" + str(type_[1])

            if type_ not in avg_exp_filtered:
                avg_exp_filtered[type_] = {}

            for gene_exp in avg_exp[type_]:
                if gene in gene_exp and rate > 10.0:
                    avg_exp_filtered[type_][gene] =  gene_exp[1]
                elif gene in gene_exp and rate <= 10.0:
                    avg_exp_filtered[type_][gene] = '0.0'


# Open nested dictionary as a pandas df
df = pd.DataFrame.from_dict(avg_exp_filtered)

# Save columns (cell-condition) and rows (genes)
columns = df.columns
rows = df.index.values

# Function to translate normalised 10X counts to log2 TPM values
def tpm_log2_converter(num):
    try:
        num=float(num)
    except:
        return 'NaN'
    if num > 0:
        tpm = float(num) * 100
        result = np.log2(tpm)
        return result

# Empty df for the log2 TPM values
df_cells_log2 = pd.DataFrame(columns=columns, index=rows)

# Filling the table
for column in columns:
    df_cells_log2[column] = df[column].apply(tpm_log2_converter)

# Empty df for filtered results
df_cells_log2_filtered = pd.DataFrame(columns=columns, index=rows)


# Going through each cell-condition column and apply z-score filter
for i in columns:
    print(i)

    fpkm = df_cells_log2[i].tolist()

    fpkm = np.array(fpkm)
    print(fpkm)
    fpkm_filtered = fpkm[np.logical_not(np.isnan(fpkm))]

    # creating the Gauss-curve
    kernel = gaussian_kde(fpkm_filtered)

    # creating the X axis -> divide the list for 100 units, xi numpy lists contains that 100 values
    # expected value = most oftest value, middle of Gaus-curve
    xi = np.linspace(fpkm_filtered.min(), fpkm_filtered.max(), 100)

    # calculate y for each x point
    yi = kernel.evaluate(xi)

    # expected value calculation, which x is by the max y value? (np.argmax(yi) = position)
    mu = xi[np.argmax(yi)]

    # fpkm > mu  = list of boolean values; mean of values right from the expected value
    U = fpkm_filtered[fpkm_filtered > mu].mean()

    # calculation of standard deviation
    sigma = (U - mu) * np.sqrt(np.pi / 2)

    # new score: deviancy from the mean divided by sigma (standard deviation)
    # z-value: relative value - deviation from the mean in the st.dev in the data -> Gaus-curve:  0.001 - 3 deviations
    zFPKM = (fpkm - mu) / sigma

    score_list = [fpkm[list(zFPKM).index(x)] if x > -3 else 'NaN' for x in zFPKM]

    s = pd.Series(score_list, index=rows)
    print(s)
    df_cells_log2_filtered[i] = s


# Writing out results
df_cells_log2_filtered.to_csv("zscore_10p_avg_exp")



