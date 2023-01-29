### scrublet
import scrublet as scr
import pandas as pd 
import numpy as np
import os
import matplotlib.pyplot as plt

txt_name = []
sample_name = []
### read data
for file in os.listdir("./count_matrix_filter/"):
    file_name = file.replace("_filter.txt","")
    txt_name.append(file_name)
    sample_name.append(file_name)

i=0
for file in os.listdir("./count_matrix_filter/"):
    txt_name[i] = pd.read_table(file,sep="\t",index_col=0,header=0).T
    i +=1
# restore predicted doublets
doublet_cell = []

### run scrublet one by one
count_matrix = txt_name[12]
scrub = scr.Scrublet(count_matrix, expected_doublet_rate=0.025)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
pd.value_counts(predicted_doublets)
scrub.plot_histogram()
plt.savefig("./tmp/his.png")
# adjust threshold 
predicted_doublets = scrub.call_doublets(threshold=0.22)
pd.value_counts(predicted_doublets)
scrub.plot_histogram()
plt.savefig("./tmp/his.png")

data_tmp = count_matrix[predicted_doublets]

# all doublets collect
doublet_cell.extend(data_tmp.index.tolist())

# save singlets
predicted_singlets = list(np.array(1)-predicted_doublets)
keep = [bool(i) for i in predicted_singlets]

data_keep = count_matrix[keep]
data_keep = data_keep.T
data_keep.to_csv("./count_matrix_scrublet/"+sample_name[12]+".csv")

###### save 
doublet_cell = pd.DataFrame(doublet_cell,columns="Cell")

doublet_cell.to_csv("./data_out/scrublet_doublets.csv")


