import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import gzip
import matplotlib.pyplot as plt
import pickle

input_dir = '/lustre/daystar/CellRanger/MT_MCE/outs/count/filtered_feature_bc_matrix/' #aggregate 된 파일
input_dir = '/lustre/daystar/Postech/Resources/CellRanger/WT/filtered_feature_bc_matrix/' #WT
input_dir = '/lustre/daystar/Postech/Resources/CellRanger/MCE/filtered_feature_bc_matrix/' #MCE
#counts
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
#barcodes
with gzip.open(input_dir + 'barcodes.tsv.gz', 'rt') as f:
    barcodes = [line.strip() for line in f.readlines()] 
#genes
with gzip.open(input_dir + 'features.tsv.gz', 'rt') as f:
    genes = np.array([line.strip().split('\t')[1] for line in f])

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.076)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

#threshold = 0.2
#predicted_doublets = scrub.call_doublets(threshold=threshold)

scrub.plot_histogram()
plt.show()

scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.show()


import pandas as pd
df = pd.DataFrame({
    'barcode': barcodes,
    'doublet_score': doublet_scores,
    'predicted_doublet': predicted_doublets,
    'sample': ['MCE'] * len(barcodes)
})

df['predicted_doublet'].value_counts()
print(df.head())

with open('/lustre/daystar/Postech/scrublet/scrublet_object_MCE.pkl', 'wb') as f:
    pickle.dump(scrub, f)

df.to_csv('/lustre/daystar/Postech/scrublet/scrublet_results_MCE.csv', index = False)


####################################
base_dir = "/lustre/sedolkang/coloncancer/scRNA_dataset/02_count/"
sample_dirs = sorted([d for d in os.listdir(base_dir) if not d.startswith('.')])

output_dir = "/lustre/sedolkang/coloncancer/Results/Minseo/Re-harmony/scrublet/"

sample_dirs = sorted([d for d in os.listdir(base_dir) if not d.startswith('.')])
all_results = []

for sample in sample_dirs:
    sample_path = os.path.join(base_dir, sample, "outs", "filtered_feature_bc_matrix")

    # Matrix load
    counts_matrix = scipy.io.mmread(os.path.join(sample_path, "matrix.mtx.gz")).T.tocsc()
    barcodes = pd.read_csv(os.path.join(sample_path, "barcodes.tsv.gz"), header=None)[0].values

    # Scrublet 실행
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # 결과 저장 (개별)
    result_df = pd.DataFrame({
        'barcode': barcodes,
        'doublet_score': doublet_scores,
        'predicted_doublet': predicted_doublets
    })
    result_df['sample'] = sample  # 샘플 이름 추가
    all_results.append(result_df)

    out_csv = os.path.join(output_dir, f"{sample}_scrublet_results.csv")
    result_df.to_csv(out_csv, index=False)

# --- 전체 결과 병합 저장 ---
merged_df = pd.concat(all_results, ignore_index=True)
merged_df.to_csv(os.path.join(output_dir, "scrublet_results_all_samples.csv"), index=False)

