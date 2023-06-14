import os
import glob
import pickle
import pandas as pd
import numpy as np
import anndata as ad

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature, Regulon
from frozendict import frozendict

import seaborn as sns

'''
https://resources.aertslab.org/cistarget/
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
'''

DATA_FOLDER = "~/tmp"
RESOURCES_FOLDER = "~/resources"
DATABASE_FOLDER = "~/databases/"
SCHEDULER = "123.122.8.24:8786"
DATABASES_GLOB = os.path.join(
    '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/SCENIC/', 'hg38_*genes_vs_motifs.rankings.feather')
MOTIF_ANNOTATIONS_FNAME = '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
MM_TFS_FNAME = '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/SCENIC/allTFs_hg38.txt'
# only contains genes expressed in >= 5% of the total oligo population (n=8571)
SC_EXP_FNAME = "~/SingleCellProjects/MyProjects/SexDifferences/SCENIC/Astro_celdaModules/astro_matrix.csv"
# celdaGeneWeightsToModules.csv"
CELDA_GEN_WGT = "/home/brasel/SingleCellProjects/MyProjects/SexDifferences/SCENIC/Astro_celdaModules/celdaModuleScoresByCell.csv"
REGULONS_FNAME = "/home/brasel/SingleCellProjects/MyProjects/SexDifferences/SCENIC/Astro_celdaModules/regulons_astro_celda.p"
MOTIFS_FNAME = "/home/brasel/SingleCellProjects/MyProjects/SexDifferences/SCENIC/Astro_celdaModules/motifs_astro_celda.csv"


# read in single cell data
ex_matrix = pd.read_csv(SC_EXP_FNAME,  header=0, index_col=0).T
ex_matrix.shape

# read in TFs
tf_names = load_tf_names(MM_TFS_FNAME)
tf_names2 = load_tf_names(
    '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GENIE3/allHumanTFs.txt')
tmp = set(tf_names)
tmp.update(tf_names2)
tf_names = list(tmp)
del tmp

# merge single cell and Module score data
celdaModules = pd.read_csv(CELDA_GEN_WGT,  header=0, index_col=0).T
ex_matrix = ex_matrix.loc[:, ex_matrix.columns.isin(tf_names)]
merg = pd.merge(celdaModules, ex_matrix, left_index=True, right_index=True)

# load ranking databases
db_fnames = glob.glob(DATABASES_GLOB)


def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]


dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs

# infer co-expression modules
# adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
# Add in module scores instead of genes here
adjacencies = grnboost2(merg, tf_names=tf_names, verbose=True)
with open('/home/brasel/SingleCellProjects/MyProjects/SexDifferences/SCENIC/Astro_celdaModules/adjacencies.pickle', "wb") as f:
    pickle.dump(adjacencies, f)

adjacencies.loc[adjacencies["target"] == "L17"]
#################################################################
#           TF target  importance
# 937   L3MBTL4    L17  323.763358
# 998     TSHZ2    L17  198.112018
# 753     HMGB1    L17  183.729736
# 734      RFX4    L17  137.770484
# 820      RORA    L17  127.955549
# 756     SMAD9    L17  125.335955
# 693      SOX5    L17  104.715062
# 225      ZIC1    L17   96.688727
# 542      BNC2    L17   95.078858
# 213      MYLK    L17   80.568080
# 795       FOS    L17   72.744527
# 1156    HIF3A    L17   51.342189
# 767     MBNL2    L17   49.514494
# 540      NFIB    L17   43.204033
# 48       NFIA    L17   39.546327
#################################################################

# create modules from adjacencies
# modules = list(modules_from_adjacencies(adjacencies, ex_matrix))#this might be where I can throw in the values from celda
modules = list(modules_from_adjacencies(adjacencies, merg))
with open('/home/brasel/SingleCellProjects/MyProjects/SexDifferences/SCENIC/Astro_celdaModules/modules.pickle', "wb") as f:
    pickle.dump(modules, f)

l17module = []
for m in modules:
    if 'L17' in m.genes:
        l17module.append(m)

summary = set([(m.name, m.gene2weight['L17']) for m in l17module])
summary = sorted(list(summary), key=lambda x: x[1])
#################################################################
# hits
# ('Regulon for ZNF609', 28.027300957299477),
# ('Regulon for TCF7L1', 28.525633897867333),
# ('Regulon for MSI2', 31.850817888909678),
# ('Regulon for JUNB', 34.029328657623424),
# ('Regulon for CERS6', 38.83294243400417),
# ('Regulon for NFIA', 39.54632702791812),
# ('Regulon for HIF3A', 51.342188613208364),
# ('Regulon for FOS', 72.74452661330734),
# ('Regulon for MYLK', 80.56807999564526),
# ('Regulon for BNC2', 95.0788583862785),
# ('Regulon for ZIC1', 96.6887270559364),
# ('Regulon for SMAD9', 125.33595480946666),
# ('Regulon for RFX4', 137.7704844013234),
# ('Regulon for TSHZ2', 198.11201793182562),
# ('Regulon for L3MBTL4', 323.76335843218465)
#################################################################
