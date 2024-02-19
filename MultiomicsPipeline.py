
import numpy as np
#from sklearn.cross_validation import StratifiedKFold
#from sklearn.metrics.pairwise import cosine_similarity
#from sklearn.metrics import jaccard_score
import pandas as pd
from argparse import ArgumentParser
import os
import scipy.stats
from DistanceClosure import *
from utils import *
from PairNodesTogether import *
#from Size_resolution_log import *
import subprocess
import sys

args = ArgumentParser()
args.add_argument('--pheno', help='Input csv file with the clinical data', dest='pheno',  default=None)
args.add_argument('--exp_mrna', help='Expression matrix for mrna data type', dest='exp_mrna',  default=None)
args.add_argument('--exp_micro', help='Expression matrix for micro data type', dest='exp_micro',  default=None)
args.add_argument('--exp_methy', help='Expression matrix for methy data type', dest='exp_methy',  default=None)
#args.add_argument('--clust_colname', help='Name of the column that corresponds to the communities', dest='clust_colname',  default=None)
args.add_argument('--path', help='Path where to create the directories ', dest='path',  default=None)
args.add_argument('--bootstrap', help= 'Performs the hierarchical clustering bootstrapping \
    in a resolution range to pick optimal resolution cutoff', dest= 'bootstrap', default= False)
args.add_argument('--remove_clust', help= 'Cluster name to be removed from clinical analysis', dest= 'remove_clust', default= None)
#args.add_argument('--backbone', help='Input txt file with the backbone edges', dest='backbone',  default=None)
#args.add_argument('--out', help='Output txt file with the network in the suitable ncol format for Molti software',dest='out', default=None)

op = args.parse_args()

def coef_var_selection(mat):
    '''
    PARAMETERS:
    -----------
    mat: expression matrix with probes in the rows and individuls in columns.
    
    '''
    coef_var=pd.Series(scipy.stats.variation(mat, axis=1))
    coef_var.index= mat.index.values
    #print('var',len(coef_var))
    #print('index var',coef_var.index)
    #coef_var= coef_var.sort_values(ascending=False)
    sel_genes=coef_var.index.values[coef_var > np.quantile(a=coef_var, q=0.75)]
    #print('max sel genes', max(coef_var[coef_var > np.quantile(a=coef_var, q=0.75)]))
    #print('sel_genes',sel_genes)
    #sel_genes=coef_var.index.values[0:100]
    mat= mat.loc[sel_genes,:]
    print(mat.head)
    return [mat,sel_genes]

def create_dir_and_save_file(dirName, file=None, mat=None):
    '''
    PARAMETERS:
    -----------
    dirName: name of the directory to be created.
    file: name of the file to be saved.
    mat: variable that will be saved in the file.

    '''

    if not os.path.exists(dirName):
        os.mkdir(dirName)
        if mat is not None:
            if isinstance(mat, str):
                fd= open(os.path.join(dirName,file), 'w')
                fd.write(mat)
                fd.close()      
            elif isinstance(mat, pd.DataFrame):
                mat.to_csv(os.path.join(dirName,file), sep=' ', header=True, index=True)
            print("Directory " , dirName ,  " Created ")
    else: 
        if mat is not None:
            if isinstance(mat,str):   
                fd= open(os.path.join(dirName,file), 'w')
                fd.write(mat)
                fd.close()
            elif isinstance(mat,pd.DataFrame):    
                mat.to_csv(os.path.join(dirName, file), sep=' ', header=True, index=True)
            print("Directory " , dirName ,  " already exists")

def ClusteringMatrix(Dict, ids_test, ids_train):
    '''
    PARAMETERS:
    -----------
    Dict: dictionary with id as key and list with the ids in the community as value.
    ids_test: ids of the test samples
    '''
    Dict_filt= {k: Dict[k] for k in ids_test}
    #print(ids_test)
    #print(ids_train)
    mat= pd.DataFrame(data= np.zeros((len(ids_test), len(ids_train))), index= ids_test, columns= ids_train)    
    #print(mat)
    for i,id1 in enumerate(mat.index.values):
        mat.iloc[i, mat.columns.isin(Dict_filt[id1])]= 1
        mat.iloc[i, ~mat.columns.isin(Dict_filt[id1])]= 0        
    return mat
def FilteringTestNodes(backbone, test_nodes):
    '''
    PARAMETERS:
    ----------
    backbone: file with the backbone edges.
    test_nodes: list with the test edges in the cross validation.

    '''
    fd= open(backbone, 'r')
    t=""
    for line in fd:
        Line= line.rstrip().split(',')
        if Line[0] in test_nodes and Line[1] in test_nodes:
            if t == "":
                t= ','.join(Line) + '\n'
            else:
                t= t + ','.join(Line) + '\n'
    fd.close()
    return t
    
# Pheno data
pheno_copd_inter=pd.read_csv(op.pheno)
llista=[]
llista_ids=[]

if op.exp_mrna is not None:
    # Expression data
    mrna= pd.read_csv(op.exp_mrna, sep= ' ')
    # Probes names as index
    mrna= mrna.set_index('Unnamed: 0')
    print('Dim',mrna.shape)
    print('intersect:', np.intersect1d(mrna.columns, pheno_copd_inter['ID'].values))
    mrna_cross= mrna.loc[:, np.intersect1d(mrna.columns, pheno_copd_inter['ID'].values)]
    print('Dim matrix mRNA:',mrna_cross.shape)
    ids_mrna= pheno_copd_inter.loc[(pheno_copd_inter['mRNA_teixit'] == 'SI'),'ID'].values
    llista_ids.append(ids_mrna)
    print('Length ids mRNA:', len(ids_mrna))
    print('Filtering of mRNA IDs done')
    mrna_cv= coef_var_selection(mrna_cross)[0]
    print('Filtering of mRNA probes data done')
    corr_mrna= mrna_cv.corr(method= 'pearson')
    print('Calculation of mRNA correlations done')
    # Distance closure
    create_dir_and_save_file(dirName= os.path.join(op.path,'Cor_data'), \
        file= 'Cor_mrna', mat= corr_mrna)
    backbone_mrna=write_backbone_network(os.path.join(op.path,'Cor_data/Cor_mrna'))
    create_dir_and_save_file(dirName=os.path.join(op.path,'Backbone_data'),\
        file= 'Backbone_mrna', mat= backbone_mrna)

    # To ncol format
    backbone_mrna_ncol= write_ncol_file(file= os.path.join(op.path,'Cor_data/Cor_mrna'), \
        backbone= os.path.join(op.path,'Backbone_data/Backbone_mrna'))
    create_dir_and_save_file(dirName=os.path.join(op.path,'Backbone_data'),\
        file= 'Backbone_mrna.ncol', mat= backbone_mrna_ncol)
    print('Transformation to ncol format done in mRNA')
    llista.append('Backbone_data/Backbone_mrna.ncol')

if op.exp_micro is not None:
    # Expression data
    micro= pd.read_csv(op.exp_micro, sep= ' ')
    # Probes names as index
    micro= micro.set_index('Unnamed: 0')
    # Filter individuals in the dataframe
    micro_cross= micro.loc[:,np.intersect1d(micro.columns, pheno_copd_inter['ID'].values)]
    print('Dim matrix miRNA:',micro_cross.shape)
    ids_micro= pheno_copd_inter.loc[(pheno_copd_inter['small_teixit'] == 'SI'),'ID'].values
    llista_ids.append(ids_micro)
    print('Length ids miRNA:', ids_micro)
    # Selection of the most variable probes
    micro_cv= coef_var_selection(micro_cross)[0]
    print('Filtering of miRNA probes data done')
    # Correlations
    corr_micro= micro_cv.corr(method= 'pearson')
    # Create dir and save files
    create_dir_and_save_file(dirName= os.path.join(op.path,'Cor_data'), \
        file= 'Cor_micro', mat= corr_micro)
    # Distance closure
    backbone_micro= write_backbone_network(os.path.join(op.path,'Cor_data/Cor_micro'))
    create_dir_and_save_file(dirName=os.path.join(op.path,'Backbone_data'),\
        file= 'Backbone_micro', mat= backbone_micro)
  
    print('Calculation of miRNA backbone networks done')
    backbone_micro_ncol= write_ncol_file(file= os.path.join(op.path,'Cor_data/Cor_micro'),
        backbone= os.path.join(op.path,'Backbone_data/Backbone_micro'))
    create_dir_and_save_file(dirName=os.path.join(op.path,'Backbone_data'),\
        file= 'Backbone_micro.ncol', mat= backbone_micro_ncol)
    print('Transformation to ncol format done in miRNA')
    llista.append('Backbone_data/Backbone_micro.ncol')

if op.exp_methy is not None:
    # Expression data
    methy= pd.read_csv(op.exp_methy, sep=' ')
    # Probes names as index
    methy= methy.set_index('Unnamed: 0')
    # Filter individuals in the dataframe
    methy_cross= methy.loc[:,np.intersect1d(methy.columns, pheno_copd_inter['ID'].values)]
    print('Dim matrix methylome:',methy_cross.shape)
    print('Filtering of IDs done')
    ids_methy= pheno_copd_inter.loc[(pheno_copd_inter['Metilacio_teixit']== 'SI'),'ID'].values
    llista_ids.append(ids_methy)
    print('Length ids methy:', ids_methy)
    # Selection of the most variable probes
    methy_cv=coef_var_selection(methy_cross)[0]
    print('Filtering of methylome probes data done')
    # Calculate correlation
    corr_methy= methy_cv.corr(method= 'pearson')
    #corr_methy.to_csv("/home/nuria/Desktop/ProjectsResources/multiomics/multiomicspipeline_test_tosendtoreview/bla_methy.csv")
    #corr_methy= abs(corr_methy)
    print('Calculation of correlations done')
    # Creater dir and save file
    create_dir_and_save_file(dirName= os.path.join(op.path,'Cor_data'), \
        file= 'Cor_methy', mat= corr_methy)
    # Distance closure
    backbone_methy= write_backbone_network(os.path.join(op.path,'Cor_data/Cor_methy'))
    create_dir_and_save_file(dirName=os.path.join(op.path,'Backbone_data'),\
        file= 'Backbone_methy', mat= backbone_methy)
   
    print('Calculation of methylation backbone networks done')
    # Transform network to ncol format
    backbone_methy_ncol= write_ncol_file(file= os.path.join(op.path,'Cor_data/Cor_methy'),
        backbone= os.path.join(op.path,'Backbone_data/Backbone_methy'))
    create_dir_and_save_file(dirName=os.path.join(op.path,'Backbone_data'),\
        file= 'Backbone_methy.ncol', mat= backbone_methy_ncol)
    print('Transformation to ncol format done in methylome')
    llista.append('Backbone_data/Backbone_methy.ncol')

# Run molti in all the resolutions
create_dir_and_save_file(dirName=os.path.join(op.path,'Resolution'))
if len(llista) == 3:
    #tuning_resolution(reso_min= 5, reso_max= 2000, fileA= os.path.join(op.path,llista[0]),\
        #fileB= os.path.join(op.path,llista[1]), fileC= os.path.join(op.path,llista[2]),\
            #output= os.path.join(op.path,'Resolution/Communities_with_only_COPD_p_'))
    tuning_resolution(reso_min= 5, reso_max= 210, fileA= os.path.join(op.path,llista[0]),\
        fileB= os.path.join(op.path,llista[1]), fileC= os.path.join(op.path,llista[2]),\
            output= os.path.join(op.path,'Resolution/Communities_with_only_COPD_p_'))

elif len(llista) == 2:
    tuning_resolution(reso_min= 5, reso_max= 210, fileA= os.path.join(op.path,llista[0]),\
        fileB= os.path.join(op.path,llista[1]), output= os.path.join(op.path,'Resolution/Communities_with_only_COPD_p_'))
elif len(llista) == 1:
    tuning_resolution(reso_min= 5, reso_max= 210, fileA= os.path.join(op.path,llista[0]), \
        output= os.path.join(op.path,'Resolution/Communities_with_only_COPD_p_'))

# Run rscript
flat_list = set([item for sublist in llista_ids for item in sublist])
print('Flat_list:',len(flat_list))
ids= ','.join(flat_list)
print(ids)
if op.bootstrap:
    subprocess.call(['Rscript', '--vanilla', 'Optimal_resolution_cutoff.r', '--directory',
    op.path, '--id', ids, '--pheno', op.pheno])
    print('Bootstrapping in all the resolution ranges done')


reso = input("Choose the modularity resolution limit: ")

if len(reso) == 1:
    print(reso)
    regex= "^Com.*[^\d](([0-{a}]\\.[0-9])|({b}\\.[0]))$".format(a=(int(reso)-1), b=int(reso))

elif len(reso) == 2:
    if int(reso[1])-1 > 0:
        regex= "^Com.*[^\d](([0-9]\\.[0-9])|([0-1][0-{a}]\\.[0-9])|({b}\\.[0]))$".format(a=int(reso[1])-1, b=int(reso))
    else:
        regex= "^Com.*[^\d](([0-9]\\.[0-9])|([0-1][0-9]\\.[0-9])|({b}\\.[0]))$".format(a=int(reso[1])-1, b=int(reso))


purge(dir= os.path.join(op.path, 'Resolution'), pattern=regex)

if op.remove_clust is None:
    subprocess.call(['Rscript', '--vanilla', 'MultiomicsPipeline.r', '--directory',
    op.path, '--id', ids, '--pheno', op.pheno])
else:
    subprocess.call(['Rscript', '--vanilla', 'MultiomicsPipeline.r', '--directory',
    op.path, '--id', ids, '--pheno', op.pheno, '--remove_clust', op.remove_clust])

