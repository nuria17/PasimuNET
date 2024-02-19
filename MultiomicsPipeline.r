
library(fpc)
library(plotrix)
library(WGCNA)
library(dplyr)
library(tibble)
library(purrr)
library(pvclust)
library(aricode)
library(argparse)
library(stringr)
library(doParallel)
library(foreach)
library(parallel)
library(tidyverse)
library(ggnewscale)


parser=ArgumentParser()

parser$add_argument("--directory", default=".", 
    help = "Path to set the working directory")
parser$add_argument("--id", type= 'character',
    help= "Ids to keep in a string format separated by comma")
parser$add_argument("--pheno", type= 'character',
    help= "Phenotype file")
parser$add_argument("--remove_clust", type= 'character',
    help= "Remove the specified cluster from clinical analysis", 
    default= NULL)

args= parser$parse_args()

ClusteringParsing <- function(file, phenoall){
    phenoall$Clust=NULL
    tmp= read.table(file)
    clust=grep('ClusterID:',unlist(tmp))
    llista_clust_all=list()
    for (i in 1:(length(clust))){ if (i < (length(clust))) {llista_clust_all=append(llista_clust_all,list(tmp[(clust[i]+1):(clust[i+1]-1),]))} else {llista_clust_all=append(llista_clust_all,list(tmp[(clust[i]+1):dim(tmp)[[1]],]))}}
    for (e in 1:length(llista_clust_all)) {phenoall$Clust[phenoall$ID %in% llista_clust_all[[e]]]= e}
    #for (e in 1:length(llista_clust_all)) {phenoall$Clust[phenoall$Sample_ID %in% llista_clust_all[[e]]]= paste0('Cluster', e)}
    phenoall$Clust= as.factor(phenoall$Clust)
    return(phenoall)
}

mfun4=function(x){
    #library(jaccard)
    # Copy from library jaccard, the jaccard function. 
    custom_jac= function (x, y, center = FALSE, px = NULL, py = NULL) {
        if (length(x) != length(y)) {
            stop("Two fingerprints (x and y) must be of the same length.")
        }
        if (is.null(px) | is.null(py)) {
            px <- mean(x)
            py <- mean(y)
        }
        sumxy <- sum(x & y)
        unionxy <- sum(x) + sum(y) - sumxy
        if (unionxy == 0) {
            j <- (px * py)/(px + py - px * py)
        }
        else {
            j <- sumxy/unionxy
        }
        if (center == FALSE) {
            return(j)
        }
        else {
            return(j - (px * py)/(px + py - px * py))
        }
    }

    x <- t(as.matrix(x))
    vfun <- function(x, y) {1 - custom_jac(x, y, center=FALSE)}
    m=(outer(split(x, row(x)), split(x, row(x)),
                    Vectorize(vfun)))
    res= as.matrix(m)
    colnames(res)= rownames(x)
    rownames(res)= rownames(x)
    res= as.dist(res)
    attr(res, "method") <- "jaccard"
    return(res)   

}

phenoall= read.csv(args$pheno)

# Set working directory
setwd(args$directory)

print(dim(phenoall))
ids_to_keep= str_split(args$id, ',') %>% unlist()
print(length(ids_to_keep))
phenoall= phenoall[phenoall$ID %in% ids_to_keep,]
print(phenoall$ID)
reso_dir_mydata= 'Resolution/'
files_reso= list.files(reso_dir_mydata, pattern= '[^csv]$')
files_reso= files_reso[order(as.numeric(gsub('Comm.*p_', '', files_reso)))]
print(length(files_reso))
print(files_reso)

# Patient x Resolution matrix
for(i in 1:length(files_reso)){
    phenoall1= ClusteringParsing(file= paste0(reso_dir_mydata,files_reso[[i]]), phenoall= phenoall)
    clust= phenoall1[,colnames(phenoall1) %in% c('ID', 'Clust')]
    colnames(clust)[2]= gsub("Com.*p_",'',files_reso[[i]])
    if (i == 1){
        mat = clust
    } else {
        mat= mat %>% left_join(clust, by= 'ID')
    }

}

mat2= mat %>% column_to_rownames(var='ID')
mat3=apply(mat2, 2, function(x) as.numeric(x))
rownames(mat3)= rownames(mat2)
# mat3 %>% head()


# Create new binary matrix with Patient x Resolution underscore cluster in columns
From_communities_to_binary= function(mat3){
    llista_comb_reso=purrr::map(seq_along(colnames(mat3)),
    ~paste0(colnames(mat3)[[.x]], '_', names(table(mat3[,.x])))) %>% unlist()
    mat_reso=matrix(nrow=nrow(mat3),ncol=length(llista_comb_reso))
    colnames(mat_reso)= llista_comb_reso
    rownames(mat_reso)= rownames(mat3)

    for(i in 1:ncol(mat_reso)){
        reso= unlist(strsplit(colnames(mat_reso)[i], '_'))[[1]]
        com= unlist(strsplit(colnames(mat_reso)[i], '_'))[[2]]
        mat_reso[,i]= (mat3[,reso] == com)+0
    }
    return(mat_reso)
}

mat_reso= From_communities_to_binary(mat3)

#print(head(mat_reso))
print(diag(as.matrix(mfun4(t(mat_reso)))))
boots=10
#boots=10000
fit_jac_4reso_binary= pvclust(t(mat_reso), method.hclust="ward.D2", method.dist=mfun4, parallel=TRUE,nboot=boots,iseed=1234)
#fit_jac_4reso_binary= pvclust(t(mat_reso), method.hclust="single", method.dist=mfun4, parallel=TRUE,nboot=boots,iseed=1234)

elem=unlist(strsplit(colnames(mat3)[dim(mat3)[[2]]], '_'))
reso=elem[length(elem)]
png(paste0('Hierachical_clustering_cutoff_', reso, '_nboot_',boots,'.png'), height=800, width=1500)
plot(fit_jac_4reso_binary)
dev.off()
cat('Num com: ')
con= file("stdin")
num_com = as.numeric(readLines(con,1))

print(num_com)
close(con)

print(cutree(fit_jac_4reso_binary$hclust,k= num_com))
#q()
pheno_nuria_copd= phenoall
pheno_nuria_copd$Cluster_multiomicspipeline=paste0('Cluster',cutree(fit_jac_4reso_binary$hclust,k= num_com))

pheno_nuria_copd %>% write.csv('Clinical_data.csv')

# pheno_nuria_copd= pheno_nuria_copd %>% 
# dplyr::mutate(Cluster_multiomicspipeline= ifelse(Cluster_multiomicspipeline == 'Cluster5', 'Cluster3.2',
#                                    ifelse(Cluster_multiomicspipeline == 'Cluster3', 'Cluster3.1',Cluster_multiomicspipeline)))

if(!is.null(args$remove_clust)){
    print(args$remove_clust)
    clust_to_remove=as.vector(strsplit(args$remove_clust, ',')) %>% unlist()
    clust_to_remove2=as.vector(strsplit(clust_to_remove, ' '))
    print((clust_to_remove2))
    print(length(unlist(clust_to_remove2)))
    pheno_nuria_copd= pheno_nuria_copd[!pheno_nuria_copd$Cluster_multiomicspipeline %in% clust_to_remove2,]
    print(pheno_nuria_copd)
}
## Figures

# list categorical variables
categ= pheno_nuria_copd[,c(46,47,50,77:98,100:115, 125:130, 137)] %>% as.list()
categ$GOLD1= ifelse(categ$GOLD == '1', 'SI','NO')
categ$GOLD2= ifelse(categ$GOLD == '2', 'SI','NO')
categ$GOLD3= ifelse(categ$GOLD == '3', 'SI','NO')
categ$GOLD4= ifelse(categ$GOLD == '4', 'SI','NO')
categ$Gender_Female= ifelse(categ$Genero == 'MUJER', 'SI','NO')
categ= categ[!names(categ) %in% c('Enf..Intersticial.Difusa.Pulmonar', 'placa', 'Carcinoma.microcítico', 'GOLD', 'Genero')]

# list continuous variables
cont= pheno_nuria_copd[,c(43:45,48,52,53,54,61,116:123)] %>% as.list()
cont= cont[!names(cont) %in% c('Anos_Sin_Fumar')]
cont= map(cont, function(x) as.numeric(x))

ClinicalProfiling_dendo_no_overlap = function(phenoall,colname_cluster, llista_categ, llista_cont){
    phenoall$Clust= phenoall[[colname_cluster]] %>% as.character() %>% as.factor()
    llista_categ= map(llista_categ, function(x) x[!is.na(phenoall$Clust)])
    llista_cont= map(llista_cont, function(x) x[!is.na(phenoall$Clust)])
    phenoall= phenoall[!is.na(phenoall$Clust),]
    llista_tab= list()
    llista_mrna=list()
    llista_jac_mrna=list()
    llista_methy=list()
    llista_jac_methy=list()
    llista_micro=list()
    llista_jac_micro=list()
    llista_p=list()
    llista_mean=list()
    llista_table=list()
    llista_perc_categ=list()
    llista_w=list()
    d=0
    c=0
    f=0
    for (i in 1:length(levels(phenoall$Clust))){
        # Fisher
        llista_tab= base::append(llista_tab, list(purrr::map(llista_categ,
        function(x) if(length(names(table(x))) > 1) {fisher.test(table(phenoall$Clust == levels(phenoall$Clust)[[i]],x), alternative = 'two.sided')$p.value} else {1})))
        llista_table=append(llista_table, list(purrr::map(llista_categ, function(x) table(phenoall$Clust == levels(phenoall$Clust)[[i]],x))))
        llista_perc_categ=append(llista_perc_categ, list(purrr::map(llista_categ, 
        function(x) (length(phenoall$ID[phenoall$Clust == levels(phenoall$Clust)[[i]] & x == 'SI' & !is.na(x)])/length(phenoall$ID[phenoall$Clust == levels(phenoall$Clust)[[i]]]))*100)))
        names(llista_tab)[[i]]= levels(phenoall$Clust)[[i]]
        names(llista_table)[[i]]= levels(phenoall$Clust)[[i]]
        names(llista_perc_categ)[[i]]= levels(phenoall$Clust)[[i]]
        # Wilcoxon
        llista_w= append(llista_w, list(purrr::map(llista_cont,function(x) 
        {if((length(as.numeric(x)[phenoall$Clust == levels(phenoall$Clust)[[i]] & !is.na(x)])) > 0){
            wilcox.test(as.numeric(x) ~ phenoall$Clust == levels(phenoall$Clust)[[i]])$p.value} else {1}})))
        llista_mean=append(llista_mean, list(purrr::map(llista_cont, 
        function(x) mean(x[phenoall$Clust == levels(phenoall$Clust)[[i]] & !is.na(x) & !is.na(phenoall$Clust)]))))
        names(llista_w)[[i]]= levels(phenoall$Clust)[[i]]
        names(llista_mean)[[i]]= levels(phenoall$Clust)[[i]]
        # Overlap with the communities from individual layers
        
        # Percentages of the clinical groups
        p= table(phenoall$Clust == levels(phenoall$Clust)[[i]], phenoall$group5)
        #llista_p= append(llista_p, list(list((p[2]/(p[2] + p[4] + p[6]))*100,(p[4]/(p[2] + p[4] + p[6]))*100, (p[6]/(p[2] + p[4] + p[6]))*100)))
        llista_p= append(llista_p, list(list(p[2]/(p[2]+p[4]), p[4]/(p[4]+p[2]), p[2] + p[4])))

    }
    # Clinical variables
    test=purrr::map2(llista_tab, llista_w, function(x,y) c(x,y))
    # Overlap with mRNA communities
    #llista_comb_mrna= append(llista_comb_mrna, list(llista_mrna))
    #names(llista_comb_mrna)[[a]]= combinations[[a]]
    # Overlap with methylation communities
    #llista_comb_methy= append(llista_comb_methy, list(llista_methy))
    #names(llista_comb_methy)[[a]]= combinations[[a]]
    # Overlap with miRNA communities
    #llista_comb_micro= append(llista_comb_micro, list(llista_micro))
    #names(llista_comb_micro)[[a]]= combinations[[a]]
    # Percentages of the clinical groups
    tab=llista_p[[1]] %>% unlist()
    for (i in 2:length(llista_p)){tab= rbind(tab, llista_p[[i]] %>% unlist())}
    rownames(tab)= levels(phenoall$Clust)
    #colnames(tab)= c('GOLD12', 'GOLD34', 'NON-SMOKERS')
    colnames(tab)= c('GOLD12', 'GOLD34', 'SIZE')
    

    return (list(clinical=test,percet=tab, mean=llista_mean, table=llista_table, categ_perc= llista_perc_categ))

}

cli_profiling= ClinicalProfiling_dendo_no_overlap(pheno_nuria_copd, 'Cluster_multiomicspipeline', categ, cont)

FromListToMatrix= function(vars, llista, pheno, colname_clust, type= NULL){
    mat= matrix(ncol= length(vars), nrow= length(levels(pheno[[colname_clust]] %>% as.factor())))
    colnames(mat)= c(names(vars))
    rownames(mat)= levels(pheno[[colname_clust]] %>% as.factor()) 

    for(i in 1:nrow(mat)){
        for(j in 1:ncol(mat)){
            if (type == 'cont'){
                mat[i,j]= llista$mean[[i]][[j]]
            } else if (type == 'categ'){
                mat[i,j]= llista$categ_perc[[i]][[j]]
            } else if (type == 'clinical') {
                mat[i,j]= llista$clinical[[i]][[j]]
            }           
        }
    }
    return(mat)
}

mat_cont= FromListToMatrix(cont, cli_profiling, pheno_nuria_copd,'Cluster_multiomicspipeline', type= 'cont')
colnames(mat_cont)= c('FEV1 % ref.', 'FVC % ref.', 'FEV1/FVC', 'BMI', 'Pack-year', 'Age', 'DLCO % ref.',
'White blood cells count', 'Neutrophils (%)', 'Lymphocytes (%)', 'Monocytes (%)', 'Eosinophils (%)', 'Basophils (%)',
'Hematocrit', 'Platelet count')
mat_categ= FromListToMatrix(categ, cli_profiling, pheno_nuria_copd,'Cluster_multiomicspipeline', type= 'categ')
colnames(mat_categ)= c('Emphysema', 'Asthma', 'Bronchiectasis', 'OSAHS', 'Tuberculosis', 'Cardiovascular risk',
'Diabetes', 'Dyslipidemia', 'Arterial hypertension', 'Cardiopathy', 'Arrhythmia', 'Coronary heart disease',
'Pulmonar hypertension', 'Heart failure', 'Heart valve disease', 'Renal disease', 'Liver disease', 'Immune disease',
'Atopy','Chronic infections', 'Lung cancer', 'Other cancers (no lung cancer)', 'Antibiotics', 'Bronchodilator',
'Beta-adrenergic agonists', 'Anticholinergics', 'Corticosteroids', 'Inhaled corticosteroids', 'Systemic corticosteroids',
'Cardiovascular medication', 'Amiodarone','Antiplatelet medication', 'Antiarrhythmic medication', 'Anticoagulant medication', 
'Diuretics', 'Statins', 'Vasodilators', 'Other anti-inflammatory drugs', 'Adenocarcinoma', 'Large cell carcinoma',
'Squamous cell carcinoma', 'Unclassified carcinoma', 'Carcinoid', 'GOLD1', 'GOLD2', 'GOLD3', 'GOLD4', 'Gender Female')

scale_fun= function(x) {
    return((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)))
}
# categ_na= do.call(cbind,categ) %>% as.data.frame()
# #categ_na= na.omit(categ_na)
# categ_binary=(categ_na == 'SI')+0
# mat_all=as.data.frame(cbind(do.call(cbind, cont),categ_binary))
# mat_all=na.omit(mat_all)
# print(cor(mat_all))
# clust_clinical=hclust(as.dist(1-cor(mat_all)),'ward.D2')

# Heatmap of the continous variables
#output= 'Figures/heatmap_clinical_cont_intersection.png'
mat_df= mat_cont %>%
    as.data.frame() %>%
    rownames_to_column(var='ID') %>%
    gather(Clinical_variables, value, -ID) %>%
    mutate(Communities= ID)
#mat_df$Clinical_variables= factor(mat_df$Clinical_variables, levels= unique(mat_df$Clinical_variable))
scale= mat_df %>% dplyr::group_by(.groups=Clinical_variables) %>% 
    dplyr::mutate(scaled= scale_fun(value)) %>%
    ungroup() %>% 
    as.data.frame()

# Heatmap of the categorical variables
#output= 'Figures/heatmap_clinical_categ_intersection.png'
mat_df2= mat_categ %>%
        as.data.frame() %>%
        rownames_to_column(var='ID') %>%
        gather(Clinical_variables, value, -ID) %>%
        mutate(Communities= ID)

# mat_all=cbind(mat_categ,mat_cont)
# mat_all=apply(mat_all,1,function(x) scale_fun(x))
# print(mat_all)
# mat_all=t(mat_all)
# mat_all=na.omit(mat_all)
# print(1-cor(mat_all))
# print(hclust(as.dist(1-cor(mat_all)),'ward.D2'))
# clust_clinical= hclust(as.dist(1-cor(mat_all)),'ward.D2')

#mat_df$Clinical_variables= factor(mat_df$Clinical_variables, levels= unique(mat_df$Clinical_variable))
scale2= mat_df2 %>% dplyr::group_by(.groups=Clinical_variables) %>% 
    dplyr::mutate(scaled= scale_fun(value)) %>%
    ungroup() %>% 
    as.data.frame()

scale_all= rbind(scale, scale2)
print(scale_all)
scale_all$scaled[is.na(scale_all$scaled)]=0 
print(scale_all)
#print(scale_all)
#test= spread(scale_all,key= Communities, value= scaled)
library(data.table)
test=dcast(scale_all, Clinical_variables ~ Communities, value.var = "scaled")
test= t(test)
colnames(test)= test[1,]
test= test[-1,]
print(test)
test= test[,unique(scale_all$Clinical_variables)]
test= apply(test,2,function(x) as.numeric(x))
print(test)
cor= cor(test)
cor[is.na(cor)]=0
#test=t(test)
# test=na.omit(t(test))
# test= t(test)
clust_clinical=hclust(as.dist(1-cor), 'ward.D2')
png('ttt.png')
plot(clust_clinical)
dev.off()
print(clust_clinical)
print(clust_clinical$order)
# Sort the columns according to the most severe cluster --> 4reso

#llista_cli=scale_all %>% 
#group_by(ID) %>% 
#nest() %>% 
#pull('data') %>% 
#map(~pull(.x,'scaled'))

# 4 resolution
order_c4_c3= unique(scale_all$Clinical_variables)[clust_clinical$order]
print(clust_clinical$order)
## Adding p-values ##
# Continuous variables
pheno_nuria_copd_contrast= pheno_nuria_copd[,colnames(pheno_nuria_copd) %in% c(names(cont))]
colnames(pheno_nuria_copd_contrast)= c('FEV1 % ref.', 'FVC % ref.', 'FEV1/FVC', 'BMI', 'Pack-year', 'Age', 'DLCO % ref.',
'White blood cells count', 'Neutrophils (%)', 'Lymphocytes (%)', 'Monocytes (%)', 'Eosinophils (%)', 'Basophils (%)',
'Hematocrit', 'Platelet count')
pheno_nuria_copd_contrast= pheno_nuria_copd_contrast[,colnames(pheno_nuria_copd_contrast) %in% colnames(test)]
print(pheno_nuria_copd_contrast)

kruskal=pheno_nuria_copd_contrast %>% 
gather(Clinical_variables,value) %>% 
group_by(Clinical_variables) %>%
nest() %>%
spread(key= Clinical_variables, data) %>% 
map(~kruskal.test(unlist(.x[[1]]), 
pheno_nuria_copd$Cluster_multiomicspipeline)$p.value)

bb=scale %>% group_by(Clinical_variables) %>%
nest() %>%
spread(key= Clinical_variables, data) %>% 
map(~.x[[1]])

var= scale$Clinical_variables %>% table() %>% names()
bb2= purrr::map(names(bb), function(x){bb[[x]] %>%
    rbind(list('p_value', kruskal[[x]], 'p_value', x, 0))})
scale_pval_cont=do.call(rbind, bb2)

# Categorical variables 
pheno_nuria_copd_contrast= pheno_nuria_copd[,colnames(pheno_nuria_copd) %in% c(names(categ), 'Genero', 'GOLD')]
colnames(pheno_nuria_copd_contrast)= c('Emphysema','Gender_Female','GOLD', 'Asthma', 'Bronchiectasis', 'OSAHS', 'Tuberculosis', 'Cardiovascular risk',
'Diabetes', 'Dyslipidemia', 'Arterial hypertension', 'Cardiopathy', 'Arrhythmia', 'Coronary heart disease',
'Pulmonar hypertension', 'Heart failure', 'Heart valve disease', 'Renal disease', 'Liver disease', 'Immune disease',
'Atopy','Chronic infections', 'Lung cancer', 'Other cancers (no lung cancer)', 'Antibiotics', 'Bronchodilator',
'Beta-adrenergic agonists', 'Anticholinergics', 'Corticosteroids', 'Inhaled corticosteroids', 'Systemic corticosteroids',
'Cardiovascular medication', 'Amiodarone','Antiplatelet medication', 'Antiarrhythmic medication', 'Anticoagulant medication', 
'Diuretics', 'Statins', 'Vasodilators', 'Other anti-inflammatory drugs', 'Adenocarcinoma', 'Large cell carcinoma',
'Squamous cell carcinoma', 'Unclassified carcinoma', 'Carcinoid')
colnames(pheno_nuria_copd_contrast)[colnames(pheno_nuria_copd_contrast) == 'Gender_Female']= 'Gender Female'
pheno_nuria_copd_contrast= pheno_nuria_copd_contrast[,colnames(pheno_nuria_copd_contrast) %in% c(colnames(test), 'GOLD', 'Gender Female')]


fisher=pheno_nuria_copd_contrast %>% 
gather(Clinical_variables,value) %>% 
group_by(Clinical_variables) %>%
nest() %>%
spread(key= Clinical_variables, data) %>%
map(function(x) if(dim(table(x[[1]])) > 1){fisher.test(table(unlist(x[[1]]),
as.factor(pheno_nuria_copd$Cluster_multiomicspipeline)),workspace=2e8)$p.value} else {1})

print(fisher)
fisher= append(fisher, list(GOLD1= fisher$GOLD, GOLD2= fisher$GOLD, GOLD3= fisher$GOLD, GOLD4= fisher$GOLD))
fisher$GOLD=NULL
print(fisher)


bb3=scale2 %>% 
filter(Clinical_variables %in% names(fisher)) %>%
group_by(Clinical_variables) %>%
nest() %>%
spread(key= Clinical_variables, data) %>% 
map(~.x[[1]])

bb4= map(names(bb3), function(x){bb3[[x]] %>%
    rbind(list('p_value', fisher[[x]], 'p_value', x, 0))})

scale_pval_categ=do.call(rbind, bb4)

scale_all= rbind(scale_pval_cont, scale_pval_categ)

print(scale_all$.groups)


scale_all$.groups= factor(scale_all$.groups, 
levels=order_c4_c3)

print(scale_all$.groups)


scale_all_filt= scale_all %>% filter(ID != 'p_value')
scale_all_pval= scale_all %>% filter(ID == 'p_value')
scale_all_pval= scale_all_pval %>% mutate(sig= ifelse(value < 0.05, 'p-value < 0.05',
                                            ifelse(value > 0.05 & value < 0.1, '0.05 < p-value < 0.1',
                                             'p-value > 0.1')))
scale_all_pval$sig= scale_all_pval$sig %>% as.factor()
hcl(h = c(15,135,255), c = 100, l = 65)
palette(hcl(h = c(15,135,255), c = 100, l = 65))

#output= 'Clinical_heatmap.png'
p=ggplot(scale_all_filt, mapping=aes(x= .groups,y=Communities)) +
geom_tile(aes(fill = scaled), color='black') +     
geom_text(aes(label = round(value,1)), size= 5) +
ggplot2::scale_fill_gradient('bb',low = "#eafaf1", high = "#186a3b") +
theme(text= element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

q= p + new_scale("fill") + geom_tile(data= scale_all_pval, aes(x=.groups,y= Communities, fill= sig),
color='white', size = 2) +    
geom_text(data= scale_all_pval, aes(label = round(value,2)), size= 5) +
ggplot2::scale_fill_brewer(direction = -1, palette = 10)

getwd()
png('Clinical_heatmap.png', width= 2500, height=1000)
q
dev.off()






