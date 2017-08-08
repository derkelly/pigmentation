setwd("/local1/derek/data/pigmentation/revision1/eqtl")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(qvalue)
library(ggbeeswarm)
library(pwr)
library(parallel)

# read in data and functions
source("paper_eqtl.fnc.R")

######################################################
######################## DATA ########################
######################################################

# master data matrix
data = get_melanocyte_data(getwd())
expression_norm = read.delim("expression/master_table.txt",sep="\t",header=T)

# egenes for mfsd12 region
mfsd12_egenes = read.delim("./flanking_genes/mfsd12_GeneList.txt",sep="\t",header=T,stringsAsFactors=F)
mfsd12_egenes$Gene_Symbol = fix_strings(mfsd12_egenes$Gene_Symbol)

# egenes for ddb1 region
ddb1_egenes = read.delim("./flanking_genes/ddb1_GeneList.txt",sep="\t",header=T,stringsAsFactors=F)
ddb1_egenes$Gene_Symbol = fix_strings(ddb1_egenes$Gene_Symbol)

# egenes for ddb1 region
herc2_egenes = read.delim("./flanking_genes/herc2_GeneList.txt",sep="\t",header=T,stringsAsFactors=F)
herc2_egenes$Gene_Symbol = fix_strings(herc2_egenes$Gene_Symbol)

# all egenes
egenes = rbind(mfsd12_egenes,ddb1_egenes,herc2_egenes)
posit_rs = read.delim("./genotypes/posit_rs_maps.txt",sep="\t",header=F,row.names=1,stringsAsFactors=F)
rs_posit = row.names(unique(posit_rs))
names(rs_posit) = unique(posit_rs)[,1]
                        
## variants determined important by caviar
## mfsd12_caviar = c("19:3565253","19:3544892","19:3545022","19:3566513","19:3547685")
## ddb1_caviar = c("11:61141259","11:61137147","11:61106525","11:61076372","11:61115821","11:61152630")

# variants to test for ddb1 and mfsd12
ddb1_testvars = c("rs1108769","rs10897150","rs12275843",
                  "rs11230664","rs7120594","rs12289370",
                  "rs7934735","rs2512809","rs2513329",
                  "rs2260655","rs2305465","rs148172827",
                  "rs2868351","rs7951574","rs7948623",
                  "rs11230678","rs10897155","rs57265008",
                  "rs4939519","rs1377457","rs1377458",
                  "rs12791961","rs4453253","rs4939520","rs3017597")

mfsd12_testvars = c("rs56203814","rs10424065","rs73527942","rs142317543","rs10414812","rs6510760","rs112332856","rs7254463","rs111317445")

herc2_testvars = c("rs1800404","rs1868333","rs6497271",
                   "rs6416602","rs11074323","rs4778248",
                   "rs7496305","rs6416603","rs6497277",
                   "rs10627923","rs4778252","rs6497282",
                   "rs12901047","rs12915877","rs4932618",
                   "rs4932620","rs1667395","rs1667393",
                   "rs1635167","rs2905952")

## get roadmap epigenomics data
epi_dir = "/local1/derek/data/epi_roadmap/all"

epi_rna.data = read.delim(file.path(epi_dir,"rna_seq/57epigenomes.RPKM.pc"), sep="\t", header=T, row.names=1)
epi_names = read.delim(file.path(epi_dir,"rna_seq/EG.name.txt"), sep="\t", header=F, row.names=1)


mfsd12_epirna = epi_rna.data[get_split(mfsd12_egenes$Ensemble_ID,"\\.",1),]
mfsd12_egenes.epi = mfsd12_egenes[!is.na(mfsd12_epirna[,1]),] # expressed in epigenomics roadmap?
mfsd12_epirna = mfsd12_epirna[!is.na(mfsd12_epirna[,1]),]

## get genes with melanocyte expression higher than the median across
## tissues from roadmap epigenomics
mfsd12_epirna.mel = mfsd12_epirna[(mfsd12_epirna$E059 > apply(mfsd12_epirna,1,median)) |
                                  (mfsd12_epirna$E061 > apply(mfsd12_epirna,1,median)),]

mfsd12_egenes.mel = mfsd12_egenes.epi[(mfsd12_epirna$E059 > apply(mfsd12_epirna,1,median)) |
                                      (mfsd12_epirna$E061 > apply(mfsd12_epirna,1,median)),]



ddb1_epirna = epi_rna.data[get_split(ddb1_egenes$Ensemble_ID,"\\.",1),]
ddb1_egenes.epi = ddb1_egenes[!is.na(ddb1_epirna[,1]),] # expressed in epigenomics roadmap?
ddb1_epirna = ddb1_epirna[!is.na(ddb1_epirna[,1]),]

## get genes with melanocyte expression higher than the median across
## tissues from roadmap epigenomics
ddb1_epirna.mel = ddb1_epirna[(ddb1_epirna$E059 > apply(ddb1_epirna,1,median)) |
                              (ddb1_epirna$E061 > apply(ddb1_epirna,1,median)),]

ddb1_egenes.mel = ddb1_egenes.epi[(ddb1_epirna$E059 > apply(ddb1_epirna,1,median)) |
                                  (ddb1_epirna$E061 > apply(ddb1_epirna,1,median)),]



herc2_epirna = epi_rna.data[get_split(herc2_egenes$Ensemble_ID,"\\.",1),]
herc2_egenes.epi = herc2_egenes[!is.na(herc2_epirna[,1]),] # expressed in epigenomics roadmap?
herc2_epirna = herc2_epirna[!is.na(herc2_epirna[,1]),]

## get genes with melanocyte expression higher than the median across
## tissues from roadmap epigenomics
herc2_epirna.mel = herc2_epirna[(herc2_epirna$E059 > apply(herc2_epirna,1,median)) |
                              (herc2_epirna$E061 > apply(herc2_epirna,1,median)),]

herc2_egenes.mel = herc2_egenes.epi[(herc2_epirna$E059 > apply(herc2_epirna,1,median)) |
                                  (herc2_epirna$E061 > apply(herc2_epirna,1,median)),]


egenes.mel = rbind(mfsd12_egenes.mel,ddb1_egenes.mel,herc2_egenes.mel)

#######################################################
################# REGRESSION ##########################
#######################################################

## covariates excluding percent ancestry and including percent ancestry
## covars = c("frozen.by","confluency","final.pass","days.in.culture","Growth.rate","Shape")
covars = c()

pigment.YRI.cor = cor.test(data$Pigmentation.scale,log2(data$YRI))


## mfsd12 genes and ancestry
mfsd12_ancestry = data.frame(gene = mfsd12_egenes[,2],
                             PC1_pval = sapply(mfsd12_egenes[,2],function(gene) cor.test(data$PC1,data[,gene],method="spearman")$p.value),
                             YRI_pval = sapply(mfsd12_egenes[,2],function(gene) cor.test(data$YRI,data[,gene],method="spearman")$p.value))
mfsd12_ancestry$PC1_qval = qvalue(mfsd12_ancestry$PC1_pval)$qvalues
mfsd12_ancestry$YRI_qval = qvalue(mfsd12_ancestry$YRI_pval)$qvalues
write.table(mfsd12_ancestry[order(mfsd12_ancestry$PC1_qval),],"results/mfsd12_ancestry.txt",sep="\t",quote=F,row.names=F,col.names=T)

## ddb1 genes and ancestry
ddb1_ancestry = data.frame(gene = ddb1_egenes[,2],
                             PC1_pval = sapply(ddb1_egenes[,2],function(gene) cor.test(data$PC1,data[,gene],method="spearman")$p.value),
                             YRI_pval = sapply(ddb1_egenes[,2],function(gene) cor.test(data$YRI,data[,gene],method="spearman")$p.value))
ddb1_ancestry$PC1_qval = qvalue(ddb1_ancestry$PC1_pval)$qvalues
ddb1_ancestry$YRI_qval = qvalue(ddb1_ancestry$YRI_pval)$qvalues
write.table(ddb1_ancestry[order(ddb1_ancestry$PC1_qval),],"results/ddb1_ancestry.txt",sep="\t",quote=F,row.names=F,col.names=T)

## herc2 genes and ancestry
herc2_ancestry = data.frame(gene = herc2_egenes[,2],
                             PC1_pval = sapply(herc2_egenes[,2],function(gene) cor.test(data$PC1,data[,gene],method="spearman")$p.value),
                             YRI_pval = sapply(herc2_egenes[,2],function(gene) cor.test(data$YRI,data[,gene],method="spearman")$p.value))
herc2_ancestry$PC1_bon = p.adjust(herc2_ancestry$PC1_pval,method="bonferroni")
herc2_ancestry$YRI_bon = p.adjust(herc2_ancestry$YRI_pval,method="bonferroni")
write.table(herc2_ancestry[order(herc2_ancestry$PC1_bon),],"results/herc2_ancestry.txt",sep="\t",quote=F,row.names=F,col.names=T)


## pigment.lm = lm(build_lm("Pigmentation.scale",covars), data)
## pigment.YRI.lm = lm(build_lm("Pigmentation.scale","log2(YRI)"), data)

## get all combinations of egenes and testvars
mfsd12_egenes_eqtl = data.frame(expand.grid(mfsd12_testvars,mfsd12_egenes[,2],stringsAsFactors=F))
ddb1_egenes_eqtl   = data.frame(expand.grid(ddb1_testvars,ddb1_egenes[,2],stringsAsFactors=F))
herc2_egenes_eqtl  = data.frame(expand.grid(herc2_testvars,herc2_egenes[,2],stringsAsFactors=F))


## get all beta values and p-values for mfsd12 eqtl analyses
temp.lm = mapply(function(x,y) summary(lm(data[,y] ~ data[,x] + data$PC1 + data$PC2)),mfsd12_egenes_eqtl[,1],mfsd12_egenes_eqtl[,2],SIMPLIFY=F)
mfsd12_egenes_eqtl$beta = lapply(temp.lm,function(x) x$coefficients[2,1])
mfsd12_egenes_eqtl$p.value = lapply(temp.lm,function(x) x$coefficients[2,4])


## get all beta values and p-values for ddb1 eqtl analyses
temp.lm = mapply(function(x,y) summary(lm(data[,y] ~ data[,x] + data$PC1 + data$PC2)),ddb1_egenes_eqtl[,1],ddb1_egenes_eqtl[,2],SIMPLIFY=F)
ddb1_egenes_eqtl$beta = lapply(temp.lm,function(x) x$coefficients[2,1])
ddb1_egenes_eqtl$p.value = lapply(temp.lm,function(x) x$coefficients[2,4])


## get all beta values and p-values for herc2 eqtl analyses
temp.lm = mapply(function(x,y) summary(lm(data[,y] ~ data[,x] + data$PC1 + data$PC2)),herc2_egenes_eqtl[,1],herc2_egenes_eqtl[,2],SIMPLIFY=F)
herc2_egenes_eqtl$beta = lapply(temp.lm,function(x) x$coefficients[2,1])
herc2_egenes_eqtl$p.value = lapply(temp.lm,function(x) x$coefficients[2,4])





## calculate eQTL with only the top variants (2-3) for each locus
mfsd12_topvars = c("rs10424065","rs6510760")
ddb1_topvars = c("rs7948623","rs11230664")
herc2_topvars = c("rs1800404","rs4932620","rs6497271")

mfsd12_topvar_eqtl = data.frame(expand.grid(mfsd12_topvars,mfsd12_egenes[,2],stringsAsFactors=F))
ddb1_topvar_eqtl = data.frame(expand.grid(ddb1_topvars,ddb1_egenes[,2],stringsAsFactors=F))
herc2_topvar_eqtl = data.frame(expand.grid(herc2_topvars,herc2_egenes[,2],stringsAsFactors=F))


## get all beta values and p-values for mfsd12 eqtl analyses
temp.lm = mapply(function(x,y) summary(lm(data[,y] ~ data[,x] + data$PC1 + data$PC2)),mfsd12_topvar_eqtl[,1],mfsd12_topvar_eqtl[,2],SIMPLIFY=F)
mfsd12_topvar_eqtl$beta = sapply(temp.lm,function(x) x$coefficients[2,1])
mfsd12_topvar_eqtl$p.value = sapply(temp.lm,function(x) x$coefficients[2,4])
mfsd12_topvar_eqtl$q.value = qvalue(mfsd12_topvar_eqtl$p.value)$qvalues
write.table(mfsd12_topvar_eqtl[order(mfsd12_topvar_eqtl$q.value),],"results/mfsd12_topvar_eqtl.txt",sep="\t",quote=F,row.names=F,col.names=T)


## get all beta values and p-values for ddb1 eqtl analyses
temp.lm = mapply(function(x,y) summary(lm(data[,y] ~ data[,x] + data$PC1 + data$PC2)),ddb1_topvar_eqtl[,1],ddb1_topvar_eqtl[,2],SIMPLIFY=F)
ddb1_topvar_eqtl$beta = sapply(temp.lm,function(x) x$coefficients[2,1])
ddb1_topvar_eqtl$p.value = sapply(temp.lm,function(x) x$coefficients[2,4])
ddb1_topvar_eqtl$q.value = qvalue(unlist(ddb1_topvar_eqtl$p.value))$qvalues
write.table(ddb1_topvar_eqtl[order(ddb1_topvar_eqtl$q.value),],"results/ddb1_topvar_eqtl.txt",sep="\t",quote=F,row.names=F,col.names=T)


## get all beta values and p-values for herc2 eqtl analyses
temp.lm = mapply(function(x,y) summary(lm(data[,y] ~ data[,x] + data$PC1 + data$PC2)),herc2_topvar_eqtl[,1],herc2_topvar_eqtl[,2],SIMPLIFY=F)
herc2_topvar_eqtl$beta = sapply(temp.lm,function(x) x$coefficients[2,1])
herc2_topvar_eqtl$p.value = sapply(temp.lm,function(x) x$coefficients[2,4])
herc2_topvar_eqtl$q.value = qvalue(unlist(herc2_topvar_eqtl$p.value))$qvalues
write.table(herc2_topvar_eqtl[order(herc2_topvar_eqtl$q.value),],"results/herc2_topvar_eqtl.txt",sep="\t",quote=F,row.names=F,col.names=T)



ddb1_rs2512809_vps37c.plot = ggplot(data,aes(x=jitter(rs2512809),y=vps37c)) + geom_point() + geom_smooth(method=lm) + theme_classic(base_size=20) + scale_x_continuous(breaks=c(0,1,2), label=c("C/C","C/A","A/A")) + xlab("rs2512809")
ddb1_rs2512809_rab3il1.plot = ggplot(data,aes(x=jitter(rs2512809),y=rab3il1)) + geom_point() + geom_smooth(method=lm) + theme_classic(base_size=20) + scale_x_continuous(breaks=c(0,1,2), label=c("C/C","C/A","A/A")) + xlab("rs2512809")
ddb1_rs2512809_ddb1.plot = ggplot(data,aes(x=jitter(rs2512809),y=ddb1)) + geom_point() + geom_smooth(method=lm) + theme_classic(base_size=20) + scale_x_continuous(breaks=c(0,1,2), label=c("C/C","C/A","A/A")) + xlab("rs2512809")



## ddb1_eqtl_signif = ddb1_eqtl_pval[ddb1_eqtl_pval[,3]*39 < 0.05,]
## mfsd12_eqtl_signif = mfsd12_eqtl_pval[p.adjust(mfsd12_eqtl_pval[,3],method="BH") < 0.05,]

## ddb1_YRI_eqtl_signif = ddb1_YRI_eqtl_pval[ddb1_YRI_eqtl_pval[,3]*39 < 0.05,]
## mfsd12_YRI_eqtl_signif = mfsd12_YRI_eqtl_pval[p.adjust(mfsd12_YRI_eqtl_pval[,3],method="BH") < 0.05,]





######################################
########## Permutation Null ##########
######################################


min_p = function(df,y_col,x_cols){

    y_perm = sample(df[,y_col])
    YRI = df$YRI
    return(min(sapply(x_cols, function(x_col) summary(lm(y_perm ~ df[,x_col] + df$PC1 + df$PC2))$coefficients[2,4])))
}



## parallelize
cl = makeCluster(detectCores()-2)
clusterExport(cl,c("min_p",
                   "data",
                   "mfsd12_egenes",
                   "mfsd12_testvars",
                   "ddb1_egenes",
                   "ddb1_testvars",
                   "herc2_egenes",
                   "herc2_testvars"))


## get distribution of empirical p-values for MFSD12 genes
mfsd12_emp_p = mat.or.vec(length(mfsd12_egenes[,2]),1)
mfsd12_perms = list()
names(mfsd12_emp_p) = mfsd12_egenes[,2]
for (gene in mfsd12_egenes[,2]){

    obs_minp = min(sapply(mfsd12_testvars, function(var) summary(lm(data[,gene] ~ data[,var] + data$PC1 + data$PC2))$coefficients[2,4])) # observed min_p
    perm_minp = parSapply(cl, 1:10000, function(i,df,y_col,x_cols) min_p(df,y_col,x_cols),
                          df=data[!is.na(data[,gene]),],
                          y_col=gene,
                          x_cols=mfsd12_testvars) # null min_p
    mfsd12_perms[[gene]] = perm_minp
    emp_p = mean(perm_minp < obs_minp,na.rm=T) # empirical p value: fraction of null min_p less than observed min_p
    print(c(gene,emp_p))
    if (emp_p == 0){
        mfsd12_emp_p[gene] = 1/sum(!(is.na(perm_minp)))
    } else {
        mfsd12_emp_p[gene] = mean(perm_minp < obs_minp,na.rm=T)
    }
}
mfsd12_emp_thresh = max(mfsd12_emp_p[qvalue(mfsd12_emp_p)$qvalues < 0.05])
mfsd12_snp_egenes.signif = lapply(mfsd12_egenes[,2],function(gene) mfsd12_egenes_eqtl %>%
                                                                   filter(mfsd12_egenes_eqtl[,2] == gene) %>%
                                                                   filter(p.value < quantile(mfsd12_perms[[gene]],mfsd12_emp_thresh)))
mfsd12_snp_egenes.signif = do.call("rbind", mfsd12_snp_egenes.signif[unlist(lapply(mfsd12_snp_egenes.signif,dim))[c(TRUE,FALSE)] != 0])
write.table(apply(mfsd12_snp_egenes.signif,2,as.character),"results/mfsd12_snp_egenes.txt",sep="\t",quote=F,row.names=F,col.names=T)


# write p-value and q-value
mfsd12_emp_pq = data.frame(gene=names(mfsd12_emp_p),pval=mfsd12_emp_p,qval=qvalue(mfsd12_emp_p)$qvalues)
write.table(mfsd12_emp_pq[order(mfsd12_emp_pq$qval),],"results/mfsd12_egene.txt",sep="\t",quote=F,row.names=F,col.names=T)


## get distribution of empirical p-values for MFSD12 genes
ddb1_emp_p = mat.or.vec(length(ddb1_egenes[,2]),1)
ddb1_perms = list()
names(ddb1_emp_p) = ddb1_egenes[,2]
for (gene in ddb1_egenes[,2]){

    obs_minp = min(sapply(ddb1_testvars, function(var) summary(lm(data[,gene] ~ data[,var] + data$PC1 + data$PC2))$coefficients[2,4])) # observed min_p
    perm_minp = parSapply(cl, 1:10000, function(i,df,y_col,x_cols) min_p(df,y_col,x_cols),
                          df=data[!is.na(data[,gene]),],
                          y_col=gene,
                          x_cols=ddb1_testvars) # null min_p
    ddb1_perms[[gene]] = perm_minp
    emp_p = mean(perm_minp < obs_minp,na.rm=T)
    print(c(gene,emp_p))
    if (emp_p == 0){
        ddb1_emp_p[gene] = 1/sum(!(is.na(perm_minp)))
    } else {
        ddb1_emp_p[gene] = mean(perm_minp < obs_minp,na.rm=T)
    }
}
ddb1_emp_thresh = max(ddb1_emp_p[p.adjust(ddb1_emp_p,method="BH") < 0.05])
ddb1_emp_thresh = max(ddb1_emp_p[qvalue(ddb1_emp_p)$qvalues < 0.05])
ddb1_snp_egenes.signif = lapply(ddb1_egenes[,2],function(gene) ddb1_egenes_eqtl %>%
                                                                   filter(ddb1_egenes_eqtl[,2] == gene) %>%
                                                                   filter(p.value < quantile(ddb1_perms[[gene]],ddb1_emp_thresh)))
ddb1_snp_egenes.signif = do.call("rbind", ddb1_snp_egenes.signif[unlist(lapply(ddb1_snp_egenes.signif,dim))[c(TRUE,FALSE)] != 0])
write.table(apply(ddb1_snp_egenes.signif,2,as.character),"results/ddb1_snp_egenes.txt",sep="\t",quote=F,row.names=F,col.names=T)


# write p-value and q-value
ddb1_emp_pq = data.frame(gene=names(ddb1_emp_p),pval=ddb1_emp_p,qval=qvalue(ddb1_emp_p)$qvalues)
write.table(ddb1_emp_pq[order(ddb1_emp_pq$qval),],"results/ddb1_egene.txt",sep="\t",quote=F,row.names=F,col.names=T)


## get distribution of empirical p-values for HERC2 genes
herc2_emp_p = mat.or.vec(length(herc2_egenes[,2]),1)
herc2_perms = list()
names(herc2_emp_p) = herc2_egenes[,2]
for (gene in herc2_egenes[,2]){

    obs_minp = min(sapply(herc2_testvars, function(var) summary(lm(data[,gene] ~ data[,var] + data$PC1 + data$PC2))$coefficients[2,4] # observed min_p
    perm_minp = parSapply(cl, 1:10000, function(i,df,y_col,x_cols) min_p(df,y_col,x_cols),
                          df=data[!is.na(data[,gene]),],
                          y_col=gene,
                          x_cols=herc2_testvars) # null min_p
    emp_p = mean(perm_minp < obs_minp,na.rm=T)
    herc2_perms[[gene]] = perm_minp
    print(c(gene,emp_p))
    if (emp_p == 0){
        herc2_emp_p[gene] = 1/sum(!(is.na(perm_minp)))
    } else {
        herc2_emp_p[gene] = mean(perm_minp < obs_minp,na.rm=T)
    }
}
herc2_emp_thresh = max(herc2_emp_p[p.adjust(herc2_emp_p,method="bonferroni") < 0.05])
herc2_emp_thresh = max(herc2_emp_p[qvalue(herc2_emp_p)$qvalues < 0.05])
herc2_snp_egenes.signif = lapply(herc2_egenes[,2],function(gene) herc2_egenes_eqtl %>%
                                                                   filter(herc2_egenes_eqtl[,2] == gene) %>%
                                                                   filter(p.value < quantile(herc2_perms[[gene]],herc2_emp_thresh)))
herc2_snp_egenes.signif = do.call("rbind", herc2_snp_egenes.signif[unlist(lapply(herc2_snp_egenes.signif,dim))[c(TRUE,FALSE)] != 0])
write.table(apply(herc2_snp_egenes.signif,2,as.character),"results/herc2_snp_egenes.txt",sep="\t",quote=F,row.names=F,col.names=T)


# write p-value and q-value
herc2_emp_pq = data.frame(gene=names(herc2_emp_p),pval=herc2_emp_p,qval=qvalue(herc2_emp_p)$qvalues)
write.table(herc2_emp_pq[order(herc2_emp_pq$qval),],"results/herc2_egene.txt",sep="\t",quote=F,row.names=F,col.names=T)






mfsd12_egenes = names(mfsd12_emp_p[p.adjust(mfsd12_emp_p,method="BH") < 0.05])
ddb1_egenes = names(ddb1_emp_p[p.adjust(ddb1_emp_p,method="BH") < 0.05])
herc2_egenes = names(herc2_emp_p[p.adjust(herc2_emp_p,method="bonferroni") < 0.05])

mfsd12_fdr05_gene = names(which.min(mfsd12_emp_p[p.adjust(mfsd12_emp_p,method="BH") > 0.05]))
mfsd12_fdr05_thresh = min(sapply(mfsd12_testvars, function(var) cor.test(data[,var],data[,mfsd12_fdr05_gene],method="spearman")$p.value))

ddb1_fdr05_gene = names(which.min(ddb1_emp_p[p.adjust(ddb1_emp_p,method="BH") > 0.05]))
ddb1_fdr05_thresh = min(sapply(ddb1_testvars, function(var) cor.test(data[,var],data[,ddb1_fdr05_gene],method="spearman")$p.value))

herc2_fwer05_gene = names(which.min(herc2_emp_p[p.adjust(herc2_emp_p,method="bonferroni") > 0.05]))
herc2_fwer05_thresh = min(sapply(herc2_testvars, function(var) cor.test(data[,var],data[,herc2_fdr05_gene],method="spearman")$p.value))


mfsd12_snp_egenes = mfsd12_eqtl_mel_pval[(mfsd12_eqtl_mel_pval[,1] %in% mfsd12_egenes) & (mfsd12_eqtl_mel_pval[,3] < mfsd12_fdr05_thresh),]
ddb1_snp_egenes = ddb1_eqtl_mel_pval[(ddb1_eqtl_mel_pval[,1] %in% ddb1_egenes) & (ddb1_eqtl_mel_pval[,3] < ddb1_fdr05_thresh),]
herc2_snp_egenes = herc2_eqtl_mel_pval[(herc2_eqtl_mel_pval[,1] %in% herc2_egenes) & (herc2_eqtl_mel_pval[,3] < herc2_fwer_thresh),]


















#############################
############ ASE ############
#############################

## power.table = mat.or.vec(56,56)

## reps = 1000
## ae = 0.2

## set.seed(0)
## for (i in 1:56){
##     for (j in 1:56){

##         power.table[i,j] = mean(replicate(reps, wilcox.test(abs(0.5-rbinom(i, 100, 0.5)/100),
##                                                             abs(0.5-rbinom(j, 100, 0.5-ae)/100))$p.value) < 0.05)
##     }
## }




ase_data = read.delim("ase/ASE_sample_data.txt",sep="\t",header=T,stringsAsFactors=F)
ase_data$SYMBOL = fix_strings(ase_data$SYMBOL)

## mfsd12_caviar = list(rs10424065 = c("rs56203814","rs10424065","rs73527942","rs142317543","rs10414812"),
##                      rs6510760 = c("rs111317445","rs6510760","rs112332856"))

## ddb1_caviar = list(rs7948623 = c("rs4453253","rs7948623"),
##                    rs57265008 = c("rs12289370","rs7934735","rs57265008"))

## test_vars.all = c(unlist(mfsd12_caviar),unlist(ddb1_caviar))
test_vars.all = c(unique(as.character(mfsd12_snp_egenes[,2])),unique(as.character(ddb1_snp_egenes[,2])))

ase_coding = grepl("missense",ase_data$Consequence) | grepl("synonymous",ase_data$Consequence) | grepl("splice",ase_data$Consequence)
ase_utr    = grepl("3_prime",ase_data$Consequence) & !(ase_data$SYMBOL %in% ase_data[ase_coding,]$SYMBOL)

ase_data.cu = ase_data[ase_coding | ase_utr,]

## get the variants associated with each gene
gene_variants = list()
for (gene in unique(ase_data.cu$SYMBOL)){

    gene_variants[[gene]] = unique(ase_data.cu[ase_data.cu$SYMBOL==gene,]$variantID)
}


                             

## get a wide format table with rows labeled by Sample and columns labeled by variant.
ase_table.het.wide = empty_df(list(),cols=unique(ase_data.cu$variantID),rows=unique(ase_data.cu$Sample))
ase_table.AE.wide = empty_df(list(),cols=unique(ase_data.cu$variantID),rows=unique(ase_data.cu$Sample))



for (i in 1:dim(ase_data.cu)[1]){
    ase_table.het.wide[ase_data.cu[i,]$Sample,ase_data.cu[i,]$variantID] = 1
    ase_table.AE.wide[ase_data.cu[i,]$Sample,ase_data.cu[i,]$variantID] = ase_data.cu[i,]$AE
}

tag_hets = as.matrix(ase_table.het.wide)
tag_hets[is.na(tag_hets)] = 0



## individuals heterozygous for test variants, making sure the
## individuals match across test and tag tables
test_hets = (data[data$Sample %in% row.names(ase_table.het.wide),test_vars.all]==1)*1
test_homs = !(data[data$Sample %in% row.names(ase_table.het.wide),test_vars.all]==1)*1

row.names(test_hets) = data$Sample[data$Sample %in% row.names(ase_table.het.wide)]
row.names(test_homs) = data$Sample[data$Sample %in% row.names(ase_table.het.wide)]

test_hets = as.matrix(test_hets[row.names(ase_table.het.wide),])
test_homs = as.matrix(test_homs[row.names(ase_table.het.wide),])

test_hets[is.na(test_hets)] = 0
test_homs[is.na(test_homs)] = 0



## number of individuals heterozygous for the test variant and tag variant
test_by_tag.hets = t(test_hets) %*% tag_hets
test_by_tag.homs = t(test_homs) %*% tag_hets

test_by_tag.all.long = merge(melt(t(test_hets) %*% tag_hets),melt(t(test_homs) %*% tag_hets),by=c(1,2),stringsAsFactors=F)
colnames(test_by_tag.all.long) = c("test_var","tag_var","hets","homs")



## carry out ASE for MFSD12
mfsd12_ase_results.expand = expand.grid(unique(as.character(mfsd12_snp_egenes[,2])),names(gene_variants)[names(gene_variants) %in% mfsd12_egenes],stringsAsFactors=F)
mfsd12_ase_results = data.frame(test_var = mfsd12_ase_results.expand[,1],
                                gene = mfsd12_ase_results.expand[,2],
                                tag_var = 0,
                                n_test_het = 0,
                                n_test_hom = 0,
                                wilcox_p = 0,
                                stringsAsFactors=F)



## loop over the CAVIAR-gene pairs
for (i in 1:dim(mfsd12_ase_results)[1]){

    ## get the caviar variant and the gene
    test_var = mfsd12_ase_results[i,1]
    gene   = mfsd12_ase_results[i,2]

    ## get the list of test and tag vars
    ## test_vars = mfsd12_caviar[[caviar]]
    tag_vars = gene_variants[[gene]]

    if (length(tag_vars) > 1){
        ## get a long format list of heterozygotes and homozygotes
        temp.test_by_tag.all.long = cbind(test_var,merge(melt(test_by_tag.hets[test_var,tag_vars],stringsAsFactors=F),
                                                         melt(test_by_tag.homs[test_var,tag_vars],stringsAsFactors=F),by="row.names",stringsAsFactors=F))

        ## get pairs with at least 2 individuals for hets and homs
        temp.test_by_tag.all.long.valid = temp.test_by_tag.all.long[temp.test_by_tag.all.long[,3]>1 &
                                                                    temp.test_by_tag.all.long[,4]>1,]

        temp.test_by_tag.all.long.valid[,1] = as.character(temp.test_by_tag.all.long.valid[,1])
    } else{
        temp.test_by_tag.all.long = data.frame(test_var = test_var,
                                               tag_vars = tag_vars,
                                               n_test_het = melt(test_by_tag.hets[test_var,tag_vars],stringsAsFactors=F)[,1],
                                               n_test_hom = melt(test_by_tag.homs[test_var,tag_vars],stringsAsFactors=F)[,1])

        temp.test_by_tag.all.long.valid = temp.test_by_tag.all.long[temp.test_by_tag.all.long[,3]>1 &
                                                                    temp.test_by_tag.all.long[,4]>1,]

        temp.test_by_tag.all.long.valid[,1] = as.character(temp.test_by_tag.all.long.valid[,1])
    }
        
    if (dim(temp.test_by_tag.all.long.valid)[1] > 0){
        ## get the best test-tag pair
        best.test_tag = temp.test_by_tag.all.long.valid[which.max(pwr.t2n.test(temp.test_by_tag.all.long.valid[,3],
                                                                               temp.test_by_tag.all.long.valid[,4],
                                                                               d=0.2,
                                                                               sig.level=0.05)$power),]
        
        het_indivs = row.names(test_hets)[test_hets[,test_var]==1 &
                                          tag_hets[,best.test_tag[,2]]==1]
        
        hom_indivs = row.names(test_hets)[test_homs[,test_var]==1 &
                                          tag_hets[,best.test_tag[,2]]==1]
        
        mfsd12_ase_results[i,]$test_var = test_var
        mfsd12_ase_results[i,]$tag_var  = as.character(best.test_tag[,2])
        mfsd12_ase_results[i,]$n_test_het = best.test_tag[,3]
        mfsd12_ase_results[i,]$n_test_hom = best.test_tag[,4]
        mfsd12_ase_results[i,]$wilcox_p = wilcox.test(ase_table.AE.wide[het_indivs,best.test_tag[,2]],
                                                      ase_table.AE.wide[hom_indivs,best.test_tag[,2]])$p.value
    } else{

        mfsd12_ase_results[i,3:6] = c(NA,NA,NA,NA)
    }
}






## carry out ASE for DDB1
ddb1_ase_results.expand = expand.grid(unique(as.character(ddb1_snp_egenes[,2])),names(gene_variants)[names(gene_variants) %in% ddb1_egenes],stringsAsFactors=F)
ddb1_ase_results = data.frame(test_var = ddb1_ase_results.expand[,1],
                                gene = ddb1_ase_results.expand[,2],
                                tag_var = 0,
                                n_test_het = 0,
                                n_test_hom = 0,
                                wilcox_p = 0,
                                stringsAsFactors=F)



## loop over the CAVIAR-gene pairs
for (i in 1:dim(ddb1_ase_results)[1]){

    ## get the caviar variant and the gene
    test_var = ddb1_ase_results[i,1]
    gene   = ddb1_ase_results[i,2]

    ## get the list of test and tag vars
    ## test_vars = ddb1_caviar[[caviar]]
    tag_vars = gene_variants[[gene]]

    if (length(tag_vars) > 1){
        ## get a long format list of heterozygotes and homozygotes
        temp.test_by_tag.all.long = cbind(test_var,merge(melt(test_by_tag.hets[test_var,tag_vars],stringsAsFactors=F),
                                                         melt(test_by_tag.homs[test_var,tag_vars],stringsAsFactors=F),by="row.names",stringsAsFactors=F))

        ## get pairs with at least 2 individuals for hets and homs
        temp.test_by_tag.all.long.valid = temp.test_by_tag.all.long[temp.test_by_tag.all.long[,3]>1 &
                                                                    temp.test_by_tag.all.long[,4]>1,]

        temp.test_by_tag.all.long.valid[,1] = as.character(temp.test_by_tag.all.long.valid[,1])
    } else{
        temp.test_by_tag.all.long = data.frame(test_var = test_var,
                                               tag_vars = tag_vars,
                                               n_test_het = melt(test_by_tag.hets[test_var,tag_vars],stringsAsFactors=F)[,1],
                                               n_test_hom = melt(test_by_tag.homs[test_var,tag_vars],stringsAsFactors=F)[,1])

        temp.test_by_tag.all.long.valid = temp.test_by_tag.all.long[temp.test_by_tag.all.long[,3]>1 &
                                                                    temp.test_by_tag.all.long[,4]>1,]

        temp.test_by_tag.all.long.valid[,1] = as.character(temp.test_by_tag.all.long.valid[,1])
    }
        
    if (dim(temp.test_by_tag.all.long.valid)[1] > 0){
        ## get the best test-tag pair
        best.test_tag = temp.test_by_tag.all.long.valid[which.max(pwr.t2n.test(temp.test_by_tag.all.long.valid[,3],
                                                                               temp.test_by_tag.all.long.valid[,4],
                                                                               d=0.2,
                                                                               sig.level=0.05)$power),]
        
        het_indivs = row.names(test_hets)[test_hets[,test_var]==1 &
                                          tag_hets[,best.test_tag[,2]]==1]
        
        hom_indivs = row.names(test_hets)[test_homs[,test_var]==1 &
                                          tag_hets[,best.test_tag[,2]]==1]
        
        ddb1_ase_results[i,]$test_var = test_var
        ddb1_ase_results[i,]$tag_var  = as.character(best.test_tag[,2])
        ddb1_ase_results[i,]$n_test_het = best.test_tag[,3]
        ddb1_ase_results[i,]$n_test_hom = best.test_tag[,4]
        ddb1_ase_results[i,]$wilcox_p = wilcox.test(ase_table.AE.wide[het_indivs,best.test_tag[,2]],
                                                      ase_table.AE.wide[hom_indivs,best.test_tag[,2]])$p.value
    } else{

        ddb1_ase_results[i,3:6] = c(NA,NA,NA,NA)
    }
}









## perform eQTL analysis
for (i in dim(mfsd12_candidates){

    tag_vars = gene_variants[[mfsd12_candidates[i,2]]]

    het_tag = list()
    for (tag_var in tag_vars){
        het_tag ## START HERE
    }

    het_test = (data[,mfsd12_caviar[[test_var]]]==1)*1
    row.names(het_test) = data$Sample

        for 
        tag_hets = row.names(ase_table.het.wide)[apply(ase_table.het.wide[,tag_vars],1,sum)>1]
    }

    for (test_var in names(ddb1_caviar)){

        
    }
}




## get best variant to test per gene
ase_genes = unique(ase_data$SYMBOL)

gene_topvar = list()
for (gene in ase_genes){

    var_temp = list()
    for (var in c(ddb1_vars,mfsd12_vars)){
        ## get the coding variant with the most individuals that are
        ## heterozygous for var
        var_temp[[var]]$coding_count = sum(ase_data$SYMBOL == gene &
                                           ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                           (grepl("missense",ase_data$Consequence) |
                                            grepl("synonymous",ase_data$Consequence)))

        var_temp[[var]]$coding = ase_data[ase_data$SYMBOL == gene &
                                              ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                              (grepl("missense",ase_data$Consequence) |
                                               grepl("synonymous",ase_data$Consequence)),]$variantID %>%
                                     table %>% which.max %>% names

        ## get the utr variant with the most individuals that are
        ## heterozygous for var
        var_temp[[var]]$utr_count = sum(ase_data$SYMBOL == gene &
                                        ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                        grepl("UTR",ase_data$Consequence))


        var_temp[[var]]$utr = ase_data[ase_data$SYMBOL == gene &
                                       ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                       grepl("UTR",ase_data$Consequence),]$variantID %>%
                              table %>% which.max %>% names

        ## get the intron variant with the most individuals that are
        ## heterozygous for var
        var_temp[[var]]$intron_count = sum(ase_data$SYMBOL == gene &
                                           ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                           grepl("intron",ase_data$Consequence))


        var_temp[[var]]$intron = ase_data[ase_data$SYMBOL == gene &
                                          ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                          grepl("intron",ase_data$Consequence),]$variantID %>%
                                 table %>% which.max %>% names

        ## get the splice variant with the most individuals that are
        ## heterozygous for var
        var_temp[[var]]$splice_count = sum(ase_data$SYMBOL == gene &
                                           ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                           grepl("splice",ase_data$Consequence))


        var_temp[[var]]$splice = ase_data[ase_data$SYMBOL == gene &
                                          ase_data$Sample %in% data[data[,var]==1,]$Sample &
                                          grepl("splice",ase_data$Consequence),]$variantID %>%
                                 table %>% which.max %>% names
    }

    gene_topvar[[gene]] = var_temp
}



# run Mann-Whitney for heterozygotes and homozygotes
ase_coding = as.data.frame(matrix(nrow=length(ase_genes), ncol=length(c(mfsd12_vars,ddb1_vars)))); row.names(ase_coding) = ase_genes; colnames(ase_coding) = c(mfsd12_vars,ddb1_vars); ase_coding$gene = ase_genes
ase_utr = as.data.frame(matrix(nrow=length(ase_genes), ncol=length(c(mfsd12_vars,ddb1_vars)))); row.names(ase_utr) = ase_genes; colnames(ase_utr) = c(mfsd12_vars,ddb1_vars); ase_utr$gene = ase_genes
ase_splice = as.data.frame(matrix(nrow=length(ase_genes), ncol=length(c(mfsd12_vars,ddb1_vars)))); row.names(ase_splice) = ase_genes; colnames(ase_splice) = c(mfsd12_vars,ddb1_vars); ase_splice$gene = ase_genes
for (gene in ase_genes){

    for (var in c(mfsd12_vars,ddb1_vars)){

        coding_var = gene_topvar[[gene]][[var]]$coding
        utr_var = gene_topvar[[gene]][[var]]$utr
        splice_var = gene_topvar[[gene]][[var]]$splice
        
        ## Samples that have AE information and have the coding variant
        ## described and are heterozygous for the variant
        het_coding = ase_data[ase_data$SYMBOL==gene &
                              ase_data$variantID==coding_var &
                              ase_data$Sample %in% data[data[,var]==1,]$Sample,]$AE

        ## Samples that have AE information and have the coding variant
        ## described and are heterozygous for the variant
        hom_coding = ase_data[ase_data$SYMBOL==gene & 
                              ase_data$variantID==coding_var &
                              ase_data$Sample %in% data[!(data[,var]==1),]$Sample,]$AE


        ## Samples that have AE information and have the utr variant
        ## described and are heterozygous for the variant
        het_utr = ase_data[ase_data$SYMBOL==gene & 
                           ase_data$variantID==utr_var &
                           ase_data$Sample %in% data[data[,var]==1,]$Sample,]$AE

        ## Samples that have AE information and have the utr variant
        ## described and are heterozygous for the variant
        hom_utr = ase_data[ase_data$SYMBOL==gene & 
                           ase_data$variantID==utr_var &
                           ase_data$Sample %in% data[!(data[,var]==1),]$Sample,]$AE

        
        ## Samples that have AE information and have the splice variant
        ## described and are heterozygous for the variant
        het_splice = ase_data[ase_data$SYMBOL==gene & 
                              ase_data$variantID==splice_var &
                              ase_data$Sample %in% data[data[,var]==1,]$Sample,]$AE

        ## Samples that have AE information and have the splice variant
        ## described and are heterozygous for the variant
        hom_splice = ase_data[ase_data$SYMBOL==gene & 
                              ase_data$variantID==splice_var &
                              ase_data$Sample %in% data[!(data[,var]==1),]$Sample,]$AE


        
        if (length(het_coding) > 0 & length(hom_coding) > 0){
            
            ase_coding[gene,var] = wilcox.test(abs(hom_coding),abs(het_coding))$p.value
        }
        else {
            ase_coding[gene,var] = NA
        }

        if (length(het_utr) > 0 & length(hom_utr) > 0){
            
            ase_utr[gene,var] = wilcox.test(abs(hom_utr),abs(het_utr))$p.value
        }
        else {
            ase_utr[gene,var] = NA
        }

        if (length(het_splice) > 0 & length(hom_splice) > 0){
            
            ase_splice[gene,var] = wilcox.test(abs(hom_splice),abs(het_splice))$p.value
        }
        else {
            ase_splice[gene,var] = NA
        }
    }
}













## Plots with error bars
barplot_data = data.frame(Sample = data$Sample,
                          YRI = data$YRI,
                          mfsd12 = data$mfsd12,
                          ddb1 = data$ddb1,
                          vps37c = data$vps37c,
                          rab3il1 = data$rab3il1,
                          pigment = data$Pigmentation.scale,
                          ancestry = "med_YRI", stringsAsFactors=F)

barplot_data$ancestry[barplot_data$YRI>=0.7] = "high_YRI"
barplot_data$ancestry[barplot_data$YRI<=0.1] = "low_YRI"

barplot_data$pigment_level = "medium"
barplot_data$pigment_level[barplot_data$pigment < 3] = "low"
barplot_data$pigment_level[barplot_data$pigment > 6] = "high"

barplot_data$ancestry = as.factor(barplot_data$ancestry)
barplot_data$pigment_level = as.factor(barplot_data$pigment_level)

YRI_v_mfsd12.df = barplot_data[barplot_data$ancestry=="high_YRI" |
                               barplot_data$ancestry=="low_YRI",]

YRI_v_mfsd12.plot = ggplot(YRI_v_mfsd12.df,aes(x=ancestry, y=mfsd12)) +
    geom_boxplot() + annotate("text",x=1,y=1,label="p=0.01") + theme_classic()## +
    ## geom_errorbar(aes(ymin=len-ci, ymax=len+ci),
                  ## width=.2,                    # Width of the error bars
                  ## position=position_dodge(.9))

pigment_v_mfsd12.df = barplot_data[barplot_data$pigment_level=="high" |
                               barplot_data$pigment_level=="low",]
pigment_v_mfsd12.plot = ggplot(pigment_v_mfsd12.df,aes(x=pigment_level, y=mfsd12)) +
    geom_boxplot() + annotate("text",x=1.5,y=1.5,label="p=0.70") + theme_classic()## +



YRI_v_ddb1.df = barplot_data[barplot_data$ancestry=="high_YRI" |
                               barplot_data$ancestry=="low_YRI",]

YRI_v_ddb1.plot = ggplot(YRI_v_ddb1.df,aes(x=ancestry, y=ddb1)) +
    geom_boxplot() + annotate("text",x=1,y=-1,label="p=0.028") + theme_classic()## +
    ## geom_errorbar(aes(ymin=len-ci, ymax=len+ci),
                  ## width=.2,                    # Width of the error bars
                  ## position=position_dodge(.9))

pigment_v_ddb1.df = barplot_data[barplot_data$pigment_level=="high" |
                               barplot_data$pigment_level=="low",]
pigment_v_ddb1.plot = ggplot(pigment_v_ddb1.df,aes(x=pigment_level, y=ddb1)) +
    geom_boxplot() + annotate("text",x=2,y=1,label="p=0.001") + theme_classic(base_size=20)## +


pigment_v_vps37c.df = barplot_data[barplot_data$pigment_level=="high" |
                                   barplot_data$pigment_level=="low",]
pigment_v_vps37c.plot = ggplot(pigment_v_vps37c.df,aes(x=pigment_level, y=vps37c)) +
    geom_boxplot() + annotate("text",x=2,y=1.5,label="p=0.005") + theme_classic(base_size=20)## +


pigment_v_rab3il1.df = barplot_data[barplot_data$pigment_level=="high" |
                                   barplot_data$pigment_level=="low",]
pigment_v_rab3il1.plot = ggplot(pigment_v_rab3il1.df,aes(x=pigment_level, y=rab3il1)) +
    geom_boxplot() + annotate("text",x=2,y=2,label="p=0.08") + theme_classic(base_size=20)## +


ggsave("figures/pigment_v_ddb1.png",pigment_v_ddb1.plot,width=5,height=5)
ggsave("figures/pigment_v_vps37c.png",pigment_v_vps37c.plot,width=5,height=5)
ggsave("figures/pigment_v_rab3il1.png",pigment_v_rab3il1.plot,width=5,height=5)


rs6510760_mfsd12.plot = ggplot(data[!is.na(data$mfsd12),]) +
    geom_boxplot(aes(rs6510760,mfsd12)) +
    geom_point(aes(jitter(rs6510760),mfsd12),color="gray") +
    theme_classic()





ddb1_rs57265008.lm = summary(lm(pigment ~ rs57265008 + vps37c + rab3il1, data))
pigment_ddb1.df = data.frame(ddb1 = data[-ddb1_rs57265008.lm$na.action,]$ddb1, pigment.resid = ddb1_rs57265008.lm$residuals)


vps37c_rs57265008.lm = summary(lm(pigment ~ rs57265008 + ddb1 + rab3il1, data))
pigment_vps37c.df = data.frame(vps37c = data[-vps37c_rs57265008.lm$na.action,]$vps37c, pigment.resid = vps37c_rs57265008.lm$residuals)


rab3il1_rs57265008.lm = summary(lm(pigment ~ rs57265008 + ddb1 + vps37c, data))
pigment_rab3il1.df = data.frame(rab3il1 = data[-rab3il1_rs57265008.lm$na.action,]$rab3il1, pigment.resid = rab3il1_rs57265008.lm$residuals)

rs6510760_hmg20b.plot











########################################
############### SPLICING ###############
########################################





splice_list = list()
for (var in mfsd12_testvars){

    temp_splice = list()
    for (site in splice_sites){

        logodds = log(data[,site]/(1-data[,site]))
        temp_splice[[site]] = lm(build_lm(paste("log(",site,"/(1-",site,"))",sep=""),c(covars,"log(YRI)",var)),data)
    }
    splice_list[[var]] = temp_splice
}
    









    


#######################################################
######################## PLOTS ########################
#######################################################

ddb1_YRI_eqtl_signif.plot.df = data.frame()
ddb1_eqtl_signif.plot.df = data.frame()
for (gene in unique(ddb1_eqtl_signif[,1])){
    for (posit in unique(ddb1_eqtl_signif[,2])){

        temp_lm = gene.norm[[gene]]
        temp_YRI_lm = gene.YRI.norm[[gene]]

        temp_df = data.frame(gene=gene, posit=posit, genotype=data[-temp_lm$na.action,posit], resid_exp=temp_lm$residuals)
        temp_YRI_df = data.frame(gene=gene, posit=posit, genotype=data[-temp_lm$na.action,posit], resid_exp=temp_lm$residuals)
        ddb1_eqtl_signif.plot.df = rbind(ddb1_eqtl_signif.plot.df,temp_df)
    }
}

## ddb1_eqtl_signif.bee.df = cbind(ddb1_eqtl_signif.plot.df[,c("gene","posit")],beeswarm(resid_exp ~ genotype,ddb1_eqtl_signif.plot.df,priority="none"))

ddb1_eqtl_signif.plot = ggplot(ddb1_eqtl_signif.plot.df,aes(factor(genotype),resid_exp)) + geom_boxplot() + geom_smooth(method=lm,se=F,fullrange=T) + theme_bw() + facet_grid(gene ~ posit)


YRI_pigment.df = data.frame(YRI=data[-summary(pigment.lm)$na.action,]$YRI, pigment_resid=summary(pigment.lm)$residuals)
YRI_pigment.label = data.frame(x=-2.5, y=-2, label = "p-value = 2.97e-05")
YRI_pigment.plot = ggplot(YRI_pigment.df,aes(log2(YRI),pigment_resid)) + geom_point() + geom_text(data=YRI_pigment.label,aes(x=x,y=y,label=label)) + geom_smooth(method=lm,se=F,fullrange=T) + theme_bw() + labs(x="log2(Fraction YRI Ancestry)", y="Residualized Pigmentation Score", title="Pigmentation vs. Ancestry")



egenes_combos = expand.grid(egenes[,2],egenes[,2],stringsAsFactors=F)
egenes_combos = egenes_combos[!(egenes_combos[,1]==egenes_combos[,2]),]
cor_gene.df = data.frame(gene1=egenes_combos[,1], gene2=egenes_combos[,2], rho=mat.or.vec(dim(egenes_combos)[1],1), p.value=mat.or.vec(dim(egenes_combos)[1],1))
for (i in 1:dim(egenes_combos)[1]){

    gene1 = egenes_combos[i,1]; gene2 = egenes_combos[i,2]
    gene1_exp = data[complete,gene1]; gene2_exp = data[complete,gene2]

    cor_results = cor.test(gene1_exp,gene2_exp,method="spearman")
    
    cor_gene.df[i,]$rho = cor_results$estimate
    cor_gene.df[i,]$p.value = cor_results$p.value
}

cor_gene.plot = ggplot(cor_gene.df,aes(gene1_exp,gene2_exp)) + geom_point() + geom_smooth(method=lm,se=F,fullrange=T) + theme_bw() + facet_grid(gene1 ~ gene2)

######################## THEMES ########################
## scale_colour_ptol



################# DATA #################


# Ancestry and pigmentation
# pigmentation linear model

pigment.lm = lm(pigmentation_scale ~ Jiyeon_frozen + final_pass + days_cult + confluency + Mixed_shape + Starburst_shape + growth,covariates.all.df)
pigment.YRI.lm = lm(pigmentation_scale ~ Jiyeon_frozen + final_pass + days_cult + confluency + Mixed_shape + Starburst_shape + growth + log(YRI),covariates.all.df)


# get p-value for YRI-pigmentation
YRI_pigmentation_pval = coef(summary(pigment.YRI.lm))[,4]["log(YRI)"]


# form dataframe
pigment.lm.resid = data.frame(pigment.lm$residuals)
pigment.lm.resid$Sample = row.names(pigment.lm.resid)
pigment.YRI.lm.resid = data.frame(pigment.YRI.lm$residuals)
pigment.YRI.lm.resid$Sample = row.names(pigment.YRI.lm.resid)


# ancestry and gene expression
expression_norm.long = melt(expression_norm[,c(1,6:111)],id.vars="Gene.Symbol")
colnames(expression_norm.long) = c("gene","Sample","expression_norm")

expression_norm.resid = data.frame()
expression_norm.YRI.resid = data.frame()

expression_norm.lm = list()
expression_norm.YRI.lm = list()
for (gene in egenes){

    gene_expression = data.frame(t(expression_norm[expression_norm$Gene.Symbol==gene,6:111])) # get gene expression
    colnames(gene_expression) = "expression"
    gene_expression$Sample = row.names(gene_expression)

    # data for linear models
    gene_lm.df = merge(gene_expression,covariates.all.df,by="Sample",all=T)
    row.names(gene_lm.df) = gene_lm.df$Sample # keep track of sample names

    # calculate linear models
    gene_lm.lm = lm(expression ~ Jiyeon_frozen + final_pass + days_cult + confluency + Mixed_shape + Starburst_shape + growth, gene_lm.df) # linear model w/o YRI
    gene_lm.YRI.lm = lm(expression ~ Jiyeon_frozen + final_pass + days_cult + confluency + Mixed_shape + Starburst_shape + growth + log(YRI), gene_lm.df) # linear model w/ YRI

    expression_norm.resid = rbind(expression_norm.resid, cbind(gene,names(gene_lm.lm$residuals),gene_lm.lm$residuals))
    expression_norm.YRI.resid = rbind(expression_norm.YRI.resid, cbind(gene,names(gene_lm.YRI.lm$residuals),gene_lm.YRI.lm$residuals))

    expression_norm.lm[[gene]] = gene_lm.lm
    expression_norm.YRI.lm[[gene]] = gene_lm.YRI.lm
}

colnames(expression_norm.resid) = c("gene","Sample","expression_norm.resid")
colnames(expression_norm.YRI.resid) = c("gene","Sample","expression_norm.YRI.resid")

expression_norm.resid = factor2numeric(expression_norm.resid,"expression_norm.resid")
expression_norm.YRI.resid = factor2numeric(expression_norm.YRI.resid,"expression_norm.YRI.resid")

expression_norm.resid.wide = cast(expression_norm.resid, Sample~gene)
expression_norm.YRI.resid.wide = cast(expression_norm.YRI.resid, Sample~gene)





############## PUBLICATION PLOTS ##############

# Epigenomics Roadmap gene expression
epi_dir = "/local1/derek/data/epi_roadmap/all"

epi_rna.data = read.delim(file.path(epi_dir,"rna_seq/57epigenomes.RPKM.pc"), sep="\t", header=T, row.names=1)
epi_names = read.delim(file.path(epi_dir,"rna_seq/EG.name.txt"), sep="\t", header=F, row.names=1)

epi_meta = read.delim(file.path(epi_dir,"EID_metadata.tab"), sep="\t", header=T, row.names=1)

epi_color = epi_meta$COLOR; names(epi_color) = row.names(epi_meta)


## epi_rna.data.long = melt(epi_rna.data, id.vars="gene_id");

# get gene names
gene_ids = unlist(strsplit(as.character(expression_norm$gene_id),split="\\."))[seq(from=1,to=dim(expression_norm)[1]*2,by=2)]
names(gene_ids) = expression_norm$Gene.Symbol

# plot gene expression of genes across tissues
mfsd12_tiss_exp.data = melt(epi_rna.data[gene_ids["mfsd12"],])
colnames(mfsd12_tiss_exp.data) = c("Sample","RPKM")
mfsd12_tiss_exp.data$Tissue = epi_names[as.character(mfsd12_tiss_exp.data[,1]),1]
mfsd12_tiss_exp.data$Color = epi_meta[as.character(mfsd12_tiss_exp.data$Sample),]$COLOR


mfsd12_tiss_exp.plot = ggplot(mfsd12_tiss_exp.data,aes(x=Tissue,y=log10(RPKM),fill=Sample)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=as.character(mfsd12_tiss_exp.data$Color)) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(margin = margin(r=20)), axis.title.x = element_text(margin = margin(t=20)), axis.text.x=element_text(angle = -90, hjust = 0)) +
    ggtitle("MFSD12 expression in Epigenomics Roadmap")




important_genes = c("ddb1","cyb561a3","fzr1","tmem138","mfsd12","hmg20b","gipc3")

all_tiss_exp.data = melt(cbind(as.factor(important_genes),epi_rna.data[gene_ids[important_genes],]))
colnames(all_tiss_exp.data) = c("egene","Sample","RPKM")
all_tiss_exp.data$Tissue = factor(epi_names[as.character(all_tiss_exp.data$Sample),1], levels=epi_names[as.character(all_tiss_exp.data$Sample),1])
epi_roadmap_col = factor(epi_meta[as.character(unique(all_tiss_exp.data$Sample)),]$COLOR, levels=epi_meta[as.character(unique(all_tiss_exp.data$Sample)),]$COLOR)

all_tiss_exp.plot = ggplot(all_tiss_exp.data,aes(x=Tissue,y=RPKM,fill=Sample)) + geom_bar(stat="identity") + facet_grid(egene ~ .,scale="free_y") + scale_fill_manual(values=as.character(epi_roadmap_col)) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(margin = margin(r=20)), axis.title.x = element_text(margin = margin(t=20)), axis.text.x=element_text(angle = -90, hjust = 0)) + ggtitle("All gene expression in Epigenomics Roadmap")
b











#################################
########## EPI ROADMAP ##########
#################################


# Epigenomics Roadmap gene expression
epi_dir = "/local1/derek/data/epi_roadmap/all"

epi_rna.data = read.delim(file.path(epi_dir,"rna_seq/57epigenomes.RPKM.pc"), sep="\t", header=T, row.names=1)
epi_names = read.delim(file.path(epi_dir,"rna_seq/EG.name.txt"), sep="\t", header=F, row.names=1)

epi_meta = read.delim(file.path(epi_dir,"EID_metadata.tab"), sep="\t", header=T, row.names=1)

epi_color = as.character(epi_meta$COLOR); names(epi_color) = row.names(epi_meta)


## get gene names
gene_ids = get_split(as.character(ddb1_egenes.mel$Ensemble_ID),"\\.",1)
names(gene_ids) = ddb1_egenes.mel$Gene_Symbol
## gene_ids["vps37c"] = "ENSG00000167987"


important_genes = c("ddb1","cyb561a3","fzr1","tmem138","mfsd12","hmg20b","gipc3","ppp1r32","sdhaf2","cpsf7")## ,"vps37c")

important_genes = ddb1_egenes.mel[,2]

all_tiss_exp.data = melt(cbind(as.factor(important_genes),epi_rna.data[gene_ids[important_genes],]))
colnames(all_tiss_exp.data) = c("egene","Sample","RPKM")
all_tiss_exp.data$Tissue = factor(epi_names[as.character(all_tiss_exp.data$Sample),1], levels=unique(epi_names[as.character(all_tiss_exp.data$Sample),1]))
epi_roadmap_col = as.character(epi_meta[as.character(unique(all_tiss_exp.data$Sample)),]$COLOR)
names(epi_roadmap_col) = factor(unique(all_tiss_exp.data$Sample))
epi_roadmap_col["E000"] = "#000000"

all_tiss_exp.plot = ggplot(all_tiss_exp.data,aes(x=Sample,y=RPKM,fill=Sample)) +
    geom_bar(stat="identity") + facet_grid(egene ~ .,scale="free_y") +
    scale_fill_manual(values=epi_roadmap_col) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(margin = margin(r=20)), axis.title.x = element_text(margin = margin(t=20)), axis.text.x=element_text(angle = -90, hjust = 0)) +
    ggtitle("All gene expression in Epigenomics Roadmap")

ggsave("figures/ddb1_egenes_mel_epiroadmap_expression.pdf",all_tiss_exp.plot,width=8,height=15)

ggsave("figures/ddb1_egenes_mel_epiroadmap_expression.png",all_tiss_exp.plot,width=8,height=15)







## perform with ddb1
ddb1_ase_results.expand = expand.grid(names(ddb1_caviar),names(gene_variants)[names(gene_variants) %in% ddb1_egenes.mel[,2]],stringsAsFactors=F)
ddb1_ase_results = data.frame(caviar_var = ddb1_ase_results.expand[,1],
                                gene = ddb1_ase_results.expand[,2],
                                test_var = 0,
                                tag_var = 0,
                                n_test_het = 0,
                                n_test_hom = 0,
                                wilcox_p = 0,
                                stringsAsFactors=F)

## loop over the CAVIAR-gene pairs
for (i in 1:dim(ddb1_ase_results)[1]){

    ## get the caviar variant and the gene
    caviar = ddb1_ase_results[i,1]
    gene   = ddb1_ase_results[i,2]

    ## get the list of test and tag vars
    test_vars = ddb1_caviar[[caviar]]
    tag_vars = gene_variants[[gene]]

    if (length(tag_vars) > 1){
        ## get a long format list of heterozygotes and homozygotes
        temp.test_by_tag.all.long = merge(melt(test_by_tag.hets[test_vars,tag_vars],stringsAsFactors=F),
                                          melt(test_by_tag.homs[test_vars,tag_vars],stringsAsFactors=F),by=c(1,2),stringsAsFactors=F)

        ## get pairs with at least 2 individuals for hets and homs
        temp.test_by_tag.all.long.valid = temp.test_by_tag.all.long[temp.test_by_tag.all.long[,3]>1 &
                                                                    temp.test_by_tag.all.long[,4]>1,]

        temp.test_by_tag.all.long.valid[,1] = as.character(temp.test_by_tag.all.long.valid[,1])
        temp.test_by_tag.all.long.valid[,2] = as.character(temp.test_by_tag.all.long.valid[,2])
    } else{
        temp.test_by_tag.all.long = data.frame(test_vars = names(test_by_tag.hets[test_vars,tag_vars]),
                                               tag_vars = tag_vars,
                                               n_test_het = melt(test_by_tag.hets[test_vars,tag_vars],stringsAsFactors=F)[,1],
                                               n_test_hom = melt(test_by_tag.homs[test_vars,tag_vars],stringsAsFactors=F)[,1])

        temp.test_by_tag.all.long.valid = temp.test_by_tag.all.long[temp.test_by_tag.all.long[,3]>1 &
                                                                    temp.test_by_tag.all.long[,4]>1,]

        temp.test_by_tag.all.long.valid[,1] = as.character(temp.test_by_tag.all.long.valid[,1])
        temp.test_by_tag.all.long.valid[,2] = as.character(temp.test_by_tag.all.long.valid[,2])
    }
        
    if (dim(temp.test_by_tag.all.long.valid)[1] > 0){
        ## get the best test-tag pair
        best.test_tag = temp.test_by_tag.all.long.valid[which.max(pwr.t2n.test(temp.test_by_tag.all.long.valid[,3],
                                                                               temp.test_by_tag.all.long.valid[,4],
                                                                               d=0.2,
                                                                               sig.level=0.05)$power),]
        
        het_indivs = row.names(test_hets)[test_hets[,best.test_tag[,1]]==1 &
                                          tag_hets[,best.test_tag[,2]]==1]
        
        hom_indivs = row.names(test_hets)[test_homs[,best.test_tag[,1]]==1 &
                                          tag_hets[,best.test_tag[,2]]==1]
        
        ddb1_ase_results[i,]$test_var = best.test_tag[,1]
        ddb1_ase_results[i,]$tag_var  = best.test_tag[,2]
        ddb1_ase_results[i,]$n_test_het = best.test_tag[,3]
        ddb1_ase_results[i,]$n_test_hom = best.test_tag[,4]
        ddb1_ase_results[i,]$wilcox_p = wilcox.test(ase_table.AE.wide[het_indivs,best.test_tag[,2]],
                                                      ase_table.AE.wide[hom_indivs,best.test_tag[,2]])$p.value
    } else{

        ddb1_ase_results[i,3:7] = c(NA,NA,NA,NA,NA)
    }
}

ddb1_ase_BH_rs7948623_pval = data.frame(gene = ddb1_ase_results[ddb1_ase_results[,1]=="rs7948623",][order(p.adjust(ddb1_ase_results[ddb1_ase_results[,1]=="rs7948623",7],method="BH")),2],
                                        BH_pval = sort(p.adjust(ddb1_ase_results[ddb1_ase_results[,1]=="rs7948623",7],method="BH"),na.last=T))

ddb1_ase_BH_rs57265008_pval = data.frame(gene = ddb1_ase_results[ddb1_ase_results[,1]=="rs57265008",][order(p.adjust(ddb1_ase_results[ddb1_ase_results[,1]=="rs57265008",7],method="BH")),2],
                                        BH_pval = sort(p.adjust(ddb1_ase_results[ddb1_ase_results[,1]=="rs57265008",7],method="BH"),na.last=T))











fig_data = data.frame(hmg20b=data$hmg20b,mfsd12=data$mfsd12,vps37c=data$vps37c,rab3il1=data$rab3il1,ddb1=data$ddb1, rs6510760=c("G/G","G/A","A/A")[data$rs6510760+1],rs4939519=c("C/C","C/A","A/A")[data$rs4939519+1])



fig_data = data.frame(hmg20b=data$hmg20b,
                      mfsd12=data$mfsd12,
                      vps37c=data$vps37c,
                      rab3il1=data$rab3il1,
                      ddb1=data$ddb1,
                      herc2=data$herc2,
                      oca2=data$oca2,
                      rs6510760=factor(c("G/G","G/A","A/A")[data$rs6510760+1],levels=c("G/G","G/A","A/A")),
                      rs148172827=factor(c("ATCAA/ATCAA","ATCAA/-","-/-")[data$rs148172827+1],levels=c("ATCAA/ATCAA","ATCAA/-","-/-")),
                      rs6497271=factor(c("A/A","A/G","G/G")[data$rs6497271+1],levels=c("A/A","A/G","G/G")),
                      rs4932620=factor(c("T/T","T/C","C/C")[data$rs4932620+1],levels=c("T/T","T/C","C/C")))


mfsd12_egenes=c("mfsd12","hmg20b")
ddb1_egenes=c("ddb1","vps37c","rab3il1")
herc2_egenes=c("herc2","oca2")

fig_data.ddb1 = melt(fig_data[,c(ddb1_egenes,"rs148172827")], id="rs148172827", value="gene")
fig_data.mfsd12 = melt(fig_data[,c(mfsd12_egenes,"rs6510760")], id="rs6510760", value="gene")
fig_data.region2.herc2 = melt(fig_data[,c(herc2_egenes,"rs6497271")], id="rs6497271", value="gene")
fig_data.region3.herc2 = melt(fig_data[,c(herc2_egenes,"rs4932620")], id="rs4932620", value="gene")


colnames(fig_data.ddb1) = c("genotype","gene","normalized_expression")
colnames(fig_data.mfsd12) = c("genotype","gene","normalized_expression")
colnames(fig_data.region2.herc2) = c("genotype","gene","normalized_expression")
colnames(fig_data.region3.herc2) = c("genotype","gene","normalized_expression")


mfsd12.plot = ggplot(na.omit(fig_data.mfsd12),aes(x=genotype,y=normalized_expression)) + geom_beeswarm() + geom_smooth(method='lm') + facet_grid(. ~ gene) + theme_classic()
ddb1.plot = ggplot(na.omit(fig_data.ddb1),aes(x=genotype,y=normalized_expression)) + geom_beeswarm() + geom_smooth(method='lm') + facet_grid(. ~ gene) + theme_classic()
herc2.region2.plot = ggplot(na.omit(fig_data.region2.herc2),aes(x=genotype,y=normalized_expression)) + geom_beeswarm() + geom_smooth(method='lm') + facet_grid(. ~ gene) + theme_classic() + labs(title="rs6497271")


herc2.region3.plot = ggplot(na.omit(fig_data.region3.herc2[fig_data.region3.herc2$gene=="oca2",]),aes(x=genotype,y=normalized_expression)) + geom_beeswarm() + geom_smooth(method='lm') + theme_classic() + labs(title="rs4932620")


oca2.YRI.plot = ggplot(data,aes(x=log2(YRI),y=oca2)) + geom_beeswarm() + geom_smooth(method='lm') + theme_classic() + labs(title="YRI ancestry vs. OCA2")
