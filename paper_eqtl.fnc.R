#######################################################
###################### FUNCTIONS ######################
#######################################################


# converts factor rows to numeric
factor2numeric = function(dframe,cols){
    for (col in cols){
        dframe[,col] = as.numeric(levels(dframe[,col]))[dframe[,col]]
    }
    return(dframe)
}


# estimates genotypes based on output from bam-readcount
get_genotype = function(counts_mat){

    vars = colnames(counts_mat)
    indivs = row.names(counts_mat)
    genos_df = data.frame(mat.or.vec(length(indivs),length(vars)))
    colnames(genos_df) = vars
    row.names(genos_df) = indivs

        for (i in 1:length(vars)){

            var_list = strsplit(as.character(counts_mat[,vars[i]]),split=":")

            var_mat = t(matrix(as.numeric(unlist(var_list)),nrow=4,ncol=length(var_list)))
            var_sums = order(apply(var_mat,2,sum))

            all_miss = apply(var_mat,1,sum)<2
            
            alt_alleles = (var_mat[,var_sums[3]] > var_mat[,var_sums[4]]/4) & (var_mat[,var_sums[3]] > 1)
            ref_alleles = (var_mat[,var_sums[4]] > var_mat[,var_sums[3]]/4)
            genos_df[alt_alleles,vars[i]] = 1
            genos_df[(!ref_alleles) & alt_alleles,vars[i]] = 2
            genos_df[all_miss,vars[i]] = NA
        }

    return(genos_df)
}


# estimates allelic ratio based on output of bam-readcount
get_allelic_ratio = function(counts_mat){

    vars = colnames(counts_mat)
    indivs = row.names(counts_mat)
    allrat.df = data.frame(mat.or.vec(length(indivs),length(vars)))
    colnames(allrat.df) = vars
    row.names(allrat.df) = indivs

        for (i in 1:length(vars)){

            var_list = strsplit(as.character(counts_mat[,vars[i]]),split=":")

            var_mat = t(matrix(as.numeric(unlist(var_list)),nrow=4,ncol=length(var_list)))
            var_sums = order(apply(var_mat,2,sum))

            all_miss = apply(var_mat,1,sum)<2
            
            allelic_ratio = var_mat[,var_sums[3]]/apply(var_mat[,var_sums[3:4]],1,sum)
            allrat.df[,vars[i]] = allelic_ratio
        }

    return(allrat.df)
}




## convert genotypes between text format and numeric
convert_geno = function(df,geno_cols){

    other_col = colnames(df)[!(colnames(df) %in% geno_cols)]

    temp_geno = data.frame(mat.or.vec(dim(df)[1],length(geno_cols)))
    colnames(temp_geno) = geno_cols
    
    for (col in geno_cols){

        na_elements = is.na(df[,col])
        
        temp_geno[!na_elements & grepl("0\\1",df[,col]), col] = 1
        temp_geno[!na_elements & grepl("1\\0",df[,col]), col] = 1
        temp_geno[!na_elements & grepl("1\\1",df[,col]), col] = 2
        temp_geno[!na_elements & df[,col]=="H", col] = 1
        temp_geno[na_elements, col] = NA
    }

    return(cbind(df[,other_col],temp_geno))
}


## builds and returns an lm function
build_lm = function(y,x){

    if (length(x) == 1){
        RHS = x
    }

    else{
        RHS = paste(x, collapse=" + ")
    }
    
    return( paste(y, RHS, sep=" ~ ") )
}


# this makes strings lowercase and replaces troublesome symbols
fix_strings = function(strings){

    return(gsub("-","_",gsub(":",".",sapply(strings,tolower))))
}


# gets the nth element after splitting strings. assumes all strings
# have the same number of split_char
get_split = function(strings,split_char,element_num){

    # split the strings
    split_list = strsplit(strings,split_char)

    # get number of elements
    num_elements = length(split_list[[1]])

    # get vector logical vector to extract elements
    extract_vec = mat.or.vec(num_elements,1)
    extract_vec[element_num] = 1
    extract_vec = as.logical(extract_vec)

    return( unlist(split_list)[extract_vec] )
}


get_melanocyte_data = function(wd){

    setwd(wd)
    
    ## filenames
    ## f_expression = "expression/master_table.txt"
    f_exp_mfsd12 = "expression/DDB1_MFSD12_RegionGenes_RSEM_normalized_values_5.15.17.MFSD12.txt"
    f_exp_ddb1 = "expression/DDB1_MFSD12_RegionGenes_RSEM_normalized_values_5.15.17.DDB1.txt"
    f_exp_herc2 = "expression/herc2_peerfactor_normalized_expression.txt"
    f_spl_mfsd12 = "../splicing/mfsd12/mfsd12_splice_frac.txt"
    ## f_geno_mfsd12 = "genotypes/5.18.15.17_genotypes_NCI_sequencing_and_SKLdata_MFSD12_and_DDB1.MFSD12.txt"
    ## f_geno_ddb1 = "genotypes/5.18.15.17_genotypes_NCI_sequencing_and_SKLdata_MFSD12_and_DDB1.DDB1.txt"
    f_geno_mfsd12 = "genotypes/6.26.17_chr11.61115821.genotype_NCI_sequencing_and_SKLdata_MFSD12.txt"
    f_geno_ddb1 = "genotypes/6.26.17_chr11.61115821.genotype_NCI_sequencing_and_SKLdata_DDB1.txt"
    f_geno_herc2 = "genotypes/OCA2_VCF_for_STishkoff"
    f_covariates = "covariates/Melanocytes_020316_MFSD12_RNASEQ_genotypes.covariates.txt"
    f_ancestry = "covariates/Melanocytes_ADMIXTURE_ancestry_table_v1.txt"
    f_gencode = "gencode.v19.annotation.gene_id.map.txt"
    f_pca = "covariates/datafinal_pca3_covariates.txt"
    
    
    ## read in expression, genotype, covariate, and ancestry data
    ## expression  = read.delim(f_expression,header=T,sep="\t",stringsAsFactors=F)
    exp_mfsd12  = read.delim(f_exp_mfsd12,header=T,sep="\t",stringsAsFactors=F)
    exp_ddb1    = read.delim(f_exp_ddb1,header=T,sep="\t",stringsAsFactors=F)
    exp_herc2   = read.delim(f_exp_herc2,header=T,sep="\t",stringsAsFactors=F)
    spl_mfsd12  = read.delim(f_spl_mfsd12,header=T,sep="\t",stringsAsFactors=F)
    geno_mfsd12 = read.delim(f_geno_mfsd12,header=T,sep="\t",skip=5,stringsAsFactors=F)
    geno_ddb1   = read.delim(f_geno_ddb1,header=T,sep="\t",skip=5,stringsAsFactors=F)#; geno_ddb1 = merge(geno_ddb1,geno_mfsd12[,c("Mel_ID","Sample")])
    geno_herc2  = read.delim(f_geno_herc2,header=T,sep="\t",skip=45,stringsAsFactors=F)#; geno_ddb1 = merge(geno_ddb1,geno_mfsd12[,c("Mel_ID","Sample")])
    covariates  = read.delim(f_covariates,header=T,sep="\t")
    ancestry    = read.delim(f_ancestry,header=T,sep="\t",stringsAsFactors=F)
    gencode     = read.delim(f_gencode,header=F,sep="\t",row.names=1,stringsAsFactors=F)
    pca         = read.delim(f_pca,header=T,sep="\t",row.names=1,stringsAsFactors=F)
    
    ## format peculiarities in datasets
    colnames(covariates)=c("Sample","mRNA.sample.ID","microRNA.sample.ID","frozen.by","confluency",
                           "pigment","morphology","final.pass","days.in.culture","Pigmentation.scale",
                           "Pigmentation.category","Growth.rate","Shape")
    
    covariates = mutate(covariates, confluency = factor(confluency, levels=levels(confluency)[c(2,3,4,5,6,7,1)]),
                        pigment = as.numeric(ordered(pigment, levels=levels(pigment)[c(9,2,6,7,4,3,5,1,8)])),
                        final.pass = factor(final.pass, levels=levels(final.pass)[c(2,3,4,5,6,1)]),
                        Pigmentation.category = as.numeric(ordered(Pigmentation.category, levels=levels(Pigmentation.category)[c(2,3,1)])),
                        Growth.rate = factor(Growth.rate, levels=levels(Growth.rate)[c(2,3,1)]))

    pca = data.frame(t(pca))
    pca$Sample = row.names(pca)
    
    ## format genotype matrices
    ## ddb1_vars = grep("^rs",colnames(geno_ddb1),value=T)
    ## mfsd12_vars = grep("^rs",colnames(geno_mfsd12),value=T)
    ## herc2_vars = grep("^rs",colnames(geno_herc2),value=T)

    ## geno_ddb1[geno_ddb1 == "0/0"] = 0
    ## geno_ddb1[geno_ddb1 == "0/1"] = 1
    ## geno_ddb1[geno_ddb1 == "1/0"] = 1
    ## geno_ddb1[geno_ddb1 == "1/1"] = 2
    ## geno_ddb1[geno_ddb1 == "???"] = NA
    ## geno_ddb1[geno_ddb1 == ""] = NA
    
    ## geno_ddb1 = mutate_at(geno_ddb1,starts_with("rs"),as.numeric)
    geno_ddb1.repl = geno_ddb1 %>%
        select(starts_with("rs")) %>%
        mutate_all(funs(replace(.,startsWith(.,"0/0"),0)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"0/1"),1)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"1/0"),1)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"1/1"),2)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"???"),NA))) %>%
        mutate_all(funs(as.numeric)) %>% as.data.frame

    geno_ddb1.repl$Sample = geno_ddb1$Sample

    
    geno_mfsd12.repl = geno_mfsd12 %>%
        select(starts_with("rs")) %>%
        mutate_all(funs(replace(.,startsWith(.,"0/0"),0)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"0/1"),1)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"1/0"),1)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"1/1"),2)))  %>%
        mutate_all(funs(replace(.,startsWith(.,"???"),NA))) %>%
        mutate_all(funs(as.numeric)) %>% as.data.frame

    geno_mfsd12.repl$Sample = geno_mfsd12$Sample
## geno_mfsd12[geno_mfsd12 == "0/0"] = 0
    ## geno_mfsd12[geno_mfsd12 == "0/1"] = 1
    ## geno_mfsd12[geno_mfsd12 == "1/0"] = 1
    ## geno_mfsd12[geno_mfsd12 == "1/1"] = 2
    ## geno_mfsd12[geno_mfsd12 == "???"] = NA
    ## geno_mfsd12[geno_mfsd12 == ""] = NA

    ## geno_mfsd12 = mutate_at(geno_mfsd12,startsWith("rs"),as.numeric)


    geno_herc2.repl = geno_herc2 %>%
        select(starts_with("C")) %>%
        mutate_all(funs(replace(.,startsWith(.,"0/0"),0))) %>%
        mutate_all(funs(replace(.,startsWith(.,"0/1"),1))) %>%
        mutate_all(funs(replace(.,startsWith(.,"1/0"),1))) %>%
        mutate_all(funs(replace(.,startsWith(.,"1/1"),2))) %>%
        mutate_all(funs(as.numeric)) %>%
        t %>% as.data.frame

    
    colnames(geno_herc2.repl) = geno_herc2$ID
    geno_herc2.repl$Sample = row.names(geno_herc2.repl)
    ## geno_ddb1 = convert_geno(geno_ddb1,indiv_ddb1)
    ## geno_mfsd12 = convert_geno(geno_mfsd12,indiv_mfsd12)
    
    ## geno_ddb1.sub = cbind(indiv_ddb1,data.frame(t(geno_ddb1[,indiv_ddb1])))
    ## colnames(geno_ddb1.sub) = gsub(":",".",c("Sample",paste(paste("chr",paste(as.character(geno_ddb1$ID),sep=""),sep=""),".",sep="")))
    
    ## colnames(geno_mfsd12)[colnames(geno_mfsd12)=="Cell.line.ID"] = "Sample"
    ## posit_posits = grep("^X",colnames(geno_mfsd12))
    ## colnames(geno_mfsd12)[posit_posits] = paste(colnames(geno_mfsd12)[posit_posits],".",sep="")

    ## quantile_normalize HERC2 data
    exp_herc2.rank = exp_herc2
    
    
    ## format expression matrices
    exp_mfsd12_indiv = grep("^C",colnames(exp_mfsd12),value=T)
    exp_ddb1_indiv   = grep("^C",colnames(exp_ddb1),value=T)
    exp_herc2_indiv  = grep("^C",colnames(exp_herc2),value=T)

    exp_mfsd12_gene.id = unlist(strsplit(exp_mfsd12$gene_id,"\\."))[c(TRUE,FALSE)]
    exp_ddb1_gene.id   = unlist(strsplit(exp_ddb1$gene_id,"\\."))[c(TRUE,FALSE)]
    exp_herc2_gene.id  = unlist(strsplit(exp_herc2$gene_id,"\\."))[c(TRUE,FALSE)]
    
    exp_mfsd12.sub = cbind(exp_mfsd12_indiv,data.frame(t(exp_mfsd12[,exp_mfsd12_indiv])))
    colnames(exp_mfsd12.sub) = c("Sample",fix_strings(gencode[exp_mfsd12_gene.id,2]))

    exp_ddb1.sub = cbind(exp_ddb1_indiv,data.frame(t(exp_ddb1[,exp_ddb1_indiv])))
    colnames(exp_ddb1.sub) = c("Sample",fix_strings(gencode[exp_ddb1_gene.id,2]))

    exp_herc2.sub = cbind(exp_herc2_indiv,data.frame(t(exp_herc2[,exp_herc2_indiv])))
    colnames(exp_herc2.sub) = c("Sample",fix_strings(gencode[exp_herc2_gene.id,2]))

    
    
    expression = Reduce(function(x,y) merge(x, y, all=TRUE, by="Sample"), list(exp_mfsd12.sub,exp_ddb1.sub,exp_herc2.sub))

    
    ## create master matrix
    master_data.df = Reduce(function(x,y) merge(x, y, all=TRUE, by="Sample"), list(covariates,ancestry,pca,geno_ddb1.repl,geno_mfsd12.repl,geno_herc2.repl,expression,spl_mfsd12))

    ## master_data.df = mutate(master_data.df, Plate. = factor(Plate.))
    
    return(master_data.df)
}


## permute the y's to get empirical p-values
permute_pval = function(data,y,cov,x,n){

    temp.lm = lm(build_lm(y,cov),data)
    y_resid = summary(temp.lm)$residuals
    temp_x = data[-temp.lm$na.action,x]
    temp_df = data.frame(y=y_resid, x=temp_x)
    emp_p = mat.or.vec(n,1)
    for (i in 1:n){
        temp_df$y = sample(y_resid)
        emp_p[i] = summary(lm(y ~ x, temp_df))$coefficients["x",4]
    }

    return(emp_p)
}



## get an empty dataframe with the data specified as well as the given
## column-names and row-names
empty_df = function(data,cols,rows){
    
    df = data.frame(mat.or.vec(length(rows),length(cols)))
    if (length(data) > 0){
        for (i in length(data)){
            df[,i] = data[[i]]
        }
    }
    
    colnames(df) = cols
    row.names(df) = rows

    return(df)
}
