## the function max_pwr finds the tag variant with the top power to
## detect alternative splicing for each test_var gene combination
max_pwr = function(df,genotypes,testvar,gene){

    n_hethom = df %>%
    group_by(variantID) %>%  ## group by each tag variant
    summarise(n_het = sum(genotypes[Sample]==1,na.rm=T),     ## calculate the number of het_test het_tag and
              n_hom = sum(genotypes[Sample]!=1,na.rm=T)) %>% ## the number of hom_test het_tag
    group_by(variantID) %>%
        filter(n_het > 1 & n_hom > 1)

    if (nrow(n_hethom)>0){

        n_hethom$pwr = mapply(function(het,hom) pwr.t2n.test(het,hom,d=0.2,sig.level=0.05,alternative="two.sided")[["power"]],
                              n_hethom$n_het,n_hethom$n_hom) ## calculate the power for each variant
        best_pwr = n_hethom[which.max(n_hethom$pwr),]
        return(best_pwr)
        
    } else {
        return(NA)
    }
}



test_ase = function(df,genos,testvar,gene){

    genotypes = genos[,testvar]
    names(genotypes) = genos$Sample

    best_pwr = max_pwr(df,genotypes,testvar,gene)

    if (is.na(best_pwr)){
        return(NA)
    } else{

        het_ae = df %>%
            filter(variantID %in% best_pwr$variantID) %>%
            filter(genotypes[Sample]==1) %>%
            select(AE)

        hom_ae = df %>%
            filter(variantID %in% best_pwr$variantID) %>%
            filter(genotypes[Sample]!=1) %>%
            select(AE)
        return(wilcox.test(het_ae$AE,hom_ae$AE,alternative="two.sided")$p.value)
    }
}
