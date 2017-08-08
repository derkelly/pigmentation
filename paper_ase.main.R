ase_data = read.delim("ase/ASE_sample_data.txt",sep="\t",header=T,stringsAsFactors=F)
ase_data$SYMBOL = fix_strings(ase_data$SYMBOL)

test_vars.all = c(unique(as.character(mfsd12_snp_egenes[,2])),unique(as.character(ddb1_snp_egenes[,2])))

ase_coding = grepl("missense",ase_data$Consequence) | grepl("synonymous",ase_data$Consequence) | grepl("splice",ase_data$Consequence)
ase_utr    = grepl("3_prime",ase_data$Consequence) & !(ase_data$SYMBOL %in% ase_data[ase_coding,]$SYMBOL)

ase_data.cu = ase_data[ase_coding | ase_utr,]


## test ase MFSD12
mfsd12_ase = data.frame(expand.grid(mfsd12_testvars,mfsd12_egenes.mel[,2],stringsAsFactors=F))
colnames(mfsd12_ase) = c("testvar","gene")

mfsd12_ase$p.value = mapply(function(testvar,gene) test_ase(ase_data.cu[ase_data.cu$SYMBOL==gene,],data[,c("Sample",testvar)],testvar,gene),
                            mfsd12_ase$testvar,mfsd12_ase$gene)


## test ase DDB1
ddb1_ase = data.frame(expand.grid(ddb1_testvars,ddb1_egenes.mel[,2],stringsAsFactors=F))
colnames(ddb1_ase) = c("testvar","gene")

ddb1_ase$p.value = mapply(function(testvar,gene) test_ase(ase_data[ase_data$SYMBOL==gene,],data[,c("Sample",testvar)],testvar,gene),
                            ddb1_ase$testvar,ddb1_ase$gene)


## test ase HERC2
herc2_ase = data.frame(expand.grid(herc2_testvars,herc2_egenes.mel[,2],stringsAsFactors=F))
colnames(herc2_ase) = c("testvar","gene")

herc2_ase$p.value = mapply(function(testvar,gene) test_ase(ase_data.cu[ase_data.cu$SYMBOL==gene,],data[,c("Sample",testvar)],testvar,gene),
                            herc2_ase$testvar,herc2_ase$gene)
