###################################################
## Analysis of variants/mutations
## data from haplotypecaller (gatk) and varscan
## depth >= 10
## include deletion + insertion
## set thresholds to 25% alternate allele
## write overlap and specific genes for each sample
## include strand bias calculcated by Fisher's test (<0.05)
## QUAL: phred-scaled quality score => >30 for all snv in vcf
## filter for missense mutations
###################################################


# vcf
# ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
# 0.05 => ~13
# 0.01 => 20
# the lower p value the more likely a strand bias


dirElem  = strsplit(getwd(),'/')[[1]]
expname  = paste('_',dirElem[length(dirElem)],sep='')
minCount = 10 # minimal depth


## gatk/varscan vcf + vep file
files = list.files("snp", pattern=".vcf$", full.names=T,include.dirs=T)
files = files[grep(".varscan.",files, invert=T)]
files_vs = list.files("snp", pattern=".ohneChr.varscan.vcf$", full.names=T,include.dirs=T)
files_vep = list.files("snp", pattern=".ohneChr.varscan.vcf.txt$", full.names=T,include.dirs=T)
name = gsub(".vcf", "", gsub("snp/","", files))


## all variants of the cell lines in one table
for (i in 1:length(name)) {
	print(files[i])
	## gatk results
	gk = read.csv(files[i],header=T, stringsAsFactors=F, skip=3393, sep="\t")[,-c(3,7)]
	names(gk)[1] = "CHROM"
	gk$GT = sapply(gk[,8],function(x) strsplit(x, ":")[[1]][1])
	gk$AD = gsub(",","/",sapply(gk[,8],function(x) strsplit(x, ":")[[1]][2]))
	gk$DP = as.numeric(sapply(gk[,8],function(x) strsplit(x, ":")[[1]][3]))
	gk$FS = sapply(gk[,6],function(x) strsplit(x, ";FS=")[[1]][2])
	gk$FS = as.numeric(sapply(as.vector(gk$FS),function(x) strsplit(x, ";")[[1]][1]))
	colnames(gk)[(ncol(gk)-3):ncol(gk)] = paste(name[i],c("GT","AD","DP", "FS"), sep="_")
	gk$id = paste(gk$CHROM, gk$POS, gk$REF, gk$ALT, sep="_")
	gk = gk[gk[,grep("_DP$",colnames(gk))]>=minCount,]
	gk$POS = as.numeric(gk$POS)
	gk = gk[,-5:-8]

	## varscan variant vcf data
	vs = read.csv(files_vs[i],header=T, stringsAsFactors=F, skip=23, sep="\t")[,c(1,2,4,5,10)]
	names(vs)[1] = "CHROM"
	vs$id = paste("chr", paste(vs$CHROM, vs$POS, vs$REF, vs$ALT, sep="_"),sep="")
	j=5
	vs$GT = sapply(vs[,j],function(x) strsplit(x, ":")[[1]][1])
	vs$AD = paste(sapply(vs[,j],function(x) strsplit(x, ":")[[1]][5]), sapply(vs[,j],function(x) strsplit(x, ":")[[1]][6]), sep="/")
	vs$DP = as.numeric(sapply(vs[,j],function(x) strsplit(x, ":")[[1]][4]))
	colnames(vs)[(ncol(vs)-2):ncol(vs)] = paste(name[i],c("GT","AD","DP"), sep="_")
	vs = vs[which(vs[,grep("_DP$",colnames(vs))]>=minCount),-j]
	vs$POS = as.numeric(vs$POS)

	## overlap varscan and haplotypecaller
	vs = vs[vs$id%in%gk$id,]
	gk = gk[gk$id%in%vs$id,]

	## load variant effect information (vep)
	vep = read.csv(files_vep[i],header=T, stringsAsFactors=F,skip=56, sep="\t")
	vep[vep=="-"] = NA
	colnames(vep)[1] = c("Uploaded_variation")
	vep$Uploaded_variation = paste(gsub(":", "_",vep$Location), vep$Allele,sep="_")
	vep$Uploaded_variation = paste("chr",vep$Uploaded_variation, sep="")

	## write to file
	gk$id = paste(gk$CHROM,gk$POS,gk$ALT, sep="_")
	gk_vep_h = merge(gk,vep[,-c(2:5)], by.x="id",by.y="Uploaded_variation")
	gk_vep_h = gk_vep_h[unique(c(grep("missense",gk_vep_h$Consequence),grep("stop_gained",gk_vep_h$Consequence),grep("frameshift",gk_vep_h$Consequence))),]

	if (i==1) {
		all = gk_vep_h[,c(1:8)]
		all_vep = vep
	} else {
		all = merge(all, gk_vep_h[c(1:8)], by=c("id","CHROM","POS","REF","ALT"), all=T)
		all_vep = rbind(all_vep,vep[!vep$Uploaded_variation%in% all_vep$Uploaded_variation,])
	}
}


all = merge(all,all_vep, by.x="id",by.y="Uploaded_variation")
write.table(all,paste('tables/mut_hc_vs_vep',expname,'.tsv',sep=''), row.names=F,quote=F,sep='\t', na="")





sessionInfo()
