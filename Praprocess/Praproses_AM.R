library("scrime")
missing_symbol <- function(matriks)
{
	jenis = matriks[,1]
	geno = matriks[,2:1537]
	feno = matriks[,1538:1549]
	for(i in 1:dim(geno)[1])
	{
		id_m = which(geno[i,] == "--")
		if(length(id_m) > 0)
		{
			geno[i,id_m] = NA
		}
	}
	return(list(jenis, geno, feno))
}

missing_call <- function(genotipo, th)#th bisa 5%, 0.05
{
	id_baris = 0
	id_kolom = 0
	baris_buang = c()
	kolom_buang = c()
	ukuran = dim(genotipo)
	for(baris in 1:ukuran[1])
	{
		tot_mis = length(which(is.na(genotipo[baris,])))
		if((tot_mis/ukuran[2]) > th)
		{
			id_baris = id_baris+1
			baris_buang[id_baris] = baris
		}
	}
	genotipo2 = genotipo[-baris_buang,]
	ukuran = dim(genotipo2)
	if(ukuran[1] < 1)
	{
		print("Tidak ada sampel yg lulus QC, toleransi maf terlalu kecil")
		break;
	}
	for(kolom in 1:ukuran[2])
	{
		tot_mis = length(which(is.na(genotipo2[,kolom])))
		if((tot_mis/ukuran[1]) > th)
		{
			id_kolom = id_kolom+1
			kolom_buang[id_kolom] = kolom
		}
	}
	genotipo3 = genotipo2[,-kolom_buang]
	return(list(genotipo3, baris_buang, kolom_buang))
}

mcr = 0.05
maf = 0.05 #dapet 81 1072

print("Uploading data")
padi <- read.csv("dataFull.csv", header=FALSE, sep=",")
posisi <- read.csv("posisi.csv", header=FALSE, sep=",")
colsnp <- paste("SNP", 1:1536, sep="")
colfeno <- paste("Feno", 1:12, sep="")
colnames(padi) <- c("Jenis", colsnp, colfeno)

print("Changing missing symbol")
padi_ms <- missing_symbol(padi)
padi.jenis <- padi_ms[[1]]
padi.genotipo <- padi_ms[[2]]
padi.fenotipo <- padi_ms[[3]]

print("Arranging SNPs base on position")
padi.arrg1 <- cbind(posisi, t(padi.genotipo))
padi.arrg1 <- padi.arrg1[order(padi.arrg1[,1], padi.arrg1[,2]),]
padi.arrg <- t(padi.arrg1[,3:dim(padi.arrg1)[2]])

print("Checking missing call rate")
padi_mc <- missing_call(padi.arrg, mcr)
padi.genotipo <- padi_mc[[1]]
padi.fenotipo <- padi.fenotipo[-padi_mc[[2]],]
padi.jenis <- padi.jenis[-padi_mc[[2]]]

print("Recode SNP")
padi.genotipo <- recodeSNPs(t(padi.genotipo), first.ref = FALSE, geno = 1:3, snp.in.col = FALSE)
padi.genotipo <- t(padi.genotipo)

print("Checking minor allel frequency")
maf.genotipo <- rowMAFs(t(padi.genotipo), check = TRUE) #cek maf
id.maf.2 <- which(maf.genotipo > maf)
padi.genotipo <- padi.genotipo[,id.maf.2]

padi.genotipo.qc <- padi.genotipo
#padi.genotipo.qc[is.na(padi.genotipo.qc)] <- 0

#padi.genotipo.qc <- matrix(unlist(padi.genotipo.qc), ncol = dim(padi.genotipo.qc)[2], byrow = TRUE)
#padi.fenotipo <- matrix(unlist(padi.fenotipo), ncol = dim(padi.fenotipo)[2], byrow = TRUE)
padi.genotipo.qc <- matrix(unlist(padi.genotipo.qc), ncol = dim(padi.genotipo.qc)[2])
padi.fenotipo <- matrix(unlist(padi.fenotipo), ncol = dim(padi.fenotipo)[2])
print("Praproses is finished") #, file=padi.genotipo.qc")

#simpan jenis padi sebagai matrix
jenis_padi <- t(padi.jenis)
jenis_padi <- matrix(unlist(jenis_padi), ncol = dim(jenis_padi)[2])

#simpan data full
dataFullQC <- cbind(jenis_padi, padi.genotipo.qc, padi.fenotipo)

#export as csv file
write.table(dataFullQC, file = "dataFullQC", sep = ",", col.names = NA, qmethod = "double")
write.table(padi.genotipo.qc, file = "dataGenoQC - 005 005.csv", sep = ",", col.names = NA, qmethod = "double")
write.tabke(padi.fenotipo, file = "dataGenoQC", sep = ",", col.names = NA, qmethod = "double")
