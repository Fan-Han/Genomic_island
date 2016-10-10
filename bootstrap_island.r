# randomization for genomic islands
rm(list=ls())
library("ggplot2")
library("gridExtra")
# choose cutoff (fst, Dxy or p value)
cutoff = 4
# replicates
k = 10000

#myset <- list.files(pattern="*.p_value")
# sympatry
myset <- c("MG-G_CG-G_Autosome.island.std", "MG-G_DG-G_Autosome.island.std", "PL-Z_STF-Z_Autosome.island.std", "PL-Z_PAL-Z_Autosome.island.std")
dataset = NULL
for (i in myset){
	mydata <- read.table(i, header=T)
	name <- unlist(strsplit(i, "\\."))
	name <- name[1]
	print(name)
	
	# island
	island <- mydata[mydata$ZFST>=cutoff,]
	print(nrow(island))
	island$ID <- rep(name, nrow(island))
	island$Type <- rep("island", nrow(island))
	
	island <- subset(island, select=c(ID, Type, Dxy))
	
	dataset <- rbind(dataset, island)
	
	bk_resample = c()
	# background bootstrapping
	for (n in 1:k){
	resample <- sample(mydata$Dxy, size = nrow(island), replace = FALSE)
	ave_resample = round(mean(resample, na.rm=T), digits=4)
	
	bk_row <- c(name, "background", ave_resample)
	dataset <- rbind(dataset, bk_row)
	
	bk_resample <- c(bk_resample, ave_resample)
	}
	if(nrow(island) != 0){
	bootstrap_P = sum(bk_resample == round(mean(island$Dxy, na.rm=T), digits=4))/k
	print(bootstrap_P)
	}
}

dataset.df <- data.frame(ID=factor(dataset$ID), Type=factor(dataset$Type),Dxy=as.numeric(dataset$Dxy)) 
g_sym <- ggplot(data=dataset.df, aes(x=ID, y=Dxy, fill=Type)) + geom_boxplot(outlier.shape=NA)+ coord_flip() + scale_x_discrete(limits=c("PL-Z_PAL-Z_Autosome","PL-Z_STF-Z_Autosome","MG-G_DG-G_Autosome","MG-G_CG-G_Autosome")) + theme_bw() + scale_fill_manual(values=c("white","grey")) + xlab("") + ylab("")+ theme(legend.position="none", axis.text.y = element_blank(), strip.text.y=element_blank()) + scale_y_continuous(expand=c(0,0), limits=c(0,0.0075))

# allopatry galapagos
rm(dataset,dataset.df,mydata)
myset <- c("MTF-F_LTF-P_Autosome.island.std", "SGF-S_STF-Z_Autosome.island.std", "MG-M_STF-Z_Autosome.island.std", "MG-M_PAL-Z_Autosome.island.std")
dataset = NULL
for (i in myset){
	mydata <- read.table(i, header=T)
	name <- unlist(strsplit(i, "\\."))
	name <- name[1]
	print(name)
	
	# island
	island <- mydata[mydata$ZFST>=cutoff,]
	print(nrow(island))
	island$ID <- rep(name, nrow(island))
	island$Type <- rep("island", nrow(island))
	
	island <- subset(island, select=c(ID, Type, Dxy))
	
	dataset <- rbind(dataset, island)
	
	bk_resample = c()
	# background bootstrapping
	for (n in 1:k){
	resample <- sample(mydata$Dxy, size = nrow(island), replace = FALSE)
	ave_resample = round(mean(resample, na.rm=T), digits=4)
	
	bk_row <- c(name, "background", ave_resample)
	dataset <- rbind(dataset, bk_row)
	
	bk_resample <- c(bk_resample, ave_resample)
	}
	if(nrow(island) != 0){
	bootstrap_P = sum(bk_resample == round(mean(island$Dxy, na.rm=T), digits=4))/k
	print(bootstrap_P)
	}
}

dataset.df <- data.frame(ID=factor(dataset$ID), Type=factor(dataset$Type),Dxy=as.numeric(dataset$Dxy)) 
g_allog <- ggplot(data=dataset.df, aes(x=ID, y=Dxy, fill=Type)) + geom_boxplot(outlier.shape=NA)+ coord_flip() + scale_x_discrete(limits=c("MG-M_PAL-Z_Autosome","MG-M_STF-Z_Autosome","SGF-S_STF-Z_Autosome","MTF-F_LTF-P_Autosome")) + theme_bw() + scale_fill_manual(values=c("white","grey")) + xlab("") + ylab("")+ theme(legend.position="none", axis.text.y = element_blank(), strip.text.y=element_blank()) + scale_y_continuous(expand=c(0,0), limits=c(0,0.0075))


#allopatry cocos
rm(dataset,dataset.df,mydata)
myset <- c("PIN_CG-G_Autosome.island.std", "PIN_MG-M_Autosome.island.std", "PIN_DP-P_Autosome.island.std")
dataset = NULL
for (i in myset){
	mydata <- read.table(i, header=T)
	name <- unlist(strsplit(i, "\\."))
	name <- name[1]
	print(name)
	
	# island
	island <- mydata[mydata$ZFST>=cutoff,]
	print(nrow(island))
	island$ID <- rep(name, nrow(island))
	island$Type <- rep("island", nrow(island))
	
	island <- subset(island, select=c(ID, Type, Dxy))
	
	dataset <- rbind(dataset, island)
	
	bk_resample = c()
	# background bootstrapping
	for (n in 1:k){
	resample <- sample(mydata$Dxy, size = nrow(island), replace = FALSE)
	ave_resample =round(mean(resample, na.rm=T), digits=4)
	
	bk_row <- c(name, "background", ave_resample)
	dataset <- rbind(dataset, bk_row)
	
	bk_resample <- c(bk_resample, ave_resample)
	}
	if(nrow(island) != 0){
	bootstrap_P = sum(bk_resample == round(mean(island$Dxy, na.rm=T), digits=4))/k
	print(bootstrap_P)
	}
}

dataset.df <- data.frame(ID=factor(dataset$ID), Type=factor(dataset$Type),Dxy=as.numeric(dataset$Dxy)) 
g_alloc <- ggplot(data=dataset.df, aes(x=ID, y=Dxy, fill=Type)) + geom_boxplot(outlier.shape=NA)+ coord_flip() + scale_x_discrete(limits=c("PIN_DP-P_Autosome","PIN_MG-M_Autosome","PIN_CG-G_Autosome")) + theme_bw() + scale_fill_manual(values=c("white","grey")) + xlab("") + ylab("") + theme(legend.position="none", axis.text.y = element_blank(), strip.text.y=element_blank()) + scale_y_continuous(expand=c(0,0), limits=c(0,0.0075))


# plot together
pdf("Figure3.pdf",onefile=TRUE)
grid.arrange(g_sym,g_allog,g_alloc, ncol=1)
dev.off()
