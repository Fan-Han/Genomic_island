rm(list=ls())
library(ggplot2)
library(plotrix)
library(boot)


rm(list=ls())
k=10000


#alldata <- data.frame(type = factor(), rate = as.numeric(), se = as.numeric())
options(stringsAsFactors = FALSE)
mydata <- read.table("50kb_window.rate", header=F)
Bmean <- function(data, indices) {
    d <- data[indices] # allows boot to select sample 
    return(mean(d))
}
test <- mydata$V6[!is.na(mydata$V6)]
results <- boot(data=test, statistic=Bmean, R=k)
ci <- boot.ci(results, type="perc")


alldata <- data.frame(type = as.character("background"), rate = as.numeric(mean(mydata$V6,na.rm=TRUE)), lower = ci[[4]][4], upper=ci[[4]][5] )
# ALX1 and HMGA2 regions
ALX1 <- mydata[mydata$V4=="JH739921"&mydata$V5>=320000&(mydata$V5+50000)<=560000,]
ALX1.df <- data.frame(type=as.character("ALX1"), rate = as.numeric(mean(ALX1$V6,na.rm=TRUE)), lower = 0, upper=0)
alldata <- rbind(alldata, ALX1.df)
HMGA2 <- mydata[mydata$V4=="JH739900"&mydata$V5>=6945000&(mydata$V5+50000)<=7470000,]
HMGA2.df <- data.frame(type=as.character("HMGA2"), rate = as.numeric(mean(HMGA2$V6,na.rm=TRUE)), lower = 0, upper=0)
alldata <- rbind(alldata, HMGA2.df)

rm(results, ci, test, mydata)
file_path="/Project/b2012111_DarwinsFinch/SpeciationIslands/20160923_recalculation/Correlation_dxy_rr_fst_island/"
files <- list.files(path=file_path, pattern="\\.island")


for (i in files){
	full_path <- paste(file_path, i, sep="")
	mydata <- read.table(full_path, header=F)
	file <- unlist(strsplit(i, "\\."))
	name <- file[1]
	print(name)

	Bmean <- function(data, indices) {
    d <- data[indices] # allows boot to select sample 
    return(mean(d))
}
	test <- mydata$V17[!is.na(mydata$V17)]
	results <- boot(data=test, statistic=Bmean, R=k)
	ci <- boot.ci(results, type="perc")

	newline <- c(name, mean(mydata$V17,na.rm=TRUE), ci[[4]][4],ci[[4]][5])

	alldata <- rbind(alldata, newline)
}



alldata <- data.frame(type = factor(alldata$type), rate = as.numeric(alldata$rate), lower = as.numeric(alldata$lower), upper=as.numeric(alldata$upper))

alldata$type <- factor(alldata$type, c("background","MG-G_CG-G","MG-G_DG-G","PL-Z_STF-Z","PL-Z_PAL-Z","MTF-F_LTF-P","SGF-S_STF-Z","MG-M_STF-Z","MG-M_PAL-Z","PIN_CG-G","PIN_MG-M","PIN_DP-P","ALX1","HMGA2"))
limits <- aes(ymax = upper, ymin = lower)
pdf("island_background_RR_bootstrap.pdf")
g <- ggplot(data=alldata, aes(x=type,y=rate)) + geom_point(size=1) +geom_errorbar(aes(ymin=lower, ymax=upper), width=0.25) + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ scale_y_continuous(expand = c(0,0))
print(g)
dev.off()
