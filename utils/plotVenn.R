#usage R --vanilla < plotVenn.R --args file.txt name svtype min

cmd_args <- commandArgs(trailingOnly = TRUE)


img<-function(file,name,type,min){

require(gplots)

pdfName<-paste(c(name, type, format(Sys.time(), "%a%b%d_%H_%M_%S.pdf")), collapse="_")

dat<-read.table(file, header=TRUE)

dat<-dat[abs(dat$length) > min,]

cols<-1:(length(dat) -2 )

fdat<-dat[dat$type == type,cols]

pdf(pdfName)
venn(fdat)
dev.off()

}

img(cmd_args[1], cmd_args[2], cmd_args[3], cmd_args[4])