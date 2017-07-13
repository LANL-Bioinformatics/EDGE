#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
print(args)
cov_table_file <-args[1]
classified_result_table<-args[2]
out_prefix<-args[3]
top<-30  # only plot top30 if more than 30 unique taxas

# check input
if ( !file.exists(classified_result_table) ){
	cat( " ./classification_plot.R  readsToContigs_coverage.table ctg_class.csv out_prefix \n")
	quit(save="no")
}

# function to cap first letter 
.simpleCap <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1, 1)), substring(s, 2),
	sep = "", collapse = " ")
}

pdf(file=paste(out_prefix,".contigsClassification.pdf",sep=""),width=10,height=8)
def.par <- par(no.readonly = TRUE) # default 
if (file.exists(cov_table_file)){
	cov_table<-read.table(file=cov_table_file,header=TRUE)
}
taxa_table<-read.table(file=classified_result_table,header=TRUE,sep="\t",comment.char="",quote="")

if (length(which(taxa_table$RANK=="RANK"))>0)
{
	# Remove extra header
	taxa_table<-taxa_table[-c(which(taxa_table$RANK=="RANK")),]

	# re factor of the column
	taxa_table$ACC_COV_LEN<-factor(taxa_table$ACC_COV_LEN)

	# Factor to numeric conversion
	taxa_table$ACC_COV_LEN<-as.numeric(levels(taxa_table$ACC_COV_LEN))[taxa_table$ACC_COV_LEN]
}

if (file.exists(cov_table_file)){
	# Remove 0 coverage contigs
	cov_0_above_table<-subset(cov_table,cov_table$Base_Coverage.>0)
	#scale<-sqrt(a$Length/sd(a$Length))
	#scale[scale<0.2]=0.2
	scale<-log10(cov_0_above_table$Length) - 1.5
}
ranklist<-c("superkingdom","phylum","class","order","family","genus","species","strain")
#b$RANK<-factor(b$RANK,levels=list)
#c<-table(b$ORGANISM, b$RANK)
#color<-rainbow(nrow(c))
#barplot(c/colSums(c),col=color,border=color,xaxt='n')
#text(cex=1, x=x-.25, y=-0.1, labs, xpd=TRUE, srt=45)

margin_adj<-0
unique_col<-sample(rainbow(top))
#output summary for top 5
write(c("RANK","Top1","Top2","Top3","Top4","Top5"),file="summary_by_hitAccLength.txt",sep="\t",ncolumns=6)
write(c("RANK","Top1","Top2","Top3","Top4","Top5"),file="summary_by_topHitCount.txt",sep="\t",ncolumns=6)
for ( i in 1:7) {
	unique_c<-unique_col
	# adjust species level margin space for long species name
	if(i>6){margin_adj<-2}
  	# select subset by ran
	data.rank <- subset(taxa_table, (RANK == ranklist[i] | RANK == "unclassified") & ACC_COV_LEN > 0)
  	
  	# sum acc_cov_len by organism group
  	bp_plot_table<-aggregate(ACC_COV_LEN ~ ORGANISM, data=data.rank, FUN=sum)
  	# sort by descreasing order 
  	bp_plot_table<-bp_plot_table[order(bp_plot_table$ACC_COV_LEN,decreasing = TRUE),]
  	
  	# get contig assignemnt by top hit 
	data.rank.forTopHit<-subset(data.rank, (RANK == ranklist[i] | (RANK == "unclassified" & LENGTH == ACC_COV_LEN)))
  	top_hit_table<-merge(aggregate(ACC_COV_LEN ~ X..SEQ, data=data.rank.forTopHit, FUN=max),data.rank.forTopHit)
 
	if (file.exists(cov_table_file)){
  		# get depth and GC info by merging table with contig id 
  		top_hit_table_with_cov<-merge(cov_0_above_table,top_hit_table,by.x="ID",by.y="X..SEQ")
  	
  		# re factor of the subset
  		top_hit_table_with_cov$ORGANISM<-factor(top_hit_table_with_cov$ORGANISM)

  		#count by taxanomy assignments
  		count<-table(top_hit_table_with_cov$ORGANISM)
  		# unique taxa number
  		total_tax<-length(count)
	}else{
  		#count by taxanomy assignments
  		count<-table(top_hit_table$ORGANISM)
  		# unique taxa number
  		total_tax<-length(count)
	}
  	
  	# Only use top 30 to plot for the proper resolution
  	if (total_tax>top){
  		for_barplot_count<-sort(count,decreasing=TRUE)[1:top]
  		bp_plot_table<-head(bp_plot_table,n=top)
  	}else{
  		for_barplot_count<-sort(count,decreasing=TRUE)
  	}
  	
  	unique_org<-names(for_barplot_count)
  	# assign colors for each taxa
  	# unique_c<-rainbow(length(unique_org))
	if (file.exists(cov_table_file)){
  		top_hit_table_with_cov$COLOR<-rep("darkgrey",length(top_hit_table_with_cov$ORGANISM))
	}
  	bp_plot_table$COLOR<-rep("darkgrey",length(bp_plot_table$ORGANISM))
  	for (j in 1:length(unique_org))
  	{
		if (unique_org[j]=="unclassified"){unique_c[j]<-"lightgrey";}
		if (file.exists(cov_table_file)){
  			top_hit_table_with_cov$COLOR[which(top_hit_table_with_cov$ORGANISM==unique_org[j])]<-unique_c[j] 
		}
  		bp_plot_table$COLOR[which(bp_plot_table$ORGANISM==unique_org[j])]<-unique_c[j] 
  	}
  	
	 ## color not used in bp_plot_table, the different taxa top list between count and length
  	diff_col<-unique_c[!unique_c %in% bp_plot_table$COLOR]
  	if (length(which(bp_plot_table$COLOR=='darkgrey'))>0){
       		bp_plot_table$COLOR[which(bp_plot_table$COLOR=='darkgrey')]<-diff_col[1:length(which(bp_plot_table$COLOR=='darkgrey'))]
  	}
  
 
  	
  	## barplot # contigs vs organism
  	par(mar=c(10+margin_adj,6,4,2))
  	
  	## barplot bp contigs vs organism
  	barplot(bp_plot_table$ACC_COV_LEN,names.arg=bp_plot_table$ORGANISM,col=bp_plot_table$COLOR,las=3,cex.names=0.8,ylab="Length Of Contigs (bp)")
  	if (total_tax>top){
  	 	mtext(paste("Rank:", .simpleCap(ranklist[i]),"(",total_tax,")","Only Show Top",top),3,adj=0,cex=0.8,line = 1)
  	}else{
  		mtext(paste("Rank:", .simpleCap(ranklist[i]),"(",total_tax,")"),3,adj=0,cex=0.8,line = 1)
  	}
  	
  	barplot(for_barplot_count,col=unique_c,ylab="Number Of Contigs",las=3,cex.names=0.8)
  	if (total_tax>top){
  		mtext(paste("Rank:", .simpleCap(ranklist[i]),"(",total_tax,")","Only Show Top",top),3,adj=0,cex=0.8,line = 1)
  	}else{
  		mtext(paste("Rank:", .simpleCap(ranklist[i]),"(",total_tax,")"),3,adj=0,cex=0.8,line = 1)
  	}
  	#labs <- colnames(count)
  	#text(cex=1, x=x-.25, y=-0.1, labs, xpd=TRUE, srt=45)
  	
	if (file.exists(cov_table_file)){
  		## scatter plot: Depth Coverage vs GC. dot size propotion to (log) contig size 
  		par(mar=c(5,6,4,11+margin_adj))
  		par(xpd=FALSE)
  		plot(top_hit_table_with_cov$Avg_fold,top_hit_table_with_cov$GC.,ylab="GC (%)",xlab="Average Coverage Fold (x)",main="Contig Average fold Coverage vs. %GC",pch=20,col=top_hit_table_with_cov$COLOR,cex=scale,log="x",bty="l")
  		grid(col="grey")
  		mtext(paste("Rank:", .simpleCap(ranklist[i]),"(",total_tax,")"),3,adj=0,cex=0.8)
  		par(xpd=TRUE)
  		coord<-par("usr")
  		if (total_tax>top){
  			legend(10^coord[2],coord[4],legend=c(unique_org,paste("Not Top",top)),col=c(unique_c,"darkgrey"),bty='n',pch=20,cex=0.8)
  		}else{
  			legend(10^coord[2],coord[4],legend=unique_org,col=unique_c,bty='n',pch=20,cex=0.8)
  		}
	}
	write(c(ranklist[i],as.character(bp_plot_table$ORGANISM[1:5])),file="summary_by_hitAccLength.txt",sep="\t",ncolumns=6,append=TRUE)
  	write(c(ranklist[i],unique_org[1:5]),file="summary_by_topHitCount.txt",sep="\t",ncolumns=6,append=TRUE)
}

par(def.par)
proc.time()

tmp<-dev.off()
quit(save="no")
