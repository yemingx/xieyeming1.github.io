################## read and write ##################
df<-read.table('gene_expression.xls', sep='\t',header = TRUE)
df<-fread(paste0("PA_str_stack.xls.gz"),header=TRUE)
write.table(new_df,file = 'TAIR_gene_expression.xls',quote = FALSE,sep="\t",row.names = FALSE)

library(data.table)
args = commandArgs(trailingOnly=TRUE)
in_path=args[1]
out_path_to_file=args[2]
HiC<- read.table(in_path,sep='')
HiC_merge<- merge(HiC,HiC,by='V2',allow.cartesian = TRUE)
dt2 <- CJ(bin1 = unique(dt1$bin1), bin2 = unique(dt1$bin2))

write.table(HiC_merge, file = out_path_to_file,sep = "\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

################# sample batch operation #################
sample_list<-read.table('../sample.list')
# V1 raw sample id, V2 sample name for plotting
suffix<-".mcool"
file_name<-paste(sample_list$V1,suffix , sep="")
hicf_raw<-as.list(file_name)
names(hicf_raw)<-as.list(sample_list$V2)
hicf<-hicf_raw[1:4]

gsub('.mcool','',fileName(hic_compts))

for (x in 1:10) {
  print(x)
}
################# plot template #################
library(egg)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

+theme_presentation() +scale_fill_manual(values=c(mypal[1],mypal[2],mypal[3],mypal[4]))

# histogram 2 sample
ggplot(df,aes(x=order,color=sample))+geom_histogram(fill="white", alpha=0.5, position="identity",binwidth=1)+
  theme_classic()+xlim(0,30)+ theme(legend.position="top")

p3<-ggplot(target_sort_melt, aes(relative_pos_to_site, site, fill= m_ratio)) + geom_tile() + scale_fill_gradient(low="#F7F7F7", high=color3) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
                     axis.text=element_text(),
                     axis.title=element_text(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks.y=element_blank()) +scale_color_manual(values=color3) + 
  theme(legend.position="top") 
grid.arrange(p1, p2, p3, nrow = 1,top=(paste0(label,'    chr7, ',bin,' bp bin',", depth: ",depth)), widths=c(2,1,1)) 

# dot plot multi sample
library(ggsci)
pdf('test_xiang.pdf',width=10,height=5)
ggplot(cafl3, aes(x=time, y=conc, color=Rep))+geom_point(size=0.2)+facet_wrap(~Group)+
  xlab('Time(s)') + ylab('[Ca2+]cyt')+ theme_bw()+scale_color_jco()
dev.off()

# flip y and scientific notation y axis
ggplot(stack_df, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line()+ geom_point()+
  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd),
                   alpha=0.2, size=0.5)+
  #scale_x_continuous(breaks = c(1e5,1e6,1e7))+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4]))+xlim(0,13000000)+
 scale_y_continuous(trans = "reverse", labels = scales::scientific_format())+
  labs(subtitle='gradient: 25000 50000 100000 200000 400000 
       800000 1600000 3200000 6400000 12800000\nnagano2017 515 cells \ndipC 1071 cells\nhichew 1178 cells')

########### cor hex plot ###########
ggplot(sca_hic, aes(x=sca_hic$E1.x,y=sca_hic$E1.y)) + geom_hex(bins = 100)+ geom_abline(intercept = 0, slope = 1) +#geom_smooth(method=lm)+ 
  scico::scale_fill_scico(palette = "hawaii")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2) + ylim(-2,2) +
  labs(title="eigenvector, x:hic, y:sca_40M, hg19, bin 100k, cor 0.88")

ggplot(eigen_merge, aes(x=eigen_merge$E1.x,y=eigen_merge$E1.y)) + geom_hex(bins = 50)+ geom_abline(intercept = 0,alpha=0.2, slope = 1,linetype = "dashed") +#geom_smooth(method=lm)+ 
  scico::scale_fill_scico(palette = "lajolla",direction=-1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2) + ylim(-2,2) +
  labs(title=paste0("eigen    x:",sample_1,", \ny:",sample_2,", hg19, bin, \npearson: ",p_cor," spearman:",s_cor)) + 
  theme(legend.direction = "horizontal",legend.position = "top", legend.box = "top")

########### heatmap ###########
Heatmap(dat_mtx_pause_active, name = "log2_gro_signal", 
       cluster_columns = FALSE,col=colorRamp2(c(0,2,4), c("darkgreen", "white", "darkred")),
       show_row_names = FALSE,show_column_names = FALSE,cluster_rows = FALSE,
       column_title=paste0('HEK GRO signal From tss to tes. row_num: ',length(rownames(dat_mtx_pause_active))))

########### violin plot ###########
library(ggplot2)
library(ggsci)
library(httpgd)
library(dplyr)
library(reshape2)
npg_colors <- pal_npg()(10)

# combine all data
all_data <- rbind(dnase_non_overlap, dnase_overlap)
all_data$log2_macs_qvalue <- ifelse(all_data$V9 > 0, log2(all_data$V9), NA)

wilcox_result <- wilcox.test(dnase_overlap$V9, dnase_non_overlap$V9, 
                           alternative = "two.sided")
head(all_data)
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(log2_macs_qvalue, na.rm = TRUE),2),
    max_val = round(max(log2_macs_qvalue, na.rm = TRUE),2),
    sd = sd(log2_macs_qvalue, na.rm = TRUE))
summ
levels(factor(all_data$subgroup))
ggplot(all_data, aes(x=subgroup, y=log2_macs_qvalue, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.3, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(max_val, na.rm = TRUE)), 
            data = summ, size=4, vjust = 2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("lightblue","indianred")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$log2_macs_qvalue, na.rm = TRUE) * 1.2) +
  labs(title = 'DNAse Accessibility Comparison',
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = "log2(MACS q-value)") +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")

########### line dot plot ###########


########### line compare plot ###########
p3<-ggplot(eigen_line,aes(x=locus_id,y=E1,group=sample,color=sample)) + geom_line(stat="identity") +
  scale_color_manual(values = c('#E64B35B2','#4DBBD5B2')) +theme_bw()+ geom_hline(yintercept = 0,linetype="dashed")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ theme(legend.position="top")
pdf(paste0('eigen_',sample_1,'_',sample_2,'_',eigen_bin,'.pdf'),width=5,height=8)
grid.arrange(p3,p1,nrow=2)
cowplot::plot_grid(p1,p3,nrow=2)
dev.off()


################# loop #################
raw_target<-fread(paste0("PA_str_stack.xls.gz"),header=TRUE)
sample_group<-c('PA','2ab-PA-user','2ab-user','PA-user')
metrics<-data.frame()
plot_list <- list()

for(i in 1:length(sample_group)){
label=sample_group[i]
p<-ggplot(singleton_adaptor_raw, aes(x=merge_time))+ 
    geom_histogram(fill="white",col='black', alpha=0.5, position="identity",binwidth=1)+ 
    labs(title=label,x="merge time", y = "count")+ theme_classic()
plot_list[[i]]<-p
metrics<-rbind.data.frame(metrics,cbind.data.frame(sample=label,aligned_read=aligned_read_num)}

pdf('singleton_adaptor_merge_time.pdf',height = 5,width = 5)
cowplot::plot_grid(plotlist = plot_list, nrow = 2)
library(gridExtra)
grid.arrange(grobs=plot_list,nrow=2,ncol=2)
dev.off()
write.table(metrics,file = 'adaptor_metrics.xls',quote = FALSE,sep="\t",row.names = FALSE)

################## fill na ##################
library(zoo)
#row wise fill
symbol_to_TAIR[] <- t(apply(your_data_frame, 1, zoo::na.locf))

#column wise fill
na.locf(na.locf(your_data_frame), fromLast = TRUE)
hic$E1[is.na(hic$E1)]<-0
# remove df NA row
atac_vs_sca<-atac_vs_sca_raw[-which(is.na(atac_vs_sca_raw$GpC_m_ratio)),]
final[complete.cases(final), ]

################## ggplot2 ##################
ggplot(dat) + theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5))

myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
g+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +geom_vline()
abline(v = 16, col = "darkgreen")

pdf(paste0(attribute_name,'_methty_ratio_global_',range,'.pdf'))
attribute_raw_colsum<-rbind.data.frame(cbind.data.frame(attribute_raw_colsum_GpC,sample='GpC'),cbind.data.frame(attribute_raw_colsum_CpG,sample='CpG'))
ggplot(data=attribute_raw_colsum, aes(x=relative_pos_to_site, y=site_methy_ratio, group=sample)) +geom_line(aes(color=sample)) + geom_vline(xintercept=c(-10,10),linetype="dotted") +
theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                   axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +scale_color_manual(values=c("olivedrab","royalblue"))+
  labs(title=paste0(attribute_name," chr7    ",range," range    bin_size = 10"))
dev.off()

gradient: color=colorRampPalette(c("#F7F7F7","olivedrab"))(50)

ggplot(df, aes(x=weight, color=sex)) +
  geom_density()
# Add mean lines
p<-ggplot(df, aes(x=weight, color=sex)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")
################## plot ##################
#set index levels
df_norm$sample<-factor(df_norm$sample, levels = c('poreC_only_gpc','nome_15min_gpc','nome_30min_gpc','nome_3h_gpc'))

# Sort the data frame by the new factor levels
data_sorted <- data[order(data$Category), ]

#pdf('chr7_read_id_order.pdf')
#plot(density(nome75_read_id_order$V1), main=paste0('chr7_read_id_order'),col='red',xlim=c(0,20))
#lines(density(nome15_read_id_order$V1),col="blue",xlim=c(0,20))
#legend("topright",c("nome75","nome15"), col=c("red","blue"), cex=1.1,pch=15,pt.cex = 0.5)
#dev.off()
DFplotlong <- read.table(text='subject iq      condition RT
1       98      A         312
1       98      B         354
1       98      C         432
2       102     A         134
2       102     B         542
2       102     C         621', header=TRUE)
#
ggplot(DFplotlong, aes(iq, RT, colour = condition, linetype = condition)) +
  geom_point() +
  geom_smooth(method = lm, fullrange = TRUE, alpha = .15) +
  theme_bw() +
  labs(x = "iq", y = "reaction times") +
  scale_colour_manual(values=c("#999999","#000000", "#900009"),
                     name="condition", 
                     breaks=c("A", "B", "C"), 
                     labels = c("easy", "medium", "hard")) +
  scale_linetype_discrete(name="condition", 
                          breaks=c("A", "B", "C"), 
                          labels = c("easy", "medium", "hard"))

# line-dot single molecule graph
ggplot(seg_bin_GpC_ratio_order_2_sort,aes(x=pos,y=read_id)) + geom_point(alpha=0.9,aes(size=GpC_ratio,shape = order,color = GpC_ratio)) +
  scale_size_continuous(range = c(2,4))+scale_color_continuous(c(-0.2,0.2))+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "lightgrey"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_color_gradient2(low = "grey30", mid = "white", high = "red")+
  geom_line(stat="smooth",method = "lm",
            alpha = 0.1)+
  labs(title="SCA_WT, ")

# more than 6 shapes
gp <- ggplot(df,aes(x=t, y=y, group=sn,color=sn, shape=sn)) +
             scale_shape_manual(values=1:nlevels(df$sn)) +
             labs(title = "Demo more than 6 shapes", x="Theat (deg)", y="Magnitude") +
             geom_line() + 
             geom_point(size=3)

##################aesthetic##############
# rotate x axis label
ggplot(data = valid_report_ratio, aes(y = V1,x = sample,group=cell_lib,color=sample)) + geom_line() +
  geom_point()+labs(title="valid / report ratio")+ theme_classic()+
  scale_color_manual(values=c('darksalmon','brown3','darkred','cyan3','darkcyan','lightblue4'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################## stats ###################
quantile(x,c(0.05))
ecdf(x)(2)
################## matrix ###################
# cartesian join
HiC_merge<- merge(HiC,HiC,by=’V2’,allow.cartesian = TRUE)

# contact table to matrix
dat_mtx<-dcast(dat_dedup_fill,V2~V3,value.var = "V1.x",fun.aggregate = sum)

# matrix to contact table
dat_mtx_df2_mtx_ICE_melt<-melt(as.matrix(dat_mtx_df2_mtx_ICE))

# remove lower triangle signal
dat_mtx[lower.tri(dat_mtx)] <- NA
################## dataframe ###################
newdf<-rbind(df, data.frame(hello="hola", goodbye="ciao"))

# replace value
df$Marks[df$Names == "Sita"] <- 25

# subset dataframe
dt[dt$fct %in% vc,]
dt[!(dt$fct %in% vc),]
subset_df <- merge_df[grepl("2i", merge_df[,2]), ]

# check column value frequency
count(tss_10x, 'subtype')

# Left join
df2 <- merge(x=emp_df,y=dept_df, by="dept_id", all.x=TRUE)
################## dedup ###################
df = df[!duplicated(df$Date),]
unique(v)
v[duplicated(v)]
seq_bin_hic_dedup<-seq_bin_hic %>% group_by(pos.x,pos.y) %>% summarize(Total=sum(count))

################## datatable ###################
# looping
m_seg_nome3h_GpC[,'p_val' := mapply(function(x,y) binom.test(x=x,n=y,p=0.1,alternative='greater')$p.value,GpC_methy,GpC_all)]

# merge range
library(data.table)
b<- fread(
  "chr, Start, End
  chr7, 0, 15
  chr7, 19, 200
  chr7, 200, 300"
)

a<- fread(
  "chr, Start, End
  chr7, 0, 10
  chr7, 10, 20
  chr7, 20, 30"
)
setkey(a, chr,Start, End)
> foverlaps(b, a)
    chr Start End i.Start i.End
1: chr7     0  10       0    15
2: chr7    10  20       0    15
3: chr7    10  20      19   200
4: chr7    20  30      19   200
5: chr7    NA  NA     200   300

################## gene name convert ##################
library(org.At.tair.db)
genes<-df$gene_symbol
#genes <- c("AT2G14610","AT4G23700","AT3G26830","AT3G15950","AT3G54830","AT5G24105")
keytypes(org.At.tair.db)
select(org.At.tair.db, keys = genes,column = c('TAIR'), keytype = 'SYMBOL')
mapIds(org.At.tair.db, keys = genes,column = c('TAIR'), keytype = 'SYMBOL')

################## merge table ##################
library(data.table)
args = commandArgs(trailingOnly=TRUE)
in_1=args[1]
in_2=args[2]
out_path_to_file=args[3]
t1<- read.table(in_1,sep='')
t2<- read.table(in_2,sep='')
HiC_merge<- merge(t1,t2,by='V2',allow.cartesian = TRUE)
head(HiC_merge)
write.table(HiC_merge, file = out_path_to_file,sep = "\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.csv(mtcars, file=gzfile("mtcars.csv.gz"))

#Natural join: To keep only rows that match from the data frames, specify the argument all=FALSE.
#Full outer join: To keep all rows from both data frames, specify all=TRUE.
#Left outer join: To include all the rows of your data frame x and only those from y that match, specify all.x=TRUE.
#Right outer join: To include all the rows of your data frame y and only those from x that match, specify all.y=TRUE.

###### merge range
ranges <- merge(rangesA,rangesB,by="chrom",suffixes=c("A","B"))
ranges[with(ranges, startB <= startA & stopB >= stopA),]
#  chrom startA stopA startB stopB
#1     1    200   250    200   265
#2     5    100   105     99   106

################# aggregate ###################
# row value mean, sum
methy_merge_hic_sum<-aggregate(methy_merge_hic, by=list(methy_merge_hic$bin.x,methy_merge_hic$bin.y), sum)
aggregate(number ~ year, data=df1, mean)
dat_dedup_cis<-ddply(dat,.(V2,V3),summarise,mean(V3))

# row count sum
ddply(m_seg,.(read_id),nrow)

dt2 <- dt[,list(sumamount = sum(amount), freq = .N), by = c("id","date")]

################# scatter_plot ###################
library(ggplot2)
data("cars")
ggplot(cars, aes(x=speed, y=dist)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_smooth(method=lm, color='#2C3E50')

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}
p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)

################# statistic tests ###################
#multiple test
m_seg_nome3h_GpC$binom_padj_um<-p.adjust(m_seg_nome3h_GpC$p_val,method="BH")

########### remove character ##########
# Remove Special Characters
address_str <- "127 Anton Blvd, Apt #7 - Wilmington, DE"
new_str <- gsub('[^[:alnum:] ]','',address_str)

# Remove Multiple Characters
address_str <- "127 Anton Blvd, Apt #7 - Wilmington, DE"
new_str <- gsub('[AntonApt]','',address_str)

# Remove Single Character
address_str <- "127 Anton Blvd, Apt #7 - Wilmington, DE"
new_str <- gsub(',','',address_str)

# remove prefix and suffix
new_string_tmp <- sub(".*_", "", in_path_to_file)
new_string <- sub("*.mcool", "", new_string_tmp)

# create vector with repeated values
 x1<-rep(c(1,2,3,4,5),each=10)

########## install ##########
install.packages("path_to_source_package", repos = NULL, type = "source")

install.packages("/research/xieyeming1/proj_2023/hicVelocity_20231027/BSgenome.Mmusculus.UCSC.mm10_1.4.3.tar.gz", repos = NULL, type = "source")
sudo yum install libxml2
sudo yum install libxml2-devel
conda install -c r r-xml2

install.packages('XML', 'restfulr','rtracklayer','BSgenome')

########## colors ##########
https://github.com/thomasp85/scico

########## track ##########
https://ivanek.github.io/Gviz/reference/plotTracks.html
