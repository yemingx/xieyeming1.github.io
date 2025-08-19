violin_2compare_test
================
<yemingxie@gmail.com>
Mon Aug 18 18:15:35 CST 2025

``` r
knitr::opts_chunk$set(echo = TRUE)

library(ggplot2)
library(ggsci)
library(httpgd)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(ggsignif)

args = commandArgs(trailingOnly=TRUE)
in_path='/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/overlap_metrics'
plot_title='hek_atac_overlap_pval'
y_axis='log2_MACS_q_val'
in_file1='hek_atac_non_overlap_pval.bdg'
in_file2='hek_atac_overlap_pval.bdg'

npg_colors <- pal_npg()(10)

setwd(in_path)

in_file_1 <- read.table(in_file1, sep="\t", header=F)
in_file_1$subgroup<-in_file1
in_file_2 <- read.table(in_file2, sep="\t", header=F)
in_file_2$subgroup<-in_file2

# combine all data
all_data <- rbind(in_file_1, in_file_2)
all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)

wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4, 
                           alternative = "two.sided")
head(all_data)
```

    ##     V1    V2    V3       V4                      subgroup feature_val
    ## 1 chr1 10003 10582 78.60654 hek_atac_non_overlap_pval.bdg    6.296577
    ## 2 chr1 15593 15698  7.47299 hek_atac_non_overlap_pval.bdg    2.901686
    ## 3 chr1 32572 32667  4.47711 hek_atac_non_overlap_pval.bdg    2.162568
    ## 4 chr1 34721 34900  5.11197 hek_atac_non_overlap_pval.bdg    2.353879
    ## 5 chr1 79262 79336  3.60646 hek_atac_non_overlap_pval.bdg    1.850583
    ## 6 chr1 91365 91519  9.68647 hek_atac_non_overlap_pval.bdg    3.275971

``` r
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE))
summ
```

    ## # A tibble: 2 × 5
    ##   subgroup                          n  mean max_val    sd
    ##   <chr>                         <int> <dbl>   <dbl> <dbl>
    ## 1 hek_atac_non_overlap_pval.bdg 58400  2.73    11.5  1.06
    ## 2 hek_atac_overlap_pval.bdg     33013  3.71    11.6  1.40

``` r
levels(factor(all_data$subgroup))
```

    ## [1] "hek_atac_non_overlap_pval.bdg" "hek_atac_overlap_pval.bdg"

``` r
# pdf(paste0(plot_title,'.pdf'))
ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.3, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(max_val, na.rm = TRUE)), 
            data = summ, size=4, vjust = 2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("lightblue","indianred")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")
```

![](violin_2compare_test_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# dev.off()
```
