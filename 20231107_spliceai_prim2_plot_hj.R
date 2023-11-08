###################
# spliceai plotting
# hojinlee
# 2023.11.07
###################

library("dplyr")
library("ggplot2")
library("stringr")
library("reshape2")
library(plyr)

setwd("/Volumes/hjdrive/thy_n25/thy_recurrent/20231107/spliceai/PRIM2/")

# 1. WT ------------------
tb <- read.table("spliceai_wt_PRIM2_hj.txt", header = T, sep = "\t")
tb <- tb %>% dplyr::mutate(neither_prob = 1-acceptor_prob-donor_prob)
tb_melt <- melt(tb, id.vars = c("position", "sequence"))
tb_melt <- tb_melt %>% dplyr::group_by(position) %>% dplyr::arrange(desc(value)) %>% dplyr::slice(1)
tb_melt <- tb_melt %>% dplyr::select(c("position", "variable"))

tb <- merge(tb, tb_melt, by.x = "position", by.y = "position", all.x = T)
tb <- tb %>% dplyr::rename("consequence" = variable)

tb_splice_acceptor <- tb %>% dplyr::filter(consequence == "acceptor_prob")
tb_splice_donor <- tb %>% dplyr::filter(consequence == "donor_prob")
tb_splice <- tb %>% dplyr::filter(consequence != "neither_prob")

write.table(tb_splice, "20231107_wt_PRIM2_splicing_table_hj.txt", sep = "\t", col.names = T, row.names = T, quote = F)

# 2. MT ------------------
tb <- read.table("spliceai_mt_PRIM2_hj.txt", header = T, sep = "\t")
tb <- tb %>% dplyr::mutate(neither_prob = 1-acceptor_prob-donor_prob)
tb_melt <- melt(tb, id.vars = c("position", "sequence"))
tb_melt <- tb_melt %>% dplyr::group_by(position) %>% dplyr::arrange(desc(value)) %>% dplyr::slice(1)
tb_melt <- tb_melt %>% dplyr::select(c("position", "variable"))

tb <- merge(tb, tb_melt, by.x = "position", by.y = "position", all.x = T)
tb <- tb %>% dplyr::rename("consequence" = variable)

tb_splice_acceptor <- tb %>% dplyr::filter(consequence == "acceptor_prob")
tb_splice_donor <- tb %>% dplyr::filter(consequence == "donor_prob")
tb_splice <- tb %>% dplyr::filter(consequence != "neither_prob")

write.table(tb_splice, "20231107_mt_PRIM2_splicing_table_hj.txt", sep = "\t", col.names = T, row.names = T, quote = F)

# 3. subset transcripts ------------------
# this step is for select transcript that is used in dowstream analysis. 

tx_db <- read.table("../20231107_ucsc_refseq_all_genome_hj.txt", header = T, sep = "\t", check.names = F, comment.char = "") # comment.char = "" is for considering "#" as character. 
tx_db <- tx_db %>% dplyr::filter(name == "NM_000947.5")
tx_db$name

tx_list <- list()
for (tx in tx_db$name) {
  temp <- tx_db %>% dplyr::filter(name == tx) %>% dplyr::select(exonStarts) %>% pull 
  exonstart <- as.numeric(strsplit(temp, ",")[[1]])
  temp <- tx_db %>% dplyr::filter(name == tx) %>% dplyr::select(exonEnds) %>% pull
  exonend <- as.numeric(strsplit(temp, ",")[[1]])
  tx_list[[tx]] <- data.frame(tx = rep(tx, length(exonstart)), exonstart = exonstart, exonend = exonend)
}
tx_df <- do.call(rbind, tx_list)

tx_df$region <- paste("exon", seq(1:nrow(tx_df)), sep = "_")
tx_df$exonstart <- tx_df$exonstart + 1

write.table(tx_df, "20231107_PRIM2_tx_table_hj.txt", sep = "\t", col.names = T, row.names = T, quote = F)

# 4.SpliceAI results ------------------
splicing_wt_df <- read.table("20231107_wt_PRIM2_splicing_table_hj.txt", header = T, sep = "\t")
splicing_wt_df <- splicing_wt_df %>% dplyr::mutate(info = case_when(consequence == "acceptor_prob" ~ acceptor_prob,
                                                              consequence == "donor_prob" ~ donor_prob,
                                                              TRUE ~ NA))

splicing_wt_df <- splicing_wt_df %>% dplyr::select(c("consequence", "position", "info")) %>% dplyr::rename("region" = consequence)
splicing_wt_df <- splicing_wt_df %>% dplyr::filter(info > 0.8) # region with high score means high probability that the position is function as splice site.  

splicing_mt_df <- read.table("20231107_mt_PRIM2_splicing_table_hj.txt", header = T, sep = "\t")
splicing_mt_df <- splicing_mt_df %>% dplyr::mutate(info = case_when(consequence == "acceptor_prob" ~ acceptor_prob,
                                                                    consequence == "donor_prob" ~ donor_prob,
                                                                    TRUE ~ NA))

splicing_mt_df <- splicing_mt_df %>% dplyr::select(c("consequence", "position", "info")) %>% dplyr::rename("region" = consequence)
splicing_mt_df <- splicing_mt_df %>% dplyr::filter(info > 0.8) # region with high score means high probability that the position is function as splice site.  

# 5.combine exon and spliceAI results ------------------
library(plyr)
tx_df <- read.table("20231107_PRIM2_tx_table_hj.txt", header = T , sep = "\t")
tx_df <- tx_df %>% dplyr::select(-c("tx"))
tx_melt <- melt(tx_df, id.vars = "region")
tx_melt <- tx_melt %>% dplyr::arrange(value) %>% dplyr::rename("info" = variable, "position" = value)

tx_splicing_df <- rbind.fill(tx_melt, splicing_wt_df) %>% arrange(position) # this dataframe will be used in step 8. 

pos_min <- min(tx_splicing_df$position - 1000)
pos_max <- max(tx_splicing_df$position + 1000)

interval <- (pos_max - pos_min)/(2*nrow(tx_melt)) # interval for making mold. 

2*nrow(tx_melt) # then we need to put exons to 140 intervals.

# 6.exon_intron range (real position) ------------------
colnames(tx_df) <- c("start", "end", "region")
intron_start <- tx_df$end + 1
intron_end <- tx_df$start -1
intron_start <- intron_start[-length(intron_start)]
intron_end <- intron_end[-1]
intron_df <- data.frame(start = intron_start,
                        end = intron_end,
                        region = paste("intron", seq(1:length(intron_end)), sep = "_"))

exon_intron <- rbind(tx_df, intron_df)
exon_intron <- exon_intron %>% dplyr::mutate(length = end - start)
rownames(exon_intron) <- NULL

exon_length_sum <- exon_intron$length[str_detect(string = exon_intron$region, pattern = "exon")] %>% sum()
intron_length_sum <- exon_intron$length[str_detect(string = exon_intron$region, pattern = "intron")] %>% sum()

exon_intron <- exon_intron %>%
  dplyr::mutate(freq_exon_intron = case_when(str_detect(string = region, pattern = "exon") ~ 100*length/exon_length_sum,
                                             str_detect(string = region, pattern = "intron") ~ 100*length/intron_length_sum,
                                             TRUE ~ NA)) %>% dplyr::arrange(start)

exon_intron %>% dplyr::select(region, freq_exon_intron)
write.table(exon_intron, "20231107_exon_intron_PRIM2_hj.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# 7.mold df for visualization ------------------

## 7.1 make mold ------------------
start_pos <- cumsum(c(0, rep(interval, 2*nrow(tx_melt)-1))) %>% round()
start_end <- cumsum(c(rep(interval, 2*nrow(tx_melt)))) %>% round()

interval_mold <- data.frame(interval_start = start_pos,
                            interval_end = start_end)


write.table(interval_mold, "20231107_interval_mold_PRIM2_hj.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## 7.2 read mold df including manually assigned exons  --------
interval_mold <- read.table("20231107_interval_mold_PRIM2_hj_edit.txt", header = T, sep = "\t", fill = T)
colnames(interval_mold) <- c("interval_start", "interval_end", "exon")

exon_mold <- interval_mold %>% dplyr::filter(str_detect(exon, pattern = "exon"))

exon_mold_final <- data.frame()
for (exon_tmp in unique(exon_mold$exon)) {
  exon_mold_tmp <- exon_mold %>% dplyr::filter(exon == exon_tmp)
  interval_start_min <- min(exon_mold_tmp$interval_start)
  interval_end_min <- max(exon_mold_tmp$interval_end)
  
  exon_mold_final_tmp <- data.frame(interval_start = interval_start_min,
                                    interval_end = interval_end_min,
                                    exon = exon_tmp)
  
  exon_mold_final <- rbind.fill(exon_mold_final, exon_mold_final_tmp)
}


splicing_df_exon_edge <- data.frame()
for (i in 1:nrow(splicing_wt_df)) {
  splicing_tmp <- splicing_wt_df[i,]
  pos_tmp <- splicing_tmp$position
  exon_intron_tmp <- exon_intron %>% dplyr::filter((start == pos_tmp) | (end == pos_tmp)) 
  if (nrow(exon_intron_tmp) != 0) {
    splicing_tmp$exon_edge <- TRUE
  } else {
    splicing_tmp$exon_edge <- FALSE 
  }
  splicing_df_exon_edge <- rbind.fill(splicing_df_exon_edge, splicing_tmp)
}

## 7.3 splicing site located in exon edge --------
splicing_df_exon_edge_TRUE <- splicing_df_exon_edge %>% dplyr::filter(exon_edge == T) # this splicing site is located in exon edge.

arrow_df <- data.frame()
for (i in 1:nrow(splicing_df_exon_edge_TRUE)) {
  tmp <- splicing_df_exon_edge_TRUE[i,]
  tmp2 <- tx_melt %>% dplyr::filter(position == tmp$position)
  
  if (str_detect(tmp2$info, pattern = "start")) {
    #print("st")
    arrow_pos <- exon_mold %>% dplyr::filter(exon == tmp2$region) %>% dplyr::pull(interval_start) %>% min()
  } else if (str_detect(tmp2$info, pattern = "end")) {
    #print("en")
    arrow_pos <- exon_mold %>% dplyr::filter(exon == tmp2$region) %>% dplyr::pull(interval_end) %>% max()
  }
  
  arrow_name <- tmp$region
  arrow_df_tmp <- data.frame(arrow_pos = arrow_pos, arrow_name = arrow_name, exon = tmp2$region)
  arrow_df <- rbind.fill(arrow_df, arrow_df_tmp)
}

write.table(arrow_df, "20231107_arrow_wt_table_hj.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## 7.4 splicing site not located in exon edge --------
splicing_df_exon_edge_FALSE <- splicing_df_exon_edge %>% dplyr::filter(exon_edge == F) # this splicing site is not located in exon edge.

# manually add the position using *_arrow_table_hj.txt file.

## 7.5 read mold df including manually assigned exons (mt) -------
interval_mold <- read.table("20231107_interval_mold_PRIM2_hj_edit.txt", header = T, sep = "\t", fill = T)
colnames(interval_mold) <- c("interval_start", "interval_end", "exon")

exon_mold <- interval_mold %>% dplyr::filter(str_detect(exon, pattern = "exon"))

splicing_df_exon_edge <- data.frame()
for (i in 1:nrow(splicing_mt_df)) {
  splicing_tmp <- splicing_mt_df[i,]
  pos_tmp <- splicing_tmp$position
  exon_intron_tmp <- exon_intron %>% dplyr::filter((start == pos_tmp) | (end == pos_tmp)) 
  if (nrow(exon_intron_tmp) != 0) {
    splicing_tmp$exon_edge <- TRUE
  } else {
    splicing_tmp$exon_edge <- FALSE 
  }
  splicing_df_exon_edge <- rbind.fill(splicing_df_exon_edge, splicing_tmp)
}

## 7.6 splicing site located in exon edge (mt) --------
splicing_df_exon_edge_TRUE <- splicing_df_exon_edge %>% dplyr::filter(exon_edge == T) # this splicing site is located in exon edge.

arrow_df <- data.frame()
for (i in 1:nrow(splicing_df_exon_edge_TRUE)) {
  tmp <- splicing_df_exon_edge_TRUE[i,]
  tmp2 <- tx_melt %>% dplyr::filter(position == tmp$position)
  
  if (str_detect(tmp2$info, pattern = "start")) {
    #print("st")
    arrow_pos <- exon_mold %>% dplyr::filter(exon == tmp2$region) %>% dplyr::pull(interval_start) %>% min()
  } else if (str_detect(tmp2$info, pattern = "end")) {
    #print("en")
    arrow_pos <- exon_mold %>% dplyr::filter(exon == tmp2$region) %>% dplyr::pull(interval_end) %>% max()
  }
  
  arrow_name <- tmp$region
  arrow_df_tmp <- data.frame(arrow_pos = arrow_pos, arrow_name = arrow_name, exon = tmp2$region)
  arrow_df <- rbind.fill(arrow_df, arrow_df_tmp)
}

write.table(arrow_df, "20231107_arrow_mt_table_hj.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## 7.7 splicing site not located in exon edge (mt) --------
splicing_df_exon_edge_FALSE <- splicing_df_exon_edge %>% dplyr::filter(exon_edge == F) # this splicing site is not located in exon edge.

## 7.8 plotting --------
p <- exon_mold_final %>% ggplot2::ggplot() +
  ### make frame ###
  xlim(min(exon_mold_final$interval_start) - 2000, max(exon_mold_final$interval_end) + 2000) + 
  ylim(0, 20) + # ylim
  theme_bw() +
  theme(panel.border = element_rect(size = 2),
        aspect.ratio = 0.2, 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())


p <- p + ggplot2::geom_segment(aes(x = min(interval_start), y = 14,
                              xend = max(interval_end), yend = 14),
                          color="#808080", linewidth=1) +
  ggplot2::geom_segment(aes(x = min(interval_start), y = 4,
                            xend = max(interval_end), yend = 4),
                        color="#808080", linewidth=1)

p <- p + ggplot2::geom_rect(aes(xmin = interval_start, xmax = interval_end,
                           ymin = rep(12, length(interval_start)), ymax = rep(16, length(interval_end))),
                       fill = "#808080") +
  ggplot2::geom_rect(aes(xmin = interval_start, xmax = interval_end,
                         ymin = rep(2, length(interval_start)), ymax = rep(6, length(interval_end))),
                     fill = "#808080")

# read wt arrow
arrow_df <- read.table("20231107_arrow_wt_table_hj.txt", header = T, sep = "\t")

p <- p + ggplot2::geom_point(data = arrow_df[arrow_df$arrow_name == "acceptor_prob",],
                          aes(x = arrow_pos, y = 17),
                          shape = 25, size = 2, fill = "#BD5CE6", colour = "#BD5CE6") + 
  ggplot2::geom_point(data = arrow_df[arrow_df$arrow_name == "donor_prob",],
                      aes(x = arrow_pos, y = 11),
                      shape = 24, size = 2, fill = "#FFBC00", colour = "#FFBC00")

# read wt arrow
arrow_df <- read.table("20231107_arrow_mt_table_hj.txt", header = T, sep = "\t")

p <- p + ggplot2::geom_point(data = arrow_df[arrow_df$arrow_name == "acceptor_prob",],
                        aes(x = arrow_pos, y = 7),
                        shape = 25, size = 2, fill = "#BD5CE6", colour = "#BD5CE6") + 
  ggplot2::geom_point(data = arrow_df[arrow_df$arrow_name == "donor_prob",],
                      aes(x = arrow_pos, y = 1),
                      shape = 24, size = 2, fill = "#FFBC00", colour = "#FFBC00")


ggsave(filename = "20231107_spliceai_plot_PRIM2_hj.pdf", plot = p, width = 8, height = 4)


# 8. exon plot for score -------
wt_df <- read.table("spliceai_wt_PRIM2_hj.txt", header = T, sep = "\t")
mt_df <- read.table("spliceai_mt_PRIM2_hj.txt", header = T, sep = "\t")

arrow_wt_df <- read.table("20231107_arrow_wt_table_hj.txt", header = T, sep = "\t")
arrow_mt_df <- read.table("20231107_arrow_mt_table_hj.txt", header = T, sep = "\t")

# exon 8 is removed out. 
tx_splicing_df # 57372288, 57372355

p <- ggplot2::ggplot() +
  ### make frame ###
  xlim(57372288-20, 57372355+20) + 
  ylim(0, 80) + # ylim
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        axis.line = element_blank(),
        aspect.ratio = 0.3, 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

p <- p + ggplot2::geom_segment(aes(x = 57372288-20, y = c(10, 50),
                                   xend = 57372355+20, yend = c(10, 50)),
                               color="#808080", linewidth = 1)

p <- p + ggplot2::geom_rect(aes(xmin = 57372288, xmax = 57372355,
                                ymin = c(10, 50) - 10, ymax = c(10, 50) + 10),
                            fill = "#808080") 

wt_df %>% dplyr::filter(position == 57372288) # 0.9583246
mt_df %>% dplyr::filter(position == 57372288) # 0.02802367
wt_df %>% dplyr::filter(position == 57372355) # 0.9788005
mt_df %>% dplyr::filter(position == 57372355) # 1.244037e-06

p <- p + ggplot2::geom_rect(aes(xmin = 57372288, xmax = 57372288 + 3,
                                ymin = 60, ymax = 60 + 15*0.9583246),
                            fill = "#BD5CE6") 
p <- p + ggplot2::geom_rect(aes(xmin = 57372288, xmax = 57372288 + 3,
                                ymin = 20, ymax = 20 + 15*0.02802367),
                            fill = "#BD5CE6") 
p <- p + ggplot2::geom_rect(aes(xmin = 57372355-3, xmax = 57372355,
                                ymin = 60, ymax = 60 + 15*0.9788005),
                            fill = "#FFBC00") 
p <- p + ggplot2::geom_rect(aes(xmin = 57372355-3, xmax = 57372355,
                                ymin = 20, ymax = 20 + 15*1.244037e-06),
                            fill = "#FFBC00") 

ggsave(filename = "20231107_exon8_PRIM2_hj.pdf", plot = p, width = 6, height = 3)

