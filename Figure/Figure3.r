#simulation2 CA MSE ablation
library(ggplot2)

data <- data.frame(
  Method = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"),
  MSE = c(1.6967, 1.9636, 4.1565, 4.2362),
  sd = c(0.0208 ,0.0422,0.0316 ,0.0432)
)
data$Method <- factor(data$Method, levels = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"))
p <- ggplot(data, aes(x = Method, y = MSE, fill = Method))+ 
geom_col(width = 0.45) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CA") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768"))+ 
ylim(0, 5)+
geom_errorbar(data=data, aes(ymin = MSE - sd, ymax = MSE + sd), width = 0.2)
ggsave("Figure/simulation2_CA_MSE_ablation.pdf", p, device = "pdf")
p



#simulation2 CRC MSE ablation
library(ggplot2)

data <- data.frame(
  Method = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"),
  MSE = c(1.7354, 2.1448, 4.3059, 4.4256),
  sd = c(0.0208 ,0.0422,0.0316 ,0.0432)
)
data$Method <- factor(data$Method, levels = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"))
p <- ggplot(data, aes(x = Method, y = MSE, fill = Method))+ 
geom_col(width = 0.45) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CRC") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768"))+ 
ylim(0, 5)+
geom_errorbar(data=data, aes(ymin = MSE - sd, ymax = MSE + sd), width = 0.2)
ggsave("Figure/simulation2_CRC_MSE_ablation.pdf", p, device = "pdf")
p


#simulation2 CA Pearson ablation
library(ggplot2)

data <- data.frame(
  Method = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"),
  Pearson = c(0.2777, 0.2591, 0.2161, 0.2145),
  sd = c(0.0085, 0.0095, 0.0090 ,0.0103)
)
data$Method <- factor(data$Method, levels = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"))
p <- ggplot(data, aes(x = Method, y = Pearson, fill = Method))+ 
geom_col(width = 0.45) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CA") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768"))+ 
ylim(0, 0.3)+
geom_errorbar(data=data, aes(ymin = Pearson - sd, ymax = Pearson + sd), width = 0.2)
ggsave("Figure/simulation2_CA_Pearson_ablation.pdf", p, device = "pdf")
p


#simulation2 CRC Pearson ablation
library(ggplot2)

data <- data.frame(
  Method = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"),
  Pearson = c(0.2648, 0.2311, 0.2102, 0.2172),
  sd = c(0.0085, 0.0095, 0.0090 ,0.0103)
)
data$Method <- factor(data$Method, levels = c("mbSparse", "No AE", "No CVAE", "No AE And CVAE"))
p <- ggplot(data, aes(x = Method, y = Pearson, fill = Method))+ 
geom_col(width = 0.45) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CRC") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768"))+ 
ylim(0, 0.3)+
geom_errorbar(data=data, aes(ymin = Pearson - sd, ymax = Pearson + sd), width = 0.2)
ggsave("Figure/simulation2_CRC_Pearson_ablation.pdf", p, device = "pdf")
p


library(ggplot2)
taxa <- c("k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Micrococcales|f__Micrococcaceae|g__Rothia|s__Rothia_mucilaginosa", "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_hadrus")
taxa <- rbind(taxa, c("k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_dentium", "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii"),
                    c("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_ovatus", "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium|s__Eubacterium_hallii"),
                    c("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii", "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_finegoldii"),
                    c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_haemolysans", "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_finegoldii"),
                    c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_haemolysans", "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_stercoris"))
                      


OTU <- read.csv("PRJEB7774_data/PRJEB7774.csv", row.names="X")
OTU <- t(OTU)
OTU <- OTU / rowSums(OTU) * 10^6
OTU <- log10(OTU+1)
meta_data_case <- read.csv("PRJEB7774_data/PRJEB7774_metadata_case.csv")
case_ori <- OTU[meta_data_case$Run, ]
meta_data_control <- read.csv("PRJEB7774_data/PRJEB7774_metadata_control.csv")
control_ori <- OTU[meta_data_control$Run, ]

OTU_imputed <- read.csv("PRJEB7774_data/PRJEB7774_impute_normalize.csv", row.names="X")
OTU_imputed <- t(OTU_imputed)
case_imputed <- OTU_imputed[meta_data_case$Run, ]
control_imputed <- OTU_imputed[meta_data_control$Run, ]

for(i in 1:dim(taxa)[1]){
  x_taxa_name <- sub(".*s__(.*)", "\\1", taxa[i,1])
  y_taxa_name <- sub(".*s__(.*)", "\\1", taxa[i,2])
  file_path <- paste("Figure/PRJEB7774_", x_taxa_name, "_", y_taxa_name, ".pdf", sep = "")
  pdf(file_path , width = 12, height =12)
  x_taxa_name <- gsub("_+", " ", x_taxa_name)
  y_taxa_name <- gsub("_+", " ", y_taxa_name)
  par(mfrow = c(2,2))
  par(mar = c(1, 1, 1, 1) + 5 )
  plot(y = control_ori[, taxa[i,1]], x = control_ori[, taxa[i,2]],
       xlab = '', ylab = '', main = "control data",cex.main = 3,
       xlim = c(0, 6), ylim = c(0, 6), pch=16, col = "#BBBABB", cex = 3, xaxt="n", yaxt="n", font.main = 1)
  axis(side = 1, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5, mgp = c(1.5, 1.5, 0))
  axis(side = 2, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5,las=2)
  box(col = "gray")

  nz_idx <- intersect(which(control_ori[,taxa[i,1]] > log10(1)), which(control_ori[,taxa[i,2]] > log10(1)))
  abline(coef = lm(control_ori[,taxa[i,1]] ~ control_ori[,taxa[i,2]])$coefficients, col = "black", lwd = 4)
  abline(coef = lm(control_ori[nz_idx,taxa[i,1]] ~ control_ori[nz_idx,taxa[i,2]])$coefficients, col = alpha("#1686C5", 0.5), lwd = 4)
  legend(x='topright', legend= c("Cor:",round(cor(control_ori[,taxa[i,1]], control_ori[,taxa[i,2]]), digits = 2),             
                                 round(cor(control_ori[nz_idx,taxa[i,1]], control_ori[nz_idx,taxa[i,2]]), digits = 2)),
                                 lty = c(0,1,1), cex = 1.5, col=c("white","black", alpha("#1686C5", 0.5)), lwd = 4, bty="n")
  par(mar = c(1, 1, 1, 1) + 5 )  
  plot(y = control_imputed[, taxa[i,1]], x = control_imputed[, taxa[i,2]],
       xlab = '', ylab = '', main = "impute + control data",cex.main = 3,
      xlim = c(0, 6), ylim = c(0, 6),pch=16, col = "gray", cex = 3, xaxt="n", yaxt="n",font.main = 1)
  axis(side = 1, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5, mgp = c(1.5, 1.5, 0))
  axis(side = 2, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5,las=2)
  box(col = "gray")
  abline(coef = lm(control_imputed[,taxa[i,1]] ~ control_imputed[,taxa[i,2]])$coefficients, col = alpha("#1686C5", 0.5), lwd = 4)
  legend(x='topright', legend= c("Cor:",round(cor(control_imputed[,taxa[i,1]], control_imputed[,taxa[i,2]]), digits = 2)), 
         lty = c(0,1), cex = 1.5, col=c("white", alpha("#1686C5", 0.5)), lwd = 4, bty="n")
  par(mar = c(1, 1, 1, 1) + 5 )
  plot(y = case_ori[, taxa[i,1]], x = case_ori[, taxa[i,2]],
       xlab = '',, ylab = '',, main = "case data ", cex.main =3,
       xlim = c(0, 6), ylim = c(0, 6),pch=16, col = "gray", cex = 3, xaxt="n", yaxt="n",font.main = 1)
  axis(side = 1, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5, mgp = c(1.5, 1.5, 0))
  axis(side = 2, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5,las=2)     
  box(col = "gray")
  nz_idx <- intersect(which(case_ori[,taxa[i,1]] > log10(1)), which(case_ori[,taxa[i,2]] > log10(1)))
  abline(coef = lm(case_ori[,taxa[i,1]] ~ case_ori[,taxa[i,2]])$coefficients, col = "black", lwd = 4)
  abline(coef = lm(case_ori[nz_idx,taxa[i,1]] ~ case_ori[nz_idx,taxa[i,2]])$coefficients, col = alpha("#1686C5", 0.5), lwd = 4)
  legend(x='topright', legend= c("Cor:",round(cor(case_ori[,taxa[i,1]], case_ori[,taxa[i,2]]), digits = 2),                                
                                 round(cor(case_ori[nz_idx,taxa[i,1]], case_ori[nz_idx,taxa[i,2]]), digits = 2)), 
                                 lty = c(0,1,1), cex = 1.5, col=c("white","black", alpha("#1686C5", 0.5)), lwd = 4, bty="n")
  par(mar = c(1, 1, 1, 1) + 5)
  plot(y = case_imputed[, taxa[i,1]], x = case_imputed[, taxa[i,2]],
       xlab = '',, ylab = '', main = "impute + case data",cex.main = 3,
       xlim = c(0, 6), ylim = c(0, 6),pch=16, col = "gray", cex = 3, xaxt="n", yaxt="n",font.main = 1)
  axis(side = 1, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5, mgp = c(1.5, 1.5, 0))
  axis(side = 2, at = seq(0, 6), labels = TRUE, tck = -0.01, cex.axis=2.5,las=2)      
  box(col = "gray")
  abline(coef = lm(case_imputed[,taxa[i,1]] ~ case_imputed[,taxa[i,2]])$coefficients, col = alpha("#1686C5", 0.5), lwd = 4)
  legend(x='topright', legend= c("Cor:",round(cor(case_imputed[,taxa[i,1]], case_imputed[,taxa[i,2]]), digits = 2)),
         lty = c(0,1), cex = 1.5, col=c("white", alpha("#1686C5", 0.5)), lwd = 4, bty="n")
  mtext(x_taxa_name, side = 1, line = -2, outer = TRUE, cex = 3 ,font = 1, padj = 0.1)
  mtext(y_taxa_name, side = 2, line = -3, outer = TRUE, cex = 3 ,font = 1, padj = 0.1)
  dev.off()
}