library(ggplot2)
data <- data.frame(
  x = rep(1:3, each = 2),
  Method = rep(c("Non-imputation", "mbSparse"), times = 3),
  mean = c(0.5833333,	0.5714286,	
           0.28,	0.32,
           0.3783784,	0.4102564)
)

data$Method <- factor(data$Method, levels = c("Non-imputation", "mbSparse"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
  geom_col(position = "dodge") +
  scale_x_discrete(labels = c("Precision", "Recall", "F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "Wilcoxon") +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
  theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 18) ,axis.title.x = element_blank(), legend.title =  element_text(size = 18)) + 
  scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+
  ylab("Rate(%)")+
  scale_y_continuous(limits = c(0, 0.75)) 
ggsave("Figure/wilcoxon.tiff", p, device = "tiff")
# ggsave("Figure/difference_analisis_wilcoxon.pdf", p, device = "pdf")
p



library(ggplot2)
data <- data.frame(
  x = rep(1:3, each = 2),
  Method = rep(c("Non-imputation", "mbSparse"), times = 3),
  mean = c(0.583333,	0.6521739,	
           0.28,	0.3,
           0.3783784,	0.4109589)
)

data$Method <- factor(data$Method, levels = c("Non-imputation", "mbSparse"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
  geom_col(position = "dodge") +
  scale_x_discrete(labels = c("Precision", "Recall", "F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "LEfSe") +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
  theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 18) ,axis.title.x = element_blank(), legend.title =  element_text(size = 18)) + 
  scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+
  ylab("Rate(%)")+
  scale_y_continuous(limits = c(0, 0.75)) 
ggsave("Figure/LEfSe.tiff", p, device = "tiff")
# ggsave("Figure/difference_analisis_wilcoxon.pdf", p, device = "pdf")
p



library(ggplot2)
data <- data.frame(
  x = rep(1:3, each = 2),
  Method = rep(c("Non-imputation", "mbSparse"), times = 3),
  mean = c(0.001,	0.1492537,	
           0.001,	0.2,
           0.001,	0.1709402)
)

data$Method <- factor(data$Method, levels = c("Non-imputation", "mbSparse"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
  geom_col(position = "dodge") +
  scale_x_discrete(labels = c("Precision", "Recall", "F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "edgeR") +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
  theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 18) ,axis.title.x = element_blank(), legend.title =  element_text(size = 18)) + 
  scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+
  ylab("Rate(%)")+
  scale_y_continuous(limits = c(0, 0.75)) 
ggsave("Figure/edgeR.tiff", p, device = "tiff")
# ggsave("Figure/difference_analisis_wilcoxon.pdf", p, device = "pdf")
p


library(ggplot2)
data <- data.frame(
  x = rep(1:3, each = 2),
  Method = rep(c("Non-imputation", "mbSparse"), times = 3),
  mean = c(0.4545455,	0.5675676,	
           0.1,	0.42,
           0.1639344,	0.4827586)
)

data$Method <- factor(data$Method, levels = c("Non-imputation", "mbSparse"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
  geom_col(position = "dodge") +
  scale_x_discrete(labels = c("Precision", "Recall", "F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "ALDEx2") +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
  theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 18) ,axis.title.x = element_blank(), legend.title =  element_text(size = 18)) + 
  scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+
  ylab("Rate(%)")+
  scale_y_continuous(limits = c(0, 0.75)) 
ggsave("Figure/ALDEx2.tiff", p, device = "tiff")
# ggsave("Figure/difference_analisis_wilcoxon.pdf", p, device = "pdf")
p


library(ggplot2)
data <- data.frame(
  x = rep(1:3, each = 2),
  Method = rep(c("Non-imputation", "mbSparse"), times = 3),
  mean = c(0.001,	0.2352941,	
           0.001,	0.4,
           0.001,	0.2962963)
)

data$Method <- factor(data$Method, levels = c("Non-imputation", "mbSparse"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
  geom_col(position = "dodge") +
  scale_x_discrete(labels = c("Precision", "Recall", "F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "DESeq2 phyloseq") +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
  theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 18) ,axis.title.x = element_blank(), legend.title =  element_text(size = 18)) + 
  scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+
  ylab("Rate(%)")+
  scale_y_continuous(limits = c(0, 0.75)) 
ggsave("Figure/DESeq2_phyloseq.tiff", p, device = "tiff")
# ggsave("Figure/difference_analisis_wilcoxon.pdf", p, device = "pdf")
p