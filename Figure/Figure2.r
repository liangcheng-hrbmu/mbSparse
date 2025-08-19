#simulation2 CA MSE
library(ggplot2)

data <- data.frame(
  Method = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"),
  MSE = c(5.3451,1.6974,2.7413 ,5.1816 ,2.5589 ,4.9892 ,4.9947 ,2.7288),
  sd = c(0.0831,0.0249 ,0.0502 ,0.1211 ,0.1952 ,0.1039 ,0.07 ,0.08)
)
data$Method <- factor(data$Method, levels = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"))
p <- ggplot(data, aes(x = Method, y = MSE, fill = Method))+ 
geom_col(width = 0.8) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CA") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+ 
ylim(0, 6)+
geom_errorbar(data=data, aes(ymin = MSE - sd, ymax = MSE + sd), width = 0.2)
ggsave("Figure/simulation2_CA_MSE.pdf", p, device = "pdf")
p


#simulation2 CRC MSE
library(ggplot2)

data <- data.frame(
  Method = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"),
  MSE = c(5.4562 ,1.7411,2.7147 ,5.2634 ,2.4954 ,5.1629,5.0872, 3.6933  ),
  sd = c(0.0656 ,0.0208,0.0516 ,0.0822 ,0.1560 ,0.0754, 0.1135 ,0.1223)
)
data$Method <- factor(data$Method, levels = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"))
p <- ggplot(data, aes(x = Method, y = MSE, fill = Method))+ 
geom_col(width = 0.8) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CRC") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12)) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+ 
ylim(0, 6)+
geom_errorbar(data=data, aes(ymin = MSE - sd, ymax = MSE + sd), width = 0.2)
ggsave("Figure/simulation2_CRC_MSE.pdf", p, device = "pdf")
p


#simulation2 CA Pearson
library(ggplot2)

data <- data.frame(
  Method = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"),
  Pearson = c(0.1818,0.2672 ,0.2385 ,0.1642 ,0.1354 ,0.1619, 0.1352, 0.2769),
  sd = c(0.0065,0.0107,0.0123 ,0.0101 ,0.0072 ,0.0052, 0.0103 ,0.0112)
)
data$Method <- factor(data$Method, levels = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"))
p <- ggplot(data, aes(x = Method, y = Pearson, fill = Method))+ 
geom_col(width = 0.8) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CA") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+ 
ylim(0, 0.3)+
geom_errorbar(data=data, aes(ymin = Pearson - sd, ymax = Pearson + sd), width = 0.2)
ggsave("Figure/simulation2_CA_Pearson.pdf", p, device = "pdf")
p


#simulation2 CRC Pearson
library(ggplot2)

data <- data.frame(
  Method = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"),
  Pearson = c(0.1775 ,0.2404 ,0.1551 ,0.1485 ,0.1573 ,0.2548, 0.1435, 0.2023  ),
  sd = c(0.0059,0.0085, 0.0090 ,0.0037 ,0.0118 ,0.0071, 0.0101, 0.0092)
)
data$Method <- factor(data$Method, levels = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"))
p <- ggplot(data, aes(x = Method, y = Pearson, fill = Method))+ 
geom_col(width = 0.8) +
theme(panel.background = element_rect(fill = "white")) +
labs(title = "CRC") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_line())+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ) + 
theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth=1)) +
    scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+ 
ylim(0, 0.3)+
geom_errorbar(data=data, aes(ymin = Pearson - sd, ymax = Pearson + sd), width = 0.2)
ggsave("Figure/simulation2_CRC_Pearson.pdf", p, device = "pdf")
p


#simulation1 CA MSE
library(ggplot2)
df <- data.frame(
  x = c(40, 50, 60),
  orginal = c(8.4718, 6.9595, 5.9055),
  softImpute = c(5.2259, 5.1917, 5.2644),
  SAVER = c(8.1608, 6.7355, 5.7307),
  mbimpute = c(5.1782, 5.6688, 5.7609),
  GE_impute = c(8.2393, 6.8599, 5.8173),
  mbSparse = c(4.3539, 3.8291, 3.5912),
  Mbdenoise_pn = c(8.0258, 6.2119, 5.6321),
  Mbdenoise_lmn = c(8.4936,6.9954, 5.8793) 
)

names(df) <- c("x", "No imputation", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbSparse", "mbDenoise-pn", "mbDenoise-lmn")
df_long <- tidyr::pivot_longer(df, cols = -x, names_to = "Method", values_to = "y")
method_order <- c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn")
# colors <- c("#D4D4D4", "#FF8063", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#80718E")
colors <- c("#D4D4D4", "#C83528", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#FF8063")
plot <- ggplot(df_long, aes(x = x, y = y, color = Method)) +
  geom_line(linetype = "solid", linewidth = 1) + 
  geom_point() +  
  labs(title = "CA", x = "Rate", y = "MSE") + 
  theme_minimal() +
  scale_color_manual(values = colors, limits = method_order)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm"))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.title.x = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  scale_x_continuous(breaks = seq(40, 60, by = 10))+
  ylim(0, 9)
ggsave("Figure/simulation1_CA_MSE.pdf", plot, device = "pdf")
plot


#simulation1 CRC MSE
library(ggplot2)
df <- data.frame(
  x = c(40, 50, 60),
  orginal = c(8.1817, 6.7367, 5.7445),
  softImpute = c(5.2282, 5.3273, 5.2901),
  SAVER = c(7.9345, 6.5724, 5.6016),
  mbimpute = c(6.3207, 6.0547,5.8182),
  GE_impute = c(7.9704, 6.6738, 5.6766),
  mbSparse = c(4.4385, 3.7857, 3.5049),
  Mbdenoise_pn = c(7.6747, 6.2530, 5.5314),
  Mbdenoise_lmn = c(8.2904, 6.8598, 5.7748)
)

names(df) <- c("x", "No imputation", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbSparse", "mbDenoise-pn", "mbDenoise-lmn")
df_long <- tidyr::pivot_longer(df, cols = -x, names_to = "Method", values_to = "y")
method_order <- c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn")
# colors <- c("#D4D4D4", "#FF8063", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#80718E")
colors <- c("#D4D4D4", "#C83528", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#FF8063")
plot <- ggplot(df_long, aes(x = x, y = y, color = Method)) +
  geom_line(linetype = "solid", linewidth = 1) + 
  geom_point() +  
  labs(title = "CRC", x = "Rate", y = "MSE") + 
  theme_minimal() +
  scale_color_manual(values = colors, limits = method_order)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm"))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.title.x = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  scale_x_continuous(breaks = seq(40, 60, by = 10))+
  ylim(0, 9)
ggsave("Figure/simulation1_CRC_MSE.pdf", plot, device = "pdf")
plot



#simulation1 CA Pearson
library(ggplot2)
df <- data.frame(
  x = c(40, 50, 60),
  orginal = c(0.2319 ,0.2893 , 0.3246),
  softImpute = c(0.1015, 0.1455, 0.1528),
  SAVER = c(0.2375, 0.2938, 0.3269),
  mbimpute = c(0.0690, 0.1416, 0.1761),
  GE_impute = c(0.2228, 0.2781, 0.3191),
  mbSparse = c(0.2706, 0.3382, 0.3783), 
  Mbdenoise_pn = c(0.1985, 0.2330 ,0.2911),
  Mbdenoise_lmn = c(0.1976, 0.2577, 0.2892)
)

names(df) <- c("x", "No imputation", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbSparse", "mbDenoise-pn", "mbDenoise-lmn")
df_long <- tidyr::pivot_longer(df, cols = -x, names_to = "Method", values_to = "y")
method_order <- c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn")
# colors <- c("#D4D4D4", "#FF8063", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#80718E")
# colors <- c("#D4D4D4", "#C83528", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#FF8063")
colors <- c("#D4D4D4", "#C83528", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#FF8063")
plot <- ggplot(df_long, aes(x = x, y = y, color = Method)) +
  geom_line(linetype = "solid", linewidth = 1) + 
  geom_point() +  
  labs(title = "CA", x = "Rate", y = "Pearson") + 
  theme_minimal() +
  scale_color_manual(values = colors, limits = method_order)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm"))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.title.x = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  scale_x_continuous(breaks = seq(40, 60, by = 10)) +
  ylim(0, 0.4)
ggsave("Figure/simulation1_CA_Pearson.pdf", plot, device = "pdf")
plot


#simulation1 CRC Pearson
library(ggplot2)
df <- data.frame(
  x = c(40, 50, 60),
  orginal = c(0.2836, 0.3249 , 0.3384),
  softImpute = c(0.1438, 0.1337,0.1301),
  SAVER = c(0.2856,0.3268,0.3410),
  mbimpute = c( 0.0443, 0.1760,0.2057),
  GE_impute = c(0.2697, 0.3121 ,0.3323 ),
  mbSparse = c(0.2998, 0.3476,0.3911),
  Mbdenoise_pn = c(0.2642, 0.2934,0.2835),
  Mbdenoise_lmn = c(0.2645, 0.2965, 0.3007)
)

names(df) <- c("x", "No imputation", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbSparse", "mbDenoise-pn", "mbDenoise-lmn")
df_long <- tidyr::pivot_longer(df, cols = -x, names_to = "Method", values_to = "y")
method_order <- c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn")
colors <- c("#D4D4D4", "#C83528", "#76A5CB", "#488768", "#EAA5C4", "#1AFEEA", "#AC929F", "#FF8063")
# colors <- c("#D4D4D4", "#C83528", "#4288ED", "#86EB68", "#AA5ADD", "#1AFEEA", "#F5C548", "#E95162")
plot <- ggplot(df_long, aes(x = x, y = y, color = Method)) +
  geom_line(linetype = "solid", linewidth = 1) + 
  geom_point() +  
  labs(title = "CRC", x = "Rate", y = "Pearson") + 
  theme_minimal() +
  scale_color_manual(values = colors, limits = method_order)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm"))+
  theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.title.x = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  scale_x_continuous(breaks = seq(40, 60, by = 10)) +
  ylim(0, 0.4)
ggsave("Figure/simulation1_CRC_Pearson.pdf", plot, device = "pdf")
plot

#1000
library(ggplot2)
df <- cbind( c(0.1693,0.1684,0.1630,0.1645,	0.1724,0.1264,0.1271,0.1267,0.1278,	0.1270 ,0.1381 	,	0.1368 	,	0.1351 	,	0.1311 	,	0.1414 , 0.1351 	,	0.1359 	,	0.1292 	,	0.1337 	,	0.1417 , 0.1410 	,	0.1416 	,	0.1396 	,	0.1402 	,	0.1429 , 0.1410 	,	0.1416 	,	0.1396 	,	0.1402 	,	0.1429 ,  0.1491 	,	0.1543 	,	0.1606 	,	0.1565 	,	0.1552 ),
      c( rep("original", 5), rep("0%", 5), rep("20%", 5),  rep("40%", 5), rep("60%", 5), rep("80%", 5), rep("100%", 5)))
colnames(df) <- c("MSE", "Status")
df <- data.frame(df)
df$MSE <- as.numeric(as.character(df$MSE))
df$Status <- factor(df$Status, levels = c("original", "0%", "20%", "40%","60%", "80%", "100%"))
plot <- ggplot(df, aes(x=Status, y=MSE)) +
  geom_boxplot() +
  labs(title = "Sequencing depth 1000") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 4, 1), "cm"))+
  theme(axis.line.y = element_blank(),axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 17)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  xlab("Status") +
  ylab("MSE")+
  ylim(0.12, 0.18)
plot
ggsave("Figure/Difference_rate_ouliter_Sequencing_depth_1000.pdf", plot, device = "pdf")

#2000
library(ggplot2)
df <- cbind( c( 0.1658, 0.1687, 0.1654, 0.1677, 0.1649,0.1132, 0.1164, 0.1107, 0.1158, 0.1139 , 0.1242  , 0.1295  , 0.1250  , 0.1278  , 0.1244 , 0.1227  , 0.1264  , 0.1264  , 0.1254  , 0.1237 , 0.1328  , 0.1384  , 0.1388  , 0.1380  , 0.1334 ,0.1440  , 0.1443  , 0.1442  , 0.1451  , 0.1428 , 0.1350  , 0.1383  , 0.1472  , 0.1456  , 0.1432),
      c( rep("original", 5), rep("0%", 5), rep("20%", 5),  rep("40%", 5), rep("60%", 5), rep("80%", 5), rep("100%", 5)))
colnames(df) <- c("MSE", "Status")
df <- data.frame(df)
df$MSE <- as.numeric(as.character(df$MSE))
df$Status <- factor(df$Status, levels = c("original", "0%", "20%", "40%","60%", "80%", "100%"))
plot <- ggplot(df, aes(x=Status, y=MSE)) +
  geom_boxplot() +
  labs(title = "Sequencing depth 2000") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 4, 1), "cm"))+
  theme(axis.line.y = element_blank(),axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 17)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  xlab("Status") +
  ylab("MSE") +
  ylim(0.11, 0.18)
plot
ggsave("Figure/Difference_rate_ouliter_Sequencing_depth_2000.pdf", plot, device = "pdf")

#5000
library(ggplot2)
df <- cbind( c(0.1611, 0.1666, 0.1661, 0.1722, 0.1723, 0.0988, 0.1065, 0.1033, 0.1042, 0.1085, 0.1062  , 0.1169  , 0.1191  , 0.1166  , 0.1215 , 0.1123  , 0.1188  , 0.1226  , 0.1226  , 0.1237 , 0.1126  , 0.1247  , 0.1297  , 0.1284  , 0.1282 ,  0.1309  , 0.1321  , 0.1349  , 0.1341  , 0.1460 ,  0.1167  , 0.1225  , 0.1217  , 0.1248  , 0.1328 ),
      c(rep("original", 5), rep("0%", 5), rep("20%", 5),  rep("40%", 5), rep("60%", 5), rep("80%", 5), rep("100%", 5)))
colnames(df) <- c("MSE", "Status")
df <- data.frame(df)
df$MSE <- as.numeric(as.character(df$MSE))
df$Status <- factor(df$Status, levels = c("original", "0%", "20%", "40%","60%", "80%", "100%"))
plot <- ggplot(df, aes(x=Status, y=MSE)) +
  geom_boxplot() +
  labs(title = "Sequencing depth 5000") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 4, 1), "cm"))+
  theme(axis.line.y = element_blank(),axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 17)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  xlab("Status") +
  ylab("MSE") +
  ylim(0.09, 0.18)
plot
ggsave("Figure/Difference_rate_ouliter_Sequencing_depth_5000.pdf", plot, device = "pdf")

#10000
library(ggplot2)
df <- cbind( c(0.1643 ,0.1657, 0.1630, 0.1597, 0.1607,0.1033, 0.1014, 0.1031, 0.0993, 0.1018, 0.1136  , 0.1117  , 0.1097  , 0.1112  , 0.1090 , 0.1121  , 0.1138  , 0.1181  , 0.1144  , 0.1151 , 0.1250  , 0.1212  , 0.1200  , 0.1152  , 0.1165 , 0.1188  , 0.1257  , 0.1221  , 0.1156  , 0.1155 ,  0.1166  , 0.1185  , 0.1178  , 0.1154  , 0.1148  ),
      c( rep("original", 5), rep("0%", 5), rep("20%", 5),  rep("40%", 5), rep("60%", 5), rep("80%", 5), rep("100%", 5)))
colnames(df) <- c("MSE", "Status")
df <- data.frame(df)
df$MSE <- as.numeric(as.character(df$MSE))
df$Status <- factor(df$Status, levels = c("original", "0%", "20%", "40%","60%", "80%", "100%"))
plot <- ggplot(df, aes(x=Status, y=MSE)) +
  geom_boxplot() +
  labs(title = "Sequencing depth 10000") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 4, 1), "cm"))+
  theme(axis.line.y = element_blank(),axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 17), axis.title.x = element_text(size = 17)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  xlab("Status") +
  ylab("MSE") +
  ylim(0.09, 0.18)
plot
ggsave("Figure/Difference_rate_ouliter_Sequencing_depth_10000.pdf", plot, device = "pdf")

library(ggplot2)
data <- data.frame(
  x = rep(1:4, each = 8),
  Method = rep(c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"), times = 4),
  mean = c(0.1675,0.1270, 0.2046, 0.1684, 0.2044, 0.1597, 0.1577, 0.1900,
           0.1665,0.1136, 0.1777, 0.1649, 0.1708, 0.1565, 0.1458, 0.1579,
           0.1689,0.1043, 0.1608, 0.1609, 0.1554, 0.1554, 0.1420, 0.1535,
          0.1620, 0.1018, 0.1548, 0.1561, 0.1450, 0.1496, 0.1411, 0.1499),
 sd = c(0.0038, 0.0010, 0.0018, 0.0103, 0.0023, 0.0043, 0.0017, 0.0089,
        0.0016, 0.0021, 0.0049, 0.0010, 0.0014, 0.0017, 0.0035, 0.0021,
        0.0050, 0.0037, 0.0061, 0.0038, 0.0026, 0.0055,0.0039, 0.0049,
        0.0028, 0.0016, 0.0026, 0.0015, 0.0016, 0.0027, 0.0033, 0.0053)
)

data$Method <- factor(data$Method, levels = c("No imputation", "mbSparse", "softImpute", "Saver", "mbImpute", "GE-Impute", "mbDenoise-pn", "mbDenoise-lmn"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
geom_col(position = "dodge") +
scale_x_discrete(labels = c("1000", "2000", "5000", "10000")) +
theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
theme(panel.background = element_rect(fill = "white")) +
labs(title = "Difference Sequencing depth") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24))+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ,axis.title.x = element_blank()) + 
scale_fill_manual(values = c("#C0C0C0", "#29315F", "#AFBF82", "#488768", "#B7CBDA", "#45739F", "#DADCAC", "#87B1C3"))+
ylab("MSE")+
xlab("Sequencing depth")+
geom_errorbar(data=data, aes(ymin = mean - sd, ymax = mean + sd), width = 0.9, position = "dodge")
ggsave("Figure/Difference_Sequencing_depth_all.pdf", p, device = "pdf")
p


library(ggplot2)
df <- cbind( c( 0.1693, 0.1684, 0.1630, 0.1645, 0.1724, 0.1264, 0.1271,	0.1267,0.1278, 0.1270, 0.1658,	0.1687,	0.1654, 0.1677, 0.1649, 0.1132, 0.1164, 0.1107, 0.1139, 0.1139, 0.1611, 0.1666, 0.1661, 0.1677, 0.1723, 0.0988, 0.1065, 0.1033, 0.1042, 0.1085, 0.1643, 0.1657, 0.1630, 0.1596, 0.1607, 0.1033,	0.1014,	0.1031, 0.0993, 0.1018),
       rep( c( rep("No imputation", 5), rep("mbSparse", 5) ), 4),
       c( rep(1000, 10), rep(2000, 10), rep(5000, 10), rep(10000, 10) ) )
colnames(df) <- c("mse", "status", "seq_depth")
df <- data.frame(df)
df$mse <- as.numeric(as.character(df$mse))
df$status <- factor(df$status, levels = c("No imputation", "mbSparse"))
df$sequencing_depth <- factor(df$seq_depth, levels = c("1000", "2000", "5000", "10000"))
plot <- ggplot(df, aes(x=seq_depth, y=mse, fill=status)) +
  geom_boxplot() + 
  scale_x_discrete(limits = c("1000", "2000", "5000", "10000")) +
  xlab("Sequence Depth") +
  ylab("MSE") +
  labs(fill = "Method")+
  labs(title = "Difference Sequencing depth") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2.8, 1), "cm"))+
  theme(axis.line.y = element_blank(),axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.title.x = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),legend.title = element_text(size = 12) )+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_line(), axis.ticks.y = element_line())+
  scale_fill_manual(values = c("#C0C0C0", "#29315F")) 
ggsave("Figure/Difference_Sequencing_depth.pdf", plot, device = "pdf")
plot


library(ggplot2)
data <- data.frame(
  x = rep(1:3, each = 2),
  Method = rep(c("No imputation", "mbSparse"), times = 3),
  mean = c(0.895, 0.934,
           0.875,0.918,
           0.85,0.93)
)

data$Method <- factor(data$Method, levels = c("No imputation", "mbSparse"))
p <- ggplot(data, aes(x = factor(x), y = mean, fill = Method))+
geom_col(position = "dodge", width=0.5) +
scale_x_discrete(labels = c("Feng et al.", "Wirbel et al.", "Hale et al.")) +
theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
theme(panel.background = element_rect(fill = "white")) +
labs(title = "PR-AUC Result") +
theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), plot.subtitle = element_text(hjust = 0.5), plot.margin = unit(c(2, 2, 2, 1), "cm")) + 
theme(axis.ticks.x = element_blank(), axis.line.x = element_line(), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20))+
theme(axis.line.y = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24)) +
theme(legend.position = "bottom", legend.text = element_text(size = 12) ,axis.title.x = element_blank()) + 
scale_fill_manual(values = c("#C0C0C0", "#29315F"))+
ylab("PR-AUC")+
xlab("Datas")+
ylim(0, 1)
ggsave("Figure/PR_AUC_Result.pdf", p, device = "pdf")
p