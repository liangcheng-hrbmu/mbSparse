CRC_name <- read.csv("Table1/species_associated_with_D011236.csv")$scientific.name
CRC_name <- gsub(" ", "_", CRC_name)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_Deseq2_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_Deseq2_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_Deseq2_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                c("identified as DA", "not identified as DA")))
print(table)

result <- fisher.test(table)
print(result)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_MBAE_Deseq2_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_MBAE_Deseq2_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_MBAE_Deseq2_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

impute_table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                       dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                       c("identified as DA", "not identified as DA")))
print(impute_table)

result <- fisher.test(impute_table)
print(result)




CRC_name <- read.csv("Table1/species_associated_with_D011236.csv")$scientific.name
CRC_name <- gsub(" ", "_", CRC_name)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_wilcoxon_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_wilcoxon_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_wilcoxon_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                c("identified as DA", "not identified as DA")))
print(table)

# 运行Fisher确切性检验
result <- fisher.test(table)
print(result)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_MBAE_wilcoxon_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_MBAE_wilcoxon_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_MBAE_wilcoxon_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

impute_table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                       dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                       c("identified as DA", "not identified as DA")))
print(impute_table)

result <- fisher.test(impute_table)
print(result)



CRC_name <- read.csv("Table1/species_associated_with_D011236.csv")$scientific.name
CRC_name <- gsub(" ", "_", CRC_name)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_ALDEx2_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_ALDEx2_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_ALDEx2_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                c("identified as DA", "not identified as DA")))
print(table)

result <- fisher.test(table)
print(result)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_MBAE_ALDEx2_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_MBAE_ALDEx2_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_MBAE_ALDEx2_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

impute_table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                       dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                       c("identified as DA", "not identified as DA")))
print(impute_table)

result <- fisher.test(impute_table)
print(result)

CRC_name <- read.csv("Table1/species_associated_with_D011236.csv")$scientific.name
CRC_name <- gsub(" ", "_", CRC_name)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_lefse_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_lefse_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_lefse_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                c("identified as DA", "not identified as DA")))
print(table)

result <- fisher.test(table)
print(result)

da_name_PRJNA397219 <- read.delim("PRJNA397219_data/PRJNA397219_MBAE_lefse_taxa.txt", header = FALSE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,2])
da_name_PRJDB4176 <- read.delim("PRJDB4176_data/PRJDB4176_MBAE_lefse_taxa.txt", header = FALSE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,2])
da_name_PRJEB7774 <-read.delim("PRJEB7774_data/PRJEB7774_MBAE_lefse_taxa.txt", header = FALSE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,2])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("PRJNA397219_data/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("PRJDB4176_data/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("PRJEB7774_data/PRJEB7774.csv")$clade_name
all_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", all_name_PRJEB7774)
all_name <- union(all_name_PRJNA397219, all_name_PRJDB4176)
all_name <- union(all_name, all_name_PRJEB7774)

matches <- match(da_name, CRC_name)
counts <- table(matches)
a = dim(counts)[1]

matches <- match(all_name, CRC_name)
counts <- table(matches)
b = dim(counts)[1] - a

c = length(da_name) - a 

d =  length(all_name) - a - b - c

impute_table <- matrix(c(a, b, c, d), ncol = 2, byrow = TRUE,
                       dimnames = list(c("annotated by the term", "not annotated by the term"), 
                                       c("identified as DA", "not identified as DA")))
print(impute_table)

result <- fisher.test(impute_table)
print(result)