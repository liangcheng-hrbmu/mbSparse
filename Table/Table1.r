CRC_name <- read.csv("../data/species_associated_with_D011236.csv")$scientific.name
CRC_name <- gsub(" ", "_", data)

da_name_PRJNA397219 <- read.csv("../data/PRJNA397219/edgeR_DA.csv", header=TRUE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,1])
da_name_PRJDB4176 <- read.csv("../data/PRJDB4176/edgeR_DA.csv", header=TRUE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,1])
da_name_PRJEB7774 <- read.csv("../data/PRJEB7774/edgeR_DA.csv", header=TRUE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,1])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("../data/PRJNA397219/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("../data/PRJDB4176/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("../data/PRJEB7774/PRJEB7774.csv")$clade_name
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

da_name_PRJNA397219 <- read.csv("../data/PRJNA397219/edgeR_mbSparse_DA.csv", header=TRUE)
da_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", da_name_PRJNA397219[,1])
da_name_PRJDB4176 <- read.csv("../data/PRJDB4176/edgeR_mbSparse_DA.csv", header=TRUE)
da_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", da_name_PRJDB4176[,1])
da_name_PRJEB7774 <- read.csv("../data/PRJEB7774/edgeR_mbSparse_DA.csv", header=TRUE)
da_name_PRJEB7774 <- sub(".*s__(.*)", "\\1", da_name_PRJEB7774[,1])
da_name <- union(da_name_PRJNA397219, da_name_PRJDB4176)
da_name <- union(da_name, da_name_PRJEB7774)

all_name_PRJNA397219 <- read.csv("../data/PRJNA397219/PRJNA397219.csv")$clade_name
all_name_PRJNA397219 <- sub(".*s__(.*)", "\\1", all_name_PRJNA397219)
all_name_PRJDB4176 <- read.csv("../data/PRJDB4176/PRJDB4176.csv")$clade_name
all_name_PRJDB4176 <- sub(".*s__(.*)", "\\1", all_name_PRJDB4176)
all_name_PRJEB7774 <- read.csv("../data/PRJEB7774/PRJEB7774.csv")$clade_name
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