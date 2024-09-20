#set remove_rate
remove_rate = 0.7
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("../data/PRJEB7774/PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')



#set remove_rate
remove_rate = 0.4
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("../data/PRJEB7774/PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')


#set remove_rate
remove_rate = 0.1
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("../data/PRJEB7774/PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')


#set remove_rate
remove_rate = 0.2
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("../data/PRJEB7774/PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')



#set remove_rate
remove_rate = 0.3
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')


#set remove_rate
remove_rate = 0.5
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("../data/PRJEB7774/PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')


#set remove_rate
remove_rate = 0.6
cat("remove rate:", remove_rate * 100, '\n')
zeros_identified_array <- c()
impute_before_pearson_array <- c()
impute_after_pearson_array <- c()
for(j in 1 : 10) {
    OTU <- read.csv("../data/PRJEB7774/PRJEB7774.csv", row.names="X")
    OTU <- OTU / colSums(OTU) * 10^6 
    OTU <- log10(OTU+1)
    OTU_ori_bool <- OTU != 0
    OTU_downsample <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100, "_", j, ".csv"), row.names="X")
    OTU_impute <- read.csv(paste0("../data/PRJEB7774/PRJEB7774_remove_", remove_rate * 100,"_impute_normalize_", j, ".csv"), row.names="X")
    
    OTU_impute[OTU_impute < 0.001] <- 0

    OTU_downsample_bool <- !(OTU_downsample != 0)
    OTU_impute_bool <- OTU_impute != 0
    OTU_impute_remove_downsample_bool <- OTU_impute_bool & OTU_downsample_bool 
    OTU_ori_remove_downsample_bool <- OTU_ori_bool & OTU_downsample_bool
    OTU_impute_correct_bool <- OTU_impute_remove_downsample_bool & OTU_ori_remove_downsample_bool
    correct <- sum(OTU_impute_correct_bool, na.rm = TRUE) / sum(OTU_ori_remove_downsample_bool, na.rm = TRUE)  
    zeros_identified_array <- c(zeros_identified_array, correct)
    cat("downsampling zeros identified:", correct, '\n')

    OTU <- t(OTU)
    OTU_downsample <- t(OTU_downsample)
    OTU_impute <- t(OTU_impute)
    downsample_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU)[1]) {
        downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean + cor(OTU[i,], OTU_downsample[i,])
    }
    downsample_pearson_coefficient_mean <- downsample_pearson_coefficient_mean / dim(OTU)[1]
    impute_before_pearson_array <- c(impute_before_pearson_array, downsample_pearson_coefficient_mean)
    cat("impute_before_pearson:", downsample_pearson_coefficient_mean, '\n')

    OTU_ori_impute <- t(read.csv("../data/PRJEB7774/PRJEB7774_impute_normalize.csv", row.names="X"))
    impute_pearson_coefficient_mean = 0
    for (i in 1:dim(OTU_ori_impute)[1]) {
        impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean + cor(OTU_ori_impute[i,], OTU_impute[i,])
    }
    impute_pearson_coefficient_mean <- impute_pearson_coefficient_mean / dim(OTU_ori_impute)[1]
    impute_after_pearson_array <- c(impute_after_pearson_array, impute_pearson_coefficient_mean)
    cat("impute_after_pearson:", impute_pearson_coefficient_mean, '\n')
}
cat("zeros_identified_mean:",mean(zeros_identified_array), '\n')
cat("zeros_identified_sd:",sd(zeros_identified_array), '\n')
cat("impute_before_pearson_mean:", mean(impute_before_pearson_array), '\n')
cat("impute_before_pearson_sd:", sd(impute_before_pearson_array), '\n')
cat("impute_after_pearson_mean:", mean(impute_after_pearson_array), '\n')
cat("impute_after_pearson_sd:", sd(impute_after_pearson_array), '\n')
