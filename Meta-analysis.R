#setwd("src/Input/")

library(tools)
library(data.table)
library(plyr)
library(xlsx)

#define function to automatically import input files, assign filename as variable name
sep_importer <- function(filename) {
  assign(gsub(" ", "_", file_path_sans_ext(filename)), read.csv(filename, check.names = FALSE, stringsAsFactors = FALSE), envir = .GlobalEnv)
}

#import files
sapply(dir("src/Input/"), sep_importer)


#define study merger
study_merger <- function(main, add) {
  merge(main, add, by = "Gene names", all.x = TRUE)
}
#merge data.frames per study for normalization
rogers <- study_merger(`2007_Rogers_IgG_Early`, `2007_Rogers_IgG_Late`)
campbell <- study_merger(`2012_Campbell_Latex_Early`, `2012_Campbell_Latex_Late`)
goyette <- study_merger(`2012_Goyette_Latex_Early`, `2012_Goyette_Latex_Late`)
dill <- study_merger(`2015_Dill_Avidin_Early`, `2015_Dill_Avidin_Late`)
dill <- study_merger(dill, `2015_Dill_Calreticulin_Early`)
dill <- study_merger(dill, `2015_Dill_Calreticulin_Late`)
dill <- study_merger(dill, `2015_Dill_Complement_Early`)
dill <- study_merger(dill, `2015_Dill_Complement_Late`)
dill <- study_merger(dill, `2015_Dill_IgG_Early`)
dill <- study_merger(dill, `2015_Dill_IgG_Late`)
dill <- study_merger(dill, `2015_Dill_LPS_Early`)
dill <- study_merger(dill, `2015_Dill_LPS_Late`)
dill <- study_merger(dill, `2015_Dill_Mannan_Early`)
dill <- study_merger(dill, `2015_Dill_Mannan_Late`)
dill <- study_merger(dill, `2015_Dill_Phosphatidylserine_Early`)
dill <- study_merger(dill, `2015_Dill_Phosphatidylserine_Late`)

#set NAs to 0
rogers[is.na(rogers)] <- 0
campbell[is.na(campbell)] <- 0
goyette[is.na(goyette)] <- 0
dill[is.na(dill)] <- 0

#assign first column of each table to become row names, this makes processing easier
rogers2 <- rogers[,-1]
rownames(rogers2) <- rogers[,1]
campbell2 <- campbell[,-1]
rownames(campbell2) <- campbell[,1]
goyette2 <- goyette[,-1]
rownames(goyette2) <- goyette[,1]
dill2 <- dill[,-1]
rownames(dill2) <- dill[,1]

#start by making sure each table runs from 0 to max
rogers_zero <- rogers2 + abs(min(rogers2))
campbell_zero <- campbell2 + abs(min(campbell2))
goyette_zero <- goyette2 + abs(min(goyette2))
dill_zero <- dill2 + abs(min(dill2))

#define function for normalization
normalizer <- function(x) {
  x / max(x)
}

#normalize all tables from 0 to 1
rogers_scaled <- as.data.frame(t(apply(rogers_zero, 1, normalizer)))
goyette_scaled <- as.data.frame(t(apply(goyette_zero, 1, normalizer)))
campbell_scaled <- as.data.frame(t(apply(campbell_zero, 1, normalizer)))
dill_scaled <- as.data.frame(t(apply(dill_zero, 1, normalizer)))

#set NAs to 0
rogers_scaled[is.na(rogers_scaled)] <- 0
campbell_scaled[is.na(campbell_scaled)] <- 0
goyette_scaled[is.na(goyette_scaled)] <- 0
dill_scaled[is.na(dill_scaled)] <- 0


#place row names back into first column to facilitate merging
setDT(rogers_scaled, keep.rownames = TRUE)[]
rogers_scaled <- rename(rogers_scaled, c("rn" = "Gene names"))
setDT(goyette_scaled, keep.rownames = TRUE)[]
goyette_scaled <- rename(goyette_scaled, c("rn" = "Gene names"))
setDT(campbell_scaled, keep.rownames = TRUE)[]
campbell_scaled <- rename(campbell_scaled, c("rn" = "Gene names"))
setDT(dill_scaled, keep.rownames = TRUE)[]
dill_scaled <- rename(dill_scaled, c("rn" = "Gene names"))

#merge tables to produce final table, longest table first
pre_final <- study_merger(dill_scaled, rogers_scaled)
pre_final <- study_merger(pre_final, goyette_scaled)
pre_final <- study_merger(pre_final, campbell_scaled)

#final set NAs to 0
pre_final[is.na(pre_final)] <- 0
pre_final <- as.data.frame(pre_final)

#assign first column to be row names
final <- pre_final[,-1]
rownames(final) <- pre_final[,1]

#export as xlsx, further analysis possible in Excel
write.xlsx(final, "Final meta-analysis scaled by row.xlsx")
