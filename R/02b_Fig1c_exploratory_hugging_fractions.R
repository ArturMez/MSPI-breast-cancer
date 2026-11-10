########################
## Figure 1c: Exploratory cohort hugging fraction analysis
########################
## Loads raw proximity data, classifies cell types, and analyzes
## the proportion of cells hugging specific macrophage subtypes (M1, M2)

library(data.table)
library(dplyr)
#options(digits=8)


## ── Configure data paths ────────────────────────────────────────────────
## Adjust these to point to your local exploratory cohort data directory
data_dir <- "./data/exploratory/"


## Load raw hugging and T cell data
data_huggingmfs <- data.table::fread(file.path(data_dir, "hmf.tsv"))
data_cd3ki67 <- data.table::fread(file.path(data_dir, "cd3ki67.tsv"))


colnames(data_cd3ki67)
colnames(data_huggingmfs)

head(data_cd3ki67)
head(data_huggingmfs)

data_cd3ki67$Name
data_cd3ki67$`object info (ALL) - Envelope`

data_cd3ki67$MergerID <- paste(data_cd3ki67$Name, "___", data_cd3ki67$`object info (ALL) - Envelope`)
data_huggingmfs$MergerID <- paste(data_huggingmfs$Name, "___", data_huggingmfs$`object info - Envelope`)

unique(data_cd3ki67$`object info (ALL) - LabelName`)
unique(data_huggingmfs$`object info (ALL) - LabelName`)
table(data_huggingmfs$`object info (ALL) - LabelName`)

Merged <- inner_join(data_cd3ki67, data_huggingmfs, by="MergerID")
unique(Merged$`object info (ALL) - LabelName`)

head(Merged)
colnames(Merged)
unique(Merged$`object info (ALL) - ROIType`)


Merged$Cell_classes_by_markers <- Merged$`object info (ALL) - LabelName`
Merged <- Merged[,c(-11)]
table(Merged$Cell_classes_by_markers)
Merged$Cells_hugging_classes <- Merged$`object info - LabelName`
Merged <- Merged[,c(-27)]
table(Merged$Cells_hugging_classes)


colnames(Merged) 
Merged <- Merged[,c(-17, -18, -19, -20, -21, -34,-35)]




data <- Merged
head(data) 
colnames(data) 
str(data) 
## convert to OD
data[ , c(25:29)] <- -log10(data[, c(25:29)]/255)

names(data)[6] <- "Sample.Name"
names(data)[14] <- "Cell X Position"
names(data)[15] <- "Cell Y Position"   





###############################################################################################
####        analysis
##########################


#######################
#### re-shape

d<- data

colnames(d)
unique(d$Sample.Name)

XXX<- d %>% group_by(Cells_hugging_classes) %>% summarise(Percentage=n()/nrow(.))
pie(XXX$Percentage, labels = XXX$Cells_hugging_classes, cex=0.8)


colnames(d)
unique(d$Cell_classes_by_markers)
unique(d$Cells_hugging_classes)

d$Cell_class <- gsub('HUGGED CD8S', 'CD8', d$Cells_hugging_classes)
d$Cell_class <- gsub('HUGGING MFS', 'MFS', d$Cell_class)
unique(d$Cell_class)
table(d$Cells_hugging_classes)

unique(d$Cell_classes_by_markers)

d$Cell_classes_by_markers   <- gsub("CD8\\+", "CD8", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub(" KI67\\+", "_Ki67", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub(" KI6BY", "_Ki67 BY", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub(" KIBY", "_Ki67 BY", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub("KI67\\+", "_Ki67", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub("Ki67BY", "Ki67 BY", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub("KI67BY", "Ki67 BY", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub("CD136", "CD163", d$Cell_classes_by_markers)
d$Cell_classes_by_markers   <- gsub("BY", "by", d$Cell_classes_by_markers)

d$cell_groups1 <- d$Cell_classes_by_markers

unique(d$cell_groups1)


d$cell_groups1   <- gsub("CD163\\+CD206\\+CD68\\-", "Myeloid", d$cell_groups1)
d$cell_groups1   <- gsub("CD163\\+CD206\\-CD68\\-", "Myeloid", d$cell_groups1)
d$cell_groups1   <- gsub("CD163\\-CD206\\+CD68\\-", "Myeloid", d$cell_groups1)


d$cell_groups1   <- gsub("CD68LOW", "_l", d$cell_groups1)
d$cell_groups1   <- gsub("CD68\\+", "_h", d$cell_groups1)

d$cell_groups1   <- gsub("_l", "", d$cell_groups1)
d$cell_groups1   <- gsub("_h", "", d$cell_groups1)

unique(d$cell_groups1)
d$cell_groups1   <- gsub("CD163\\-CD206\\-", "M1", d$cell_groups1)


d$cell_groups1   <- gsub("CD163\\+CD206\\-", "OTHERS", d$cell_groups1)
d$cell_groups1   <- gsub("CD163\\+CD206\\+", "M2", d$cell_groups1)

d$cell_groups1   <- gsub("CD163\\-CD206\\+", "OTHERS", d$cell_groups1)



d$cell_groups1   <- gsub("OTHERS.*", "OTHERS", d$cell_groups1)
#d$cell_groups1   <- gsub(".*OTHERS", "OTHERS", d$cell_groups1)



unique(d$cell_groups1)

#d <- d[ which( d$cell_groups1 !="OTHERS" ), ]

## make cell subset
d$cell_subsets <- d$cell_groups1
unique(d$cell_subsets)
d$cell_subsets   <- gsub("_Ki67.*", "", d$cell_subsets)
d$cell_subsets   <- gsub(" by.*", "", d$cell_subsets)

unique(d$cell_subsets)
table(d$cell_groups1)
unique(d$cell_groups1)

        ## make Ki67 groups
        d$cell_groups2<- d$cell_groups1
        d$cell_groups2   <- gsub("_Ki67.*", "", d$cell_groups2)
        unique(d$cell_groups2)


          

        
        
        
        
        
        # M perxpective
        
        dat2 <- d
        dat2 <- dat2[ which( dat2$cell_subsets !="OTHERS" ), ]
        dat2 <- dat2[ which( dat2$cell_subsets !="Myeloid" ), ]
        unique(dat2$cell_subsets)
        dat2 <- dat2[ which( dat2$Cells_hugging_classes !="CD8" ), ]
        dat2 <- dat2[ which( dat2$Cells_hugging_classes !="HUGGED CD8S" ), ]
       
        # M1
        XXX <- dat2[ which( dat2$cell_subsets =="M1" ), ]

        XXX<- XXX %>% group_by(Cells_hugging_classes) %>% summarise(Percentage=n()/nrow(.))
        pie(XXX$Percentage, labels = XXX$Cells_hugging_classes, cex=0.8)
        XXX
        
        
        # M2
        XXX <- dat2[ which( dat2$cell_subsets =="M2" ), ]

        XXX<- XXX %>% group_by(Cells_hugging_classes) %>% summarise(Percentage=n()/nrow(.))
        pie(XXX$Percentage, labels = XXX$Cells_hugging_classes, cex=0.8)
        XXX
        
        
        
        
        # T cell perxpective
        
        dat2 <- d
        dat2 <- dat2[ which( dat2$cell_subsets !="OTHERS" ), ]
        dat2 <- dat2[ which( dat2$cell_subsets !="Myeloid" ), ]
        unique(dat2$cell_subsets)
        dat2 <- dat2[ which( dat2$Cells_hugging_classes !="MFS" ), ]
        dat2 <- dat2[ which( dat2$Cells_hugging_classes !="HUGGING MFS" ), ]
        unique(dat2$Cells_hugging_classes)
        unique(dat2$cell_groups2)
        
        # M1
        XXX <- dat2
        XXX <- XXX[ which(XXX$cell_groups2 !="CD8 by Myeloid" ), ]
        XXX <- XXX[ which(XXX$cell_groups2 !="CD8 by OTHERS" ), ]
        XXX <- XXX[ which(XXX$cell_groups2 !="CD8 by M2" ), ]
        
        
        XXX<- XXX %>% group_by(cell_groups2) %>% summarise(Percentage=n()/nrow(.))
        pie(XXX$Percentage, labels = XXX$Cells_hugging_classes, cex=0.8)
        XXX
        
        # M2
        XXX <- dat2
        XXX <- XXX[ which(XXX$cell_groups2 !="CD8 by Myeloid" ), ]
        XXX <- XXX[ which(XXX$cell_groups2 !="CD8 by OTHERS" ), ]
        XXX <- XXX[ which(XXX$cell_groups2 !="CD8 by M1" ), ]
        
        
        XXX<- XXX %>% group_by(cell_groups2) %>% summarise(Percentage=n()/nrow(.))
        pie(XXX$Percentage, labels = XXX$Cells_hugging_classes, cex=0.8)
        XXX
        
        
        