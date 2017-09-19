#Load your libraries
    library(readr)
    library(data.table)
    library(plyr)
    library(flexclust)
    library(reshape2)

##Set working directory to location of your data#############################################################################
    setwd("C:/User/.../Your_Data_Here")

#Import Data
    Pop1 = read_delim("Experiment_ID_N - Pop1 Spots in tracks statistics.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    Pop2 <- read_delim("Experiment_ID_N - Pop2 Spots in tracks statistics.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
  

#Create a subset of data with only the frame data
    Pop1_FRAME = Pop1$FRAME
    Pop2_FRAME = Pop2$FRAME

#Load the appropriate package for use in the following code
    require(flexclust)
    require(reshape2)

#Create an empty list that we will fill with the appropriate vectors
    dist_list = list()

#The meat of the code which generates distances between Pop1 and Pop2
    for (i in 0:max(Pop1_FRAME))
        {X1 = as.matrix(Pop1$POSITION_X[which(Pop1$FRAME == i)])                                                   
         Y1 = as.matrix(Pop1$POSITION_Y[which(Pop1$FRAME == i)])                                                   
            Pop1_XY = cbind(X1, Y1)                                                                                   

         X2 = as.matrix(Pop2$POSITION_X[which(Pop2$FRAME == i)])                                                   
         Y2 = as.matrix(Pop2$POSITION_Y[which(Pop2$FRAME == i)])                                                   
            Pop2_XY = cbind(X2, Y2)                                                                                   

         mat_dist = dist2(Pop2_XY, Pop1_XY)                                                                        
            rownames(mat_dist) = Pop2$TRACK_ID[which(Pop2$FRAME == i)]                                               
            colnames(mat_dist) = Pop1$TRACK_ID[which(Pop1$FRAME == i)]                                                

         mat_dist = setNames(melt(mat_dist), c("Pop1_ID", "Pop2_ID","Distance"))                                   
         mat_dist$Frame = Pop1$FRAME[which(Pop1$FRAME == i)]                                                       

         dist_list[[i+1]] <- mat_dist                                                                              
        }

#Combine and export the list of distances data frames into a single data frame
    require(data.table)
    Distances = rbindlist(dist_list)                                                                               

    write.csv(Distances, file = "Experiment_ID_N - Pop2-Pop1 Distances.csv")
  
    rm(dist_list)