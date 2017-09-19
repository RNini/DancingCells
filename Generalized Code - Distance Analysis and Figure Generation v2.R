#Loads your libraries ########################################################################
  library(readr)
  library(data.table)
  library(plyr)
  library(ggplot2)
  library(flexclust)
  library(reshape2)

#Sets the working directory ##################################################################
  setwd("C:/User/Your_Data_Here/")

#Imports Data previously generated in R ######################################################
  {require(readr)
    Distances = read.csv("./Experiment_ID_N - Pop1-Pop2 Distances.csv")
    Distances = Distances[, c(2:5)]

    Tum_Tracks = read_delim("Experiment_ID_N - Pop1 Spots in tracks statistics.txt",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
      Tum_Tracks = Tum_Tracks[, c(2:ncol(Tum_Tracks))]
    
    Pop2_Tracks = read_delim("Experiment_ID_N - Pop2 Spots in tracks statistics.txt",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
      Pop2_Tracks = Pop2_Tracks[, c(2:ncol(Pop2_Tracks))]
    
    print(max(Distances$Pop1_ID))
    print(max(Distances$Pop2_ID))
  }

#Subsets distance data to include only data points below a determined threshold ##############
  Distances = Distances[which(Distances$Distance < 35),]

#Isolates nearest neighbors (Pop2 Cells with respect to Pop1 Cells) and calculates average distances ##############################################################

  #Create an empty list that we will fill with the appropriate vectors
    Neighbors_list = list()
  
  #Calculate the number of "nearest neighbors" for each Pop1 cell in each frame
    for (i in 0:max(Distances$Pop1_ID))
      {Tum = Distances[which(Distances$Pop1_ID == i),]                             
       Tum = Tum[which(Tum$Distance < 35),]                                         
      
        if (nrow(Tum) > 0) {
      
          Persistance = count(Tum, "Pop2_ID")                                           
            colnames(Persistance) = c("Pop2_ID", "Freq_FRAME")                            
          
          Avg_Dist = as.data.frame(aggregate(Tum$Distance, list(Tum$Pop2_ID), mean))    
            colnames(Avg_Dist) = c("Pop2_ID", "Mean_Distance")                            
          
          SD_Dist = as.data.frame(aggregate(Tum$Distance, list(Tum$Pop2_ID), sd))       
            colnames(SD_Dist) = c("Pop2_ID", "SD_Distance")                               
          
          A = merge(Persistance, Avg_Dist, by = c("Pop2_ID"))                           
          Persistance = merge(A, SD_Dist, by = c("Pop2_ID"))
            Persistance$Pop1_ID = rep(mean(Tum$Pop1_ID), nrow(Persistance))            
          
          Neighbors_list[[i+1]] <- Persistance                                         
          
          rm(A, Avg_Dist, Persistance, SD_Dist, Tum)                                   
        }
      
        else {}
      
      rm(i)
      }
  
  #Concatenate the data frames into a single data frame
    {require(data.table)
    Neighbors = rbindlist(Neighbors_list)                                         
    
    Neighbors = Neighbors[,c(5,1,2,3,4)]                                          
    
    Neighbors = transform(Neighbors, Pop2_ID = as.numeric(Pop2_ID),
                          Pop1_ID = as.numeric(Pop1_ID))                        
    
    rm(A, Avg_Dist, Persistance, SD_Dist, Tum)                                    
    
    }

#Counts the frequency of Pop1-Pop2 Interactions #################

  #Create an empty list that we will fill with the appropriate vectors
    N_Pop2_list = list()

  #Generate a second data frame containing the number of Pop2 Cells  
  #associated with each Pop1 cell, as a function of time (FRAME)
    for (i in 0:max(Distances$Pop1_ID))
      {Tum = Distances[which(Distances$Pop1_ID == i),]                             
      Tum = Tum[which(Tum$Distance < 35),]                                         
      
      if (nrow(Tum) > 0) {
      
          N_Pop2 = as.data.frame(table(Tum$Frame))                                     
            colnames(N_Pop2) = c("Frame", "Freq_Pop2")                                   
          
          N_Pop2$Pop1_ID = rep(mean(Tum$Pop1_ID), nrow(N_Pop2))                      
          
          N_Pop2_list[[i+1]] <- N_Pop2                                                 
          }
      
      else {}
      
      rm(N_Pop2, Tum)
      }

#Concatenate the data frames into a single data frame
  {require(data.table)
   Pop2_Number = rbindlist(N_Pop2_list)                                         
    
   Pop2_Number = Pop2_Number[,c(3,2,1)]                                          
    
   rm(N_Pop2_list)                                                             
  }

#Tabulates the frames of first contact and final contact b/w Pop1 and Pop2 Cells ##################

  #Create an empty list that we will fill with the appropriate vectors
    T_Pop2_list = list()

  #Calculate the first/final frames of contact, and contact duration of each Pop1-Pop2 interaction
    for (i in 0:max(Distances$Pop1_ID))
      {Tum = Distances[which(Distances$Pop1_ID == i),]                             
       Tum = Tum[which(Tum$Distance < 35),]                                          
      
      if (nrow(Tum) > 0) {
        
        F_alpha = as.data.frame(aggregate(Tum$Frame, by = list(Tum$Pop2_ID), min))   
          colnames(F_alpha) = c("Pop2_ID", "First_Frame")                              
        
        F_omega = as.data.frame(aggregate(Tum$Frame, by = list(Tum$Pop2_ID), max))   
          colnames(F_omega) = c("Pop2_ID", "Final_Frame")                              
      
        F_ao = merge(F_alpha, F_omega, by = c("Pop2_ID"))                            
          F_ao$Pop1_ID = rep(mean(Tum$Pop1_ID), nrow(F_alpha))                      
        
        T_Pop2_list[[i+1]] <- F_ao                                                  
        
        rm(F_alpha, F_omega, Tum)                                                   
      }
      
      else {}
      
      rm(i,F_ao)
      }

  #Concatenate the data frames into a single data frame
    {require(data.table)
      T_Pop2 = rbindlist(T_Pop2_list)                                             
      
      T_Pop2 = T_Pop2[,c(4,1,2,3)]                                                
      
      rm(T_Pop2_list)                                                            
    }

  # Combines all data sets above into a single data frame for analysis ##########################
    {Neighbors = merge(Neighbors, T_Pop2, by = c("Pop2_ID", "Pop1_ID"))             
    Neighbors = Neighbors[,c(1,2,6,7,3,4,5)]                                      
    }
  
# Generates a data frame with Cell Event information ########################################## 
    
    #Imports and Concatenates your data
    
    #Set the working directory
    setwd("E:/Ryan/Data/For Cell Event Pipeline Development/Data/")
    
    #Import the appropriate data sets
    Pop1_CE = read_csv("Cell Event Tracks/Pop1 CE - Spots in tracks statistics.csv")
    Pop1_CE = Pop1_CE[,c(-1,-2,-5,-8,-11,-12,-13,-21,-22)]
    
    BG1 = read_csv("BackgroundR1.csv")
    BG2 = read_csv("BackgroundR2.csv")
    BG3 = read_csv("BackgroundR3.csv")
    BG4 = read_csv("BackgroundR4.csv")
    
    #Concatenate the background data
    Names = c("FRAME","REGION_AREA","MEAN_INT","SD_INT","MIN_INT","MAX_INT")
    colnames(BG1) = Names
    colnames(BG2) = Names
    colnames(BG3) = Names
    colnames(BG4) = Names
    
    BG1$Region = 1
    BG2$Region = 2
    BG3$Region = 3
    BG4$Region = 4
    
    BG = rbind(BG1,BG2,BG3,BG4)
    BG$VAR_INT = BG$SD_INT^2                            #We are calculating the variance here so we
    #can calculate the mean SD amongst ROIs
    rm(Names,BG1,BG2,BG3,BG4)
    
    # Determines the threshold for Cell-Event+ signals
    
    #Calculate the average background intensity for each frame
    Names = c("FRAME","MEAN_INT","SD_INT")
    BG_Int = list()
    
    for (i in min(BG$FRAME):max(BG$FRAME)) {
      
      BG_Sub = BG[which(BG$FRAME == i),]
      
      Int_AVG = mean(BG_Sub$MEAN_INT)                   #Calculates the average intensity of ROIs, per frame
      Int_SD = sqrt(mean(BG_Sub$VAR_INT))               #Calculates the SD by take the sqrt of mean variance
      
      BGs = as.data.frame(cbind(i, Int_AVG, Int_SD))
      
      colnames(BGs) = Names
      
      BG_Int[[i]] <- BGs
    }
    
    rm(BG_Sub,Int_AVG,Int_SD, BGs)
    
    require(data.table)
    Frame_BG = rbindlist(BG_Int)
    
    #Determine the threshold for a Cell-Event Positive signal
    Frame_BG$THRESHOLD = Frame_BG$MEAN_INT + (2 * Frame_BG$SD_INT)
    
    #Bind the above data frame to the Cell-Event Tracking data frame
    Pop1_CE = merge(Pop1_CE,Frame_BG, by = "FRAME")
    Pop1_CE = Pop1_CE[,c(-14,-15)]
    
    #Check TrackMate signal intensity data to identify Cell-Event+ Spots
    Pop1_CE$ADJ_THRESHOLD = Pop1_CE$THRESHOLD + 32768             #Adjusts the scale of the intensity
    #(Corrects for an ImageJ artifact)
    Pop1_CE$CELL_EVENT = 0                                         #Assigns a default value to Cell-Event column
    
    for (i in 1:nrow(Pop1_CE)) {                                   #Checks each row to see if mean intensity
      #exceeds the calculated threshold
      Tum = Pop1_CE[i,]
      
      if (Tum$MEAN_INTENSITY >= Tum$ADJ_THRESHOLD) {              #For intensities greater than the 
        Pop1_CE$CELL_EVENT[[i]] = 1                                #threshold, we replace the default
      }                                                             #value with 1
      
    } 
    
    #Summarize Cell Event Data to first frame for cell event signal 
    Tum_Cell_Event_List = list()
    
    for (i in min(Pop1_CE$TRACK_ID):max(Pop1_CE$TRACK_ID)) {
      
      Tum = Pop1_CE[which(Pop1_CE$TRACK_ID == i),c(1,3,16)]     #Reduce the data frame to include only necessary variables
      
      if (nrow(Tum) > 0) {                                        #Skips missing tracks
        CE = as.data.frame(i)
        colnames(CE) = "Pop1_ID"
        
        if (mean(Tum$CELL_EVENT) == 0) {                            #Sets the Cell Event frame to -1 (placeholder) if the
          CE$Cell_Event_Frame = -1                                    #the associated track is Cell-Event-
        }
        
        else {
          
          Tum = Tum[which(Tum$CELL_EVENT == 1),]
          
          CE$Cell_Event_Frame = min(Tum$FRAME)                      #Picks the minimum frame
        }
        
        Tum_Cell_Event_List[[i+1]] <- CE                            #Binds the frame and Pop1_ID dataframe to a list
      }
    }
    
    {require(data.table)
      Tum_Cell_Event = rbindlist(Tum_Cell_Event_List)             #Binds the list of data frames into a new data frame
    }
    
    
  # Adds above to previously generated data frames for analysis ##########################

    Neighbors = merge(Neighbors, Tum_Cell_Event, by = "Pop1_ID")                     
      Neighbors = Neighbors[,c(1,2,3,4,5,8,6,7)]
    
    Pop2_Number = merge(Pop2_Number, Tum_Cell_Event, by = "Pop1_ID")
 
#Subsets and processes the previously generated Pop1 data as a function
#of cell Event information ###################################################################
  
  #Split the data sets by Cell_Event_Frame
    Neighbors_Live = Neighbors[Neighbors$Cell_Event_Frame == -1,]                   
      Neighbors_Live = Neighbors_Live[Neighbors_Live$Mean_Distance > 0]                           
    
    Neighbors_Dead = Neighbors[Neighbors$Cell_Event_Frame > -1,]                    
      Neighbors_Dead = Neighbors_Dead[Neighbors_Dead$Mean_Distance > 0]               

    Pop2_NumberD = Pop2_Number[Pop2_Number$Cell_Event_Frame > -1,]                     
      Pop2_NumberD$Frame = as.numeric(Pop2_NumberD$Frame)                               
    
  #Calculate the time to death
    
    #Create an empty list that we will fill with the appropriate vectors
      Tum_T2Death_List = list()
    
    #Calculate the number of frames from first contact to Cell Event
      for (i in min(Pop2_NumberD$Pop1_ID):max(Pop2_NumberD$Pop1_ID)) {
        Tum = Pop2_NumberD[which(Pop2_NumberD == i),]
      
        if (nrow(Tum) > 0) {
          Tum = Tum[which.min(Tum$Frame),]                                            
            Tum$Tum_T2Death = (Tum$Cell_Event_Frame - Tum$Frame)                              
          
          Tum_T2Death_List[[i+1]] = Tum                                                   
          
          rm(Tum)
          }
      
        else{}
      }

    #Concatenate the above outputs to a single data frame
      {require(data.table)
       Tum_T2Death = as.data.frame(rbindlist(Tum_T2Death_List))                      
        
       Tum_T2Death = Tum_T2Death[which(Tum_T2Death$Tum_T2Death > 0),c(1,5)]          
        
       rm(Tum_T2Death_List)                                                          
      }

#Determines the number of Pop2 contacts with each Pop1 to Cell Event ##############################
     
    #Subset the data to only include interactions before Cell Event
      Pop2_Number = transform(Pop2_Number, Frame = as.numeric(Frame))                    
      Tum = Pop2_Number[which(Pop2_Number$Cell_Event_Frame > Pop2_Number$Frame),]        
    
    #Tabulate the number of interactions from first contact to Cell Event
      Pop2_Number_CE = aggregate(Freq_Pop2 ~ Pop1_ID, Tum, sum)                   
      
#Generates random Pop2 Cell Event information for workflow development #############################
    
    #Create an empty list that we will fill with the appropriate vectors
      Pop2_Cell_Event_List = list()
    
    #Generate random frames for Cell Event for Pop2 Cells
      for (i in 0:max(Distances$Pop2_ID)) {
        
        Event = rbinom(1,1,.3)                                                        
        
        if (Event == 1) {
          
          Pop2 = Pop2_Tracks[which(Pop2_Tracks$TRACK_ID == i),"FRAME"]                   
          
          Cell_Event_Frame = sample(min(Tum$FRAME):max(Tum$FRAME),
                                    1, replace = TRUE)                                
        }
        
        else {
          Cell_Event_Frame = -1                                                        
        }
        
        Pop2_ID = as.data.frame(Distances[which(Distances$Pop2_ID == i),"Pop2_ID"])      
          Pop2_ID = Pop2_ID[1,]
        
        CE = as.data.frame(cbind(Pop2_ID,Cell_Event_Frame))                            
        
        Pop2_Cell_Event_List[[i+1]] = CE                                               
        
        rm(CE, Pop2_ID, Event,Tum)                                                     
      }
      
      {require(data.table)
       Pop2_Cell_Event = rbindlist(Pop2_Cell_Event_List)                               
       rm(Pop2_Cell_Event_List)
      }
    
#Generates a data frame for Kaplan Meier Survival Curve analysis #################################
    
  #Combine Pop1_ID, Frame information, and Cell_Event_Frame information    
    Surv = Tum_Tracks[,c("TRACK_ID","FRAME")]                                                 
      colnames(Surv) = c("Pop1_ID","Frame")                                        
    
    Surv = merge(Surv, Tum_Cell_Event, by = "Pop1_ID")                           
    
  #Eliminate spurious information
  
    Surv$Cell_Event_Frame[which(Surv$Cell_Event_Frame == -1)] = NA                
      Surv1 = Surv[which(is.na(Surv$Cell_Event_Frame)),]                            
    
    Surv = Surv[which(Surv$Frame <= Surv$Cell_Event_Frame),]                      
    
    Surv = rbind(Surv,Surv1)                                                      
    
    rm(Surv1)                                                                     
  
  #Determine the duration of tracking
    
    #Create an empty list to fill with the appropriate data frames
      Duration_List = list()
  
    #Calculate first and final frames of tracking
      for (i in 0:max(Surv$Pop1_ID))
      {Tum = Surv[which(Surv$Pop1_ID == i),]                                         
      
      if (nrow(Tum) > 0) {
        
        F_alpha = as.data.frame(aggregate(Tum$Frame, by = list(Tum$Pop1_ID), min))   
        colnames(F_alpha) = c("Pop1_ID", "First_Frame")                              
        
        F_omega = as.data.frame(aggregate(Tum$Frame, by = list(Tum$Pop1_ID), max))   
        colnames(F_omega) = c("Pop1_ID", "Final_Frame")                              
        
        F_ao = merge(F_alpha, F_omega, by = c("Pop1_ID"))                            
        
        Duration_List[[i+1]] <- F_ao                                                   
        
        rm(F_alpha, F_omega, F_ao, Tum)                                                 
      }
      
      else {}
      
      rm(i)
      }
  
    #Aggregate the entries into a single data frame and merge with Surv
      Duration = rbindlist(Duration_List)                                             
        Duration$Duration = (Duration$Final_Frame - Duration$First_Frame)               
      
      Surv = merge(Surv, Duration, by = "Pop1_ID")                                   
        Surv = Surv[,c("Pop1_ID", "First_Frame","Final_Frame",
                      "Cell_Event_Frame","Duration")]                                  
      
      Surv = unique(Surv[,])                                                          
    
  #Create Censor column
    Status = is.na(Surv$Cell_Event_Frame)                                             
      Status = !Status                                                                  
    
    Surv = cbind(Surv, as.numeric(Status))                                            
      colnames(Surv) = c("Pop1_ID", "First_Frame","Final_Frame",
                        "Cell_Event_Frame","Duration", "Status")                       
         
#Exports the generated data sets for later analysis ###############################################

  setwd("./Analyzed")    
    
  #These files correspond to Pop1-specific data sets generated above
    #Field 1
      write.csv(Surv, file = "Experiment_ID_1 - Pop1 Survival.csv")                                                
      write.csv(Tum_T2Death, file = "Experiment_ID_1 - Pop1 Time to Cell Event from First Contact.csv")            
      write.csv(Pop2_Number_CE, file = "Experiment_ID_1 - Pop1-Pop2 Interactions to Cell Event.csv")           
      write.csv(Pop2_Number, file = "Experiment_ID_1 - Pop1-Pop2 Interactions.csv")                                  
    
    #Field 2
      write.csv(Surv, file = "Experiment_ID_2 - Pop1 Survival.csv")                                                
      write.csv(Tum_T2Death, file = "Experiment_ID_2 - Pop1 Time to Cell Event from First Contact.csv")            
      write.csv(Pop2_Number_CE, file = "Experiment_ID_2 - Pop1-Pop2 Interactions to Cell Event.csv")           
      write.csv(Pop2_Number, file = "Experiment_ID_2 - Pop1-Pop2 Interactions.csv")                                  
    
    #...ibid for additional groups...#

    #Field N
      write.csv(Surv, file = "Experiment_ID_N - Pop1 Survival.csv")                                                
      write.csv(Tum_T2Death, file = "Experiment_ID_N - Pop1 Time to Cell Event from First Contact.csv")            
      write.csv(Pop2_Number_CE, file = "Experiment_ID_N - Pop1-Pop2 Interactions to Cell Event.csv")           
      write.csv(Pop2_Number, file = "Experiment_ID_N - Pop1-Pop2 Interactions.csv")                                  