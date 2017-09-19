#Loads your libraries ######################################################################
  library(readr)
  library(data.table)
  library(plyr)
  library(ggplot2)
  library(flexclust)
  library(reshape2)
  library(survival)
  library(Rcmdr)

#Sets the working directory ################################################################
  setwd("C:/User/.../Your_Data_Here/")

#Generate survival curves  #################################################################

  #Import survival data sets
    Surv_1 = read.csv("Experiment_ID_1 - Pop1 Survival.csv")
    Surv_2 = read.csv("Experiment_ID_2 - Pop1 Survival.csv")
    #...ibid for additional groups...#
    Surv_N = read.csv("Experiment_ID_N - Pop1 Survival.csv")
  
  #Add identifiers to each survival data frame
    Surv_1$Group = "Group 1"                                       
      Surv_1 = Surv_1[,c(8, 2:7)]                                    
  
    Surv_2$Group = "Group 2"                                       
      Surv_2 = Surv_2[,c(8, 2:7)]                                    
  
    #...ibid for additional groups...#
 
    Surv_N$Group = "Group N"                                       
      Surv_N = Surv_N[,c(8, 2:7)]                                    
  
  #Combine survival data frames into a single data frame   
    Survival = rbind(Surv_1,Surv_2,...,Surv_N)                        
      FT = #Time Between Frames#
    
    Survival$Time = FT*Survival$Duration
  
  #Code for survival curve figure
    Survival$SurvObj = with(Survival, Surv(Time, Status == 1))
    km = survfit(SurvObj ~ Group, data = Survival, conf.type = "log-log")
    
    par(cex = 1.0, cex.main = 4.0, cex.axis = 1.5, cex.lab = 1.5)
    SurvCol = c("lawngreen", "grey1", "khaki")
    
    plot(km, xlab = "Time (Hours)", ylab = "Survival (%)",
      col = SurvCol, bg = "gray75", 
      main = "Pop1 Cell Survival",
      lty = 1:3, lwd = 3,
      yscale = 100, xscale = 60)
    
    legend("topright", legend = unique(Survival$Group),
      col = SurvCol, horiz=FALSE,
      bty='n', lty = 1:3,
      cex = 1.5, lwd = 3)    
  
#Generate figures for Pop1-Pop2 interaction persistance ####################################
  
  #Import survival data sets
    Pers_1 = read.csv("Experiment_ID_1 - Pop1-Pop2 Interactions.csv")
    Pers_2 = read.csv("Experiment_ID_2 - Pop1-Pop2 Interactions.csv")
    #...ibid for additional groups...#
    Pers_N = read.csv("Experiment_ID_N - Pop1-Pop2 Interactions.csv")
  
  #Add identifiers to each survival data frame
  
  Pers_1$Group = "Group 1"                                      
    Pers_1 = Pers_1[,c(5, 2:4)]                                     
    Pers_1 = aggregate(Freq_Pop2 ~ Group + Pop1_ID,
                data = Pers_1, simplify = FALSE, mean)               
  
  Pers_2$Group = "Group 2"                                      
    Pers_2 = Pers_2[,c(5, 2:4)]                                     
    Pers_2 = aggregate(Freq_Pop2 ~ Group + Pop1_ID,
                      data = Pers_2, simplify = FALSE, mean)       
  
  #...ibid for additional groups...#

  Pers_N$Group = "Group N"                                      
    Pers_N = Pers_N[,c(5, 2:4)]                                   
    Pers_N = aggregate(Freq_Pop2 ~ Group + Pop1_ID,
                        data = Pers_N, simplify = FALSE, mean)     
  
  #Combine Persistance data frames into a single data frame   
  Persistence = rbind(Pers_1,Pers_2,...,Pers_N)                      
  colnames(Persistence) = c("Group", "Pop1_ID","Persistence")      
  
  PNumb = as.numeric(Persistence$Persistence)
  Persistence$Persistence = PNumb                                   
  rm(PNumb)
  
  #Code for Persistance figure
    par(cex = 1.0, cex.main = 4.0, cex.axis = 1.5, cex.lab = 1.5)
    BoxCol = c("lawngreen", "grey75", "khaki")
    boxplot(Persistence ~ Group, data = Persistence,
          ylab = "Average Number of Interactions/Frame/Pop1 Cell",
          main = "Persistence of Pop1-Pop2 Contact",
          col = BoxCol)
  
  
#Generate figures for time to Cell Event from first Pop2 contact#############################
  
  #Import survival data sets
    T2D_1 = read.csv("Experiment_ID_1 - Pop1 Time to Cell Event from First Contact.csv")
    T2D_2 = read.csv("Experiment_ID_2 - Pop1 Time to Cell Event from First Contact.csv")
    #...ibid for additional groups...#
    T2D_N = read.csv("Experiment_ID_N - Pop1 Time to Cell Event from First Contact.csv")
  
  #Add identifiers to each survival data frame
    T2D_1$Group = "Group 1"                                        
      T2D_1 = T2D_1[,c(4, 2:3)]                                 
    
    T2D_2$Group = "Group 2"                                        
      T2D_2 = T2D_2[,c(4, 2:3)]                                 
    
    #...ibid for additional groups...# 

    T2D_N$Group = "Group N"                                        
      T2D_N = T2D_N[,c(4, 2:3)]                                 
  
  #Combine survival data frames into a single data frame   
    Time2Death = rbind(T2D_1,T2D_2,...,T2D_N)                    
  
  #Code for Time to Cell Death figure
      par(cex = 1.0, cex.main = 3.0, cex.axis = 1.5, cex.lab = 1.5)
      BoxCol = c("lawngreen", "grey75", "khaki")
    
    boxplot(Tum_T2Death ~ Group, data = Time2Death,
            ylab = "Time (Minutes)",
            main = "Time from First Contact to Pop1 Cell Death",
            col = BoxCol, xscale = 3)
  
#Generate figures for Total Number of Pop2 events to Cell Event #############################
  
  #Import survival data sets
    Xaction_1 = read.csv("Experiment_ID_1 - Pop1-Pop2 Interactions to Cell Event.csv")
    Xaction_2 = read.csv("Experiment_ID_2 - Pop1-Pop2 Interactions to Cell Event.csv")
    #...ibid for additional groups...#
    Xaction_N = read.csv("Experiment_ID_N - Pop1-Pop2 Interactions to Cell Event.csv")
  
  #Add identifiers to each survival data frame
  
    Xaction_1$Group = "Group 1"                                        
      Xaction_1 = Xaction_1[,c(4, 2:3)]                                 
    
    Xaction_2$Group = "Group 2"                                        
      Xaction_2 = Xaction_2[,c(4, 2:3)]                                 
    
    #...ibid for additional groups...#

    Xaction_N$Group = "Group N"                                        
      Xaction_N = Xaction_N[,c(4, 2:3)]                                 
  
  #Combine survival data frames into a single data frame   
  TTXaction = rbind(Xaction_1,Xaction_2,...,Xaction_N)                    
  
  #Code for Pop2 Events to Pop1 Death figure
  par(cex = 1.0, cex.main = 4.0, cex.axis = 1.5, cex.lab = 1.5)
  BoxCol = c("lawngreen", "grey75", "khaki")
  boxplot(Freq_Pop2 ~ Group, data = TTXaction,
          ylab = "Frequency of Interaction Events",
          main = "Interactions to Pop1 Cell Death",
          col = BoxCol)