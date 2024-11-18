library(readxl)
source("./mims_gas_functions[69].r")

####4-6 is offered as an alternative to 9-11 if you want to run an entire folder worth of MIMS files
#excelFolder <- "~/Dropbox/MIMS/Kate/FMP"
#excelFolder <- "../data/MIMS/Oct24/"
#rawFolder <- "../data/MIMS/raw_outputs/" # Add in a raw excel in the style of Gabi's raw excel sheet
#allfiles <- list.files(path = excelFolder, recursive = TRUE, pattern = "xlsx$", full.names = TRUE)

###Or, for a single target excel sheet instead of doing an entire folder
rawFolder <- "../data/MIMS/Oct24/"
allfiles <- "../data/MIMS/Oct24/CalculatedWaterColumn_outputs.xlsx" #Target file
#excelFolder <- "~/TestingMIMS"
#***or*** define excelFolder and then read in allfiles as I do above (lines 7 and 9)

O2.ArSat <- oxy_sat_field/Ar_sat_field
N2.ArSat <- N2_sat_field/Ar_sat_field


#Here, you define the columns of interest
tarCol <- c("O2/Ar", "N2/Ar") #Need to match your excel file names
satCol <- c("O2.ArSat", "N2.ArSat") #You can rename these

#Formats and reads in the raw files in the best manner for the gather_data function
read_raw_file <- function(filePath){
  nCols <- max(count.fields(filePath, sep = ','))
  firstRead <- read.csv(filePath, header = FALSE, col.names = paste0("V",seq_len(nCols)), 
                        fill = TRUE)
  truthVector <- cumsum(complete.cases(firstRead)) != 0
  skipRows <- length(truthVector[truthVector == FALSE])
  
  filedf <- read.csv(filePath, skip = skipRows+1)
  range <- c(min(filedf$Index), max(filedf$Index))
  return(list(filedf, range))
}

#Averages data, either with the raw or just the excel document
gather_data <- function(excelfile, folder_with_raw_files){
  name <- read_excel(excelfile)
  rundate <- format(as.Date(name$rundate[1]), "%Y_%m_%d")
  name <- name[!is.na(name$SampleID),]
  name_range <- c(min(name$Index), max(name$Index))
  
  #Read in the raw values, starting with file "*a.csv"
  #Allow for a "*b.csv" file to exist, check if the ranges match it if they fail to match the a
  #Need to build in allowance for partial matches
  techniq <- 1
  if(file.exists(paste0(folder_with_raw_files, "/", rundate, "a.csv"))){
    pars <- read_raw_file(paste0(folder_with_raw_files, "/", rundate, "a.csv"))
    rawname <- pars[[1]]
    raw_range <- pars[[2]]
    if (raw_range[1] <= name_range[1] && raw_range[2] >= name_range[2]){
      print(paste0("For ", excelfile, ", all data will come from the raw file a."))
      techniq <- 2
    }
    }else if(file.exists(paste0(folder_with_raw_files, "/", rundate, "b.csv"))){
      pars <- read_raw_file(paste0(folder_with_raw_files, "/", rundate, "b.csv"))
      rawname <- pars[[1]]
      raw_range <- pars[[2]]
      if (raw_range[1] <= name_range[1] && raw_range[2] >= name_range[2]){
        print(paste0("For ", excelfile, ", all data will come from the raw file b."))
        techniq <- 3
      }
    } else {print(paste0("For ", excelfile, ", all data will come from the excel file."))}
  
  datedf <- data.frame()
  if (techniq == 2 | techniq == 3){
    #Iterate through each sample
    for (j in unique(name$Samp)){
      sub <- name[(name$Samp == j),]
      indexes <- c(min(sub$Index), max(sub$Index))
      data <- rawname[(rawname$Index >= indexes[1]) & (rawname$Index <= indexes[2]),]
      data <- data[ , -which(names(data) %in% c("Index"))]
      data$Time <- as.numeric(as.POSIXct(data$Time, format = "%m/%d/%Y %I:%M:%S %p"))
      
      #Sneak in the inHg to mmHg correction
      newdat <- data.frame(lapply(colMeans(data), type.convert), stringsAsFactors=FALSE)
      newdat$Time <- as.POSIXct(newdat$Time, origin = "1970-01-01")
      temp <- mean(sub$Temp)
      press <- mean(sub$Pressure) * 25.4
      metadat <- data.frame("Samp" = sub[["Samp"]][1], "SampleID" = sub[["SampleID"]][1], 
                            "Pressure" = press, "Temp" = temp, "Calibnum" = sub[["Calibnum"]][1],
                            "Depth" = sub[["Depth"]][1], 
                            "WatDens" = watdens(temp), "OSat" = osat1(temp, press), 
                            "NSat" = nsat(temp, press), "ArSat" = arsat(temp, press), 
                            "O2.ArSat" = osat1(temp, press)/arsat(temp, press),
                            "N2.ArSat" = nsat(temp, press)/arsat(temp, press))
      
      metadat$Calibnum[is.na(metadat$Calibnum)] <- 
        as.numeric(as.character(sub[["Sampnum"]][1][is.na(sub[["Calibnum"]][1])]))
      
      tempdf <- merge(metadat, newdat)
      
      datedf <- rbind(datedf, tempdf)
    }
    
  } else(
    #Iterate through each sample
    for (j in unique(name$Samp)){
      sub <- name[(name$Samp == j),]
      colnames <- names(sub)
      start <- which(colnames == "Time")
      end <- which(colnames == "Sampleset") - 1
      data <- sub[,start:end]
      data$Time <- as.numeric(as.POSIXct(data$Time, format = "%m/%d/%Y %H:%M:%S"))
      
      #Sneak in the inHg to mmHg correction
      newdat <- data.frame(lapply(colMeans(data), type.convert), stringsAsFactors=FALSE)
      newdat$Time <- as.POSIXct(newdat$Time, origin = "1970-01-01")
      temp <- mean(sub$Temp)
      press <- mean(sub$Pressure) * 25.4
      metadat <- data.frame("Samp" = sub[["Samp"]][1], "SampleID" = sub[["SampleID"]][1], 
                            "Pressure" = press, "Temp" = temp, "Calibnum" = sub[["Calibnum"]][1],
                            "Depth" = sub[["Depth"]][1], 
                            "WatDens" = watdens(temp), "OSat" = osat1(temp, press), 
                            "NSat" = nsat(temp, press), "ArSat" = arsat(temp, press), 
                            "O2.ArSat" = osat1(temp, press)/arsat(temp, press),
                            "N2.ArSat" = nsat(temp, press)/arsat(temp, press))
      
      metadat$Calibnum[is.na(metadat$Calibnum)] <- 
        as.numeric(as.character(sub[["Sampnum"]][1][is.na(sub[["Calibnum"]][1])]))
      
      tempdf <- merge(metadat, newdat)
      
      datedf <- rbind(datedf, tempdf)
    }
  )
  
  return(datedf)
}

#Computes the concentration ratios of whatever target columns you've defined in line 16
getRatio <- function(df, targCol, satCol){
  newcolname <- paste0(targCol, ".Conc")
  
  #Pull out differences of each calibnum
  for (i in 1:length(unique(df$Calibnum))){
    sub <- df[df$Calibnum == i | df$Calibnum == i+1,]
    sub$TimeDiff <- sub$Time - min(sub$Time)
    
    calibs <- sub[sub$Depth == "Stdl" | sub$Depth == "Stdh",]
    #Linear between low and high temps
    yvals <- c(0, calibs[[satCol]], 0)
    xvals <- c(0, calibs[[targCol]], 0)
    tvals <- c(min(calibs$TimeDiff), calibs$TimeDiff, max(calibs$TimeDiff))
    offsetfun <- lm(yvals ~ xvals + tvals)
    coeff1 <- coef(summary(offsetfun))[1]
    coeff2 <- coef(summary(offsetfun))[2]
    coeff3 <- coef(summary(offsetfun))[3]
    
    #Multiply o2:ar of samples by offset for its temp
    sub[[newcolname]] <- (coeff1 + sub[[targCol]] * coeff2 + sub$TimeDiff * coeff3)
    
    df[[newcolname]][!is.na(base::match(df$Samp, sub$Samp))]<- sub[[newcolname]]
  }
  
  return(df)
}

#Run the above functions appropriately for the type of input you're using
find_Real_Ratios <- function(allfiles, rawFolder){
  #Allowances for running multiple excel files at once
  if (length(allfiles) > 1){
    allAveragedDfs <- list() 
    for (excelFile in allfiles){
      averaged_data <- gather_data(excelFile, rawFolder)
      for (i in 1:length(targets)){
        averaged_data <- getRatio(averaged_data, targets[i], saturations[i])
      }
      #Create a list of the output data that is named by the excel file
      allAveragedDfs[[excelFile]] <- averaged_data
    }
    return(allAveragedDfs)
  }else { #Just one excel file
    excelFile <- allfiles
    averaged_data <- gather_data(excelFile, rawFolder)
    for (i in 1:length(targets)){
      averaged_data <- getRatio(averaged_data, targets[i], saturations[i])
    }
    return(averaged_data)
  }
}

#Finally, actually run all of ^^that^^
averaged <- find_Real_Ratios(allfiles, rawFolder)
