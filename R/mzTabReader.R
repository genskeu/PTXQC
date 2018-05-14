########################################################################################################################
#'
#' Class which holds the data from a .mzTab.
#' 
#'
#' @export
#'

read.MzTab <- function(file){
  mztab <- mzTab$new()
  mztab$read.mzTab(file)
}
 
#' Memberfunction read.mzTab() takes path to .mzTab file
#' Data from .mzTab sections is stored in seperate data.frames
#'
#' Usage: foo <- mzTab$new()
#'        foo <- foo$read.mzTab("file.mztab")
#'        metadata_from_file.mzTab <- foo$MTD 
#'        
#' @import methods 
#' @importFrom rowr cbind.fill
#' @importFrom plyr rename
#' @exportClass mzTab
#' 
#' @export mzTab
#'  

mzTab = setRefClass("mzTab",
                    fields = list(ALL="data.frame",
                                  MTD="data.frame", ##
                                  PRT="data.frame", ##
                                  PEP="data.frame", ##
                                  PSM="data.frame", ##
                                  SML="data.frame", ##
                                  rawFileMapping="data.frame" ##
                    ),
                    methods = list(
                      initialize = function(){
                        .self$ALL = data.frame();
                        .self$MTD = data.frame();
                        .self$PRT = data.frame();
                        .self$PEP = data.frame();
                        .self$PSM = data.frame();
                        .self$SML = data.frame();
                        .self$rawFileMapping <- data.frame();
                      },
                      #function to read mzTab---------------------------------------------------------------
                      read.mzTab = function(file){
                        
                        #get max number of columns
                        nr <- max(count.fields(file, sep = "\t", quote = "'"))
                        cn <- paste0("c", 1:nr)
                        data <- read.table(file=file, col.names = cn , sep="\t", fill=TRUE, na.strings=c("null","NA", ""), stringsAsFactors = FALSE)
                        
                        .self$ALL <- data
                        .self$MTD <- data[which(data[,1]=="MTD"),]
                        
                        .self$PRT <- data[which(data[,1]=="PRT"),]
                        if (length(.self$PRT) != 0) {
                          colnames(.self$PRT) <- as.character(data[which(data[,1]=="PRH"),])
                        }
                        
                        .self$PEP <- data[which(data[,1]=="PEP"),]
                        if (length(.self$PEP) != 0) {
                          colnames(.self$PEP) <- as.character(data[which(data[,1]=="PEH"),])
                        }
                        
                        .self$PSM <- data[which(data[,1]=="PSM"),]
                        if (length(.self$PSM) != 0) {
                          colnames(.self$PSM) <- as.character(data[which(data[,1]=="PSH"),])
                        }
                        
                        .self$SML <- data[which(data[,1]=="SML"),]
                        if (length(.self$SML) != 0) {
                          colnames(.self$SML) <- as.character(data[which(data[,1]=="SMH"),])
                        }
                        
                      }
                    )
)
########################################################################################################################
#'
#' Memberfunction to generate data.frame similar to MaxQuants evidence.txt from an mzTab object
#'
#' Usage: foo <- mzTab$new()
#'        foo <- foo$read.mzTab("file.mztab")
#'        metadata_from_file.mzTab <- foo$MTD 
#'    
#' @importFrom plyr rename
#' @name mzTab$get_df_evd

mzTab$methods(
  get_df_evd = function(.self, add_fs_col) {
    # the base for the df_evd is the PSM section of the mzTab file
    # or rbind(.self$PSM, .self$PRT etc.)
    
    #create df with number of rows = number of rows fom mzTab 
    df_evd <- data.frame()
    
    # add columns:
          # "sequence" = "modified.sequence",
          # "accession" = "proteins",
          # "PSM_ID" = "id",
          # "raw_source_file" = "raw.file",
          # "exp_mass_to_charge" = "m.z",
          # "retention_time" = "retention.time",
          # "search_engine_score[1]" = "score",
          # "opt_Match_Time_Difference" = "match.time.difference",
          # "retention_time_window" = "retention.length"

    if("sequence" %in% colnames(.self$PSM)){
      df_evd <- cbind.fill(df_evd, .self$PSM, fill = NA)
      }else if ("sequence" %in% colnames(.self$PEP)) {
        df_evd <- cbind.fill(df_evd, .self$PEP, fill = NA)
      }

    #iff match.time.difference column isn't present in df_evd:
    #add match.time.difference column filled with NA
    
    if (!("match.time.difference" %in% colnames(df_evd))) {
      df_evd[ ,"match.time.difference"] <- NA 
    }
    
    # # helper function
    # pasteUnique <- function(vector) {
    #   return(paste(unique(vector), collapse = ";"))
    # }
    # 
    # # merge rows with same PSM_ID value
    # # if PSM_ID contains no data use sequence instead
    # # brauche ich gar nicht?
    # if ("PSM_ID" %in% colnames(.self$PSM)) {
    #   if (all(is.na(.self$PSM$PSM_ID))) {
    #     df_evd <- .self$PSM
    #     df_evd <-
    #       aggregate(. ~ PSM_ID,
    #                 pasteUnique,
    #                 na.action = na.pass,
    #                 data = df_evd[which(colnames(df_evd) != "")])
    #   } else{
    #     df_evd <-
    #       aggregate(. ~ sequence,
    #                 pasteUnique,
    #                 na.action = na.pass,
    #                 data = df_evd[which(colnames(df_evd) != "")])
    #   }
    # } else{
    #   warning("Dataframe for EVD Metrics could not be created. Column named PSM_ID was not found.")
    #}
    
    #### rename the columns (column names determined by metric functions) -------------------------
    #Mapping anpassen! spectraref != raw.file, raw_source_file soll final zu opt_raw_source_file
     df_evd <- 
       rename(df_evd,
         c(
           "sequence" = "modified.sequence",
           "accession" = "proteins",
            #"PSM_ID" = "id",
           "raw_source_file" = "raw.file",
           "exp_mass_to_charge" = "m.z",
           "retention_time" = "calibrated.retention.time",
           "original_retention_time" = "retention.time",
           "search_engine_score[1]" = "score",
           "opt_Match_Time_Difference" = "match.time.difference",
           "retention_time_window" = "retention.length",
           "opt_abundance1" = "intensity",
           "unique_id" = "id"
         )
       )
    
    #### set datatypes for cloumns ------------------------
     
     #possible values for opt_is_contaminent from PEP section: "+, -" -> PTXQC needs bool
     
     if("opt_isContaminent" %in% colnames(df_evd)){
     df_evd$opt_isContaminant = sapply(
       X = df_evd$opt_isContaminant,
       FUN = function(x){
         if(x=="+"){x = TRUE}
         else if(x=="-"){x = FALSE}
         else{x = NA}
         },
       USE.NAMES = TRUE
       )
     }
     
     try(df_evd$retention.time <-as.numeric(df_evd$retention.time))
     try(df_evd$calibrated.retention.time <-as.numeric(df_evd$calibrated.retention.time))
     try(df_evd$id <-as.numeric(df_evd$id))
     try(df_evd$m.z <-as.numeric(df_evd$m.z))
     try(df_evd$charge <-as.numeric(df_evd$charge))
     try(df_evd$score <-as.numeric(df_evd$score))
     try(df_evd$match.time.difference <-as.numeric(df_evd$match.time.difference))
     try(df_evd$retention.length <-as.numeric(df_evd$retention.length))
     try(df_evd$intensity <-as.numeric(df_evd$intensity))

    #### add "extra" columns needed for the evd metrics but not immediate present in the mzTabfile -----------
    
    #Die in openMS generierten mzTab verwenden die Protein description als accession
    #daher ist das mapping vorerst nicht möglich/nötig
     
    # if ("proteins" %in% colnames(df_evd))
    # {
    #   # protein name info is extracted from the PRT section of the mzTab file
    #   protein.names = sapply(
    #     X = .self$PSM$accession,
    #     FUN = function(accession) {
    #       return(.self$PRT$description[which(.self$PRT$accession == accession)[1]])
    #     },
    #     USE.NAMES = TRUE
    #   )
    #   # map protein names onto df_evd
    #   df_evd$protein.names = sapply(
    #     X = df_evd$proteins,
    #     FUN = function(protein) {
    #       paste0(protein.names[unlist(strsplit(protein, ";"))], collapse = ";")
    #     },
    #     USE.NAMES = FALSE
    #   )
    #   rm(protein.names)
    # }
    
    if (all(c("charge", "m.z") %in% colnames(df_evd)))
    {
      df_evd$mass = df_evd$charge * df_evd$m.z
    }
    
    if ("raw.file" %in% colnames(df_evd))
    {
      #if(df_$raw.file) falls raw.file einen Filepath entält 
      try(df_evd$raw.file <- gsub(pattern = "(.*\\/)",replacement = "",x = df_evd$raw.file))
      # df_evd$raw.file = gsub(pattern = ":.*$",replacement = "",x = df_evd$raw.file)
      # # test cases had file names with characters that cause proplems with regular expressions (e.g. [,(...)
      # df_evd$raw.file =  gsub(pattern = "\\[",replacement = " ",x = df_evd$raw.file)
      # df_evd$raw.file =  gsub(pattern = "\\]",replacement = " ",x = df_evd$raw.file)
      #
      df_evd$fc.raw.file = .self$get_file_mapping(add_fs_col,df_evd,.self)
      
    }
    
    ##### Gehacke ################################################
    #für mbrAlign mzTab, da in jeder Zeile ein charge gegeben ist. -> createReport() bricht bei qualUniform()
    df_evd <- subset(df_evd,!(is.na(df_evd["retention.time"])))
    df_evd[,"type"]<- "MULTI-MATCH"
    df_evd[,"retention.time.calibration"] <- df_evd$retention.time - df_evd$calibrated.retention.time
    #retention.time und len darf nicht NA sein -> if is.na(all(PEP$retention...)) -> rename columns
    return(df_evd)
  }
)

########################################################################################################################
#'
#' Memberfunctions
#'
#' Usage: foo <- mzTab$new()
#'        foo <- foo$read.mzTab("file.mztab")
#'        metadata_from_file.mzTab <- foo$MTD 
#'    
#' @importFrom plyr rename
#'  
#' @name mzTab$methods
mzTab$methods(
  #### function to generate summary dataframe (df_smy) from an mzTab object ---------------------
  get_df_smy = function(.self, add_fs_col) {
  },
  #### function to generate parAll dataframe (df_parAll) from an mzTab object ---------------------
  get_df_parAll = function(.self) {
  },
  #### function to generate protein groups dataframe (df_pg) from an mzTab object ---------------------
  df_pg = function(.self) {
  },
  #### function to generate msms dataframe (df_msms_s) from an mzTab object ---------------------
  df_msms_s = function(.self) {
  },
  #### function to generate msms dataframe (df_msmsScan_h) from an mzTab object ---------------------
  df_msmsScan_h = function(.self) {
  },
  
  get_file_mapping = function(add_fs_col, df, .self) {
    # code für Mapping aus dem MQReader (Leicht abgeändert)
    if (add_fs_col &
        "raw.file" %in% colnames(df))
    {
      ## check if we already have a mapping
      if (nrow(.self$rawFileMapping) == 0)
      {
        .self$rawFileMapping = getShortNames(unique(df$raw.file), add_fs_col)
        ## indicate to outside that a new table is ready
        ##.self$mapping.creation = .$getMappingCreation()['auto']
      }
      cat(paste0("Adding fc.raw.file column ..."))
      
      ## do the mapping
      df$fc.raw.file = as.factor(.self$rawFileMapping$to[match(df$raw.file, .self$rawFileMapping$from)])
      ## check for NA's
      if (any(is.na(df$fc.raw.file)))
      {
        ## if mapping is incomplete
        #missing = unique(df_evd$raw.file[is.na(df_evd$fc.raw.file)])
        # if (.$mapping.creation == .$getMappingCreation()['user'])
        # {
        #   ## the user has re-run MaxQuant with more Raw files,
        #   ## but the old _filename_sort.txt file was used to read the (now incomplete mapping)
        #   warning("Incomplete mapping file '", .$external.mapping.file, "'.\nAugmenting shortened Raw files:\n  " %+%
        #             paste(missing, collapse="\n  ", sep="") %+% ".\nEdit the table if necessary and re-run PTXQC.")
        #   ## augment
        #   addon = .$getShortNames(missing, add_fs_col, nrow(.$raw_file_mapping) + 1)
        #   .$raw_file_mapping = rbind(.$raw_file_mapping,
        #                              addon)
        # } else {
        #   stop("Hithero unknown Raw files: " %+% paste(missing, collapse=", ", sep="") %+% " occurred in file '" %+% file %+% "' which were not present in previous txt files.")
        # }
      }
      cat(paste0(" done\n"))
    }
    return(df$fc.raw.file)
  }
)
