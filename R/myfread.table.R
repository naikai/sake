#' A Cat Function
#'
#' This function allows you to read data with fread (faster)
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# plot tSNE result 
myfread.table <- function(filepath, check.platform=T, header=T, sep="\t", detect.file.ext=T){
   require(data.table)
   require(tools)
   require(magrittr)

   ext <- file_ext(filepath)

   if(detect.file.ext){
      if (ext=="csv"){
        sep=","
      }else if (ext=="out" || ext=="tsv" || ext=="txt"){
          sep="\t"
      }else{
          warning("File format doesn't support, please try again")
          return(NULL)
      }
   }

   header <- read.table(filepath, nrows = 1, header = FALSE, sep=sep, stringsAsFactors = FALSE)
   first_data <- read.table(filepath, nrows=1, sep=sep, skip=1)
   if(length(header)==length(first_data)){
      cols <- c("character", rep("numeric", length(header)-1))
   }else if(length(header)==length(first_data)-1){
      cols <- c("character", rep("numeric", length(header)))
   }
   rawdata <- fread(filepath, header=F, sep=sep, skip=1, colClasses=cols)

    ### Again. Add more checking in case there are duplicate rownames
    # read in the first column and check, if there is duplicated rownames (mostly gene names)
    # then send out a warning and then make it unique
    if(sum(duplicated(rawdata$V1))>0){
      warning("There are duplicated rownames in your data")
      warning("Please double check your gene count table")
      rawdata$V1 <- make.names(rawdata$V1, unique=TRUE)
    }

   ### Add more checking in case there are duplicated column names
   # make.names(names, unique=TRUE)
   rawdata %<>% setDF %>% set_rownames(.$V1) %>%
                  '['(colnames(.) != "V1") #%>% as.numeric

   # data doesn't have colnames for first row (rownames)
   if(length(header) == dim(rawdata)[2]){
      # colnames(rawdata) <- unlist(header)
      colnames(rawdata) <- make.names(unlist(header), unique=TRUE)
   }else if (length(header) == dim(rawdata)[2] + 1){
      # colnames(rawdata) <- unlist(header)[-1]
      colnames(rawdata) <- make.names(unlist(header)[-1], unique=TRUE)
   }

   # Add checking data platform
   if(check.platform){
      rawdata <- detect_genecnt_platform(rawdata)
   }

   return(rawdata)
}
