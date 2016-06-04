#' Get ClaNA Group info
#'
#' Extract which group each point belongs to
#' @param Build_out 
#' @keywords ClaNC
#' @export
#' @examples
#' get_ClaNC_group()
get_ClaNC_group <- function(Build_out){
	idx <- apply(Build_out$cntrds, 1, function(x) which.max(abs(x-Mode(x))))
	summary <- sapply(1:max(idx), function(x) names(idx[idx==x]))
	colnames(summary) <- paste0("Group", 1:ncol(summary))
	return(summary)
}


#' ClaNC run 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
run_ClaNC <- function(data, groups, ColSideColors, active.features=20, est.num=20, select.features=10, file.prefix="ClaNC", skip.est=F){
# data - expression data, row(genes), column(samples)
# groups - assigned groups for each column(samples) in data 
# active.features - how many features we want to randomly select to test cv.errors 
# est.num - how many times we want to run the test for cv.errors 
# select.features - how many features(genes) we want to select from each group for the final results 
# skip.est - do you want to skip estimate cv errors and go straight to extract features from each group

## Load the example data sets in the "data" directory.  These are simulated 
## microarray data sets, consisting of 10 arrays each from 4 classes; the 
## first 10 arrays are from class 1, the second 10 from class 2, etc.  And 
## each array has 1500 genes on it.  The "train" data set is intended for 
## training purposes; these are the data you would use in training and 
## building your classifier.  The "test" data are intended for testing a 
## classifier, using data that was not used in the classifier building 
## process.  You may or may not have test data.  If you intend to use all of 
## your samples to build the classifier and rely on cross-validation to 
## estimate the resulting error rates, you need not have *any* test data.  

# Y_train = read.delim("data/trainExample.txt", header = T)
# Y_test = read.delim("data/testExample.txt", header = T)
# rownames(Y_train) = Y_train[, 1]
# rownames(Y_test) = Y_test[, 1]
# Y_train = as.matrix(Y_train[, -1])
# Y_test = as.matrix(Y_test[, -1])
   Y_train <- as.matrix(data)

## Create vector for class membership.  For these example data, we just 
## need 10 ones followed by 10 twos, then 10 threes, then 10 fours.
# id = rep(1:4, each = 10)
#Levels: Basal Her2 LumA LumB Normal
      idd <- as.numeric(factor(groups))
## The gene names are just the row names from the data matrices, the same 
## for both the example train and test data sets.
# gene_names = rownames(Y_train)
      gene_names <- rownames(Y_train)
## We'll just set the class names to be the numbers 1, 2, 3, 4.
# class_names = 1:4
      class_names <- 1:length(unique(idd))
## Specify how many active features(genes) to select from each group 
      active.features <- active.features 

## Now ready to use cross-validation to estimate the error rates for 
## classifiers of different sizes (different numbers of genes used in 
## building the classifier).
## THEN View the estimated error rates associated with different feature-set 
## sizes.  By default, cvClanc() will assess the classifiers built using 1, 
## 2, 3, ..., 10 features.  In this example, an estimated 100% accuracy is 
## attained with as few as 5 features per class.
      if (skip.est!=T){
         cv_total <- NULL
         cv_out <- NULL
         for (i in 1:est.num){
               cv_out = cvClanc(Y_train, idd, active=1:active.features, prior="class")
				# plot(1:active.features, cv_out$overallErrors, type = "l", lwd = 2, col = "blue", xlab = "Number of features", ylab = "CV error")
                cv_total = rbind(cv_total, cbind(1:active.features, cv_out$overallErrors))
				# qplot(1:active.features, cv_out$overallErrors, type = "l", lwd = 2, col = "blue", xlab = "Number of features", ylab = "CV error")
            }
         colnames(cv_total) <- c("NumFeatures", "CV_Error")
            file.prefix <- paste0(file.prefix, ".estnum", est.num, ".act", active.features)
            cv_total <- as.data.frame(cv_total)
            p <- ggplot(data=cv_total, aes(x=factor(NumFeatures), y=CV_Error))
            p <- p + geom_boxplot(aes(fill = factor(NumFeatures)))
            p <- p + xlab("Number of Features") 
            pdf(paste0(file.prefix, ".pdf"), height=10, width=16)
            print(p)
            dev.off()
      }

## Assuming we're happy with using 5 features per class as the basis for 
## our model, we can now train and build our final classifier.
   file.prefix <- paste0(file.prefix, ".sel", select.features)
      train_out = trainClanc(Y_train, idd, gene_names)
# print(train_out$classCntrds)
      build_out = buildClanc(Y_train, idd, class_names, train_out, active = select.features)
      write.table(build_out$geneNames, paste0(file.prefix, ".genelist.txt"), row.names = F, quote=F)

      clanc.group <- get_ClaNC_group(build_out)
      write.csv(clanc.group, paste0(file.prefix, ".group_features.csv"), quote=F, row.names=F)

### Expression profile for genes that overlapped with PAM50 
      ClaNC.genes <- build_out$geneNames
# plot_CalNC_pam50(build_out, data, groups=groups, file.prefix=paste0(file.prefix, ".overlap.PAM50.geneExp"))
### Heatmap expression profile for extracted genes(features)
      sub.data <- extract_data_by_genelist(data, ClaNC.genes)
      myHeatmap(data=sub.data, color = color, ColSideColors = ColSideColors, file.prefix=paste0(file.prefix, ".extracted.heatmap"), save.image=T, scale="row")

      return(build_out)
## And if we have test data, we can test the classifier on them.  In our 
## example, we make one misclassification (for one of the samples in class 4.
#test_out = testClanc(Y_test, idd, gene_names, build_out)
      }


