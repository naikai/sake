library(shiny)
library(shinyBS)
library(ggplot2)
library(data.table)
library(magrittr)
library(corrplot)
library(biomaRt)
library(BiocParallel)
library(AnnotationHub)
library(RColorBrewer)
library(matrixStats)
library(tools)
library(snowfall)
library(NMF)
library(DT)
library(gage)
library(gageData)
library(pathview)

shinyServer(function(input, output, session) {
  fileinfo <- reactive({
    if (input$selectfile == 'upload'){
      inFile <- input$file1
      req(inFile)
      validate(
        need(!is.null(inFile), "Please upload a gene count data set")
      )
      filename <- basename(inFile$name)
      filepath <- inFile$datapath
    }else if (input$selectfile == "preload"){
      req(input$inputdata)
      validate(
        need(input$inputdata!="", "Please select a gene count data set")
      )
      filename <- basename(input$inputdata)
      filepath <- file.path(raw_fd, input$inputdata)
    }else if (input$selectfile == "saved"){
      inFile <- input$rda
      req(inFile)
      validate(
        need(!is.null(inFile), "Please upload a .rda file")
      )
      filename <- basename(inFile$name)
      filepath <- inFile$datapath
    }else{
      return()
    }

    fileinfo <- list()
    fileinfo[['name']] <- filename
    fileinfo[['path']] <- filepath
    return(fileinfo)
  })

  ### Output Prefix ###
  File <- reactive({
    if (input$selectfile == 'upload'){
      File <- input$file1$name
    }else if (input$selectfile == "preload"){
      File <- input$inputdata
    }else if (input$selectfile == "saved"){
      File <- input$rda$name
    }else{
      return()
    }
    return(File)
  })
  File_ext <- reactive({
    file_ext(File())
  })
  file_prefix <- reactive({
    if(input$selectfile == "saved"){
      gsub(paste0(".", File_ext()), "", File())
    }else{
      paste( gsub(paste0(".", File_ext()), "", File()),
             paste0("k", input$num_cluster),
             paste0(input$top.num,input$algorithm),
             paste0("nrun", input$nrun),
             paste0("predict_", input$predict),
             sep=".")
    }
  })

  rda <- reactive({
    req(fileinfo())
    req(input$selectfile == "saved")
    filepath <- fileinfo()[['path']]
    validate(
      need(File_ext() == "rda" || File_ext() == "RData",
           message = paste("file extension", File_ext(), "unknown, please try again"))
    )
    rda <- local({load(filepath); environment()})

    # start checking for all the required variables
    validate(
      need(!is.null(rda$rawdata), "Variable 'rawdata' does not exist in your .rda file, please check it and upload again") %then%
      need(nrow(rda$rawdata) > 0, "Variable 'rawdata' in your .rda file is empty, please check it and upload again") %then%
      need(!is.null(rda$nmfres), "Variable 'nmfres' does not exist in your .rda file, please check it and upload again")
      # need(!is.null(rda$tsne_2d), "Variable 'tsne_2d' does not exist in your .rda file, please check it and upload again") %then%
      # need(!is.null(rda$tsne_3d), "Variable 'tsne_3d' does not exist in your .rda file, please check it and upload again") #%then%
      # need(!is.null(rda$dds), "Variable 'dds' does not exist in your .rda file, please check it and upload again")
    )
    return(rda)
  })

  observeEvent(rda(), {
    tsne_2d$data <- rda()$tsne_2d
    tsne_3d$data <- rda()$tsne_3d
  })

  rawdata <- reactive({
    req(fileinfo())
    filepath <- fileinfo()[['path']]
    n <- 2
    withProgress(message = 'Loading data', value = 0, {
      incProgress(1/n, detail = "Usually takes ~10 seconds")
      if(input$selectfile == "saved"){
        rawdata <- rda()$rawdata
      }else{
        rawdata <- myfread.table(filepath, check.platform=T, sep=input$sep, detect.file.ext=FALSE)
      }
    })
    colnames(rawdata) <- gsub("\\.", "-", colnames(rawdata))
    return (rawdata)
  })

  #' Raw data
  output$rawdatatbl = DT::renderDataTable({
    rawdata <- rawdata()
    validate(
      need(ncol(rawdata)>1,
           "Please modify the parameters on the left and try again!") %then%
        need(!any(duplicated(t(rawdata))),
             "There are columns with the exact same value in your data, Please check and try again!")
    )
    withProgress(message = 'Displying datatable', value = NULL, {
      DT::datatable(head(rawdata, n=20), rownames= TRUE,
                    options = list(scroll = TRUE,
                                   scrollX = TRUE,
                                   scrollY = TRUE,
                                   pageLength = 8
                    )
      ) %>% formatRound(1:ncol(rawdata), 3)
    })
  }, server=TRUE)

  #' Info box
  output$nameBox <- renderInfoBox({
    infoBox(
      "FILENAME", fileinfo()[['name']], icon = icon("list"),
      color = "green"
    )
  })
  output$sampleBox <- renderInfoBox({
    rawdata <- rawdata()
    infoBox(
      "SAMPLE", ncol(rawdata), icon = icon("list"),
      color = "blue"
    )
  })
  output$featureBox <- renderInfoBox({
    rawdata <- rawdata()
    infoBox(
      "FEATURE", nrow(rawdata), icon = icon("list"),
      color = "yellow"
    )
  })

  # transform_data <- reactiveValues(data=NULL, rawdata=NULL)
  transform_data <- reactiveValues(data=data.frame(a=NULL, b=NULL),
                                   rawdata=data.frame(a=NULL, b=NULL))
  run_selsamp <- reactiveValues(go=FALSE)

  run_trans <- function(rawdata){
    transdata <- rawdata

    n <- 2
    withProgress(message = 'Transforming data', value = 0, {
      incProgress(1/n, detail = "Takes around 5~10 seconds")

      if(input$normdata=="takeRPM"){
        transdata <- rpm(transdata, colSums(transdata))
      }else if(input$normdata=="takesizeNorm"){
        sf <- norm_factors(transdata)
        transdata <- t(t(transdata)/sf)
      }else if(input$normdata=="takeuq"){
        transdata <- uq(transdata)
      }

      transdata <- rmv_constant_0(transdata, pct=0.95, minimum=0)

      if(input$transformdata=="takevst"){
        transdata <- vst(round(transdata))
      }else if(input$transformdata=="takelog"){
        validate(
          need(!any(transdata<0), "Please make sure all the numbers in the data are positive values")
        )
        transdata <- log2(transdata + 1)
      }
    })
    return(transdata)
  }

  observeEvent(rawdata(), {
    transdata <- run_trans(rawdata())
    transform_data$data <- transdata
    transform_data$rawdata <- rawdata()
  })
  observeEvent(input$run_transf, {
    transdata <- run_trans(rawdata())
    transform_data$data <- transdata
    transform_data$rawdata <- rawdata()
  })

  output$dl_transformdata <- downloadHandler(
    filename <- function() {
      paste(file_prefix(), "filtered.txt", sep=".")
    },
    content = function(file) {
      data <- transform_data$data
      colnames(data) <- gsub("\\.", "-", colnames(data))
      write.table(data, file = file, quote=F, row.names=T, sep="\t")
    }
  )

  output$readDistrib <- renderPlotly({
    validate(
      need(dim(transform_data$data)[1] > 1, "Did you upload your file or select preloaded data yet?")
    )
    total <- colSums(transform_data$data)
    genecov <- colSums(transform_data$rawdata > 0)
    run_selsamp$go <- TRUE
    textlab <- paste(names(total),
                     paste("Total:", signif(total/1000000,3), "Million"),
                     paste("Gene Cov:", genecov),
                     sep="<br>")
    plot_ly(x=genecov, y=total, text=textlab, source="readdist", #type="scattergl",
            mode="markers", hoverinfo="text",
            themes="Catherine", marker = list(size=10, opacity=0.65)) %>%
            layout(yaxis = list(title = "Total transcripts counts", zeroline=TRUE),
                   xaxis = list(title = "Gene Coverage", zeroline=TRUE),
                   dragmode = "select",
                   showlegend=FALSE)
  })

  sel_samp <- reactive({
    # Get subset based on selection
    event.data <- event_data("plotly_selected", source="readdist")

    # If NULL dont do anything
    if(is.null(event.data) == T) return(NULL)
    if(length(event.data$curveNumber) == 0) return(NULL)

    # Get selected samples
    sel_idx <- subset(event.data, curveNumber == 0)$pointNumber + 1
    non_sel_idx <- subset(event.data, curveNumber == 0)$pointNumber + 1
    if(max(sel_idx) > ncol(transform_data$data)) return(NULL)

    res <- data.frame("Total_Reads"=colSums(transform_data$data)[sel_idx],
                      "Gene_Coverage"=colSums(transform_data$rawdata>0)[sel_idx])
    return(res)
  })

  output$selsamp_tb <- DT::renderDataTable({
    DT::datatable(sel_samp(),
                  rownames= TRUE,
                  options = list(scrollX = TRUE,
                                 scrollY = TRUE,
                                 dom = 'tipr',
                                 pageLength = 5
                  )
    )
  }, server=TRUE)

  output$selsamp_txt <- renderPrint({
    cat("Selected samples\n\n")
    idx <- input$selsamp_tb_rows_selected
    if(length(idx) > 0){
      cat(paste(rownames(sel_samp())[idx], sep = ', '))
    }else{
      cat(NULL)
    }
  })

  ### Filter selected samples from user
  observeEvent(input$run_filtsamp, {
    idx <- input$selsamp_tb_rows_selected
    filtsamp <- rownames(sel_samp())[idx]
    if(!is.null(filtsamp)){
      idx <- match(filtsamp, colnames(transform_data$data))
      if(length(idx) > 0){
        temp <- transform_data$data[, -idx]
        transform_data$data <- temp
        temp <- transform_data$rawdata[, -idx]
        transform_data$rawdata <- temp
      }else{
        warning("Can not find selected sample names in transform_data")
      }
    }
  })

  merged <- reactive({
    rawdata <- transform_data$data
    genelist <- NULL

    if(!is.null(input$selected_samples) && sum(input$selected_samples %in% colnames(rawdata))==length(input$selected_samples)){
      rawdata <- rawdata[, input$selected_samples]
    }

    if(input$list_type=='Top Ranks'){
      select.genes <- select_top_n(apply(rawdata,1,input$math), n=input$top.num, bottom=ifelse(input$orders=="bottom", T, F))
      genelist <- names(select.genes)
    }else if(input$list_type=='Upload gene list'){
      validate(
        need(!is.null(input$genefile), "Please upload a file with list of genes, with column name 'Gene' ")
      )
      genefile <- read.table(input$genefile$datapath, header=T, sep="\t")

      ext <- file_ext(input$genefile[1])
      n <- 2
      withProgress(message = 'Loading gene list', value = 0, {
        incProgress(1/n, detail = "Usually takes ~10 seconds")
        if (ext=="csv"){
          genefile <- read.table(input$genefile$datapath, header=T, sep =',', stringsAsFactors = FALSE)
        }else if (ext=="out" || ext=="tsv" || ext=="txt"){
          genefile <- read.table(input$genefile$datapath, header=T, sep ='\t', stringsAsFactors = FALSE)
        }else{
          print("File format doesn't support, please try again")
          return(NULL)
        }
      })
      validate(
        need(sum(colnames(genefile)=='Gene'), "Please make sure genefile contains column name 'Gene'")
      )
      validate(
        need(length(intersect(genefile$Gene, rownames(rawdata)))>=2, "We did not find enough genes that overlap with the expression data (<2). Please check again and provide a new list! ")
      )
      genelist <- genefile$Gene
    }else if(input$list_type == 'Whole transcriptome'){
      genelist <- rownames(rawdata)
    }

    idx <- unique(match(genelist, rownames(rawdata)))
    merged <- rawdata[idx[!is.na(idx)], ]
    return(merged)
  })

  output$filterdatatbl = DT::renderDataTable({
    filter_data <- merged()
    DT::datatable(filter_data, rownames= TRUE,
                  selection = 'single',
                  options = list(scrollX = TRUE,
                                 scroller = TRUE,
                                 pageLength = 5
                  )
    ) %>% formatRound(1:ncol(filter_data), 2)
  }, server=TRUE)


  callModule(feature, "sample", reactive({ merged() }))
  callModule(network, "gene", reactive({ merged() }))


  #' Scatter plot
  observeEvent(transform_data$data, {
    transform_data <- transform_data$data
    updateSelectizeInput(session, 'cor_sample1',
                         server = TRUE,
                         choices = sort(as.character(colnames(transform_data))),
                         selected = sort(as.character(colnames(transform_data)))[1]
    )
    updateSelectizeInput(session, 'cor_sample2',
                         server = TRUE,
                         choices = sort(as.character(colnames(transform_data))),
                         selected = sort(as.character(colnames(transform_data)))[2]
    )
  })
  scatterData <- reactive({
    transform_data <- transform_data$data
    validate(
      need(input$cor_sample1 != "", "Please select one sample for X axis"),
      need(input$cor_sample2 != "", "Please select one sample for Y axis")
    )
    x <- transform_data[, input$cor_sample1]
    y <- transform_data[, input$cor_sample2]
    if(input$scatter_log){
      if(input$transformdata == "takelog"){
        message("Already did log transformation")
      }else{
        x <- log2(x+1)
        y <- log2(y+1)
      }
    }
    res <- data.frame(x=x, y=y)
    return(res)
  })
  output$scatterPlot <- renderPlotly({
    transform_data <- transform_data$data
    x <- scatterData()$x %>% as.numeric
    y <- scatterData()$y %>% as.numeric
    plot_range <- c(0, min(x,y):max(x,y), max(x,y) + 0.2)
    reg = lm(y ~ x)
    modsum = summary(reg)
    r <- signif(cor(x,y), 3)
    R2 = signif(summary(reg)$r.squared, 3)
    textlab <- paste(paste0("x = ", signif(x,3)), paste0("y = ", signif(y,3)), rownames(transform_data), sep="<br>")

    n <- 2
    withProgress(message = 'Generating Scatter Plot', value = 0, {
      # p <- (x=~x, y=~y, text=~textlab, type="scattergl", #source = "scatter",
      #              mode="markers", hoverinfo="text",
      #              #themes="Catherine",
      #              marker = list(size=8, opacity=0.65)) %>%
      p <- plot_ly() %>%
       add_trace(x=~x, y=~y, text=~textlab, type="scattergl",
                   mode="markers", hoverinfo="text",
                   marker = list(size=8, opacity=0.65)) %>%
       add_trace(x=~plot_range, y=~plot_range, type="scattergl", mode = "lines", name = "Ref", line = list(width=2)) %>%
            layout(title = paste("cor:", r),
                   xaxis = list(title = input$cor_sample1, zeroline=TRUE),
                   yaxis = list(title = input$cor_sample2, zeroline=TRUE),
                   showlegend=FALSE)
    })
    # if(input$show_r2){
    #   p <- p %>% layout( annotations = list(x = x, y = y, text=paste0("R2=", R2), showarrow=FALSE) )
    # }

    if(input$show_fc){
      df <- data.frame(x = plot_range[plot_range>=1], y=plot_range[plot_range>=1]-1)
      p <- p %>%
        add_trace(data=df, x= ~x, y= ~y, type = "scattergl", mode = "lines", name = "2-FC", line = list(width=1.5, color="#fec44f")) %>%
        add_trace(data=df, x= ~y, y= ~x, type = "scattergl", mode = "lines", name = "2-FC", line = list(width=1.5, color="#fec44f"))
    }
    p #%>% toWebGL()
  })
  output$scatterdatatbl = DT::renderDataTable({
    x <- scatterData()$x
    y <- scatterData()$y
    if(input$transformdata == "takelog"){
      logFC = x - y
    }else{
      esp <- 0.001
      logFC = log2((x+esp)/(y+esp))
    }
    GeneCard <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                              rownames(transform_data$data), "'>", rownames(transform_data$data), "</a>")
    scatter_data <- data.frame(Gene=GeneCard,
                               X = x,
                               Y = y,
                               logFC = logFC,
                               meanFC = mean(x+y)*logFC) %>%
                    dplyr::arrange(desc(abs(logFC)))
    DT::datatable(scatter_data, rownames= FALSE,
                  options = list(scrollX = TRUE,
                                 pageLength = 8),
                  escape = FALSE
    ) %>% formatRound(2:5, 3)
  }, server=TRUE)


  ### Main part of the program ###
  nmf_seed <- eventReactive(input$runNMF, {
    nmf_seed <- as.numeric(input$nmf_seed)
    if(nmf_seed == 0){
      set.seed(as.integer(Sys.time()))
      sample(1:999999, 1)
    }
  })
  nmf_res <- eventReactive(input$runNMF, {
    merged <- merged()
    validate(
      need(is.numeric(input$num_cluster) & is.whole(input$num_cluster) & input$num_cluster>1, "Please provide a valid positive numeric integer value for 'Num of clusters(K)'") %then%
      need(input$num_cluster <= ncol(merged), "Number of clusters(K) must be smaller than number of samples\nPlease Try selecting another K") %then%
      need(is.numeric(input$nrun) & is.whole(input$nrun) & input$nrun>1, "Please provide a valid numeric integer value for 'nrun'")
    )

    if(input$selectfile == "saved"){
      nmfres <- rda()$nmfres
    }else{
      if(ncol(merged) >= 300){
        # toggleModal(session, "nmfModal1", toggle = "open")

        # Run NMF through YABI
        # test_yabi <- try(system("yabish login yabi buffalo! backends"))
        # Will remove yabi later, right now just don't let yabi command start and run NMF on the shiny01 server
        test_yabi <- try(system("naikai"))
        if(test_yabi == 0){
          createAlert(session, "YabiAlert", "YabiAlert1", title = "Running NMF on Amazon Cloud", style = "success",
                      content = paste0("Sample size is ", ncol(merged),
                                       ". We are running tiis job on the cloud and will notify you when the results are ready.<br>"),
                      append = FALSE)
          # pop up window asking for email
          data_path <- "/mnt/sake-uploads"
          # data_path <- "~/Desktop/sake-uploads"
          output <- file.path(data_path, paste0(file_prefix(), ".txt"))
          write.table(merged, output, sep="\t", quote=F)
          command <- paste("su - centos -c '/home/centos/yabish_NMF_email.sh",
                           "-d", output,
                           "-t", nrow(merged),
                           "-k", input$num_cluster,
                           "-r", input$nrun,
                           "-m", input$mode,
                           "-n", "FALSE",
                           "-a", input$algorithm,
                           "-s", input$nmf_seed,
                           "-c", 8,
                           "-x", "FALSE",
                           "-f", 0,
                           "-q", "FALSE'")
          # send email to the user and stop sake
          # t1 <- isolate(try(system(command, intern = TRUE)))
          print("Sending command to yabi")
          t1 <- try(system(command, wait = FALSE))
          print(paste0("Command sent:", command))
          if(t1 == 0){
            Sys.sleep(5)
            closeAlert(session, "YabiAlert1")
            Sys.sleep(5)
            print("closing YabiAlert and stopping sake")
            # stopApp(11)
            stop(paste("sent jobs to yabi and stop sake:", command))
          }else{
            createAlert(session, "YabiAlert", "YabiAlert2", title = "Error message from running Yabi", style = "danger",
                        content = t1,
                        append = FALSE)
            Sys.sleep(5)
            stop(paste("command creates erros:", command))
          }
        }else{
          createAlert(session, "NMFAlert", "NMFAlert1", title = "WARNING", style = "warning",
                      content = paste0("Sample size is ", ncol(merged), ". This might take long to run. <br>"),
                      append = FALSE)
        }
      }

      # Run NMF on local server
      ptm <- proc.time()
      withProgress(message = 'Running NMF', value = 0, {
        incProgress(1/2, detail = "Takes around 30~60 seconds for each K")
        nmfres <- myNMF(merged,
                        mode=input$mode,
                        cluster=input$num_cluster,
                        nrun=input$nrun,
                        ncores = cores(),
                        algorithm = input$algorithm,
                        seed = nmf_seed()
        )
        b <- proc.time() - ptm
        print("### Time to run NMF ### ")
        print(b)
      })

      closeAlert(session, "NMFAlert1")
      closeAlert(session, "YabiAlert2")
    }
    return(nmfres)
  })

  ### Download NMF ###
  observeEvent(input$runNMF, {
    output$dlNMF_UI <- renderUI({
      tagList(
        column(width=2, br(), downloadButton("dlNMF", "Download NMF result", class = 'dwnld'))
      )
    })
  })

  output$dlNMF <- downloadHandler(
    filename <- function() {
      paste( file_prefix(), "nmf_results.zip", sep=".")
    },
    content = function(file) {
      tmpdir <- tempdir()
      current_dir <- getwd()
      setwd(tmpdir)

      filenames <- sapply(c("nmf_Groups", "nmf_Features", "Original_plus_NMF"),
                          function(x) paste(file_prefix(), x, "csv", sep="."))
      write.csv(nmf_groups(), filenames[1], quote=F, row.names=F)
      write.csv(nmf_features(), filenames[2], quote=F, row.names=F)
      write.csv(ori_plus_nmfResult()[["Total"]], filenames[3], quote=F)

      for(i in 1:(length(ori_plus_nmfResult())-1)){
        filenames <- c(filenames, paste(file_prefix(), paste0("onlyNMF", i), "csv", sep="."))
        write.csv(ori_plus_nmfResult()[[paste("NMF", i)]], filenames[length(filenames)], quote=F)
      }

      nmfres <- nmf_res()
      rawdata <- rawdata()
      nmfres_filename <- paste(file_prefix(), "nmfres", "rda", sep=".")
      save(rawdata, nmfres, file=nmfres_filename)
      zip(zipfile=file, files=c(filenames, nmfres_filename))
      setwd(as.character(current_dir))
    },
    contentType = "application/zip"
  )

  ### Estimate K ###
  runSummary <- reactive({
    req(class(nmf_res()) == "NMF.rank")
    nmf_res <- nmf_res()
    summary <- nmf_summary(nmf_res)
    if(input$mode=="estim"){
      summary <- t(summary)
    }
    return(summary)
  })
  output$estimSummary = DT::renderDataTable({
    runSummary <- runSummary()
    DT::datatable(runSummary, rownames= TRUE,
                  extensions = 'Buttons',
                  options = list(dom = 'Bfrtip',
                                 buttons =
                                   list('copy', 'print', list(
                                     extend = 'collection',
                                     buttons = list(list(extend='csv',
                                                         filename = filename),
                                                    list(extend='excel',
                                                         filename = filename),
                                                    list(extend='pdf',
                                                         filename= filename)),
                                     text = 'Download'
                                   )),
                                 scrollX = TRUE,
                                 pageLength = 12
                  )
    ) %>% formatRound(1:ncol(runSummary), 2)
  }, server=FALSE)

  output$estimPlot <- renderPlot({
    validate(
      need(class(nmf_res()) == "NMF.rank", "Seems like you haven't run NMF 'estim run' yet\nOr the loaded .rda file is from 'real run' result")
    )
    nmf_res <- nmf_res()
    if(input$estimtype=="consensus"){
      consensusmap(nmf_res)
    }else if(input$estimtype=="quality"){
      plot(nmf_res)
    }
  }, height=777)


  ### NMF plots ###
  output$nmfplot <- renderPlot({
    validate(
      need(class(nmf_res()) == "NMFfitX1", "Seems like you haven't run NMF 'real' run yet\nOr the loaded .rda file is from 'estim run' result")
    )
    nmf_res <- nmf_res()
    nmf_plot(nmf_res, type=input$plottype, silorder=input$nmfplot_silhouette,
             subsetRow = ifelse(input$select_feature_num>0, as.numeric(input$select_feature_num), TRUE ))
  },height = 777)

  output$silhouetteplot <- renderPlot({
    validate(
      need(class(nmf_res()) == "NMFfitX1", "Seems like you haven't run NMF 'real' run yet")
    )
    nmf_res <- nmf_res()
    nmf_silhouette_plot(nmf_res, type=input$plottype)
  },height = 777)

  output$dl_nmf_estimplot <- downloadHandler(
    filename <- function() {
      paste(file_prefix(), input$mode, "nmf.pdf", sep=".")
    },
    content = function(file) {
      nmf_res <- nmf_res()
      pdf(file, width=input$estim_pdf_w, height=input$estim_pdf_h)
      if(input$mode == "estim"){
        consensusmap(nmf_res)
        print(plot(nmf_res))
      }
      dev.off()
    }
  )

  output$dl_nmf_realplot <- downloadHandler(
    filename <- function() {
      paste(file_prefix(), input$mode, "nmf.pdf", sep=".")
    },
    content = function(file) {
      nmf_res <- nmf_res()
      pdf(file, width=input$real_pdf_w, height=input$real_pdf_h)
      if(input$mode == "real"){
        nmf_plot(nmf_res, type="samples", silorder=input$nmfplot_silhouette)
        nmf_plot(nmf_res, type="features", silorder=input$nmfplot_silhouette,
                 subsetRow = ifelse(input$select_feature_num>0, as.numeric(input$select_feature_num), TRUE ))
        nmf_plot(nmf_res, type="consensus")
        if(!is.null(tsne_2d$data)){
          print(nmfplottsne())
        }
      }
      dev.off()
    }
  )

  nmf_groups <- reactive({
    validate(
      need(class(nmf_res()) == "NMFfitX1", "Seems like you haven't run NMF 'real' run yet")
    )
    withProgress(message = 'Extracting Groups info', value = NULL, {
      nmf_extract_group(nmf_res(), type=input$predict)
    })
  })
  output$nmfGroups <- DT::renderDataTable({
    nmf_groups <- nmf_groups()
    filename <- "nmf_groups"
    colnames(nmf_groups)[which(colnames(nmf_groups)=="nmf_subtypes")] <- "Group"
    # Get subset based on selection from scatter plot
    event.data <- event_data("plotly_selected", source = "nmfGroups")
    if(is.null(event.data) == TRUE){
      sel_idx <- NULL
    }else{
      sel_idx <- subset(event.data, curveNumber == 0)$pointNumber + 1
    }
    DT::datatable(nmf_groups, rownames=FALSE, escape=-1,
                  selection = list(selected = sel_idx),
                  options = list(dom = 'rtip',
                                 scrollX = TRUE,
                                 pageLength = 8,
                                 order=list(list(1,'asc'))
                  )
    ) %>% formatRound(3:ncol(nmf_groups), 3)
  }, server=FALSE)

  # Add t-SNE plot next to NMF group
  observeEvent(input$run_nmftSNE, {
    plot_data <- plot_data()[['plot_data']]
    nsamples <- ncol(plot_data()[['plot_data']])
    if(input$nmftsne_perplexity*3 > (nsamples-1)) {
           createAlert(session, "perplexityAlert", "perplexityAlert1", title = "WARNING", style = "warning",
                       content = paste("Perpleixty", input$nmftsne_perplexity,
                                       "is too large compared to num of samples", nsamples),
                       append = FALSE)
    }else{
      closeAlert(session, "perplexityAlert1")

      withProgress(message = 'Running t-SNE', value = NULL, {
        incProgress(1/3, detail = "For t-SNE 2D")
        try(
          tsne_2d$data <- run_tsne(plot_data, iter=input$tsne_iter, dims=2,
                                   initial_dims = input$tsne_pca_num, theta = input$tsne_theta,
                                   perplexity=input$nmftsne_perplexity, cores=cores())
        )
        incProgress(2/3, detail = "For t-SNE 3D")
        try(
          tsne_3d$data <- run_tsne(plot_data, iter=input$tsne_iter, dims=3,
                                   initial_dims = input$tsne_pca_num, theta = input$tsne_theta,
                                   perplexity=input$nmftsne_perplexity, cores=cores())
        )
      })
    }
  })

  nmfplottsne <- reactive({
    plot_data <- plot_data()[['plot_data']]
    nsamples <- ncol(plot_data)
    validate(
      need(input$mode=="real", "This part won't work for 'estim' module\nPlease Try NMF 'Real' run") %then%
      need(!is.null(tsne_2d$data), "Please hit 'Run t-SNE' button")
    )
    color <- NULL
    # Will add an option to select by default column names from samples
    if(length(nmf_groups()$nmf_subtypes)>0){
      idx <- match(colnames(plot_data), nmf_groups()$Sample_ID)
      nmf_subtypes <- nmf_groups()$nmf_subtypes[idx]
      # color <- create.brewer.color(nmf_subtypes, length(unique(nmf_subtypes)), "naikai")
      color <- nmf_subtypes
    }else if(length(ColSideColors()[['name']][,1])>0){
      color <- ColSideColors()[['name']][,1]
    }else{
      stop("Error: color scheme undefined")
    }
    # # default by user selection
    point_size <- input$plot_point_size - 5
    if(input$nmf_prob){
      if(input$predict=="consensus"){
        point_size <- point_size + point_size * (1-nmf_groups()$sil_width)
      }else if(input$predict == "samples"){
        point_size <- point_size + point_size * (1-nmf_groups()$prob)
      }
    }
    nmf_tsne_res <- tsne_2d$data
    naikai <- plot_tsne(nmf_tsne_res, color=color, alpha=input$plot_point_alpha,
                        add.label = input$plot_label, save.plot=F, real.plot=F,
                        add.legend = input$plot_legend,
                        point.size=point_size, label.size=input$plot_label_size)
    ### Check if user select rows (samples)
    s = input$nmfGroups_rows_selected
    if (length(s)){
      if(!input$nmf_prob){
        sub_color <- create.brewer.color(nmf_groups()$nmf_subtypes, length(unique(color)), "naikai")
        sub_tsne_res <- parse_tsne_res(nmf_tsne_res) %>% as.data.frame %>% dplyr::slice(s)
        naikai <- naikai + geom_point(data=sub_tsne_res, aes(x=x, y=y), size=point_size+2, colour=sub_color[s], alpha=input$plot_point_alpha)
      }
    }
    p <- ggplotly(naikai, source="nmfGroups") %>% layout(dragmode = "select")
    return(p)
  })
  output$nmftsneplot <- renderPlotly({
    withProgress(message = 'Genrating t-SNE plot', value = NULL, {
      nmfplottsne()
    })
  })


  output$nmfgrp_var <- renderPlotly({
    nmf_groups <- nmf_groups()
    plot_data <- plot_data()[['plot_data']]
    idx <- match(colnames(plot_data), nmf_groups$Sample_ID)

    colnames(plot_data) <- paste0("NMF", nmf_groups$nmf_subtypes[idx])
    dd <- plot_data %>% t %>% reshape2::melt(.) %>% group_by(Var2, Var1) %>% summarise(var=var(value))

    num_clus <- dd$Var1 %>% levels %>% length
    mycolor <- create.brewer.color(1:num_clus, num=num_clus, name="naikai")
    dd$Var1 <- factor(dd$Var1, levels = paste0("NMF", 1:num_clus))

    a <- ggplot(dd, aes(var, fill=Var1, colour=Var1)) +
          geom_density(alpha=input$cl_alpha) + scale_x_log10() +
          scale_colour_manual(values = mycolor) +
          scale_fill_manual(values = mycolor) +
          labs(x="Log (variance)", y="Density", colour="", fill="") +
          theme_bw()
    p <- plotly_build(a)

    # remove duplicated name
    p$data <- lapply(p$data, FUN = function(x){
      x$name <-  sub("\\((\\S+\\d),.*", "\\1", x$name, perl=TRUE)
      x$line$width <- 1.2
      return(x)
    })
    p$layout$legend$font$size <- as.numeric(input$cl_lgsize)
    p$layout$xaxis$tickfont$size <- input$cl_tickfont
    p$layout$xaxis$titlefont$size <- input$cl_axisfont
    p$layout$xaxis$titlefont$family <- "Arial"
    p$layout$yaxis$tickfont$size <- input$cl_tickfont
    p$layout$yaxis$titlefont$size <- input$cl_axisfont
    p$layout$yaxis$titlefont$family <- "Arial"

    p
  })

  output$nmfgrp_expgene <- renderPlotly({
    rawdata <- rawdata()
    nmf_groups <- nmf_groups()
    idx <- match(colnames(rawdata), nmf_groups$Sample_ID)
    rawdata <- rawdata[, !is.na(idx)]
    colnames(rawdata) <- paste0("NMF", nmf_groups$nmf_subtypes[idx[!is.na(idx)]])
    dd <- data.frame(Group=colnames(rawdata), Genes=colSums(rawdata>0))

    num_clus <- dd$Group %>% levels %>% length
    mycolor <- create.brewer.color(1:num_clus, num=num_clus, name="naikai")

    a <- ggplot(dd, aes(Genes, fill=Group, colour=Group)) +
          geom_density(alpha=input$cl_alpha) +
          scale_fill_manual(values = mycolor, guide="none") +
          scale_color_manual(values = mycolor) +
          labs(x="No. of expressed genes", y="Density", colour="", fill="") +
          theme_bw()
    p <- plotly_build(a)

    # remove duplicated name
    p$data <- lapply(p$data, FUN = function(x){
      x$name <-  sub("\\((\\S+\\d),.*", "\\1", x$name, perl=TRUE)
      x$line$width <- 1.2
      return(x)
    })
    p$layout$legend$font$size <- as.numeric(input$cl_lgsize)
    p$layout$xaxis$tickfont$size <- input$cl_tickfont
    p$layout$xaxis$titlefont$size <- input$cl_axisfont
    p$layout$xaxis$titlefont$family <- "Arial"
    p$layout$yaxis$tickfont$size <- input$cl_tickfont
    p$layout$yaxis$titlefont$size <- input$cl_axisfont
    p$layout$yaxis$titlefont$family <- "Arial"

    p
  })

  output$nmfgrp_coef <- renderPlotly({
    nmf_groups <- nmf_groups()
    plot_data <- plot_data()[['plot_data']]
    idx <- match(colnames(plot_data), nmf_groups$Sample_ID)
    nmf_groups <- nmf_groups[idx, ]

    colnames(plot_data) <- paste0("NMF", nmf_groups$nmf_subtypes)
    pair_cor <- function(x) {
      res <- WGCNA::cor(x, nThreads = 4)
      return(res[lower.tri(res)])
    }

    num_clus <- nmf_groups$nmf_subtypes %>% unique %>% length
    mycolor <- create.brewer.color(1:num_clus, num=num_clus, name="naikai")
    num_poss_cor <- nmf_groups$nmf_subtypes %>% table %>% apply(., 1, function(x) x * (x-1) / 2) %>% sum
    group_cor <- data.frame(matrix(NA_real_, nrow=num_poss_cor, ncol=2))

    j <- 1
    for(i in 1:num_clus){
      idx <- nmf_groups$nmf_subtypes == i
      res <- pair_cor(plot_data[, idx])
      group_cor[j:(j+length(res)-1), ] <- cbind(paste0("NMF", i), res)
      j <- j + length(res)
    }
    colnames(group_cor) <- c("NMF", "Cor")
    group_cor$Cor <- as.numeric(group_cor$Cor)

    a <- ggplot(data=group_cor, aes(x=NMF, y=Cor, color=NMF)) +
      geom_boxplot(aes(color = NMF)) +
      stat_summary(fun.y=mean, shape=19, col='black', geom='point') +
      theme_bw() +
      scale_colour_manual(values = mycolor) +
      labs(x="", y="Cor", colour="")
    # ggplotly(a)

    p <- plotly_build(a)
    # remove outliers for plotly boxplot
    p$data <- lapply(p$data, FUN = function(x){
      if(x$type == "box"){
        x$marker = list(opacity = 0)
      }
      return(x)
    })
    p$layout$legend$font$size <- as.numeric(input$cl_lgsize)
    p$layout$xaxis$tickfont$size <- input$cl_tickfont
    p$layout$xaxis$titlefont$size <- input$cl_axisfont
    p$layout$xaxis$titlefont$family <- "Arial"
    p$layout$yaxis$tickfont$size <- input$cl_tickfont
    p$layout$yaxis$titlefont$size <- input$cl_axisfont
    p$layout$yaxis$titlefont$family <- "Arial"

    p
  })

  nmf_features <- reactive({
    validate(
      need(input$mode=="real", "This part won't work for eastimating k\nPlease Try NMF 'Realrun'") %then%
      need(class(nmf_res())=="NMFfitX1", "Please rerun NMF with 'real run'")
    )
    nmf_res <- nmf_res()
    # add more options to select features here. then compare with ClaNC
    # Add rank by MAD expression
    if(input$select_method=="total"){
      res <- nmf_extract_feature(nmf_res, method="total") %>% unique
    }
    else if(input$select_method=="default"){
      res <- nmf_extract_feature(nmf_res, method=input$select_method, manual.num = as.numeric(input$select_feature_num)) %>%  unique
      validate(
        need(nrow(res)!=1 & !is.na(res$Gene), "None of the features passed through default NMF filtering criteria, please manually specify number of genes from each group or use 'select features by Rank'") %then%
        need(nrow(res)>1, "Not enought features are selected by default! \nMaybe try with more input features or try select features by Rank")
      )
    }else if(input$select_method=="rank"){
      res <- nmf_extract_feature(nmf_res, rawdata = merged(), method=input$select_method, manual.num = as.numeric(input$select_feature_num),
                          FScutoff=input$select_FScutoff, math=input$select_math) %>% unique
      validate(
        need(nrow(res)>1, "Not enought features! \nMaybe try lower the FeatureScore cutoff")
      )
    }
    return(res)
  })
  nmf_features_annot <- reactive({
    data <- nmf_features() %>% as.data.table
    merged <- NULL

    n <- 2
    withProgress(message = 'Extracting features info', value = 0, {
      incProgress(1/n, detail = "Takes around 5~10 seconds")

      ext.func.data <- ext_gene_function(data$Gene, ensembl = ensembl())
      setkey(data, "Gene")
      setkey(ext.func.data, "hgnc_symbol")
      merged <- ext.func.data[data, nomatch=0, allow.cartesian = TRUE]
      merged <- as.data.frame(merged) %>%
                dplyr::arrange(Group, desc(featureScore))
      merged <- merged[, c(4,2,7,6,8)]
      merged <- unique(merged)
      merged$featureScore <- round(merged$featureScore, 3)
      merged$prob <- round(merged$prob, 3)
    })
    return(merged)
  })
  output$nmfFeatures <- DT::renderDataTable({
    filename <- "nmf_features"
    DT::datatable(nmf_features_annot(), rownames=FALSE, escape=-1,
                  extensions = 'Buttons', selection = 'multiple',
                  options = list(dom = 'Bfrtip',
                                 buttons =
                                   list('copy', 'print', list(
                                     extend = 'collection',
                                     buttons = list(list(extend='csv',
                                                         filename = filename),
                                                    list(extend='excel',
                                                         filename = filename),
                                                    list(extend='pdf',
                                                         filename= filename)),
                                     text = 'Download'
                                   )),
                                 pageLength = 50,
                                 autoWidth = TRUE
                  )
    )
  }, server = FALSE)

  output$nmf_vioplot <- renderPlotly({
    req(nmf_features_annot())
    validate(
      need(!is.null(input$nmfFeatures_rows_selected), "Please select a gene")
    )
    withProgress(message = 'Generating box/violin plot', value = NULL, {
      gene <- nmf_features_annot()[input$nmfFeatures_rows_selected, "GeneCard"] %>%
              gsub(".*'>(.*)</a>", "\\1", .)
      gene.data <- transform_data$rawdata[gene, ]

      gene.data <- gene.data %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(genename = rowname) %>%
        tidyr::gather(Sample, Expr, -genename) %>%
        dplyr::mutate(NMF = factor(paste0("NMF", rep(nmf_groups()$nmf_subtypes, each=length(gene)))))
      if(input$sel_vioscale == "raw"){
        print('Use original scale')
      }else if(input$sel_vioscale == "log2"){
        gene.data$Expr <- log2(gene.data$Expr+1)
      }else if(input$sel_vioscale == "log10"){
        gene.data$Expr <- log10(gene.data$Expr+1)
      }

      # match the coloring from PCA, t-SNE, and heatmap
      num_clus <- gene.data$NMF %>% levels %>% length
      mycolor <- create.brewer.color(1:num_clus, num=num_clus, name="naikai")

      if(input$sel_grp == "Manual"){
          col.num <- as.numeric(factor(mycolor, levels = mycolor))
          col.lvl <- as.character(mycolor)
          col.lvl <- col.lvl[as.numeric(input$sel_grp_order)]
          mycolor <- col.lvl[col.num]
      }
      # add option to select either boxplot or violin plot
      if(input$sel_vioplot == "violin"){
        a <- ggplot(data=gene.data, aes(x=NMF, y=Expr, color=NMF)) +
          geom_violin(aes(color = NMF), scale="width", width=0.6, show.legend = FALSE) +
          geom_jitter(aes(color=NMF), alpha=0.6, width=0.1, show.legend = FALSE) +
          theme_bw() +
          facet_wrap(~ genename, ncol = input$sel_ncol) +
          scale_colour_manual(values = mycolor) +
          labs(x="", y="", colour="")
      }else if(input$sel_vioplot == "box"){
        a <- ggplot(data=gene.data, aes(x=NMF, y=Expr, color=NMF)) +
          geom_boxplot(aes(color = NMF), show.legend = FALSE) +
          geom_jitter(aes(color=NMF), alpha=0.6, width=0.1, show.legend = FALSE) +
          theme_bw() +
          facet_wrap(~ genename, ncol = input$sel_ncol) +
          scale_colour_manual(values = mycolor) +
          labs(x="", y="", colour="")
      }else{
        warning(paste("Unknown plot type:", input$sel_vioplot, "Please check again"))
      }

      # custom range for Y-axis is only being triggered when ymax > 0
      if(!is.na(input$sel_ymax) & input$sel_ymax>0){
        validate(
          need(input$sel_ymax > input$sel_ymin, "Please make sure 'ymax' is greater than 'ymin'")
        )
        a <- a + ylim(c(input$sel_ymin, input$sel_ymax))
      }

      p <- plotly_build(a)
      # remove outliers for plotly boxplot
      p[[1]]$data <- lapply(p[[1]]$data, FUN = function(x){
        if(x$type == "box"){
          x$marker = list(outliercolor = 0)
        }
        return(x)
      })

      p[[1]]$layout$margin$l <- p[[1]]$layout$margin$l + 5
      p[[1]]$layout$legend$font$size <- input$sel_legend_size
      p[[1]]$layout$showlegend = input$sel_show_legend

      # loop through ncol(Xaxis) and gene (Yaxis)
      num_xaxis <- grep("xaxis", names(p[[1]]$layout)) %>% length
      num_yaxis <- grep("yaxis", names(p[[1]]$layout)) %>% length
      for(i in 1:num_xaxis){
        p[[1]]$layout[[paste0("xaxis", ifelse(i>1, i, ""))]]$tickfont$size <- input$sel_xfont_size
      }
      for(i in 1:num_yaxis){
        p[[1]]$layout[[paste0("yaxis", ifelse(i>1, i, ""))]]$tickfont$size <- input$sel_yfont_size
      }
      for(i in 1:length(gene)){
        p[[1]]$layout$annotations[[i]]$font$size <- input$sel_title_size
      }

      p
    })
  })

  ori_plus_nmfResult <-reactive({
    rawdata <- rawdata()
    validate(
      need(!is.null(nmf_groups()), "Please run NMF and get the result first")
    )
    res <- list()
    idx <- match(colnames(rawdata), nmf_groups()$Sample_ID)
    newdata <- rawdata
    colnames(newdata) <- paste(paste0("NMF",nmf_groups()$nmf_subtypes[idx]), colnames(newdata), sep="_")
    res[["Total"]] <- newdata

    # add separate files for each NMF group
    for(i in 1:length(unique(nmf_groups()$nmf_subtypes))){
      subidx <- nmf_groups()$nmf_subtypes == i
      res[[paste("NMF", i)]] <- newdata[, idx[subidx]]
    }

    return(res)
  })

  ### Statistics ###
  output$runStatSummary <- renderTable({
    runSummary()
    }, sanitize.text.function = function(x) x, options = list(orderClasses = TRUE, lengthMenu = c(10,25,50,100, NROW(merged)), pageLength=25)
  )
  output$runStatPlot <- renderPlot({
    nmf_res <- nmf_res()
    plot(nmf_res)
  }, height=666)


  ensembl <- reactive({
    # Need to add a testing here when the BioMart is down
    # Or find offline annotation db
    # ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="jul2016.archive.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  })

  #' Heatmap setting
  geneListInfo <- reactive({
    rawdata <- transform_data$data
    gene.list <- NULL
    gene.groups <- NULL

    if(input$heatmapGene=='nmf'){
      nmf_features <- nmf_features()
      validate(
        need(nrow(nmf_features)>0, "NMF not run yet? \nPlease run NMF classification and then try again!")
      )
      gene.list <- nmf_features$Gene
    }
    else if(input$heatmapGene=='rank'){
      select.genes <- select_top_n(apply(rawdata,1,input$heat.math), n=as.numeric(input$heat.top.num), bottom=F)
      gene.list <- names(select.genes)
    }
    else if(input$heatmapGene=='preload'){
      validate(
        need(!is.null(input$predefined_list), "Please select at least one predefined gene list")
      )
      for(i in 1:length(input$predefined_list)){
        gene.list <- c(gene.list, pre_genelist[[ input$predefined_list[i] ]]$Gene)
        # if file contains 'Group' info in one of the columns, use that for RowColors
        if("Group" %in% colnames(pre_genelist[[ input$predefined_list[i] ]])){
          gene.groups <- c(gene.groups, pre_genelist[[ input$predefined_list[i] ]]$Group)
        }
      }
    }
    else if(input$heatmapGene=='manual'){
      validate(
        need(!is.null(input$manual_genes), "You didn't select any genes. Please select genes from the data") %then%
          need(length(input$manual_genes)>=2, "Please select at least two genes! ")
      )
      gene.list <- input$manual_genes
    }
    else if(input$heatmapGene=='upload'){
      validate(
        need(!is.null(input$heatmapGeneFile), "Please upload a file with list of genes, with column name 'Gene' ")
      )
      ext <- file_ext(input$heatmapGeneFile)[1]
      if (ext=="csv"){
        genefile <- read.csv(input$heatmapGeneFile$datapath, header=T, stringsAsFactors = FALSE)
      }else if (ext=="out" || ext=="tsv" || ext=="txt"){
        genefile <- read.table(input$heatmapGeneFile$datapath, header=T, sep="\t", stringsAsFactors = FALSE)
      }else{
        print("File format doesn't support, please use '.csv', '.txt', '.out', '.tsv' and try again")
        return(NULL)
      }
      validate(
        need(sum(colnames(genefile)=='Gene'), "Please make sure genefile contains column name 'Gene'")
      )
      gene.list <- genefile$Gene
      if("Group" %in% colnames(genefile)){
        gene.groups <- genefile$Group
      }
    }

    validate(
      need(length(intersect(gene.list, rownames(rawdata)))>=2, "We did not find enough genes that overlap with the expression data (<2). \nPlease check again and provide a new list! \nMaybe this is from species other than Human (default gene list are from Human)")
    )

    # only select unique gene name on rawdata
    idx <- unique(match(gene.list, rownames(rawdata)))
    na.genes <- gene.list[is.na(idx)]
    if(length(na.genes)>0){
      print(c("These genes cannot be found in the data:", paste(na.genes, collapse=", ")))
    }
    data.rownames <- rownames(rawdata)[idx[!is.na(idx)]]

    # match back those gene names to genelist to get idx for gene.groups
    ii <- match(data.rownames, gene.list)
    if(!is.null(gene.groups)){
      gene.groups <- gene.groups[ii]
    }
    gene.list <- data.rownames

    # finish loading or selecting gene list file
    res <- list()
    res[['list']] <- gene.list
    res[['group']] <- gene.groups
    return(res)
  })


  OrdColopts <- reactive({
    c("Default", "Hierarchical", "Filename", "GeneExpr")
  })
  ColClrByopts <- reactive({
    c("Filename")
  })
  observeEvent(rawdata(), {
    updateSelectizeInput(session, 'OrdCol',
                         server = TRUE,
                         choices = as.character(OrdColopts()),
                         selected = "Hierarchical"
    )
    updateSelectizeInput(session, 'ColClrBy',
                         server = TRUE,
                         choices = as.character(ColClrByopts()),
                         selected = "Filename"
    )
  })
  observeEvent(input$runNMF, {
    updateSelectizeInput(session, 'OrdCol',
                         server = TRUE,
                         choices = as.character(c(OrdColopts(), "NMF Group")),
                         selected = "NMF Group"
    )
    updateSelectizeInput(session, 'ColClrBy',
                         server = TRUE,
                         choices = as.character(c(ColClrByopts(), "NMF Group")),
                         selected = "NMF Group"
    )
  })

  # separate out this for PCA and t-SNE plots
  plot_data <- reactive({
    rawdata <- transform_data$data
    genelist <- geneListInfo()[['list']]

    idx <- match(genelist, rownames(rawdata))
    plot_data <- rawdata[idx,]

    res <- list()
    res[["plot_data"]] <- plot_data
    res[["genelist"]] <- genelist
    return(res)
  })

  heatmap_data <- reactive({
    rawdata <- transform_data$data
    genelist <- geneListInfo()[['list']]
    gene.groups <- geneListInfo()[['group']]

    # only select unique gene name on rawdata
    idx <- match(genelist, rownames(rawdata))
    heatmap_data <- rawdata[idx,]

    if(input$OrdCol == 'Filename'){    # Add options to sort by column groups
      column.names <- strsplit(colnames(heatmap_data), split="\\_")
      validate(
        need(input$sortcolumn_num <= length(column.names[[1]]), "Group num is too big, select a smaller number!")
      )
      idx <- order(sapply(column.names, function(x) x[[input$sortcolumn_num]]))
    }else if(input$OrdCol == 'GeneExpr'){ # Add options to sort by expression value of specific gene
      validate(
        need(input$sortgene_by!="", "Please select a gene!")
      )
      idx <- order(subset(heatmap_data, rownames(heatmap_data)==input$sortgene_by))
    }
    else if(input$OrdCol == 'NMF Group') { # Add options to sort by nmf groups
      idx <- nmf_groups()$nmf_subtypes %>% order
    }else{
      idx <- colnames(heatmap_data)
    }
    heatmap_data <- heatmap_data[, idx]

    res <- list()
    res[["heatmap_data"]] <- heatmap_data
    res[["genelist"]] <- genelist
    res[["gene.groups"]] <- gene.groups
    return(res)
  })

  color <- reactive({
    color <- switch(input$heatColorScheme,
                    "RdBu_Brewer" = colorRampPalette(rev(brewer.pal(input$color_num, "RdBu")))(51),
                    "RdYlBu_Brewer" = colorRampPalette(rev(brewer.pal(input$color_num, "RdYlBu")))(51),
                    "Red_Green" = colorRampPalette(c("red2", "white", "green"))(input$color_num),
                    "Red_Blue" = colorRampPalette(c("blue", "white", "red"))(input$color_num),
                    "Rainbow" = colorRampPalette(c("lightblue", "blue", "white", "gold2", "red2"))(input$color_num),
                    "Black_White" = colorRampPalette(c("black", "white"))(input$color_num)
    )
    return(color)
  })
  ColScheme <- reactive({
    ColScheme <- list()
    ColScheme <- c(input$ColScheme1, input$ColScheme2, input$ColScheme3, input$ColScheme4, input$ColScheme5, input$ColScheme6)
    return(ColScheme)
  })

  ColSideColors <- reactive({
    heatmap_data <- heatmap_data()[['heatmap_data']]
    res <- list()
    if(input$ColClrBy=='NMF Group'){
      req(nmf_groups)
      # rematch the sample order to the sorted order based on selection from input$OrdCol
      idx <- match(colnames(heatmap_data), nmf_groups()$Sample_ID)
      nmf_subtypes <- nmf_groups()$nmf_subtypes[idx]
      res[['color']] <- create.brewer.color(nmf_subtypes, length(unique(nmf_subtypes)), "naikai") %>% as.matrix
      res[['name']] <- paste0("NMF", nmf_subtypes) %>% as.matrix
    }else if(input$ColClrBy == 'Filename'){
      res <- name_to_color(colnames(heatmap_data), split_pattern ="\\_",
                           num_color = 6,
                           ColScheme = ColScheme()[1:6] )
                           # num_color = input$ColSideColorsNum,
                           # ColScheme = ColScheme()[1:input$ColSideColorsNum] )
      # add options to select which column to show
      if(is.null(input$SelColSideColors)){
        res[['color']] <- NULL
        res[['name']] <- NULL
      }else{
        res[['color']] <- res[['color']][, as.numeric(input$SelColSideColors)] %>% as.matrix()
        res[['name']] <- res[['name']][, as.numeric(input$SelColSideColors)] %>% as.matrix()
      }
    }
    return(res)
  })

  RowSideColors <- reactive({
    groups <- heatmap_data()[['gene.groups']]
    row.color <- NULL
    res <- list()
    res[['color']] <- NULL
    res[['name']] <- NULL

    if(!is.null(groups)){
      row.color <- rbind(row.color, create.brewer.color(groups, input$RowColScheme.num, input$RowColScheme))
      res[["color"]] <- as.matrix(row.color)
      res[["name"]] <- as.matrix(groups)
    }
    return(res)
  })

  Colv <- reactive({
    Colv <- NULL
    validate(
      need(nrow(heatmap_data()[['heatmap_data']] > 2), "Did not get enough genes that overlapped with your data\nPlease try another gene list or use rank from data set")
    )
    if(input$OrdCol=="Hierarchical"){
      Colv <- heatmap_data()[['heatmap_data']] %>% t %>% dist(method=input$distance) %>% hclust(method=input$linkage) %>% as.dendrogram
      ColMean <- colMeans(heatmap_data()[['heatmap_data']], na.rm=T)
      Colv <- reorderfun(Colv, ColMean)
    }
    return(Colv)
  })

  Rowv <- reactive({
    Rowv <- NULL
    validate(
      need(nrow(heatmap_data()[['heatmap_data']] > 2), "Did not get enough genes that overlapped with your data\nPlease try another gene list or use rank from data set")
    )
    if(input$OrdRow=="hier"){
      Rowv <- heatmap_data()[['heatmap_data']] %>% dist(method=input$distance) %>% hclust(method=input$linkage) %>% as.dendrogram
      RowMean <- rowMeans(heatmap_data()[['heatmap_data']], na.rm=T)
      Rowv <- reorderfun(Rowv, RowMean)
    }
    return(Rowv)
  })

  dendrogram <- reactive({
    if(is.null(Rowv()) & is.null(Colv())){
      return("none")
    }else if(is.null(Rowv()) & !is.null(Colv())){
      return("column")
    }else if(!is.null(Rowv()) & is.null(Colv())){
      return("row")
    }else if(!is.null(Rowv()) & !is.null(Colv())){
      return("both")
    }
  })

  ### Observe and update input selection ###
  observe({
    transform_data <- transform_data$data
    if(!is.null(transform_data)){
      if (length(rownames(transform_data)) > 1){
        # update the render function for selectize - update manual if user choose to select its own markers
        updateSelectizeInput(session, 'manual_genes',
                             server = TRUE,
                             choices = as.character(rownames(transform_data))
        )
      }
    }
  })
  observe({
    if (!is.null(geneListInfo()[['list']])){
      updateSelectizeInput(session, 'sortgene_by',
                           server = TRUE,
                           choices = as.character(geneListInfo()[['list']])
      )
    }
  })

  output$d3heatmap <- renderD3heatmap({
    heatmap_data <- heatmap_data()[['heatmap_data']]
    ColSideColors <- ColSideColors()[["color"]]
    if(!is.null(ColSideColors)){
      ColSideColors <- t(ColSideColors[, ncol(ColSideColors):1])
    }
    myd3Heatmap(
      as.matrix(heatmap_data),
      type="heatmap",
      dist.m=input$distance,
      hclust.m=input$linkage,
      scale="row",
      scale.method=input$scale_method,
      scale.first=input$scale_first,
      RowSideColors = RowSideColors()[["color"]],
      ColSideColors = ColSideColors,
      color = color(),
      Colv = Colv(),
      Rowv = Rowv(),
      cexCol=input$cexRow,
      cexRow=input$cexCol,
      na.rm = T,
      show_grid = FALSE,
      dendrogram = dendrogram(),
      dendro.ord = "manual",
      anim_duration = 0
    )
  })

  plotInput <- function(){
    heatmap_data <- heatmap_data()[['heatmap_data']]
    filename <- paste(input$type, input$math, paste0(input$orders, paste0(nrow(merged()), "genes")), sep=".")
    myHeatmap.3(as.matrix(heatmap_data),
                type="heatmap",
                color=color(),
                RowSideColors = RowSideColors()[["color"]],
                RowSideColors.name=RowSideColors()[["name"]],
                ColSideColors = ColSideColors()[["color"]],
                ColSideColorsSize=input$ColSideColorsSize,
                ColSideColors.name=ColSideColors()[["name"]],
                cor.method=input$cor.method,
                dist.m=input$distance,
                hclust.m=input$linkage,
                title=input$title,
                Colv = Colv(),
                Rowv = Rowv(),
                scale="row",
                scale.method=input$scale_method,
                scale.first=input$scale_first,
                dendro.ord="manual",
                dendrogram=dendrogram(),
                cexCol=input$cexCol,
                cexRow=input$cexRow,
                col.legend = input$col_legend,
                row.legend = input$row_legend,
                file.prefix=filename,
                save.image=F
    )
  }

  output$dl_heatmap<- downloadHandler(
    filename <- function() {
      paste(file_prefix(), "heatmap.pdf", sep="_")
    },
    content = function(file) {
      pdf(file, width = input$heat_pdf_wd, height = input$heat_pdf_ht)
      plotInput()
      dev.off()
    }
  )

  plot_colopts <- reactive({
    c("Default", "Filename", "GeneExpr", "ReadDepth", "NumExpresGenes")
  })
  observeEvent(rawdata(), {
    updateSelectizeInput(session, 'pt_col',
                         server = TRUE,
                         choices = as.character(plot_colopts()),
                         selected = "Filename"
    )
    choices <- rawdata() %>% colnames() %>% strsplit("_") %>% .[[1]] %>% length %>% seq(1, .) %>% as.character()
    updateSelectizeInput(session, 'pt_file_grp',
                         server = TRUE,
                         choices = choices,
                         selected = choices[1]
    )
  })
  observeEvent(input$pt_file_grp, {
    # this should change with respect to the selection of filename group
    choices <- rawdata() %>% colnames() %>% strsplit("_") %>%
      sapply(., function(x) x[as.numeric(input$pt_file_grp)]) %>% unique %>% length %>%
      seq(1, .) %>% as.character()
    updateSelectizeInput(session, 'pt_grp_order',
                         server = TRUE,
                         choices = choices,
                         selected = choices[1]
    )
  })
  observeEvent(transform_data$data, {
    updateSelectizeInput(session, 'pt_allgene',
                         server = TRUE,
                         choices = as.character(rownames(transform_data$data)),
                         selected = as.character(rownames(transform_data$data)[1])
    )
  })

  observeEvent(input$runNMF,{
    updateSelectizeInput(session, 'pt_col',
                         server = TRUE,
                         choices = as.character(c(plot_colopts(), "NMF Group", "NMF Feature")),
                         selected = 'NMF Group'
    )
    updateSelectizeInput(session, 'pt_nmfgene',
                         server = TRUE,
                         choices = as.character(nmf_features()$Gene),
                         selected = as.character(nmf_features()$Gene[1])
    )
  })
  observeEvent(nmf_groups()$nmf_subtypes, {
    choices <- nmf_groups()$nmf_subtypes %>% unique %>% length %>% seq(1, .) %>% as.character()
    updateSelectizeInput(session, 'sel_grp_order',
                         server = TRUE,
                         choices = choices,
                         selected = choices[1]
    )
  })

  point_col <- reactive({
    # col <- toRGB("steelblue")
    col <- "steelblue"
    plot_data <- plot_data()[['plot_data']]
    group <- NULL
    if(input$pt_col == "Default"){
      group <- "Default"
    }else if(input$pt_col == "NMF Group"){
      req(nmf_groups)
      idx <- match(colnames(plot_data), nmf_groups()$Sample_ID)
      nmf_subtypes <- nmf_groups()$nmf_subtypes[idx]
      col <- create.brewer.color(nmf_subtypes, length(unique(nmf_subtypes)), "naikai")
      group <- paste0("NMF", nmf_subtypes)
    }else if(input$pt_col == "NMF Feature"){
      req(input$pt_nmfgene)
      idx <- match(colnames(plot_data), colnames(transform_data$data))
      gene_data <- transform_data$data[input$pt_nmfgene, idx] %>% as.numeric
      col <- create.brewer.color(gene_data, num = 9, name="YlOrRd")
      group <- input$pt_nmfgene
    }else if(input$pt_col == "Filename"){
      filename <- colnames(plot_data)
      filename <- sapply(strsplit(filename, "_"), function(x) paste0(x[as.numeric(input$pt_file_grp)], collapse = "_"))
      col <- create.brewer.color(filename, length(unique(filename)), "naikai") %>% as.matrix
      group <- filename %>% as.matrix

      # add option to manual select the order of coloring
      if(input$pt_sel_grp_order == "Manual"){
        col <- col[, 1]
        col.num <- as.numeric(factor(col, levels = unique(col) %>% as.character))
        col.lvl <- unique(col)
        col.lvl <- col.lvl[as.numeric(input$pt_grp_order)]
        col <- col.lvl[col.num] %>% as.matrix
      }
    }else if(input$pt_col == "GeneExpr"){
      req(input$pt_allgene)
      idx <- match(colnames(plot_data), colnames(transform_data$data))
      gene_data <- transform_data$data[input$pt_allgene, idx] %>% as.numeric
      col <- create.brewer.color(gene_data, num = 9, name="YlOrRd")
      group <- input$pt_allgene
    }else if(input$pt_col == "ReadDepth"){
      req(rawdata())
      idx <- match(colnames(plot_data), colnames(rawdata()))
      gene_data <- rawdata()[, idx] %>% colSums() %>% as.numeric
      col <- create.brewer.color(gene_data, num = 9, name="YlOrRd")
      group <- "Read Depth"
    }else if(input$pt_col == "NumExpresGenes"){
      req(rawdata())
      idx <- match(colnames(plot_data), colnames(rawdata()))
      gene_data <- colSums(rawdata()[, idx]>0) %>% as.numeric
      col <- create.brewer.color(gene_data, num = 9, name="YlOrRd")
      group <- "Num Expressed Genes"
    }else{
      warning("Wrong point color assigning method!")
    }

    res <- list()
    res[['color']] <- col
    res[['group']] <- group
    return(res)
  })

  pca_data <- reactive({
    data <- plot_data()[['plot_data']]
    withProgress(message = 'Calculating PCs', value = NULL, {
      pc <- prcomp(t(data), center = TRUE)

      vars<- pc$sdev^2
      vars<- signif(vars/sum(vars) * 100, 2)
      title = paste(paste0("PC", 1:length(vars)), paste0("(", vars, "% Var)"))

      projection <- as.data.frame(pc$x)
    })

    res <- list()
    res[['proj']] <- projection
    res[['title']] <- title
    return(res)
  })

  #' PCA plot
  output$pcaPlot_2D <- renderPlotly({
    projection <- pca_data()[['proj']]
    title <- pca_data()[['title']]

    validate(
      need(as.numeric(input$pca_x) != as.numeric(input$pca_y), "PC on the X-axis is the same as PC on the Y-axis.\nPlease select different number for these two axes")
    )
    idx <- c(as.numeric(input$pca_x), as.numeric(input$pca_y))

    withProgress(message = 'Generating 2d PCA plot', value = NULL, {
      projection <- projection %>%
        dplyr::select(idx)  %>%
        set_colnames(c("pca_x", "pca_y")) %>%
        mutate(color = point_col()[['color']] %>% as.character(),
               group = point_col()[['group']] %>% as.character(),
               Sample = colnames(plot_data()[['plot_data']])) %>%
        arrange(group)

      if(input$plot_label){
        t <- list( size=input$plot_label_size, color=toRGB("grey50") )
        p <- plot_ly(projection, x = ~pca_x, y = ~pca_y, type="scatter", mode="markers+text",
                     color = ~group, colors = ~unique(color),
                     text = ~Sample, hoverinfo="text", textfont=t, textposition="top middle",
                     marker = list(color = ~color, size=input$plot_point_size+3, opacity=input$plot_point_alpha))
      }else{
        p <- plot_ly(projection, x = ~pca_x, y = ~pca_y, type="scatter", mode="markers",
                     color = ~group, colors = ~unique(color),
                     text = ~Sample, hoverinfo="text",
                     marker = list(color = ~color, size=input$plot_point_size+3, opacity=input$plot_point_alpha)) #%>%
        # toWebGL()
      }
      p %>% layout(showlegend = input$plot_legend,
                   xaxis = list(title = title[idx[1]], zeroline=FALSE),
                   yaxis = list(title = title[idx[2]], zeroline=FALSE))
    })
  })
  output$pcaPlot_3D <- renderPlotly({
    projection <- pca_data()[['proj']]
    title <- pca_data()[['title']]
    idx <- c(1,2,3)

    withProgress(message = 'Generating 3d PCA plot', value = NULL, {
      projection <- projection %>%
        dplyr::select(idx)  %>%
        set_colnames(c("pca_x", "pca_y", "pca_z")) %>%
        mutate(color = point_col()[['color']] %>% as.character(),
               group = point_col()[['group']] %>% as.character(),
               Sample = colnames(plot_data()[['plot_data']])) %>%
        arrange(group)

      if(input$plot_label){
        p <- plot_ly(data=projection, x = ~pca_x, y = ~pca_y, z = ~pca_z, type="scatter3d", mode="markers+text",
                      color = ~group, colors = ~unique(color), text = ~Sample, hoverinfo="text",
                      marker = list(color = ~color, size=input$plot_point_size, opacity=input$plot_point_alpha))
      }else{
        p <- plot_ly(data=projection, x = ~pca_x, y = ~pca_y, z = ~pca_z, type="scatter3d", mode="markers",
                color = ~group, colors = ~unique(color), text = ~Sample, hoverinfo="text",
                marker = list(color = ~color, size=input$plot_point_size, opacity=input$plot_point_alpha)) #%>%
        # toWebGL()
      }

      p %>% layout(showlegend = input$plot_legend,
                   scene = list(
                     xaxis = list(title = title[1]),
                     yaxis = list(title = title[2]),
                     zaxis = list(title = title[3]))
      )
    })
  })

  #' t-SNE plot
  cores <- reactive({
    avail_cores <- parallel::detectCores()
    if(avail_cores==1){
      return(1)
    }else if(as.numeric(input$ncores) <= avail_cores){
      return(as.numeric(input$ncores))
    }else{
      return(avail_cores)
    }
  })
  tsne_2d <- reactiveValues(data=NULL)
  tsne_3d <- reactiveValues(data=NULL)

  observeEvent(input$runtSNE, {
    nsamples <- ncol(plot_data()[['plot_data']])
    if(input$tsne_perplexity*3 > (nsamples-1)) {
      createAlert(session, "visualperplexityAlert", "visualperplexityAlert1", title = "WARNING", style = "warning",
                  content = paste("Perpleixty", input$tsne_perplexity,
                                  "is too large compared to num of samples", nsamples),
                  append = FALSE)
    }else {
      closeAlert(session, "visualperplexityAlert1")

      plot_data <- plot_data()[['plot_data']]
      withProgress(message = 'Running t-SNE', value = NULL, {
        incProgress(1/3, detail = "For t-SNE 2D")
        tsne_2d$data <- run_tsne(plot_data, iter=input$tsne_iter, dims=2,
                                 perplexity=input$tsne_perplexity, cores=cores())
        incProgress(2/3, detail = "For t-SNE 3D")
        tsne_3d$data <- run_tsne(plot_data, iter=input$tsne_iter, dims=3,
                                 perplexity=input$tsne_perplexity, cores=cores())
      })
    }
  })

  plot_tsne_2d <- reactive({
    nsamples <- ncol(plot_data()[['plot_data']])
    validate(
      need(!is.null(tsne_2d$data), "Please click 'Run t-SNE' button")
    )
    tsne_out <- isolate(tsne_2d$data)
    withProgress(message = 'Genrating t-SNE plot', value = 0, {
      incProgress(2/3, detail = "Usually takes 15~20 seconds")
      color <- "steelblue"

      projection <- parse_tsne_res(tsne_out) %>%
        data.frame %>%
        mutate(group = point_col()[['group']] %>% as.character(),
               color = point_col()[['color']] %>% as.character(),
               Sample = rownames(.)) %>%
        arrange(group)

      min.cost <- signif(tsne_out$itercosts[length(tsne_out$itercosts)], 2)
      if(input$plot_label){
        t <- list( size=input$plot_label_size, color=toRGB("grey50") )
        p <- plot_ly(projection, x = ~x, y = ~y, type="scatter", mode="markers+text",
                     color = ~group, colors = ~unique(color),
                     text= ~Sample, hoverinfo="text", textposition="top middle", textfont=t,
                     marker = list(color = ~color, size=input$plot_point_size+3, opacity=input$plot_point_alpha))
      }else{
        p <- plot_ly(projection, x = ~x, y = ~y, type="scatter", mode="markers", text= ~Sample, hoverinfo="text",
                     color = ~group, colors = ~unique(color),
                     marker = list(color = ~color, size=input$plot_point_size+3, opacity=input$plot_point_alpha)) #%>%
          # toWebGL()
      }
      p %>% layout(showlegend = input$plot_legend,
                   title = title,
                   xaxis = list(title = "Component 1", zeroline=FALSE),
                   yaxis = list(title = "Component 2", zeroline=FALSE))
    })
  })
  output$tsneplot_2d <- renderPlotly({
    plot_tsne_2d()
  })
  output$tsneplot_3d <- renderPlotly({
    nsamples <- ncol(plot_data()[['plot_data']])
    if(is.null(tsne_3d$data)) return ()

    projection <- isolate(as.data.frame(tsne_3d$data$Y)) %>%
        set_colnames(c("tsne_1", "tsne_2", "tsne_3")) %>%
        mutate(group = point_col()[['group']] %>% as.character(),
               color = point_col()[['color']] %>% as.character(),
               Sample = rownames(.)) %>%
        arrange(group)

    if(input$plot_label){
      p <- plot_ly(data=projection, x = ~tsne_1, y = ~tsne_2, z = ~tsne_3, type="scatter3d", mode="markers+text",
                   color = ~group, colors = ~unique(color), text= ~Sample, hoverinfo="text",
                   marker = list(color = ~color, size=input$plot_point_size, opacity=input$plot_point_alpha))
    }else{
      p <- plot_ly(data=projection, x = ~tsne_1, y = ~tsne_2, z = ~tsne_3, type="scatter3d", mode="markers",
                   color = ~group, colors = ~unique(color), text= ~Sample, hoverinfo="text",
                   marker = list(color = ~color, size=input$plot_point_size, opacity=input$plot_point_alpha)) %>%
        toWebGL()
    }

    p %>% layout(showlegend = input$plot_legend,
                 scene = list(
                   xaxis = list(title = "Component 1"),
                   yaxis = list(title = "Component 2"),
                   zaxis = list(title = "Component 3"))
    )
  })


  ### DESeq2 ###
  observeEvent(input$runNMF, {
    validate(
      need(is.numeric(input$num_cluster), "num of cluster is not a valid numeric number")
    )
    updateSelectizeInput(session, 'de_group1',
                         server = TRUE,
                         choices = as.character(paste0("NMF", sort(unique(nmf_groups()$nmf_subtypes)))),
                         selected = "NMF1"
    )
    updateSelectizeInput(session, 'de_group2',
                         server = TRUE,
                         choices = as.character(paste0("NMF", sort(unique(nmf_groups()$nmf_subtypes)))),
                         selected = "NMF2"
    )
  })

  deseq_res <- eventReactive(input$runDESeq, {
    if(input$selectfile == "saved"){
      if(!is.null(rda()$dds)){
        return(rda()$dds)
      }
    }

    rawdata <- rawdata()
    if(input$de_conv2int)
      rawdata <- round(rawdata)

    if(!all(sapply(rawdata, is.whole))){
      createAlert(session, "alert", "exampleAlert", title = "WARNING", style = "warning",
                  content = "Data contains non-integer value, Please use raw count data for DESeq2. <br>
                             Or Use the option above to 'Force convert to integer' <b>(Results are not statistically rigorous)</b>",
                  append = FALSE)
      return(NULL)
    }else{
      closeAlert(session, "exampleAlert")
    }

    # since the user might have filtered samples during the QC steps,
    # rematch sample_ID and only kept the ones that are presented in NMF res
    idx <- match(nmf_groups()$Sample_ID, colnames(rawdata))
    rawdata <- rawdata[, idx]

    register(BiocParallel::MulticoreParam(workers = cores()))
    colData <- data.frame(Group = paste0("NMF", nmf_groups()$nmf_subtypes))
    ddsfeatureCounts <- DESeq2::DESeqDataSetFromMatrix(countData = rawdata,
                                                       colData = colData,
                                                       design = ~ Group)
    ### Set betaPrior=FALSE to go with MLE LFC to get simple LFC = (avg in group2/ avg in group1)
    # dds <- DESeq(ddsfeatureCounts, parallel=T, betaPrior=FALSE)
    withProgress(message = 'Running DESeq2', value = NULL, {
      dds <- DESeq2::DESeq(ddsfeatureCounts, parallel=T)
    })
    return(dds)
  })

  filt_deseq_res <- reactive({
    req(deseq_res())
    withProgress(message = 'Summarizing DESeq2 result', value = NULL, {
      register(BiocParallel::MulticoreParam(workers = 8))
      DESeq2::results(deseq_res(), parallel = TRUE,
                      contrast = c("Group", input$de_group1, input$de_group2)) %>%
        as.data.frame %>%
        subset(padj <= input$de_alpha) %>%
        .[order(.$padj), ]
    })
  })

  output$deseq_table <- DT::renderDataTable({
    validate(
      need(nrow(filt_deseq_res()) > 0, "No significant DE results found, Maybe try lower your alpha threshold?")
    )
    filename <- fileinfo()[['name']]
    GeneCard <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                       rownames(filt_deseq_res()), "'>", rownames(filt_deseq_res()), "</a>")
    filt_res <- data.frame(Gene=GeneCard,
                           filt_deseq_res())
    DT::datatable(filt_res, selection='single', rownames=FALSE,
                  extensions = 'Buttons',
                  escape = FALSE,
                  options = list(dom = 'Bfrtip',
                                 buttons =
                                   list('colvis', 'copy', 'print', list(
                                     extend = 'collection',
                                     buttons = list(list(extend='csv',
                                                         filename = filename),
                                                    list(extend='excel',
                                                         filename = filename)),
                                     text = 'Download'
                                   )),
                                 scrollY = TRUE,
                                 pageLength = 12,
                                 autoWidth = TRUE
                  )
    ) %>% formatRound(2:ncol(filt_res), 3)
  }, server = FALSE)

  output$deseq_boxplot <- renderPlotly({
    req(deseq_res())
    validate(
      need(!is.null(input$deseq_table_rows_selected), "Please select a gene")
    )
    withProgress(message = 'Generating boxplot', value = NULL, {
      norm.data <- DESeq2::counts(deseq_res(), normalized=TRUE)
      gene <- rownames(filt_deseq_res())[input$deseq_table_rows_selected]
      gene.data <- norm.data[gene, ]
      gene.data <- data.frame(Sample = names(gene.data),
                              Expr = gene.data,
                              NMF = factor(paste0("NMF", nmf_groups()$nmf_subtypes)))

      # match the coloring from PCA, t-SNE, and heatmap
      num_clus <- gene.data$NMF %>% levels %>% length
      mycolor <- create.brewer.color(1:num_clus, num=num_clus, name="naikai")

      a <- ggplot(data=gene.data, aes(x=NMF, y=Expr, color=NMF)) +
        geom_boxplot(aes(color = NMF)) +
        geom_jitter(aes(color=NMF, text=Sample), alpha=0.6, width=0.1) +
        theme_bw() +
        scale_colour_manual(values = mycolor) +
        labs(title=gene, x="", y="Expression", colour="")
      p <- plotly_build(a)
      # remove outliers for plotly boxplot
      p$data <- lapply(p$data, FUN = function(x){
        if(x$type == "box"){
          x$marker = list(opacity = 0)
        }
        return(x)
      })
      p$layout$margin$l <- p$layout$margin$l + 15
      p$layout$annotations[[1]]$x <- -0.1
      p$layout$xaxis$tickfont$size <- 8
      p$layout$annotations[[1]]$font$size <- 11
      p$layout$yaxis$tickfont$size <- 10
      p
    })
  })


  ### Pathway analysis ###
  nmf_group_feature_rank <- reactive({
    validate(
      need(is.numeric(input$num_cluster), "num of cluster is not a numeric value") %then%
      need(!is.null(nmf_groups()), "Please run NMF first")
    )
    rawdata <- transform_data$data
    # rawdata <- rawdata[rowMeans(rawdata) >= input$min_rowMean, ]
    nmf_group <- nmf_groups()$nmf_subtypes
    num_nmf_subtypes <- nmf_group %>% unique %>% length

    if(input$rank_method=="featureScore"){
      print('Do we want to use featureScore, there are only few genes')
    }else if(input$rank_method=="logFC"){
      nmf_feature_rank <- data.frame(matrix(NA_real_, nrow=nrow(rawdata), ncol=num_nmf_subtypes))

      for(i in 1:num_nmf_subtypes){
        idx <- nmf_group==i
        nmf_feature_rank[, i] <- rawdata %>% as.data.frame %>%
                                  dplyr::mutate(logFC = log(rowMeans(.[idx])/rowMeans(.[!idx]))) %>%
                                  dplyr::select(logFC)
      }
      colnames(nmf_feature_rank) <- paste0("NMF", 1:num_nmf_subtypes)
      rownames(nmf_feature_rank) <- rownames(rawdata)
    }
    return(nmf_feature_rank)
  })

  geneset.file <- reactive({
    validate(
      need(!is.null(input$predefined_list), "Please select at least one predefined gene sets")
    )
    file <- pre_geneset[[ input$predefined_list ]]
    return (file)
  })

  ncpus <- reactive({
    cpus <- c(16,8,6,4,2,1)
    cpus[min(which(sapply(cpus, function(x) input$nPerm%%x)==0))]
  })

  observeEvent(input$runNMF, {
    updateSelectizeInput(session, 'pathway_group',
                         server = TRUE,
                         choices = as.character(paste0("NMF", sort(unique(nmf_groups()$nmf_subtypes))))
    )
  })
  pathview.species <- reactive({
    species <- input$species
    if(species == "Human" | species == "human" | species == "hg19" | species == "hsa"){
      pathview.species <- "hsa"
    }else if(species == "Mouse" | species == "mouse" | species == "mm9" | species == "mmu"){
      pathview.species <- "mmu"
    }else if(species == "Drosophila" | species == "dm3" | species == "dm6"){
      pathview.species <- "dm3"
    }
  })
  id.org <- reactive({
    species <- input$species
    if(species == "Human" | species == "human" | species == "hg19" | species == "hsa"){
      id.org <- "Hs"
    }else if(species == "Mouse" | species == "mouse" | species == "mm9" | species == "mmu"){
      id.org <- "Mm"
    }
  })
  kegg.gs <- reactive({
    species <- input$species
    withProgress(message = 'Extracting KEGG pathway information', value = NULL, {
      if(species == "Human" | species == "human" | species == "hg19" | species == "hsa"){
        kg.human<- kegg.gsets("human")
        kegg.gs<- kg.human$kg.sets[kg.human$sigmet.idx]
      }else if(species == "Mouse" | species == "mouse" | species == "mm9" | species == "mmu"){
        kg.mouse<- kegg.gsets("mouse")
        kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
      }
    })
    return(kegg.gs)
  })
  go.gs <- reactive({
    species <- input$species
    withProgress(message = 'Extracting GO term information', value = NULL, {
      if(species == "Human" | species == "human" | species == "hg19" | species == "hsa"){
        go.gs <- go.gsets(species="human", pkg.name=NULL, id.type = "eg")
      }else if(species == "Mouse" | species == "mouse" | species == "mm9" | species == "mmu"){
        go.gs <- go.gsets(species="mouse", pkg.name=NULL, id.type = "eg")
      }
    })
    return(go.gs)
  })

  gage.exp.fc <- eventReactive(input$run_enrich, {
    rankdata <- nmf_group_feature_rank()
    num_groups <- ncol(rankdata)
    emp <- matrix(NA, nrow(rankdata), 1)
    res <- rep(list(emp), num_groups)
    names(res) <- colnames(rankdata)

    withProgress(message = 'Calculating LogFC between groups', value = NULL, {
      for (i in 1:num_groups){
        # In logFC, error might come from the fact that some of them are Inf or -Inf
        idx <- is.finite(rankdata[, i])
        sub_rankdata <- rankdata[idx, i, drop=F]
        # ID conversion
        id.map.sym2eg <- id2eg(ids=rownames(sub_rankdata), category = "SYMBOL", org=id.org())
        gene.entrez <- mol.sum(mol.data = sub_rankdata, id.map = id.map.sym2eg, sum.method = "mean")
        res[[i]] <- gene.entrez[, 1]
        res[[i]] <- res[[i]][!is.na(res[[i]])]
      }
    })
    return(res)
  })

  go.Res <- reactive({
    gage.exp.fc <- gage.exp.fc()
    gsets <- NULL
    res <- rep(list(NA), length(gage.exp.fc))
    names(res) <- names(gage.exp.fc)

    if(input$EnrichType == 'GO'){
      gsets <- go.gs()
      gsets <- gsets$go.sets[ gsets$go.subs[[input$goTerm]] ]
    }else if(input$EnrichType == 'KEGG'){
      gsets <- kegg.gs()
    }

    validate(
      need(length(gsets)>0, "Did not find any Gene Sets, maybe your gene count data is not from 'Human'?")
    )

    withProgress(message = 'Running enrichment analysis', value = NULL, {
      for (i in 1:length(gage.exp.fc)){
        exp.fc <- gage.exp.fc[[i]]
        goRes <- gage(exp.fc, gsets = gsets, ref = NULL, samp = NULL, same.dir=ifelse(input$samedir, TRUE, FALSE))
        res[[i]] <- goRes
      }
    })
    return(res)
  })

  go_summary <- reactive({
    goRes <- go.Res()
    genes <- goRes[[1]]$greater %>% rownames
    res <- list()
    num_genes <- 5
    hi_res <- lo_res <- matrix(NA, num_genes * length(goRes), length(goRes))
    colnames(hi_res) <- colnames(lo_res) <- names(goRes)
    total_hi_names <- total_lo_names <- rep(NA, num_genes * length(goRes))

    for(i in 1:length(goRes)){
      range <- ((i-1) * num_genes + 1) : (i * num_genes)
      hi_genenames <- goRes[[i]]$greater %>% rownames %>% head(n = num_genes)
      hi_res[range, ] <- sapply(goRes, function(x) x$greater[hi_genenames, "p.val"])
      total_hi_names[range] <- hi_genenames

      if(!is.na(goRes[[i]]$less[1])){
        lo_genenames <- goRes[[i]]$less %>% rownames %>% head(n = num_genes)
        lo_res[range, ] <- sapply(goRes, function(x) x$less[lo_genenames, "p.val"])
        total_lo_names[range] <- lo_genenames
      }
    }
    rownames(hi_res) <- total_hi_names
    rownames(lo_res) <- total_lo_names
    # save only the unique GO:Term
    res[["hi_sub_pval"]] <- hi_res[!duplicated(rownames(hi_res)), ]
    res[["lo_sub_pval"]] <- lo_res[!duplicated(rownames(lo_res)), ]

    # save total results
    res[["hi_total_pval"]] <- sapply(goRes, function(x) x$greater[genes, "p.val"])
    res[["hi_total_qval"]] <- sapply(goRes, function(x) x$greater[genes, "q.val"])
    res[["lo_total_pval"]] <- sapply(goRes, function(x) x$less[genes, "p.val"])
    res[["lo_total_qval"]] <- sapply(goRes, function(x) x$less[genes, "q.val"])

    return(res)
  })

  output$goplot_hi <- renderPlotly({
    req(go_summary())
    validate(
      need(!is.na(go_summary()[["hi_sub_pval"]]), "Did not find any GO terms for greater, maybe you select the wrong 'species'?")
    )
    withProgress(message = 'Generating summary plot for GO term (Hi)', value = NULL, {
      aa <- isolate(
        melt(go_summary()[["hi_sub_pval"]]) %>%
        set_colnames(c("GO", "NMF", "P")) %>%
        mutate(P = replace(P, is.na(P), 1))
      )

      if(input$go_logscale){
        p <- plot_ly(aa, x = ~NMF, y = ~GO, z = ~log10(P), type = "heatmap", colors = input$go_colorscale)
      }else{
        p <- plot_ly(aa, x = ~NMF, y = ~GO, z = ~P, type = "heatmap", colors = input$go_colorscale)
      }

      m <- list( l = input$go_leftmargin, r = 10, b = 40, t = 20)
      l <- list(margin = m,
                xaxis = list(title = "", tickfont=list(family = "Arial", size = input$go_xfontsize), zeroline=FALSE),
                yaxis = list(title = "", tickfont=list(family = "Arial", size = input$go_yfontsize), zeroline=FALSE))

      colorbar_format <- list(family = "Arial", size = 8)
      naikai <- plotly_build(p)
      naikai$x$data[[1]]$colorbar <- list(title = "P", thickness = input$go_barwidth, len=0.5, titlefont=colorbar_format, tickfont=colorbar_format)
      naikai$x$layout <- l
      naikai
    })
  })

  output$goplot_lo <- renderPlotly({
    req(go_summary())
    validate(
      need(!is.na(go_summary()[["lo_sub_pval"]]), "Did not find any GO terms for less, maybe you select the wrong 'species'?")
    )
    withProgress(message = 'Generating summary plot for GO term (Low)', value = NULL, {
      aa <- isolate(
        melt(go_summary()[["lo_sub_pval"]]) %>%
        set_colnames(c("GO", "NMF", "P")) %>%
        mutate(P = replace(P, is.na(P), 1))
      )

      if(input$go_logscale){
        p <- plot_ly(aa, x = ~NMF, y = ~GO, z = ~log10(P), type = "heatmap", colors = input$go_colorscale)
      }else{
        p <- plot_ly(aa, x = ~NMF, y = ~GO, z = ~P, type = "heatmap", colors = input$go_colorscale)
      }

      m <- list( l = input$go_leftmargin, r = 10, b = 40, t = 20)
      l <- list(margin = m,
                   xaxis = list(title = "", tickfont=list(family = "Arial", size = input$go_xfontsize), zeroline=FALSE),
                   yaxis = list(title = "", tickfont=list(family = "Arial", size = input$go_yfontsize), zeroline=FALSE))

      colorbar_format <- list(family = "Arial", size = 8)
      naikai <- plotly_build(p)
      naikai$x$data[[1]]$colorbar <- list(title = "P", thickness = input$go_barwidth, len=0.5, titlefont=colorbar_format, tickfont=colorbar_format)
      naikai$x$layout <- l
      naikai
    })
  })

  go_table <- reactive({
    goRes <- go.Res()
    validate(
      need(input$pathway_group!="", "Please select the enrichment result for the group you are interest in")
    )
    goRes <- goRes[[input$pathway_group]]
    if (is.null(goRes$greater)){
      stop("No significant gene sets found using this q.val cutoff")
    }else{
      return(goRes$greater)
    }
  })

  output$go_summary <- DT::renderDataTable({
    filename <- paste(input$EnrichType, input$pathway_group, input$goTerm, sep = "_")
    idx <- go_table()[, "q.val"] <= input$qval_cutoff
    go_table <- go_table()[which(!is.na(idx) & idx), ]
    if(nrow(go_table) == 0){
      stop("No significant gene sets found using this q.val cutoff")
    }
    go_table <- go_table[, c(3,4,2,5)]
    datatable(go_table, selection = 'single',
              extensions = c('Buttons', 'Responsive'),
              options = list(dom = 'Bfrtip',
                             buttons =
                               list('colvis', 'copy', 'print', list(
                                 extend = 'collection',
                                 buttons = list(list(extend='csv',
                                                     filename = filename),
                                                list(extend='excel',
                                                     filename = filename)),
                                 text = 'Download'
                               )),
                             pageLength = 6,
                             autoWidth = TRUE
                             # columnDefs = list(list(width = '10px', targets = "_all"))
              )
    ) %>% formatRound(1:3, 3)
  }, server = FALSE)


  output$pathviewImg <- renderImage({
    go_table <- go_table()
    s.idx <- input$go_summary_rows_selected[length(input$go_summary_rows_selected)]
    validate(
      need(length(s.idx)>0, "Please select one KEGG pathway from above.")
    )
    withProgress(message = 'Rendering pathview map', value = NULL, {
      path.pid <- substr(rownames(go_table)[s.idx], 1, 8)
      out.suffix <- "nmf"
      pathview(gene.data=gage.exp.fc()[[input$pathway_group]], pathway.id=path.pid, species=pathview.species(), out.suffix=out.suffix)
    })

    filename <- path.pid
    outfile <- paste(filename, out.suffix, 'png', sep=".")
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 1000,
         height = 777,
         alt = "This is alternate text")
  }, deleteFile = FALSE)

  nmf_feature_rank_for_gsea <- reactive({
    species <- input$species
    nmf_feature_rank <- nmf_group_feature_rank()
    # check if genenames contains all CAPITAL LETTERS (after removing all numeric characters)
    cap_genes <- grepl("^[[:upper:]]+$", gsub('[[:digit:]]+', '', rownames(nmf_feature_rank))) %>% sum
    cap_pct <- cap_genes/ nrow(nmf_feature_rank)
    # some of the genes might contains '.', '-', 'orf', etc, so we use percentage of genes qualified as filtering criteria

    # if(species == "Mouse" | species == "mouse" | species == "mm9" | species == "mmu"){
    if(cap_pct < 0.1){
      n <- 2
      withProgress(message = 'Assuming the input gene symbols are from Mouse, Converting to Human..', value = 0, {
        incProgress(1/n, detail = "Takes around 10~15 seconds")
        mouse_human_convert_genes <- gene_convert_mouse_to_human(rownames(nmf_feature_rank))
        ## After conversion, some of the mouse genes might point to the same human gene
        ## So select thoide mouse genes that has highest mean expression across groups?
        idx <- duplicated(mouse_human_convert_genes$HGNC.symbol)
        if (sum(idx)>0){
          warning("Some of the mouse genes are mapped back to the same human genes..")
        }
        max_dup_mean_lfc_genes <-
          cbind(mouse_human_convert_genes, nmf_feature_rank) %>%
          mutate(avg_lfc=apply(nmf_feature_rank, 1, mean),
                 abs_lfc=apply(nmf_feature_rank, 1, function(x) max(abs(x)))
          ) %>%
          dplyr::group_by(HGNC.symbol) %>%
          dplyr::filter(min_rank(desc(avg_lfc))==1) %>%
          dplyr::slice(which.max(abs_lfc)) %>%
          ungroup    # this will give Error: corrupt 'grouped_df', workaround is to ungroup it below
        print(dim(max_dup_mean_lfc_genes))
        ## Return final de-duplicated genes only
        nmf_feature_rank <- max_dup_mean_lfc_genes %>%
          dplyr::select(starts_with("NMF")) %>%
          as.data.frame
        rownames(nmf_feature_rank) <- max_dup_mean_lfc_genes$HGNC.symbol
      })
    }
    return(nmf_feature_rank)
  })

  output$downloadPathData <- downloadHandler(
    filename <- function() {
      paste( file_prefix(), "pathway.zip", sep=".")
    },
    content = function(file) {
      tmpdir <- tempdir()
      current_dir <- getwd()
      setwd(tempdir())

      ### Save GO Term results
      go_summary <- go_summary()
      total_idx <- grep("total", names(go_summary))
      go_filenames <- paste(file_prefix(), paste0("Enrichment_", names(go_summary)[total_idx]), "txt", sep=".")
      for(i in 1:length(total_idx)){
        write.table(go_summary[[ total_idx[i] ]], go_filenames[i], quote=F, row.names = TRUE, col.names=NA, sep="\t")
      }
      ### Rank file for GSEA
      # modify it to save one file for each rank ###
      nmf_group_feature_rank <- nmf_feature_rank_for_gsea()
      num_samples <- ncol(nmf_group_feature_rank)
      filenames <- paste(file_prefix(), paste0("group_feature_rank", 1:num_samples), "rnk", sep=".")
      for (i in 1:num_samples){
        sub.data <- data.frame(Gene=rownames(nmf_group_feature_rank), LogFC=nmf_group_feature_rank[,i])
        idx <- is.finite(sub.data$LogFC)
        write.table(sub.data[idx, ], filenames[i], quote=F, row.names = F, sep="\t")
      }

      zip(zipfile=file, files=c(go_filenames, filenames))
      setwd(as.character(current_dir))
    },
    contentType = "application/zip"
  )

})
