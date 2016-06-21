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
      validate(
        need(!is.null(inFile), "Please upload a gene count data set")
      )
      filename <- basename(inFile$name)
      filepath <- inFile$datapath
    }else if (input$selectfile == "preload"){
      validate(
        need(input$inputdata!="", "Please select a gene count data set")
      )
      filename <- basename(input$inputdata)
      filepath <- file.path(raw_fd, input$inputdata)
    }

    fileinfo <- list()
    fileinfo[['name']] <- filename
    fileinfo[['path']] <- filepath
    return(fileinfo)
  })

  rawdata <- reactive({
    filepath <- fileinfo()[['path']]
    n <- 2
    withProgress(message = 'Loading data', value = 0, {
      incProgress(1/n, detail = "Usually takes ~10 seconds")
      rawdata <- myfread.table(filepath, check.platform=T, sep=input$sep, detect.file.ext=FALSE)
    })
    return (rawdata)
  })

  #' Raw data
  output$rawdatatbl = DT::renderDataTable({
    rawdata <- rawdata() %>% head(n=100)
    validate(
      need(ncol(rawdata)>1,
           "Please modify the parameters on the left and try again!") %then%
      need(!any(duplicated(t(rawdata))),
           "There are duplicated columns in your data, Please check and try again!")
    )
    DT::datatable(rawdata, rownames= TRUE,
                  options = list(scroll = TRUE,
                                 scrollX = TRUE,
                                 scrollY = TRUE,
                                 pageLength = 8
                  )
    )
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

  transform_data <- reactive({
    rawdata <- rawdata()
    transdata <- rmv_constant_0(rawdata, pct=0.9, minimum=0)

    n <- 2
    withProgress(message = 'Transforming data', value = 0, {
      incProgress(1/n, detail = "Takes around 5~10 seconds")

      if(input$transformdata=="takevst"){
        transdata <- vst(round(transdata))
      }else if(input$transformdata=="takelog"){
        validate(
          need(!any(transdata<0), "Please make sure all the numbers in the data are positive values")
        )
        transdata <- log2(transdata + 1)
      }

      if(input$normdata=="takeRPM"){
        transdata <- rpm(transdata, colSums(transdata))
      }else if(input$normdata=="takesizeNorm"){
        sf <- norm_factors(transdata)
        transdata <- t(t(transdata)/sf)
      }
    })
    return(transdata)
  })

  output$transdatatbl = DT::renderDataTable({
    transform_data <- transform_data()
    filename <- fileinfo()[['name']]
    DT::datatable(transform_data, rownames= TRUE,
                  extensions = 'Buttons',
                  options = list(dom = 'Bfrtip',
                                 buttons =
                                   list('colvis', 'copy', 'print', list(
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
                                 scrollY = TRUE,
                                 pageLength = 8,
                                 order=list(list(2,'desc'))
                  )
    ) %>% formatRound(1:ncol(transform_data), 2)
  }, server=TRUE)

  output$readDistrib = renderPlotly({
    total <- colSums(transform_data())
    samples <- names(total)
    plot_ly(x=samples, y=total, type = "bar", marker = list(color = toRGB("gray65")))
  })
  output$geneCoverage = renderPlotly({
    genecnt <- apply(transform_data(), 2, function(x) sum(x>0))
    samples <- names(genecnt)
    plot_ly(x=samples, y=genecnt, type = "bar", marker = list(color = toRGB("gray65")))
  })

  merged <- reactive({
    rawdata <- transform_data()
    genelist <- NULL

    if(!is.null(input$selected_samples) && sum(input$selected_samples %in% colnames(rawdata))==length(input$selected_samples)){
      rawdata <- rawdata[, input$selected_samples]
    }

    if(input$list_type=='Rank from data'){
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
  callModule(feature, "gene", reactive({ t(merged()) }))


  #' Sample correlation plot
  output$sampleCorPlot <- renderPlot({
    transform_data <- transform_data()
    n <- 2
    withProgress(message = 'Calculating correlation', value = 0, {
      incProgress(1/n, detail = "Takes around 10 seconds")
      M <- cor(transform_data)
      p.mat <- cor_mtest(transform_data)
    })
    # col<- colorRampPalette(c("green","white", "blue", "white", "red"))(100)
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA", "#79AEDD", "#FFFFFF", "#2E9988", "#2B4444"))

    withProgress(message = 'Plotting..', value = 0, {
      incProgress(1/n, detail = "Takes around 10 seconds")
      corrplot(M, method="color", col=rev(col(200)),
               type="upper", order="hclust",
               addCoef.col = "black", # Add coefficient of correlation
               tl.col="black", tl.srt=45, tl.cex=input$cor_sam_lab_cex, #Text label color and rotation
               number.cex = input$cor_num_lab_cex,
               p.mat = p.mat, sig.level = 0.01, insig = "blank",
               diag=FALSE
      )
    })
  }, height=700)

  #' Scatter plot
  observeEvent(transform_data(), {
    transform_data <- transform_data()
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
    transform_data <- transform_data()
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
    transform_data <- transform_data()
    x <- scatterData()$x
    y <- scatterData()$y
    reg = lm(y ~ x)
    modsum = summary(reg)
    r <- signif(cor(x,y), 3)
    R2 = signif(summary(reg)$r.squared, 3)
    textlab <- paste(paste0("x = ", signif(x,3)), paste0("y = ", signif(y,3)), rownames(transform_data), sep="<br>")
    n <- 2
    withProgress(message = 'Generating Scatter Plot', value = 0, {
      p <- plot_ly(x=x, y=y, text=textlab, type="scattergl", source = "scatter",
                   mode="markers", hoverinfo="text",
                   themes="Catherine",
                   marker = list(size=8, opacity=0.65)) %>%
            layout(title = paste("cor:", r),
                   xaxis = list(title = input$cor_sample1, zeroline=TRUE),
                   yaxis = list(title = input$cor_sample2, zeroline=TRUE),
                   showlegend=FALSE)
      p <- add_trace(p, x=plot_range, y=plot_range, type = "scattergl", mode = "lines", name = "Ref", line = list(width=2))
    })
    if(input$show_r2){
      p <- p %>% layout( annotations = list(x = x, y = y, text=paste0("R2=", R2), showarrow=FALSE) )
    }

    plot_range <- min(x,y):max(x,y)
    if(input$show_fc){
      p <- add_trace(p, x=plot_range, y=plot_range[plot_range<=(max(plot_range)-1)]+1, type = "scattergl", mode = "lines", name = "2-FC", line = list(width=1.5, color="#fec44f"))
      p <- add_trace(p, x=plot_range[plot_range>=1], y=plot_range[plot_range>=1]-1, type = "scattergl", mode = "lines", name = "2-FC", line = list(width=1.5, color="#fec44f"))
    }
    p
  })
  # output$dlScatterPlots <- downloadHandler(
  #   filename = function() { paste("test", '.png', sep='') },
  #   content = function(file) {
  #     plotly_IMAGE(patientCircleInput(), format = "png", out_file = file)
  #   }
  # )
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
                              rownames(transform_data()), "'>", rownames(transform_data()), "</a>")
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


  ### Output Prefix ###
  file_prefix <- reactive({
    if (input$selectfile == 'upload'){
      File <- input$file1$name
    }else if (input$selectfile == "preload"){
      File <- input$inputdata
    }
    ext <- file_ext(File)
    paste( gsub(paste0(".", ext), "", File),
           paste0("k", input$num_cluster),
           paste0(input$top.num,input$algorithm),
           paste0("nrun", input$nrun),
           paste0("predict_", input$predict),
           sep=".")
  })

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
    ptm <- proc.time()
    n <- 2
    withProgress(message = 'Running NMF', value = 0, {
      incProgress(1/n, detail = "Takes around 30~60 seconds")
      res <- myNMF(merged,
                   mode=input$mode,
                   cluster=input$num_cluster,
                   nrun=input$nrun,
                   algorithm = input$algorithm,
                   seed = nmf_seed()
      )
    })

    b <- proc.time() - ptm
    print("### Time to run NMF ### ")
    print(b)

    return(res)
  })

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
      need(class(nmf_res()) == "NMF.rank", "Seems like you haven't run NMF 'estim run' yet")
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
      need(class(nmf_res()) == "NMFfitX1", "Seems like you haven't run NMF 'real' run yet")
    )
    nmf_res <- nmf_res()
    nmf_plot(nmf_res, type=input$plottype, silorder=T,
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
        # nmf_estim_plot(nmf_res)
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
        nmf_plot(nmf_res, type="samples", silorder=T)
        nmf_plot(nmf_res, type="features", silorder=T,
                 subsetRow = ifelse(input$select_feature_num>0, as.numeric(input$select_feature_num), TRUE ))
        nmf_plot(nmf_res, type="consensus")
        if(!is.null(tsne_2d$data)){
          print(nmfplottsne())
        }
      }
      dev.off()
    }
  )

  ### NMF results
  nmf_groups <- reactive({
    validate(
      need(class(nmf_res()) == "NMFfitX1", "Seems like you haven't run NMF 'real' run yet")
    )
    nmf_res <- nmf_res()
    withProgress(message = 'Extracting Groups info', value = NULL, {
      nmf_extract_group(nmf_res, matchConseOrder = input$matchOrder, type=input$predict )
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
                                 pageLength = 8,
                                 order=list(list(1,'asc'))
                  )
    ) %>% formatRound(3:ncol(nmf_groups), 3)
  }, server=FALSE)

  # Add t-SNE plot next to NMF group
  observeEvent(input$run_nmftSNE, {
    heatmap_data <- heatmap_data()[['heatmap_data']]
    withProgress(message = 'Running t-SNE', value = NULL, {
      incProgress(1/3, detail = "For t-SNE 2D")
      try(
      tsne_2d$data <- run_tsne(heatmap_data, iter=input$tsne_iter, dims=2,
                               perplexity=input$nmftsne_perplexity, cores=cores())
      )
      incProgress(2/3, detail = "For t-SNE 3D")
      tsne_3d$data <- run_tsne(heatmap_data, iter=input$tsne_iter, dims=3,
                               perplexity=input$nmftsne_perplexity, cores=cores())
    })
  })
  nmfplottsne <- reactive({
    nsamples <- ncol(heatmap_data()[['heatmap_data']])
    validate(
      need(input$nmftsne_perplexity*3 < (nsamples-1),
           message = paste("Perpleixty", input$nmftsne_perplexity,
                           "is too large compared to num of samples", nsamples))
    )
    validate(
      need(input$mode=="real", "This part won't work for 'estim' module\nPlease Try NMF 'Real' run")
    )
    validate(
      need(!is.null(tsne_2d$data), "Please hit 'Run t-SNE' button")
    )
    color <- NULL
    # Will add an option to select by default column names from samples
    if(length(nmf_groups()$nmf_subtypes)>0){
      color <- nmf_groups()$nmf_subtypes
    }else if(length(ColSideColors()[['name']][,1])>0){
      color <- ColSideColors()[['name']][,1]
    }else{
      stop("Error: color scheme undefined")
    }
    # # default by user selection
    point_size <- input$plot_point_size - 6
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


  nmf_features <- reactive({
    validate(
      need(input$mode=="real", "This part won't work for eastimating k\nPlease Try NMF 'Realrun'") %then%
      need(class(nmf_res())=="NMFfitX1", "Please rerun NMF with 'real run'")
    )
    nmf_res <- nmf_res()
    # add more options to select features here. then compare with ClaNC
    # Add rank by MAD expression
    if(input$select_method=="default"){
      res <- nmf_extract_feature(nmf_res, method=input$select_method, manual.num = as.numeric(input$select_feature_num))
      validate(
        need(nrow(res)>1, "Not enought features are selected by default! \nMaybe try with more input features or try select features by Rank")
      )
    }else if(input$select_method=="rank"){
      res <- nmf_extract_feature(nmf_res, rawdata = merged(), method=input$select_method, manual.num = as.numeric(input$select_feature_num),
                          FScutoff=input$select_FScutoff, math=input$select_math)
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
      merged <- ext.func.data[data, nomatch=0]
      # Sort and drop unwanted columns
      merged <- as.data.frame(merged[order(Group, -featureScore),])
      # merged$hgnc_symbol <- NULL
      # merged <- merged[, c(1,4,2,7,6,3,5,8)]
      merged <- merged[, c(4,2,7,6,8)]
      merged <- unique(merged)
      merged$featureScore <- round(merged$featureScore, 3)
      merged$prob <- round(merged$prob, 3)
    })
    return(merged)
  })
  output$nmfFeatures <- DT::renderDataTable({
    filename <- fileinfo()[['name']]
    DT::datatable(nmf_features_annot(), rownames=FALSE, escape=-1,
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
                                 pageLength = 50,
                                 autoWidth = TRUE
                  )
    )
  }, server = FALSE)

  ori_plus_nmfResult <-reactive({
    rawdata <- rawdata()
    validate(
      need(!is.null(nmf_groups()), "Please run NMF and get the result first")
    )
    idx <- match(colnames(rawdata), nmf_groups()$Sample_ID)
    colnames(rawdata) <- paste(paste0("NMF",nmf_groups()$nmf_subtypes[idx]), colnames(rawdata), sep="_")
    return(rawdata)
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
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  })

  #' Heatmap setting
  geneListInfo <- reactive({
    rawdata <- transform_data()
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

  heatmap_data <- reactive({
    rawdata <- transform_data()
    genelist <- geneListInfo()[['list']]
    gene.groups <- geneListInfo()[['group']]

#     if(!is.null(input$selected_samples) && sum(input$selected_samples %in% colnames(rawdata))==length(input$selected_samples)){
#       rawdata <- rawdata[, input$selected_samples]
#     }
#     if(input$minExp>0){
#       if(input$minType=='Mean'){
#         rawdata <- rawdata[rowMeans(rawdata) >= input$minExp, ]
#       }else if(input$minType=="Median"){
#         rawdata <- rawdata[rowMedians(rawdata) >= input$minExp, ]
#       }else if (input$minType=="Sds"){
#         rawdata <- rawdata[rowSds(rawdata) >= input$minExp, ]
#       }else if (input$minType=="Min"){
#         rawdata <- rawdata[rowMin(rawdata) >= input$minExp, ]
#       }else if (input$minType=="Max"){
#         rawdata <- rawdata[rowMax(rawdata) >= input$minExp, ]
#       }
#     }

    # only select unique gene name on rawdata
    idx <- match(genelist, rownames(rawdata))
    heatmap_data <- rawdata[idx,]

    if(input$OrdCol == 'group'){    # Add options to sort by column groups
      column.names <- strsplit(colnames(heatmap_data), split="\\_")
      validate(
        need(input$sortcolumn_num <= length(column.names[[1]]), "Group num is too big, select a smaller number!")
      )
      idx <- order(sapply(column.names, function(x) x[[input$sortcolumn_num]]))
      heatmap_data <- heatmap_data[, idx]
    }else if(input$OrdCol == 'gene'){       # Add options to sort by expression value of specific gene
      validate(
        need(input$sortgene_by!="", "Please select a gene!")
      )
      idx <- order(subset(heatmap_data, rownames(heatmap_data)==input$sortgene_by))
      heatmap_data <- heatmap_data[, idx]
    }

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
    ColScheme <- c(input$ColScheme1, input$ColScheme2, input$ColScheme3, input$ColScheme4)
    return(ColScheme)
  })

  ColSideColors <- reactive({
    heatmap_data <- heatmap_data()[['heatmap_data']]
    res <- name_to_color(colnames(heatmap_data), split_pattern ="\\_",
                         num_color = input$ColSideColorsNum,
                         ColScheme = ColScheme()[1:input$ColSideColorsNum] )
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
    if(input$OrdCol=="hier"){
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

  ### Observe and update input selection ###
  observe({
    transform_data <- transform_data()
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
      dendrogram = "both",
      dendro.ord = "manual",
      anim_duration = 0
    )
  })

  plot_colopts <- reactive({
    c("Default", "Filename")
  })
  observe({
    updateSelectizeInput(session, 'pt_col',
                         server = TRUE,
                         choices = as.character(plot_colopts()),
                         selected = "Filename"
    )
  })
  observeEvent(input$runNMF,{
    updateSelectizeInput(session, 'pt_col',
                         server = TRUE,
                         choices = as.character(c(plot_colopts(), "NMF")),
                         selected = 'NMF'
    )
  })
  point_col <- reactive({
    col <- toRGB("steelblue")
    if(input$pt_col == "Default"){
      print("default color")
    }else if(input$pt_col == "NMF"){
      nmf_subtypes <- nmf_groups()$nmf_subtypes
      col <- create.brewer.color(nmf_subtypes, length(unique(nmf_subtypes)), "naikai")
    }else if(input$pt_col == "Filename"){
      col <- ColSideColors()[["color"]][, 1]
    }else{
      warning("Wrong point color assigning method!")
    }
    return(col)
  })

  #' PCA plot
  output$pcaPlot_2D <- renderPlotly({
    data <- heatmap_data()[['heatmap_data']]
    pc <- prcomp(t(data))
    projection <- as.data.frame(pc$x)
    pca_x <- as.numeric(input$pca_x)
    pca_y <- as.numeric(input$pca_y)

    if(input$plot_label){
      t <- list( size=input$plot_label_size, color=toRGB("grey50") )
      p <- plot_ly(x=projection[,pca_x], y=projection[,pca_y], mode="markers+text",
                   text=rownames(projection), hoverinfo="text", textfont=t, textposition="top middle",
                   # marker = list(color=toRGB("steelblue"), size=input$plot_point_size+3, opacity=input$plot_point_alpha))
                   marker = list(color=point_col(), size=input$plot_point_size+3, opacity=input$plot_point_alpha))
    }else{
      p <- plot_ly(x=projection[,pca_x], y=projection[,pca_y], mode="markers",
                   text=rownames(projection), hoverinfo="text",
                   # marker = list(color=toRGB("steelblue"), size=input$plot_point_size+3, opacity=input$plot_point_alpha))
                   marker = list(color=point_col(), size=input$plot_point_size+3, opacity=input$plot_point_alpha))
    }
    p %>% layout(xaxis = list(title = paste0("PC", pca_x), zeroline=FALSE),
                 yaxis = list(title = paste0("PC", pca_y), zeroline=FALSE))
  })
  output$pcaPlot_3D <- renderPlotly({
    data <- heatmap_data()[['heatmap_data']]
    pc <- prcomp(t(data))
    projection <- as.data.frame(pc$x)
    N <- nrow(projection)
    plot_ly(x=projection[,1], y=projection[,2], z=projection[,3], type="scatter3d", mode="markers",
                 # color=rownames(projection), text=rownames(projection), hoverinfo="text",
                 color=point_col(), text=rownames(projection), hoverinfo="text",
                 marker = list(color=point_col(), size=input$plot_point_size, opacity=input$plot_point_alpha)) %>%
        layout( scene = list(
                 xaxis = list(title = "PC1"),
                 yaxis = list(title = "PC2"),
                 zaxis = list(title = "PC3")))
  })

  #' t-SNE plot
  cores <- reactive({
    cpus <- c(16,8,6,4,2,1)
    # cpus[min(which(sapply(cpus, function(x) input$nPerm%%x)==0))]
    return(8)
  })
  tsne_2d <- reactiveValues(data=NULL)
  tsne_3d <- reactiveValues(data=NULL)
  observeEvent(input$runtSNE, {
    heatmap_data <- heatmap_data()[['heatmap_data']]
    withProgress(message = 'Running t-SNE', value = NULL, {
      incProgress(1/3, detail = "For t-SNE 2D")
      tsne_2d$data <- run_tsne(heatmap_data, iter=input$tsne_iter, dims=2,
                               perplexity=input$tsne_perplexity, cores=cores())
      incProgress(2/3, detail = "For t-SNE 3D")
      tsne_3d$data <- run_tsne(heatmap_data, iter=input$tsne_iter, dims=3,
                               perplexity=input$tsne_perplexity, cores=cores())
    })
  })
  plot_tsne_2d <- reactive({
    nsamples <- ncol(heatmap_data()[['heatmap_data']])
    validate(
      need(input$tsne_perplexity*3 < (nsamples-1),
           message = paste("Perpleixty", input$tsne_perplexity,
                            "is too large compared to num of samples", nsamples)) %then%
      need(!is.null(tsne_2d$data), "Please click 'Run t-SNE' button")
    )
    tsne_out <- isolate(tsne_2d$data)
    withProgress(message = 'Genrating t-SNE plot', value = 0, {
      incProgress(2/3, detail = "Usually takes 15~20 seconds")
      color <- "steelblue"

      projection <- parse_tsne_res(tsne_out)
      projection$color <- color
      min.cost <- signif(tsne_out$itercosts[length(tsne_out$itercosts)], 2)
      title <- paste("min.cost=", min.cost)
      colors <- create.brewer.color(projection$color, length(unique(color)), "naikai")

      if(input$plot_label){
        t <- list( size=input$plot_label_size, color=toRGB("grey50") )
        p <- plot_ly(projection, x=x, y=y, mode="markers+text",
                     text=rownames(projection), hoverinfo="text", textposition="top middle", textfont=t,
                     marker = list(color=point_col(), size=input$plot_point_size+3, opacity=input$plot_point_alpha))
      }else{
        p <- plot_ly(projection, x=x, y=y, mode="markers", text=rownames(projection), hoverinfo="text",
                     marker = list(color=point_col(), size=input$plot_point_size+3, opacity=input$plot_point_alpha))
      }
      p %>% layout(title = title,
                   xaxis = list(title = "Component 1", zeroline=FALSE),
                   yaxis = list(title = "Component 2", zeroline=FALSE))
    })
  })
  output$tsneplot_2d <- renderPlotly({
    plot_tsne_2d()
  })
  output$tsneplot_3d <- renderPlotly({
    nsamples <- ncol(heatmap_data()[['heatmap_data']])
    validate(
      need(input$tsne_perplexity*3 < (nsamples-1),
           message = paste("Perpleixty", input$tsne_perplexity,
                            "is too large compared to num of samples", nsamples))
    )
    if(is.null(tsne_3d$data)) return ()
    projection <- as.data.frame(tsne_3d$data$Y)
    labels <- rownames(projection)
    N <- nrow(projection)
    plot_ly(x=projection[,1], y=projection[,2], z=projection[,3], type="scatter3d", mode="markers",
                 color=point_col(), text=labels, hoverinfo="text",
                 # marker = list(size=input$plot_point_size, opacity=input$plot_point_alpha)) %>%
                 marker = list(color=point_col(), size=input$plot_point_size, opacity=input$plot_point_alpha)) %>%
        layout( scene = list(
                 xaxis = list(title = "Component 1"),
                 yaxis = list(title = "Component 2"),
                 zaxis = list(title = "Component 3")))
  })


  ### DESeq2 ###
  observeEvent(input$runNMF, {
    validate(
      need(is.numeric(input$num_cluster), "num of cluster is not a valid numeric number")
    )
    updateSelectizeInput(session, 'de_group1',
                         server = TRUE,
                         choices = as.character(paste0("NMF", 1:input$num_cluster)),
                         selected = "NMF1"
    )
    updateSelectizeInput(session, 'de_group2',
                         server = TRUE,
                         choices = as.character(paste0("NMF", 1:input$num_cluster)),
                         selected = "NMF2"
    )
  })

  deseq_res <- eventReactive(input$runDESeq, {
    rawdata <- rawdata()
    if(input$de_conv2int)
      rawdata <- round(rawdata)

    if(!all(sapply(rawdata, is.whole))){
      createAlert(session, "alert", "exampleAlert", title = "Oops", style = "warning",
                  content = "Data contains non-integer value, Please use raw count data for DESeq2. <br>
                             Or Use the option above to 'Force convert to integer'",
                  append = FALSE)
    }else{
      closeAlert(session, "exampleAlert")
    }

    cores <- parallel::detectCores()
    register(BiocParallel::MulticoreParam(workers = cores))
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
    DESeq2::results(deseq_res(), contrast = c("Group", input$de_group1, input$de_group2)) %>%
      as.data.frame %>%
      subset(padj <= input$de_alpha) %>%
      .[order(.$padj), ]
  })

  output$deseq_table <- DT::renderDataTable({
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
                                                         filename = filename),
                                                    list(extend='pdf',
                                                         filename= filename)),
                                     text = 'Download'
                                   )),
                                 scrollY = TRUE,
                                 pageLength = 12,
                                 autoWidth = TRUE
                  )
    ) %>% formatRound(2:ncol(filt_res), 3)
  }, server = TRUE)

  output$deseq_boxplot <- renderPlotly({
    validate(
      need(!is.null(input$deseq_table_rows_selected), "Please select a gene")
    )
    norm.data <- DESeq2::counts(deseq_res(), normalized=TRUE)
    gene <- rownames(filt_deseq_res())[input$deseq_table_rows_selected]
    gene.data <- norm.data[gene, ]
    gene.data <- data.frame(Sample = names(gene.data),
                            Expr = gene.data,
                            NMF = as.factor(paste0("NMF", nmf_groups()$nmf_subtypes)))
    gene.data$NMF <- factor(gene.data$NMF, levels = rev(levels(gene.data$NMF)))

    plot_ly(data = gene.data, x=NMF, y=Expr, type = "box", color = NMF,
            boxpoints = "all", jitter = 0.3, pointpos = 0) %>%
      layout(xaxis = list(title = "", zeroline=TRUE),
             yaxis = list(title = "Expression", zeroline=TRUE),
             showlegend=TRUE)
  })


  ### Pathway analysis ###
  nmf_group_feature_rank <- reactive({
    validate(
      need(is.numeric(input$num_cluster), "num of cluster is not a numeric value")
    )
    rawdata <- transform_data()
    rawdata <- rawdata[rowMeans(rawdata) >= input$min_rowMean, ]
    nmf.group <- nmf_groups()$nmf_subtypes
    if(input$rank_method=="featureScore"){
      print('Do we want to use featureScore, there are only few genes')
    }else if(input$rank_method=="logFC"){
      nmf_feature_rank <- data.frame(matrix(NA_real_, nrow=nrow(rawdata), ncol=input$num_cluster))

      for(i in 1:input$num_cluster){
        idx <- nmf.group==i
        nmf_feature_rank[, i] <- rawdata %>% as.data.frame %>%
                                  dplyr::mutate(logFC = log(rowMeans(.[idx])/rowMeans(.[!idx]))) %>%
                                  dplyr::select(logFC)
      }
      colnames(nmf_feature_rank) <- paste0("NMF", 1:input$num_cluster)
      rownames(nmf_feature_rank) <- rownames(rawdata)
    }
    ### Check whether they are mouse or human, then convert mouse to human if needed ###
#     if(input$species=="Mouse"){
#       n <- 2
#       withProgress(message = 'Converting Mouse gene symbols to Human..', value = 0, {
#         incProgress(1/n, detail = "Takes around 10~15 seconds")
#         mouse_human_convert_genes <- gene_convert_mouse_to_human(rownames(nmf_feature_rank))
#         ## After conversion, some of the mouse genes might point to the same human gene
#         ## So select thoide mouse genes that has highest mean expression across groups?
#         idx <- duplicated(mouse_human_convert_genes$HGNC.symbol)
#         if (sum(idx)>0){
#           warning("Some of the mouse genes are mapped back to the same human genes..")
#         }
#         print(dim(nmf_feature_rank))
#         max_dup_mean_lfc_genes <-
#           cbind(mouse_human_convert_genes, nmf_feature_rank) %>%
#           mutate(avg_lfc=apply(nmf_feature_rank, 1, mean),
#                  abs_lfc=apply(nmf_feature_rank, 1, function(x) max(abs(x)))
#           ) %>%
#           dplyr::group_by(HGNC.symbol) %>%
#           dplyr::filter(min_rank(desc(avg_lfc))==1) %>%
#           dplyr::slice(which.max(abs_lfc)) %>%
#           ungroup    # this will give Error: corrupt 'grouped_df', workaround is to ungroup it below
#         print(dim(max_dup_mean_lfc_genes))
#         ## Return final de-duplicated genes only
#         nmf_feature_rank <- max_dup_mean_lfc_genes %>%
#           dplyr::select(starts_with("NMF")) %>%
#           as.data.frame
#         rownames(nmf_feature_rank) <- max_dup_mean_lfc_genes$HGNC.symbol
#       })
#     }
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
                         choices = as.character(paste0("NMF", 1:input$num_cluster))
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
    # data(go.sets.hs)
    # data(go.subs.hs)
    # go.gs <- go.sets.hs[go.subs.hs[[input$goTerm]]]
    species <- input$species
    withProgress(message = 'Extracting GO term information', value = NULL, {
      if(species == "Human" | species == "human" | species == "hg19" | species == "hsa"){
        go.gs <- go.gsets(species="human", pkg.name=NULL, id.type = "eg")
      }else if(species == "Mouse" | species == "mouse" | species == "mm9" | species == "mmu"){
        go.gs <- go.gsets(species="mouse", pkg.name=NULL, id.type = "eg")
      }
    })
    # go.res <- go.gs$go.sets[ go.gs$go.subs[[input$goTerm]] ]
    return(go.gs)
  })
  gage.exp.fc <- reactive({
    rankdata <- nmf_group_feature_rank()
    validate(
      need(input$pathway_group!="", "Please select the enrichment result for the group you are interest in")
    )
    # In logFC, error might come from the fact that some of them are Inf or -Inf
    idx <- is.finite(rankdata[, input$pathway_group])
    rankdata <- rankdata[idx, input$pathway_group, drop=F]
    # ID conversion
    id.map.sym2eg <- id2eg(ids=rownames(rankdata), category = "SYMBOL", org=id.org())
    gene.entrez <- mol.sum(mol.data = rankdata, id.map = id.map.sym2eg, sum.method = "mean")
    deseq2.fc <- gene.entrez[, 1]
    exp.fc=deseq2.fc
    return(exp.fc)
  })
  go.Res <- reactive({
    gage.exp.fc <- gage.exp.fc()
    gsets <- NULL
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
      goRes <- gage(gage.exp.fc, gsets = gsets, ref = NULL, samp = NULL, same.dir=ifelse(input$samedir, TRUE, FALSE))
    })
    return(goRes)
  })

  go_table <- reactive({
    goRes <- go.Res()
    if (is.null(goRes$greater)){
      stop("No significant gene sets found")
    }else{
      return(goRes$greater)
    }
  })

  output$go_summary <- DT::renderDataTable({
    datatable(go_table(), selection = 'single',
              options = list(pageLength=6,
                             autoWidth = TRUE,
                             columnDefs = list(list(width = '50px', targets = "_all"))
              )
    ) %>% formatRound(1:ncol(go_table()), 3)
  }, server = TRUE)


  output$pathviewImg <- renderImage({
    go_table <- go_table()
    s.idx <- input$go_summary_rows_selected[length(input$go_summary_rows_selected)]
    validate(
      need(length(s.idx)>0, "Please select one GO-term/KEGG pathway from above.")
    )
    path.pid <- substr(rownames(go_table)[s.idx], 1, 8)
    # path.pid <- substr(s.idx, 1, 8)
    out.suffix <- "nmf"
    pathview(gene.data=gage.exp.fc(), pathway.id=path.pid, species=pathview.species(), out.suffix=out.suffix)

    filename <- path.pid
    outfile <- paste(filename, out.suffix, 'png', sep=".")
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 1000,
         height = 777,
         alt = "This is alternate text")
  }, deleteFile = FALSE)

  ### Download NMF ###

  output$downloadNMFData <- downloadHandler(
    filename <- function() {
      paste( file_prefix(), "nmf_results.zip", sep=".")
    },
    content = function(file) {
      tmpdir <- tempdir()
      current_dir <- getwd()
      setwd(tmpdir)
      # print(tempdir())
      filenames <- sapply(c("nmf_Groups", "nmf_Features", "Original_plus_NMF", "runSummary"),
                          function(x) paste(file_prefix(), x, "csv", sep="."))
      write.csv(nmf_groups(), filenames[1], quote=F, row.names=F)
      write.csv(nmf_features(), filenames[2], quote=F, row.names=F)
      write.csv(ori_plus_nmfResult(), filenames[3], quote=F)
      write.csv(runSummary(), filenames[4], quote=F)
      # print (filenames)
      zip(zipfile=file, files=filenames)
      setwd(as.character(current_dir))
    },
    contentType = "application/zip"
  )

  output$downloadPathData <- downloadHandler(
    filename <- function() {
      paste( file_prefix(), "pathway.zip", sep=".")
    },
    content = function(file) {
      tmpdir <- tempdir()
      current_dir <- getwd()
      setwd(tempdir())
      ### modify it to save one file for each rank ###
#       filenames <- sapply(c("group_feature_rank"),
#                           function(x) paste(file_prefix(), x, "csv", sep="."))
#       write.csv(nmf_group_feature_rank(), filenames[1], quote=F, row.names=T)

      ### modify it to save one file for each rank ###
      num_samples <- ncol(nmf_group_feature_rank())
      filenames <- paste(file_prefix(),
                         paste0("group_feature_rank", 1:num_samples),
                         "rnk",
                         sep=".")
      for (i in 1:num_samples){
        sub.data <- data.frame(Gene=rownames(nmf_group_feature_rank()), LogFC=nmf_group_feature_rank()[,i])
        idx <- is.finite(sub.data$LogFC)
        write.table(sub.data[idx, ], filenames[i], quote=F, row.names = F, sep="\t")
      }
      zip(zipfile=file, files=filenames)
      setwd(as.character(current_dir))
    },
    contentType = "application/zip"
  )

})
