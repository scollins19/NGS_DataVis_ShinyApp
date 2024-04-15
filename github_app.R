#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
require(stringi)
require(shiny)
require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
require("BSgenome.Hsapiens.UCSC.hg19")
require("plyranges")
require("tidyverse")
require("RColorBrewer")
require("ggpubr")
require(rstatix)
require(gridExtra)
require(cowplot)
require("shinyjs")
require(shinyWidgets)
require(EnsDb.Hsapiens.v75)
require(stringr)
require(shinythemes)
require(shinyFiles)

##################
source("deeptools_functions.r") # Path to file with functions
#################

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  theme = shinytheme("cosmo"),
  
  useShinyjs(),

    # Application title
    titlePanel("Make Heatmap / Average profile / boxplots"),
   # actionButton("reset","Clear"),
    
  actionButton("clear", "Restart"),
  
  tabsetPanel(
    tabPanel("DSB",fluid=TRUE,
      
  
  tabsetPanel(
    tabPanel("Matrix", fluid = TRUE,
             sidebarLayout(
               sidebarPanel( 
                 h3("Input Data"),
                 
                 h4("Select bigwigs"),
                 shinyFilesButton("paths", "Choose a file" ,
                                  title = "Please select a file:", multiple = T,
                                  buttonType = "default", class = NULL),
                 textOutput("txt_file"),   
                 h4(" \n "),
                 textInput('names', 'Enter labels for bigwig(s) (comma separated)', "name1,name2"),
                 h4("Select bed file"),
                 shinyFilesButton("beds", "Choose a file" ,
                                  title = "Please select a file:", multiple = T,
                                  buttonType = "default", class = NULL),
                 textOutput("bed_file"),   
                 h4(" \n "),
                             textInput('bnames', 'Enter labels for bed file(s) (comma separated)', "80_DSB"),
                             textInput('basepairs', 'Window size (in bps)', "10000"),
                             textInput('binsize', 'bin size (divide basepairs by 200 or 100)', "50"),
                 textInput('fileprefix', 'Filename (prefix only)', "ENDseq_80DSB"),
                             actionButton("matrix", "Compute Matrix"),
                             width=8),
               
               mainPanel(
                 verbatimTextOutput("matrix_status")
               )
             
              )
    ),
    
    tabPanel("Heatmap", fluid = F,
             sidebarLayout(
               sidebarPanel(h3("Plot options"),
                            textInput('plottitlehm', 'Plot Title', "ENDseq at 80 DSB and control"),
                            textInput('hm_order', 'Condition to order heatmap rows (one of the bigwig labels)', "name1"),
                            textInput('v_order', 'Bigwig order', "name1,name2"),
                            textInput('b_order', 'Bed file order', "bed1,bed2"),
                            selectInput('brewer_palettehm', 'Colour palette (look up R Brewer Palettes)', choices = c("Blues", "BrBG", "BuGn", 
                                                                                                                      "BuPu", "GnBu", "Greens", 
                                                                                                                      "Greys", "Oranges", "OrRd", 
                                                                                                                      "PiYG", "PRGn", "PuBu", "PuBuGn", 
                                                                                                                      "PuOr", "PuRd", "Purples", "RdBu", 
                                                                                                                      "RdGy", "RdPu", "RdYlBu", "RdYlGn",
                                                                                                                      "Reds", "Spectral", "YlGn", 
                                                                                                                      "YlGnBu", "YlOrBr", "YlOrRd" )),
                            selectInput('pal_dir','Switch Colour palette direction', choices=c("fwd","rev")),
                            textInput('hm_max', 'Scale max value (Default or a number)', "Default"),
                            textInput('hm_min', 'Scale min value (Default or a number)', "Default"),
                            actionButton("hm", "Make Heatmap"),
                            h3("Download Options"),
                            textInput('pwidthhm', 'Plot Width', "5"),
                            textInput('pheighthm', 'Plot Height', "8"),
                            width=3),
               
               
               mainPanel(
                 plotOutput("Heatmap", height = 700,width=500),
                 downloadButton("downloadhm", "Download Heatmap")
               
               )
             )
    ),
    
    
    tabPanel("Profile", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(h3("Plot options"),
                            textInput('plottitleap', 'Plot Title', "ENDseq at 80 DSB and control"),
                            textInput('v_orderap', 'Bigwig order', "name1,name2"),
                            textInput('b_orderap', 'Bed file order', "bed1,bed2"),
                            selectInput('brewer_paletteap', 'Colour palette (look up R Brewer Palettes)', choices = c("Blues", "BrBG", "BuGn", 
                                                                                                                      "BuPu", "GnBu", "Greens", 
                                                                                                                      "Greys", "Oranges", "OrRd", 
                                                                                                                      "PiYG", "PRGn", "PuBu", "PuBuGn", 
                                                                                                                      "PuOr", "PuRd", "Purples", "RdBu", 
                                                                                                                      "RdGy", "RdPu", "RdYlBu", "RdYlGn",
                                                                                                                      "Reds", "Spectral", "YlGn", 
                                                                                                                      "YlGnBu", "YlOrBr", "YlOrRd"  )),
                            textInput('yaxis', 'Y-axis label', "END-seq average norm read count"),
                            textInput('ap_max', 'Y-axis max value (Default or a number)', "Default"),
                            textInput('ap_min', 'Y-axis min value (Default or a number)', "Default"),
                            selectInput('facetap', 'Select Layout', choices = c("Peaks","variable")),
                            actionButton("ap", "Make Average Profile"),
                            h3("Download Options"),
                            textInput('pwidthap', 'Plot Width', "7"),
                            textInput('pheightap', 'Plot Height', "6"),
                            width=3),
               
               
               mainPanel(
                 plotOutput("AvgProf", height = 400,width=500),
                 downloadButton("downloadap", "Download Average Profile")
                 
               )
             )
    ),
    
    
    tabPanel("Boxplot", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(h3("Plot options"),
                            textInput('plottitlebp', 'Plot Title', "ENDseq at 80 DSB and control"),
                            textInput('v_orderbp', 'Bigwig order', "name1,name2"),
                            textInput('b_orderbp', 'Bed file order', "bed1,bed2"),
                            selectInput('brewer_palettebp', 'Colour palette (look up R Brewer Palettes)', choices = c("Blues", "BrBG", "BuGn", 
                                                                                                                      "BuPu", "GnBu", "Greens", 
                                                                                                                      "Greys", "Oranges", "OrRd", 
                                                                                                                      "PiYG", "PRGn", "PuBu", "PuBuGn", 
                                                                                                                      "PuOr", "PuRd", "Purples", "RdBu", 
                                                                                                                      "RdGy", "RdPu", "RdYlBu", "RdYlGn",
                                                                                                                      "Reds", "Spectral", "YlGn", 
                                                                                                                      "YlGnBu", "YlOrBr", "YlOrRd"  )),
                            textInput('yaxisbp', 'Y-axis label', "END-seq average norm read count"),
                            textInput('bp_max', 'Y-axis max value (Default or a number)', "Default"),
                            textInput('bp_min', 'Y-axis min value (Default or a number)', "Default"),
                            selectInput('facetbp', 'Switch Layout (variable or Peaks)', choices = c("Peaks","variable")),
                            actionButton("bp", "Make Boxplot"),
                            h3("Download Options"),
                            textInput('pwidthbp', 'Plot Width', "5"),
                            textInput('pheightbp', 'Plot Height', "8"),
                            width=4),
               
               
               mainPanel(
                 plotOutput("Boxplot", height = 700,width=600),
                 downloadButton("downloadbp", "Download Boxplot"),
                 downloadButton("downloadbpstats", "Download Boxplot with p-val")
                 
               )
             )
    )
    
    
    
    
  )
  
    ),
  

    tabPanel("Gene Bodies",fluid=TRUE,
             
             
             tabsetPanel(
               tabPanel("Matrix", fluid = TRUE,
                        sidebarLayout(
                          sidebarPanel( 
                            h3("Input Data"),
                            
                            h4("Select bigwigs"),
                            shinyFilesButton("pathsgb", "Choose a file" ,
                                             title = "Please select a file:", multiple = T,
                                             buttonType = "default", class = NULL),
                            textOutput("txt_filegb"),   
                            h4(" \n "),
                            textInput('namesgb', 'Enter labels for bigwig(s) (comma separated)', "name1,name2"),
                            h4("Select bed file - Gene names need to be in 4th column"),
                            shinyFilesButton("bedsgb", "Choose a file" ,
                                             title = "Please select a file:", multiple = T,
                                             buttonType = "default", class = NULL),
                            textOutput("bed_filegb"),   
                            h4(" \n "),
                            textInput('bnamesgb', 'Enter labels for bed file(s) (comma separated)', "bed1,bed2"),
                            textInput('basepairsgb', 'Window -/+ TSS/TES (in bps)', "2000"),
                            textInput('binsizegb', 'bin size (divide basepairs by 200 or 100)', "50"),
                            textInput('fileprefixgb', 'Filename (prefix only)', "ENDseq_damaged_genes"),
                            actionButton("matrixgb", "Compute Matrix"),
                            width=10),
                          
                          mainPanel(
                            verbatimTextOutput("matrix_statusgb")
                          )
                          
                        )
               ),
               
               tabPanel("Heatmap", fluid = F,
                        sidebarLayout(
                          sidebarPanel(h3("Plot options"),
                                       textInput('plottitlehmgb', 'Plot Title', "ENDseq at 80 DSB and control"),
                                       textInput('hm_ordergb', 'Condition to order heatmap rows (one of the bigwig labels)', "name1"),
                                       
                                       textInput('v_ordergb', 'Bigwig order', "name1,name2"),
                                       textInput('b_ordergb', 'Bed file order', "bed1,bed2"),
                                       
                                       selectInput('brewer_palettehmgb', 'Colour palette (look up R Brewer Palettes)', choices = c("Blues", "BrBG", "BuGn", 
                                                                                                                                   "BuPu", "GnBu", "Greens", 
                                                                                                                                   "Greys", "Oranges", "OrRd", 
                                                                                                                                   "PiYG", "PRGn", "PuBu", "PuBuGn", 
                                                                                                                                   "PuOr", "PuRd", "Purples", "RdBu", 
                                                                                                                                   "RdGy", "RdPu", "RdYlBu", "RdYlGn",
                                                                                                                                   "Reds", "Spectral", "YlGn", 
                                                                                                                                   "YlGnBu", "YlOrBr", "YlOrRd" )),
                                       selectInput('pal_dirgb','Switch Colour palette direction', choices=c("fwd","rev")),
                                       textInput('hm_maxgb', 'Scale max value (Default or a number)', "Default"),
                                       textInput('hm_mingb', 'Scale min value (Default or a number)', "Default"),
                                       actionButton("hmgb", "Make Heatmap"),
                                       h3("Download Options"),
                                       textInput('pwidthhmgb', 'Plot Width', "5"),
                                       textInput('pheighthmgb', 'Plot Height', "8"),
                                       width=4),
                          
                          
                          mainPanel(
                            plotOutput("Heatmapgb", height = 700,width=500),
                            downloadButton("downloadhmgb", "Download Heatmap")
                            
                          )
                        )
               ),
               
               
               tabPanel("Profile", fluid = TRUE,
                        sidebarLayout(
                          sidebarPanel(h3("Plot options"),
                                       textInput('plottitleapgb', 'Plot Title', "ENDseq at 80 DSB and control"),
                                       
                                       textInput('v_orderapgb', 'Bigwig order', "name1,name2"),
                                       textInput('b_orderapgb', 'Bed file order', "bed1,bed2"),
                                       
                                       selectInput('brewer_paletteapgb','Colour palette (look up R Brewer Palettes)', choices = c("Blues", "BrBG", "BuGn", 
                                                                                                                                  "BuPu", "GnBu", "Greens", 
                                                                                                                                  "Greys", "Oranges", "OrRd", 
                                                                                                                                  "PiYG", "PRGn", "PuBu", "PuBuGn", 
                                                                                                                                  "PuOr", "PuRd", "Purples", "RdBu", 
                                                                                                                                  "RdGy", "RdPu", "RdYlBu", "RdYlGn",
                                                                                                                                  "Reds", "Spectral", "YlGn", 
                                                                                                                                  "YlGnBu", "YlOrBr", "YlOrRd" )),
                                       textInput('yaxisgb', 'Y-axis label', "END-seq average norm read count"),
                                       textInput('ap_maxgb', 'Y-axis max value (Default or a number)', "Default"),
                                       textInput('ap_mingb', 'Y-axis min value (Default or a number)', "Default"),
                                       selectInput('facetapgb', 'Switch Layout', choices=c("Peaks","variable")),
                                       actionButton("apgb", "Make Average Profile"),
                                       h3("Download Options"),
                                       textInput('pwidthapgb', 'Plot Width', "7"),
                                       textInput('pheightapgb', 'Plot Height', "6"),
                                       width=4),
                          
                          
                          mainPanel(
                            plotOutput("AvgProfgb", height = 400,width=500),
                            downloadButton("downloadapgb", "Download Average Profile")
                            
                          )
                        )
               ),
               
               
               tabPanel("Boxplot", fluid = TRUE,
                        sidebarLayout(
                          sidebarPanel(h3("Plot options"),
                                       textInput('plottitlebpgb', 'Plot Title', "ENDseq at 80 DSB and control"),
                                       
                                       textInput('v_orderbpgb', 'Bigwig order', "name1,name2"),
                                       textInput('b_orderbpgb', 'Bed file order', "bed1,bed2"),
                                       
                                       selectInput('brewer_palettebpgb', 'Colour palette (look up R Brewer Palettes)', choices = c("Blues", "BrBG", "BuGn", 
                                                                                                                                   "BuPu", "GnBu", "Greens", 
                                                                                                                                   "Greys", "Oranges", "OrRd", 
                                                                                                                                   "PiYG", "PRGn", "PuBu", "PuBuGn", 
                                                                                                                                   "PuOr", "PuRd", "Purples", "RdBu", 
                                                                                                                                   "RdGy", "RdPu", "RdYlBu", "RdYlGn",
                                                                                                                                   "Reds", "Spectral", "YlGn", 
                                                                                                                                   "YlGnBu", "YlOrBr", "YlOrRd"  )),
                                       textInput('yaxisbpgb', 'Y-axis label', "END-seq average norm read count"),
                                       textInput('bp_maxgb', 'Y-axis max value (Default or a number)', "Default"),
                                       textInput('bp_mingb', 'Y-axis min value (Default or a number)', "Default"),
                                       selectInput('facetbpgb', 'Switch Layout', choices=c("Peaks","variable") ) ,
                                       actionButton("bpgb", "Make Boxplot"),
                                       h3("Download Options"),
                                       textInput('pwidthbpgb', 'Plot Width', "5"),
                                       textInput('pheightbpgb', 'Plot Height', "8"),
                                       width=4),
                          
                          
                          mainPanel(
                            plotOutput("Boxplotgb", height = 700,width=600),
                            downloadButton("downloadbpgb", "Download Boxplot"),
                            downloadButton("downloadbpstatsgb", "Download Boxplot with p-val")
                            
                          )
                        )
               )
               
               
               
               
             )
             
    )
  
  
  )
  
  
  

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  
  rbwfiles <- reactiveValues(data=0)
  rbedfiles <- reactiveValues(data=0)
  
  rResult <- reactiveValues(data=0)
  
  #volumes <- c(LEGUBE1='/media/LEGUBE1',LEGUBE2='/media/LEGUBE2/')
  
  observe({  
    shinyFileChoose(input, "paths", roots = volumes, session = session)
    
    
    if(!is.null(input$paths)){
      # browser()
      bigwig_files <-parseFilePaths(volumes, input$paths)
      output$txt_file <- renderText(as.character(bigwig_files$datapath))
      rbwfiles$data <- as.character(bigwig_files$datapath)
    }
  })
  
  observe({  
    shinyFileChoose(input, "beds", roots =  volumes, session = session)
    
    
    if(!is.null(input$beds)){
      # browser()
      bed_files<-parseFilePaths(volumes, input$beds)
      output$bed_file <- renderText(as.character(bed_files$datapath))
      rbedfiles$data <- as.character(bed_files$datapath)
    }
  })
  
  
  observeEvent(input$matrix, {
  
    showModal(modalDialog("Computing matrix...", footer=NULL))
      
      bws <- rbwfiles$data 
      bpath <- rbedfiles$data 
      
    message(bws)
    message(bpath)
    
    bw_names <- unlist(strsplit(input$names,","))
    bed_names <- unlist(strsplit(input$bnames,","))
    bps <- as.numeric(input$basepairs)
    bin <- as.numeric(input$binsize)
    
    message(bps)
    message(bin)
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    #deeptoolspath <- "/home/scollins/anaconda3/envs/deeptools/bin/computeMatrix"
    

    
    for(i in 1:length(bws)){
      
      for(j in 1:length(bpath)){
        
        message(paste0(deeptoolspath," reference-point -S ", bws[i], 
                      " -R ", bpath[j],
                      " --referencePoint center -a ", half_bps, " -b ", half_bps,
                      " --samplesLabel ", bw_names[i],
                      " --outFileName ", bw_names[i],bed_names[j],sizekb,".ma.gz",
                      " --binSize ", bin,
                      " --sortRegions keep -p 25")
        )
        
        
        
        system(paste0(deeptoolspath," reference-point -S ", bws[i], 
                      " -R ", bpath[j],
                      " --referencePoint center -a ", half_bps, " -b ", half_bps,
                      " --samplesLabel ", bw_names[i],
                      " --outFileName ", bw_names[i],bed_names[j],sizekb,".ma.gz",
                      " --binSize ", bin,
                      " --sortRegions keep -p 25")
        )
        
        
      }
      
    }
    
      
      df <- plot_group(
        FILES = c( do.call(paste0, expand.grid(bw_names,bed_names,sizekb,".ma.gz"))
        ), 
        vars = c( rep(bw_names,length(bed_names))
        ),
        regions=c( rep(bed_names,each=length(bw_names))
        ),
        strand = "same"
      )
      
      rResult$data <- df
      
      removeModal()
      output$matrix_status <- renderText({"Matrix Done!"})
      
      system(paste0("rm ",working_dir,"/*.ma.gz"))

  })
  

  
  observeEvent(input$ap, {
    
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data
    
    bw_names <- unlist(strsplit(input$names,","))
    bed_names <- unlist(strsplit(input$bnames,","))
    
    bps <- as.numeric(input$basepairs)
    bin <- as.numeric(input$binsize)
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    ## save data to individual variables for plotting
    dat.plot <- rResult$data[[2]]
    
    y_axis_title=  as.character(input$yaxis) # y axis label for average profile
    
    var_order = input$v_orderap  %>% as.character %>% strsplit(.,",") %>% unlist() # order to show samples (needs to be exact sample as "vars" input in function)
    
    regions_order = input$b_orderap %>% as.character %>% strsplit(.,",") %>% unlist() # Order to show regions if more than one used. If not, just name the one
    
    hm_order_var = as.character(input$hm_order) # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
    
    brewer_palette=  as.character(input$brewer_paletteap) # brewer color palette used for heatmap
    
    
    
    dat.plot$Peaks <- factor(dat.plot$Peaks,  levels=c(regions_order))
    dat.plot$variable <- factor(dat.plot$variable, levels=c(var_order))
    
    if(input$ap_max == "Default"){
      maxval <- max(dat.plot$Value)
    }else{
      maxval <- as.numeric(input$ap_max)
    }
    
    if(input$ap_min == "Default"){
      minval <- min(dat.plot$Value)
    }else{
      minval <- as.numeric(input$ap_min)
    }
    
    
    
    if(input$facetap == "Peaks"){
      
      colorvar <- "variable"
      facetval <- "Peaks"
      nameval <- "Treatment"
      
      if(length(bw_names)>3){
        mycolor <- rev(brewer.pal(length(bw_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
      
    }else if(input$facetap == "variable"){
      colorvar <- "Peaks"
      facetval <- "variable"
      nameval <- "Regions"
      
      if(length(bed_names)>3){
        mycolor <- rev(brewer.pal(length(bed_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
    }
    
    
  
   # mycolor <- (brewer.pal(9,brewer_palette)[c(5,7)])
     
    avgp <- reactive({ggplot( na.omit( dat.plot), aes( Window, Value, colour=get(colorvar)) ) +
      labs( list( title = "", x = "", y = " " )) +
      geom_line(size=0.7, alpha=0.7) +
      theme_classic(base_size = 8) +
      theme(axis.text.y = element_text(size=7),
            axis.text.x = element_text(size=7),
            axis.title.x = element_text(size=7),
            plot.title = element_text(size=9, face="bold", hjust=0.5),
            legend.position = "bottom") + ylab(y_axis_title) +
      
      scale_x_continuous(name = paste0("distance from center"),
                         breaks = c(1, length(levels(as.factor(dat.plot$Window))) / 2 ,length(levels(as.factor(dat.plot$Window)))),
                         labels = c(paste0("-", paste0( (half_bps/1000),"kb" ) ), 'center', paste0( (half_bps/1000),"kb" )  ) 
      ) +
      geom_vline(xintercept = length(levels(as.factor(dat.plot$Window))) / 2, linetype="dashed", color="grey") +
        scale_color_manual(values=mycolor, name=nameval) +
        ggtitle(paste0(as.character(input$plottitleap)," -/+ ",(half_bps/1000),"kb")) + facet_grid(~get(facetval)) + ylim(c(minval,maxval))
      
    })
    
  
      
      output$AvgProf <- renderPlot({
        avgp()
      })
      

        output$downloadap <- downloadHandler(
          filename = function() {
            paste0(input$fileprefix,"_",sizekb, "_AvgProf.pdf")
          },
          content = function(file) {
            pdf(file=file, width=as.numeric(input$pwidthap), height=as.numeric(input$pheightap))
            print(avgp())
            dev.off()
          }
        )
    
    })
  
  
  
  observeEvent(input$hm, {
    
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data
    
    bw_names <- unlist(strsplit(input$names,","))
    bed_names <- unlist(strsplit(input$bnames,","))
    
    bps <- as.numeric(input$basepairs)
    bin <- as.numeric(input$binsize)
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    ## save data to individual variables for plotting
    dh<- rResult$data[[1]]
    
    y_axis_title=  as.character(input$yaxis) # y axis label for average profile
    
    var_order = input$v_order  %>% as.character %>% strsplit(.,",") %>% unlist() # order to show samples (needs to be exact sample as "vars" input in function)
    
    regions_order = input$b_order %>% as.character %>% strsplit(.,",") %>% unlist() # Order to show regions if more than one used. If not, just name the one
    
    hm_order_var = as.character(input$hm_order) # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
    
    brewer_palette=  as.character(input$brewer_palettehm) # brewer color palette used for heatmap
    
    
    ## MAKE HEATMAP
    
    dh$variable <- factor(dh$variable, levels=c(var_order)) # set variable factors
    dh$Peaks <- factor(dh$Peaks, levels=c(regions_order)) # set regions factors
    
    dat.heatmap <- dh
    
    hm_cutoff = quantile(dat.heatmap$value, 0.95)
    hm_cutoff_m  = quantile(dat.heatmap$value, 0.01)
    lim <- pmax(abs(hm_cutoff), abs(hm_cutoff_m))
    
    
    if(input$hm_max == "Default"){
      dat.heatmap$value[dat.heatmap$value>hm_cutoff] <- hm_cutoff
      maxval <- lim
    }else{
      dat.heatmap$value[dat.heatmap$value>as.numeric(input$hm_max)] <- as.numeric(input$hm_max)
      maxval <- as.numeric(input$hm_max)
    }
    if(input$hm_min == "Default"){
      dat.heatmap$value[dat.heatmap$value<hm_cutoff_m] <- hm_cutoff_m
      minval <- 0
    }else{
      dat.heatmap$value[dat.heatmap$value<as.numeric(input$hm_min)] <- as.numeric(input$hm_min)
      minval <- as.numeric(input$hm_min)
    }
    
    
    if(input$pal_dir=="fwd"){
      d <- 1
    }else{
      d<- -1
    }
    
    
    ord <- subset(dat.heatmap, variable==hm_order_var) %>%
      dplyr::group_by(site) %>%
      dplyr::summarise(value=mean(value))
    ord <- ord$site[order(-ord$value)]
    
    dat.heatmap$site <- factor(dat.heatmap$site, levels=rev(ord))
    
    if(length(levels(as.factor(dat.heatmap$site)))>100){
      ytxt <- theme(axis.text.y =element_blank() )
    }else{
      ytxt <-  theme(axis.text.y =element_text(size=5) )
    }
    
    ## plot heatmap
    gghm <- reactive({ggplot(dat.heatmap, aes(window, site)) +
        geom_raster(aes(fill = value)) +
        facet_grid(Peaks~variable, scales = "free_y", space = "free_y") +
        ylab(" ") +
        scale_x_discrete(name = paste0("Distance from center (bp)"),
                         breaks = c(1, length(levels(as.factor(dat.heatmap$window))) / 2 ,length(levels(as.factor(dat.heatmap$window)))),
                         labels = c(paste0("-", paste0( (half_bps/1000),"kb" ) ), 'center', paste0( (half_bps/1000),"kb" )  )   ) +
        #geom_vline(xintercept = length(levels(as.factor(dat.heatmap$window))) / 2, linetype="dashed", color="black", size=0.3) +
        theme_classic(base_size = 7) +
        theme(
              axis.text.x = element_text(size=6, face="bold"),
              axis.title.x = element_text(size=7),
              plot.title = element_text(size=7, face="bold", hjust=0.5),
              legend.position = "right") + ytxt +
        #scale_fill_gradient2(low="blue", mid="white",high="red", midpoint = 0, limits=c(-lim,lim),na.value = "black") +
        #scale_color_gradient2(low="blue", mid="white",high="red", midpoint = 0, limits=c(-lim,lim),na.value = "black") +
        scale_fill_distiller(palette = brewer_palette, direction = d, na.value = "white", limits=c(minval,maxval)) +
        scale_color_distiller(palette = brewer_palette, direction = d, na.value = "white", limits=c(minval,maxval)) +
        ggtitle(paste0(as.character(input$plottitlehm)," -/+ ",(half_bps/1000),"kb")) + facet_grid(Peaks~variable, scales="free_y")
    })
    
    output$Heatmap <- renderPlot({
      gghm()
    })
    
    output$downloadhm <- downloadHandler(
      filename = function() {
        paste0(input$fileprefix,"_",sizekb, "_Heatmap.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthhm), height=as.numeric(input$pheighthm))
        print(gghm())
        dev.off()
      }
    )
    
  })
  
  
  
  
  observeEvent(input$bp, {
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data
    
    bw_names <- unlist(strsplit(input$names,","))
    bed_names <- unlist(strsplit(input$bnames,","))
    
    bps <- as.numeric(input$basepairs)
    bin <- as.numeric(input$binsize)
    
    half_bps <- bps/2
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    ## save data to individual variables for plotting
    dat.boxplot <- rResult$data[[3]]
    
    y_axis_title=  as.character(input$yaxisbp) # y axis label for average profile
    
    var_order = input$v_orderbp  %>% as.character %>% strsplit(.,",") %>% unlist() # order to show samples (needs to be exact sample as "vars" input in function)
    
    regions_order = input$b_orderbp %>% as.character %>% strsplit(.,",") %>% unlist() # Order to show regions if more than one used. If not, just name the one
    
    hm_order_var = as.character(input$hm_order) # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
    
    brewer_palette=  as.character(input$brewer_palettebp) # brewer color palette used for heatmap
    
    
    dat.boxplot$Peaks <- factor(dat.boxplot$Peaks, levels=c(regions_order)) # set regions factors
    dat.boxplot$variable <- factor(dat.boxplot$variable, levels=c(var_order)) # set regions factors
    
    
    if(input$bp_max == "Default"){
      maxval <- max(dat.boxplot$Value)
    }else{
      maxval <- as.numeric(input$bp_max)
    }
    
    if(input$bp_min == "Default"){
      minval <- min(dat.boxplot$Value)
    }else{
      minval <- as.numeric(input$bp_min)
    }
    
    
    
    
    
    ss<- NULL
    ss1<- NULL
    
    
    if(length(bw_names)>1){
      stat.test <- as.data.frame(dat.boxplot) %>%
        group_by(Peaks) %>%
        wilcox_test(Value ~ variable, paired=T) %>%
        add_significance()
      
      ss <- tableGrob(stat.test) 
    }else if(length(bw_names)==1){
      ss <- NULL
    }
    
    if(length(bed_names)>1){
      stat.test <- as.data.frame(dat.boxplot) %>%
        group_by(variable) %>%
        wilcox_test(Value ~ Peaks, paired=F) %>%
        add_significance()
      
      ss1 <- tableGrob(stat.test) 
    }else if(length(bed_names)==1){
      ss1 <- NULL
    }
    
    
    
    if(input$facetbp == "Peaks"){
      
      colorvarbp <- "variable"
      facetvalbp <- "Peaks"
      namevalbp <- "Treatment"
      
      if(length(bw_names)>3){
        mycolor <- rev(brewer.pal(length(bw_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
  
      
      
    }else if(input$facetbp == "variable"){
      colorvarbp <- "Peaks"
      facetvalbp <- "variable"
      namevalbp <- "Regions"
      
      if(length(bed_names)>3){
        mycolor <- rev(brewer.pal(length(bed_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
    }
  
    
    
    ggbp <- reactive({ggplot(na.omit(dat.boxplot), aes(x = get(colorvarbp), y = (Value), fill=get(colorvarbp))) +
      geom_boxplot() +
      theme_classic() +
      theme(#axis.text.x  = element_text(size=12, angle=90, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y  = element_text(size=5),
        axis.title = element_text(angle = 0, hjust = 0.5, size=15),
        legend.position = "right",
        legend.text = element_text(size=7),
        legend.title = element_text(size=9),
        plot.title = element_text(hjust = 0.5, face="bold", size=7)) +
      scale_fill_manual(values = mycolor, name=namevalbp) +
      xlab(namevalbp) + ggtitle(paste0(as.character(input$plottitlebp)," -/+ ",(half_bps/1000),"kb")) + ylab(y_axis_title) +
        facet_grid(~get(facetvalbp)) + ylim(c(minval,maxval))
    
  })
    
    
    if(length(bw_names)>1 & length(bed_names)>1){
      bpnew <- reactive({plot_grid(ggbp(),ss,ss1,nrow=3,rel_heights = c(1,0.3,0.3))})
    }else if(length(bw_names)>1 & length(bed_names)==1){
      bpnew <-  reactive({plot_grid(ggbp(),ss,nrow=2,rel_heights = c(1,0.3))})
    }else if(length(bw_names)==1 & length(bed_names)>1){
      bpnew <- reactive({plot_grid(ggbp(),ss1,nrow=2,rel_heights = c(1,0.3))})
    }
    
    
    output$Boxplot <- renderPlot({
     bpnew()
    })
    
    
    output$downloadbp <- downloadHandler(
      filename = function() {
        paste0(input$fileprefix,"_",sizekb, "_Boxplot.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthbp), height=as.numeric(input$pheightbp))
        print(ggbp())
        dev.off()
      }
    )
    
    
    output$downloadbpstats <- downloadHandler(
      filename = function() {
        paste0(input$fileprefix,"_",sizekb, "_withstats_Boxplot.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthbp), height=as.numeric(input$pheightbp))
        print(bpnew())
        dev.off()
      }
    )
    
  })
  
  
  
  observe({  
    shinyFileChoose(input, "pathsgb", roots = volumes, session = session)
    
    
    if(!is.null(input$pathsgb)){
      # browser()
      bigwig_files <-parseFilePaths(volumes, input$pathsgb)
      output$txt_filegb <- renderText(as.character(bigwig_files$datapath))
      rbwfiles$data <- as.character(bigwig_files$datapath)
    }
  })
  
  observe({  
    shinyFileChoose(input, "bedsgb", roots =  volumes, session = session)
    
    
    if(!is.null(input$bedsgb)){
      # browser()
      bed_files<-parseFilePaths(volumes, input$bedsgb)
      output$bed_filegb <- renderText(as.character(bed_files$datapath))
      rbedfiles$data <- as.character(bed_files$datapath)
    }
  })
  
  
  observeEvent(input$matrixgb, {
    
    showModal(modalDialog("Computing matrix...", footer=NULL))
    
    
    
    
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data 
    
    message(bws)
    message(bpath)
    
    bw_names <- unlist(strsplit(input$namesgb,","))
    bed_names <- unlist(strsplit(input$bnamesgb,","))
    bps <- as.numeric(input$basepairsgb)
    bin <- as.numeric(input$binsizegb)
    
    message(bps)
    message(bin)
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    #deeptoolspath <- "/home/scollins/anaconda3/envs/deeptools/bin/computeMatrix"
    

    bpathnew <- c()
    bed_namesnew <- c()
    
    Chr.V <- paste0("chr",c(1:22,"X"))
    for(j in 1:length(bpath)){
      
      allgenes <- genes(EnsDb.Hsapiens.v75)
      seqlevels(allgenes) <- paste0("chr",seqlevels(allgenes))
      allgenes <- allgenes[allgenes$gene_biotype=="protein_coding" & seqnames(allgenes) %in% Chr.V,]
      
      bedfile <- bpath[j]
      bedfile <- read_table(bedfile,col_names = F)
      
      genesbed <-as.vector(unlist( bedfile[,4]))
      newbed <- allgenes[allgenes$gene_name %in% genesbed,]
      write.table(as.data.frame(newbed)[,c(1:3,7,4,5)],paste0(bed_names[j],"_new.bed"),sep="\t",col.names = F,row.names = F,quote=F)
      
      bpathnew <- c(bpathnew,paste0(bed_names[j],"_new.bed") )

    }
    
    
    bpath <- bpathnew
   # bed_names <- bed_namesnew
    
    
    for(i in 1:length(bws)){
      
      for(j in 1:length(bpath)){
        
        
        system(paste0(deeptoolspath," scale-regions -S ", bws[i], 
                      " -R ", bpath[j],
                      " -a ", bps, " -b ", bps,
                      " --samplesLabel ", bw_names[i],
                      " --regionBodyLength 5000 ",
                      " --outFileName ", bw_names[i],bed_names[j],sizekb,".ma.gz",
                      " --binSize ", bin,
                      " --sortRegions keep -p 25")
        )
        
        
      }
      
      
      
    }
    
    
    files <- do.call(paste0, expand.grid(bw_names,bed_names,sizekb,".ma.gz"))
   
      df <- plot_group(
        FILES = c( do.call(paste0, expand.grid(bw_names,bed_names,sizekb,".ma.gz"))
        ), 
        vars = c( rep(bw_names,length(bed_names))
        ),
        regions=c( rep(bed_names,each=length(bw_names))
        ),
        strand = "same"
      )
      
      rResult$data <- df
      
      removeModal()
      output$matrix_statusgb <- renderText({"Matrix Done!"})
      
      system(paste0("rm ",working_dir,"/*.ma.gz"))
      
    
  
  })
  
  
  
  observeEvent(input$apgb, {
    
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data
    
    bw_names <- unlist(strsplit(input$namesgb,","))
    bed_names <- unlist(strsplit(input$bnamesgb,","))
    
    bps <- as.numeric(input$basepairsgb)
    bin <- as.numeric(input$binsizegb)
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    ## save data to individual variables for plotting
    dat.plot <- rResult$data[[2]]
    
    y_axis_title=  as.character(input$yaxisgb) # y axis label for average profile
    
    var_order = input$v_orderapgb  %>% as.character %>% strsplit(.,",") %>% unlist() # order to show samples (needs to be exact sample as "vars" input in function)
    
    regions_order = input$b_orderapgb %>% as.character %>% strsplit(.,",") %>% unlist() # Order to show regions if more than one used. If not, just name the one
    
    brewer_palette=  as.character(input$brewer_paletteapgb) # brewer color palette used for heatmap
    
  
    dat.plot$Peaks <- factor(dat.plot$Peaks,  levels=c(regions_order))
    dat.plot$variable <- factor(dat.plot$variable, levels=c(var_order))
    
    if(input$ap_maxgb == "Default"){
      maxval <- max(dat.plot$Value)
    }else{
      maxval <- as.numeric(input$ap_maxgb)
    }
    
    if(input$ap_mingb == "Default"){
      minval <- min(dat.plot$Value)
    }else{
      minval <- as.numeric(input$ap_mingb)
    }
    
    
    
    if(input$facetapgb == "Peaks"){
      
      colorvar <- "variable"
      facetval <- "Peaks"
      nameval <- "Treatment"
      
      if(length(bw_names)>3){
        mycolor <- rev(brewer.pal(length(bw_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
      
    }else if(input$facetapgb == "variable"){
      colorvar <- "Peaks"
      facetval <- "variable"
      nameval <- "Regions"
      
      if(length(bed_names)>3){
        mycolor <- rev(brewer.pal(length(bed_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
    }
    
    
    
    # mycolor <- (brewer.pal(9,brewer_palette)[c(5,7)])
    
    avgp <- reactive({ggplot( na.omit( dat.plot), aes( Window, Value, colour=get(colorvar)) ) +
        labs( list( title = "", x = "", y = " " )) +
        geom_line(size=0.7, alpha=0.7) +
        theme_classic(base_size = 8) +
        theme(axis.text.y = element_text(size=7),
              axis.text.x = element_text(size=7),
              axis.title.x = element_text(size=7),
              plot.title = element_text(size=9, face="bold", hjust=0.5),
              legend.position = "bottom") + ylab(y_axis_title) +
        
        scale_x_continuous(name = paste0("distance from center"),
                           breaks = c(1, as.numeric(bps/bin) ,  as.numeric((bps/bin) + (5000/bin)) , as.numeric((bps/bin) + (5000/bin) + (bps/bin)) ),
                           labels = c(paste0("-",sizekb),"TSS","TES",paste0("+",sizekb))
        ) +
        geom_vline(xintercept = as.numeric(bps/bin), linetype="dashed", color="grey") +
        geom_vline(xintercept =  as.numeric(  (bps/bin) + (5000/bin) ) , linetype="dashed", color="grey") +
        scale_color_manual(values=mycolor, name=nameval) +
        ggtitle(paste0(as.character(input$plottitleapgb)," -/+ ",sizekb)) + facet_grid(~get(facetval)) + ylim(c(minval,maxval))
      
    })
    
    
    
    output$AvgProfgb <- renderPlot({
      avgp()
    })
    
    
    output$downloadapgb <- downloadHandler(
      filename = function() {
        paste0(input$fileprefixgb,"_",sizekb, "_AvgProf.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthapgb), height=as.numeric(input$pheightapgb))
        print(avgp())
        dev.off()
      }
    )
    
  })
  
  
  
  observeEvent(input$hmgb, {
    
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data
    
    bw_names <- unlist(strsplit(input$namesgb,","))
    bed_names <- unlist(strsplit(input$bnamesgb,","))
    
    bps <- as.numeric(input$basepairsgb)
    bin <-  as.numeric(input$binsizegb)
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    ## save data to individual variables for plotting
    dh<- rResult$data[[1]]
    
    y_axis_title=  as.character(input$yaxisgb) # y axis label for average profile
    
    var_order = input$v_ordergb  %>% as.character %>% strsplit(.,",") %>% unlist() # order to show samples (needs to be exact sample as "vars" input in function)
    
    regions_order = input$b_ordergb %>% as.character %>% strsplit(.,",") %>% unlist() # Order to show regions if more than one used. If not, just name the one
    
    
    hm_order_var = as.character(input$hm_ordergb) # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
    
    brewer_palette=  as.character(input$brewer_palettehmgb) # brewer color palette used for heatmap
    
    
    ## MAKE HEATMAP
    
    dh$variable <- factor(dh$variable, levels=c(var_order)) # set variable factors
    dh$Peaks <- factor(dh$Peaks, levels=c(regions_order)) # set regions factors
    
    
    
    dat.heatmap <- dh
    
    hm_cutoff = quantile(dat.heatmap$value, 0.95)
    hm_cutoff_m  = quantile(dat.heatmap$value, 0.01)
    lim <- pmax(abs(hm_cutoff), abs(hm_cutoff_m))
    
    
    if(input$hm_maxgb == "Default"){
      dat.heatmap$value[dat.heatmap$value>hm_cutoff] <- hm_cutoff
      maxval <- lim
    }else{
      dat.heatmap$value[dat.heatmap$value>as.numeric(input$hm_maxgb)] <- as.numeric(input$hm_maxgb)
      maxval <- as.numeric(input$hm_maxgb)
    }
    if(input$hm_mingb == "Default"){
      dat.heatmap$value[dat.heatmap$value<hm_cutoff_m] <- hm_cutoff_m
      minval <- 0
    }else{
      dat.heatmap$value[dat.heatmap$value<as.numeric(input$hm_mingb)] <- as.numeric(input$hm_mingb)
      minval <- as.numeric(input$hm_mingb)
    }
    
    
    
    
    
    if(input$pal_dirgb=="fwd"){
      d <- 1
    }else{
      d<- -1
    }
    
    
    
    ord <- subset(dat.heatmap, variable==hm_order_var) %>%
      dplyr::group_by(site) %>%
      dplyr::summarise(value=mean(value))
    ord <- ord$site[order(-ord$value)]
    
    dat.heatmap$site <- factor(dat.heatmap$site, levels=rev(ord))
    
    b2 <- as.numeric(bps/bin)
    b3 <- as.numeric((bps/bin) + (5000/bin))
    b4 <- as.numeric((bps/bin) + (5000/bin) + (bps/bin))
    
    ## plot heatmap
    gghmgb <- reactive({ggplot(dat.heatmap, aes(window, site)) +
        geom_raster(aes(fill = value)) +
        facet_grid(Peaks~variable, scales = "free_y", space = "free_y") +
        ylab(" ") +
        scale_x_discrete(name = c("TSS -> TES"),
                           breaks = c(1,b2,b3,b4 ),
                           labels = c(paste0("-",sizekb) ,"TSS","TES", paste0("+",sizekb) )
        ) +
        geom_vline(xintercept = as.numeric(bps/bin), linetype="dashed", color="grey") +
        geom_vline(xintercept =  as.numeric((bps/bin) + (5000/bin)) , linetype="dashed", color="grey") +
        theme_classic(base_size = 7) +
        theme(axis.text.y = element_text(size=5),
              axis.text.x = element_text(size=6, face="bold"),
              axis.title.x = element_text(size=7),
              plot.title = element_text(size=7, face="bold", hjust=0.5),
              legend.position = "right") +
        #scale_fill_gradient2(low="blue", mid="white",high="red", midpoint = 0, limits=c(-lim,lim),na.value = "black") +
        #scale_color_gradient2(low="blue", mid="white",high="red", midpoint = 0, limits=c(-lim,lim),na.value = "black") +
        scale_fill_distiller(palette = brewer_palette, direction = d, na.value = "white", limits=c(minval,maxval)) +
        scale_color_distiller(palette = brewer_palette, direction = d, na.value = "white", limits=c(minval,maxval)) +
        ggtitle(paste0(as.character(input$plottitlehmgb)," -/+ ",sizekb)) + facet_grid(Peaks~variable, scales="free_y")
    })
    
    output$Heatmapgb <- renderPlot({
      gghmgb()
    })
    
    output$downloadhmgb <- downloadHandler(
      filename = function() {
        paste0(input$fileprefixgb,"_",sizekb, "_Heatmap.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthhmgb), height=as.numeric(input$pheighthmgb))
        print(gghmgb())
        dev.off()
      }
    )
    
  })
  
  
  
  
  observeEvent(input$bpgb, {
    
    bws <- rbwfiles$data 
    bpath <- rbedfiles$data
    
    bw_names <- unlist(strsplit(input$namesgb,","))
    bed_names <- unlist(strsplit(input$bnamesgb,","))
    
    bps <- as.numeric(input$basepairsgb)
    bin <- as.numeric(input$binsizegb)
    
    half_bps <- bps/2
    
    half_bps <- bps/2
    
    sizekb <- paste0( (bps/1000),"kb" )
    
    ## save data to individual variables for plotting
    dat.boxplot <- rResult$data[[3]]
    
    y_axis_title=  as.character(input$yaxisbpgb) # y axis label for average profile
    
    var_order = input$v_orderbpgb  %>% as.character %>% strsplit(.,",") %>% unlist() # order to show samples (needs to be exact sample as "vars" input in function)
    
    regions_order = input$b_orderbpgb %>% as.character %>% strsplit(.,",") %>% unlist()
    
    brewer_palette=  as.character(input$brewer_palettebpgb) # brewer color palette used for heatmap
    
    
    dat.boxplot$Peaks <- factor(dat.boxplot$Peaks, levels=c(regions_order)) # set regions factors
    dat.boxplot$variable <- factor(dat.boxplot$variable, levels=c(var_order)) # set regions factors
    
    
    if(input$bp_maxgb == "Default"){
      maxval <- max(dat.boxplot$Value)
    }else{
      maxval <- as.numeric(input$bp_maxgb)
    }
    
    if(input$bp_mingb == "Default"){
      minval <- min(dat.boxplot$Value)
    }else{
      minval <- as.numeric(input$bp_mingb)
    }
    
    
    
    
    
    ss<- NULL
    ss1<- NULL
    
    
    if(length(bw_names)>1){
      stat.test <- as.data.frame(dat.boxplot) %>%
        group_by(Peaks) %>%
        wilcox_test(Value ~ variable, paired=T) %>%
        add_significance()
      
      ss <- tableGrob(stat.test) 
    }else if(length(bw_names)==1){
      ss <- NULL
    }
    
    if(length(bed_names)>1){
      stat.test <- as.data.frame(dat.boxplot) %>%
        group_by(variable) %>%
        wilcox_test(Value ~ Peaks, paired=F) %>%
        add_significance()
      
      ss1 <- tableGrob(stat.test) 
    }else if(length(bed_names)==1){
      ss1 <- NULL
    }
    
    
    
    if(input$facetbpgb == "Peaks"){
      
      colorvarbp <- "variable"
      facetvalbp <- "Peaks"
      namevalbp <- "Treatment"
      
      if(length(bw_names)>3){
        mycolor <- rev(brewer.pal(length(bw_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
      
    }else if(input$facetbpgb == "variable"){
      colorvarbp <- "Peaks"
      facetvalbp <- "variable"
      namevalbp <- "Regions"
      
      if(length(bed_names)>3){
        mycolor <- rev(brewer.pal(length(bed_names),brewer_palette))
      }else{
        mycolor <- rev(brewer.pal(3,brewer_palette))
      }
      
      
    }
    
    
    
    ggbp <- reactive({ggplot(na.omit(dat.boxplot), aes(x = get(colorvarbp), y = (Value), fill=get(colorvarbp))) +
        geom_boxplot() +
        theme_classic() +
        theme(#axis.text.x  = element_text(size=12, angle=90, hjust=1, vjust=0.5),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=5),
          axis.title = element_text(angle = 0, hjust = 0.5, size=15),
          legend.position = "right",
          legend.text = element_text(size=7),
          legend.title = element_text(size=9),
          plot.title = element_text(hjust = 0.5, face="bold", size=7)) +
        scale_fill_manual(values = mycolor, name=namevalbp) +
        xlab(namevalbp) + ggtitle(paste0(as.character(input$plottitlebpgb)," -/+ ",sizekb)) + ylab(y_axis_title) +
        facet_grid(~get(facetvalbp)) + ylim(c(minval,maxval))
      
    })
    
    
    if(length(bw_names)>1 & length(bed_names)>1){
      bpnew <- reactive({plot_grid(ggbp(),ss,ss1,nrow=3,rel_heights = c(1,0.3,0.3))})
    }else if(length(bw_names)>1 & length(bed_names)==1){
      bpnew <-  reactive({plot_grid(ggbp(),ss,nrow=2,rel_heights = c(1,0.3))})
    }else if(length(bw_names)==1 & length(bed_names)>1){
      bpnew <- reactive({plot_grid(ggbp(),ss1,nrow=2,rel_heights = c(1,0.3))})
    }
    
    
    output$Boxplotgb <- renderPlot({
      bpnew()
    })
    
    
    output$downloadbpgb <- downloadHandler(
      filename = function() {
        paste0(input$fileprefixgb,"_",sizekb, "_Boxplot.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthbpgb), height=as.numeric(input$pheightbpgb))
        print(ggbp())
        dev.off()
      }
    )
    
    
    output$downloadbpstatsgb <- downloadHandler(
      filename = function() {
        paste0(input$fileprefixgb,"_",sizekb, "_withstats_Boxplot.pdf")
      },
      content = function(file) {
        pdf(file=file, width=as.numeric(input$pwidthbpgb), height=as.numeric(input$pheightbpgb))
        print(bpnew())
        dev.off()
      }
    )
    
  })
  
  
  
  
  
  
  
  observeEvent(input$clear, {
    
    output$matrix_status <- renderText({" "})
    output$Heatmap <- renderPlot({" "})
    output$AvgProf <- renderPlot({" "})
    output$Boxplot <- renderPlot({" "})
    
    output$matrix_statusgb <- renderText({" "})
    output$Heatmapgb <- renderPlot({" "})
    output$AvgProfgb <- renderPlot({" "})
    output$Boxplotgb <- renderPlot({" "})
    
  })
  
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
