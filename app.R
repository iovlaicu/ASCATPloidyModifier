library(shiny)
library(ASCAT.sc)
library(ASCAT)
library(shinybusy)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(shinyalert)
library(shinydashboard)
library(shinyFiles)
library(spatstat.geom)
library(dipsaus)
#source("plotSunrise2.R")
#source("myplotSol.R")
source("PlotASCATAdjusted.R")
source("plotSunriseASCAT.R")
source("PlotLogRASCAT.R")
source("refitProfile_ASCAT.R")
source("plotSegmentedASCAT.R")

reslocal <- NULL
bclocal <- NULL
pathrdata <- NULL
rdata <- NULL
coords <- NULL
localOpt <- NULL
logRplot <- NULL
newPlot <- NULL
highlight <- NULL
options(spinner.color="#0F829F", spinner.color.background="#ffffff", spinner.size=1)

getIndex <- function(sample){
  
  index <- which(names(reslocal$allTracks.processed)==sample)
  return (index)
  
}

getSamples <- function() {
  
  if( !is.null(reslocal)){
    
    return (names(reslocal$purity))}
  
}

#######################################################################

######################### INTERFACE ###################################

#######################################################################

ui <- navbarPage(id="nav_page",
                 title="ASCATFit",   theme = shinytheme("spacelab"),
                 tabPanel("Welcome",
                          busy_start_up(
                            loader = tags$img(
                              src = "LoaderASCAT.gif",
                              width = 700
                            ),
                            text = "Loading ...",
                            color = "#0F829F",
                            timeout = 3000,
                            background = "white",
                            mode = "auto"
                          ),
                          tags$head(
                            tags$style(HTML(
                              ".hover-effect {color: #FFFFFF ; border-radius: 30px; height:150px; width:300px; background-color: #3e0533; border-color: #0F829F; padding:5px; font-size:4vh; border-width: 5px;  opacity: 0.90;
                              }",
                              ".hover-effect:hover {color: #FFFFFF ; border-radius: 30px; background-color: #3e0533; border-color: #0F829F; padding:5px; font-size:5vh; border-width: 5px;
                        opacity: 1;}",
                              ".button-container {
             width: 300px;}",
                              ".click-effect:active {
      color: #FFFFFF ; border-radius: 30px; background-color: #011b3d; border-color: #0F829F; padding:5px; font-size:5vh; border-width: 5px;
                        opacity: 1;
    
    }",
                              "code {
                display:block;
                padding:9.5px;
                margin: auto;
                width: 1100px;
                height: 500px;
                font-size:14px;
                color: #011b3d;
                line-height:4px;
                word-break:break-all;
                word-wrap:break-word;
                white-space:pre-wrap;
                background-color:#FFFFFF;
                border:6px solid #011b3d;
                border-radius:4px; 
            }",
                              "pre {
                  display:block;
                  padding:9.5px;
                  margin: auto;
                  width: 350px;
                  height: 100px;
                  font-size:14px;
                  color: #833e03;
                    line-height:5px;
                  word-break:break-all;
                  word-wrap:break-word;
                  white-space:pre-wrap;
                  background-color:#fae6d4;
                    border:4px solid #833e03;
                  border-radius:4px; 
                }",
                              "em {
                  display:block;
                  padding:9.5px;
                  margin: auto;
                  width: 800px;
                  height: 95px;
                  font-size:14px;
                  color: #833e03;
                    line-height:5px;
                  word-break:break-all;
                  word-wrap:break-word;
                  white-space:pre-wrap;
                  background-color:#fae6d4;
                    border:4px solid #833e03;
                  border-radius:4px; 
                }"
                            ))),
                          #######################################################################
                          
                          ######################### WELCOME TAB ################################
                          
                          #######################################################################
                          img(src='ASCAT_background.jpg', align = "center", style="width:100%; max-width:100%; position: absolute; z-index:-1;"),
                          div(style = "height:70px"),  
                          # ba4a00
                          fluidRow(column(6, align="center",h1(strong("ASCAT",style={'color: #0F829F; font-family: arial black ,sans-serif; text-shadow:
  # # 0 0 7px #fff,
  # # 0 0 10px #fff,
  # # 0 0 21px #fff,
  # # 0 0 42px ##0F829F,
  # # 0 0 82px ##fff;  font-size: 80px'}))), column(6, align="center",
                                                    h1(strong("Ploidy Modifier",style={'color: #0F829F; font-family: arial black ,sans-serif; text-shadow:
  # # 0 0 7px #fff,
  # # 0 0 10px #fff,
  # # 0 0 21px #fff,
  # # 0 0 42px ##0F829F,
  # # 0 0 82px ##fff;  font-size: 80px'}))))
                          ,br(), br(), br(), div(style = "height:50px"),
                          fluidRow(column(6, align="center", dropdown(code(h4(strong("1. Choose the data", style={'font-family: Arial'}), align = "left"),
                                                                           p("Please choose the ASCAT rdata object to load. You may load objects that contain the ASCAT results for only one sample.", align = "left"),  
                                                                           p("Press Start to begin modifying the profiles. Note that it takes a few seconds to load the object.", align = "left"),
                                                                           h4(strong("2. Modify the profiles \n ", style={'font-family: Arial; align-text: left'}), align = "left"),
                                                                           p("The purity & ploidy of the whole sample can be changed by clicking on the sunrise plot.", align = "left"), 
                                                                           p("You can additionally modify the copy number of 2 different segments (the major and minor allele of one segment).", align = "left"),
                                                                           p("Modifications can be done either by numeric input or directly on the Modifier Station plot.", align = "left"),
                                                                           h4(strong("3. Save the modified profiles ", style={'font-family: Arial'}), align = "left"),
                                                                           p("Once you are happy with the results, you can save the ASCAT rdata object.", align = "left")),
                                                                      size= "lg", circle= FALSE, status = "info", label = "Help", width = "1200px", inputId="help")),
                                   column(6, align="center", shinyFilesButton("get_file", "Load data" ,
                                                                              title = "Choose the ASCAT rdata object to load", multiple = FALSE,
                                                                              buttonType = "primary", class = NULL, style = "width:60px, height: 50px")
                                          
                                   )),
                          br(), div(style = "height:100px"),
                          fluidRow(column(6, align="center", offset = 3,
                                          
                                          actionButtonStyled("start",
                                                             label="Start",
                                                             class="hover-effect click-effect",
                                                             width = "300px"
                                          ))
                                   
                          )
                          
                 ),
                 
                 #######################################################################
                 
                 ######################### MODIFIER TAB ################################
                 
                 #######################################################################
                 
                 tabPanel("Modifier",  shinyWidgets::useShinydashboard(), 
                          fluidRow(box(width=12, title="Original", status="primary", solidHeader=TRUE, 
                                       column(width=8,withSpinner(plotOutput("profile"),type=3)),
                                       column(width=4,plotOutput("sunrise1", height='450px')))
                                   
                          ),
                          useShinyjs(),
                          useShinyalert(),
                          fluidRow(
                            
                          #   column(width=2, offset = 1,
                          #                 dropdownButton(
                          #                   tags$h3("Choose Sample"),
                          #                   br(),
                          #                   selectInput(
                          #                     "samples",
                          #                     label= NULL,
                          #                     choices= getSamples(),
                          #                     selected = NULL,
                          #                     multiple = FALSE,
                          #                     selectize = FALSE
                          #                   ), size= "lg", circle= FALSE, status = "info", icon = icon("list",verify_fa = FALSE), width = "450px",
                          #                   label= "Select sample",
                          #                   actionButton("view", label = "View", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 2px")
                          #                   
                          #                 )
                          # ),
                          column(width = 1, offset = 1, dropdownButton(tags$h3("Customize plot visualization"),
                                 checkboxInput("plotRaw", "Plot raw segments", FALSE),
                                
                                 checkboxInput("hideSeg", "Hide main segments", FALSE),
                              
                                 checkboxInput("logR", "View raw LogR", FALSE),
                                 checkboxInput("BAF", "View raw BAF", FALSE),
                                 size= "lg", circle= FALSE, status = "info", label = "Settings", width = "250px"
                                 )),
                          column(width = 2, offset = 1,
                                 dropdownButton(tags$h3("Choose ploidy and purity on the sunrise graph"),
                                                br(),
                                                h5("To change the ploidy and purity of the whole sample directly on the Working station profile, please click the point on the sunrise plot corresponding to the desired ploidy/purity values. When 'Automatic optima' is enabled; the closest optimal purity/ploidy pair will be selected."),
                                                br(),
                                                checkboxInput("optima", "Automatic optima", FALSE),
                                                
                                                 #checkboxInput("sizeP", "Increase size of datapoints", FALSE),
                                                size= "lg", circle= FALSE, status = "info", label = "Modify purity & ploidy", width = "450px"
                                 )
                                 
                          ),
                          
                          column(width = 5, offset=1,
                                 dropdown(tags$h3("Choose a different way to modify the profile:"),
                                          
                                          dropdownButton(tags$h3("Modify the copy number of 1 segment"),br(),
                                                         selectInput(
                                                           "Chr",
                                                           "Choose the chromosome",
                                                           c(1:22, "X", "Y"),
                                                           selected = NULL,
                                                           multiple = FALSE,
                                                           selectize = FALSE,
                                                           width = "75%"
                                                           
                                                         ),
                                                         textInput("cn1", "Major allele copy number", value = "", placeholder = NULL, width = "75%"),
                                                        
                                                         textInput("cn2", "Minor allele copy number", value = "",  placeholder = NULL, width = "75%"),
                                                         actionButton("modify", label = "Apply", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                                         size= "lg", circle= FALSE, status = "info",right= TRUE, label = "Refit segments", width = "450px"),
                                          br(),
                                          dropdownButton(tags$h3("Modify segment on graph"),
                                                         h4("To modify the copy number of 2 different segments directly on the Working station profile, please click on the segment you wish to modify, then on its desired position. Repeat for the second segment, then click on 'Apply'"),
                                                         #actionButton("start_refit", label = "Enable clicks", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                    
                                                         actionButton("refit", label = "Apply", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #ba4a00; font-size:100%; border-width: 3px"),
                                                         size= "lg", circle= FALSE, status = "info", right= TRUE, label = "Refit segments on graph", width = "450px"
                                          ),
                                          br(),
                                          
                                          size= "lg", circle= FALSE, status = "info", label = "Additional tools", width = "800px", inputId="menu", up=TRUE))
                          
                          ), br(),
                          
                          
                          fluidRow(box(width=12, title="Modifier Working Station", status="primary", solidHeader=TRUE,column(width=8, withSpinner(plotOutput("profile2", click="profile2_click"),type=3)), column(width = 4, plotOutput("sunrise2", click = "sunrise2_click", height="460px")))),
                          fluidRow(column(width = 3,  offset=1, actionButton("discard", label = "Reset profile",  icon = icon("arrows-rotate",verify_fa = FALSE), style="color: #FFFFFF ; background-color: #ba4a00; border-color: #0F829F; font-size:100%; border-width: 3px")),
                                   #column(width = 3, offset=1,  downloadButton("savetxt", label = "Save profiles", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #0F829F; font-size:100%; border-width: 3px")),
                                   column(width = 3, offset=1, downloadButton("save", label = "Save .Rda", style="color: #FFFFFF ; background-color: #ba4a00; border-color: #0F829F; font-size:100%; border-width: 3px"))), br(), br()))


#######################################################################

######################### SERVER ######################################

#######################################################################


server <- function(input, output, session) {
  
  vals <- reactiveVal()
  volumes = getVolumes()
  coords <- reactiveValues(x=NULL,y=NULL)
  output$profile <- renderPlot(NULL)
  output$profile2 <- renderPlot(NULL)
  
  observeEvent(input$profile2_click, {
    coords$x <- c(coords$x,input$profile2_click$x)
    coords$y <- c(coords$y,input$profile2_click$y)
  })
  
  observeEvent(input$sunrise2_click, {
    coords$x <- input$sunrise2_click$x
    coords$y <- input$sunrise2_click$y
  })
  
  sampleName <- reactive({
    input$samples
  })
  optValue <- reactive({input$optima})
  result <- reactive({ reslocal })
  
  shiftv <- reactive({
    input$ploidy
  })
  
  chrs <- reactive({
    list(input$Chr,input$cn1, input$cn2)
  })
  
  plotlog <- reactiveVal(FALSE)
  plotBAF <- reactiveVal(FALSE)
  plotRawSeg <- reactiveVal(FALSE)
  hideCN <- reactiveVal(FALSE)
  
  #########################################################################
  
  ######################## FILE CHOOSER ###################################
  
  #########################################################################
  
  observe({
    shinyFileChoose(input, "get_file", roots = volumes, session = session)
    
    if(!is.null(input$get_file)){
      
      file_selected<-parseFilePaths(volumes, input$get_file)
      pathrdata <<- file_selected
      
      
    }
  })
  
  observe({ toggle(id="start", condition=(input$get_file>=1))})
  
  #########################################################################
  
  ######################## START ##########################################
  
  #########################################################################
  
  observeEvent(input$start, {
    updateNavbarPage(session=session,
                     inputId="nav_page",
                     selected="Modifier")
    tryCatch({
      
      filepath <<- as.character(pathrdata$datapath)
      
      load(filepath)
      
      reslocal <<- ascat.output
      
      bclocal <<- ascat.bc
      
      updateSelectInput(session, "samples", label=NULL,
                        choices=getSamples())
     
      
      if(!"manual"%in%names(reslocal))
      {
        original <- reslocal
        reslocal <<- NULL
        reslocal$original <<- original
        reslocal$manual <<- original
  
      }
      
      # if(!"manual"%in%names(bclocal))
      # {
      #   bclocal$manual <<- bclocal
      # 
      # }
      
    
      output$profile <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$original, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38")) })
          
      output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$manual$purity)[[1]], REF="hg38")) })
      
      #plotScGrid(reslocal$allSolutions[[1]])
      output$sunrise1 <- renderPlot({isolate(ascat.plotSunrise(reslocal$original$distance_matrix[[1]], psi_opt1 = reslocal$original$ploidy[[1]], rho_opt1 = reslocal$original$purity[[1]]))})
      #plotScGrid(reslocal$allSolutions.refitted.manual[[1]])
      output$sunrise2 <- renderPlot({localOpt <<- isolate( ascat.plotSunrise(reslocal$manual$distance_matrix[[1]], psi_opt1 = reslocal$manual$ploidy[[1]], rho_opt1 = reslocal$manual$purity[[1]], optima=TRUE))})
      
        
      removeModal()
        
    },
    error=function(e) {
      
    })
    
  })
  
  #########################################################################
  
  ######################## DISCARD/KEEP ###################################
  
  #########################################################################
  
  observeEvent(input$discard, {
    
    vals(input$samples)
    
    #if(! is.null(sampleName())){

      
      output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$original, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg = hideCN(), plot_unrounded = plotRawSeg())) })
      
      
      output$sunrise2 <- renderPlot({isolate( localOpt <<- ascat.plotSunrise(reslocal$original$distance_matrix[[1]], psi_opt1 = reslocal$manual$ploidy[[1]], rho_opt1 = reslocal$manual$purity[[1]], optima=TRUE))})
      
      
      reslocal$manual <<- reslocal$original
      #bclocal$manual <<- bclocal
    
    #}
      coords$y <<- NULL
      coords$x <<- NULL
      highlight <<- NULL
    #else{
      #return (NULL)
    #}
    
  })
  
  #########################################################################
  
  ######################## DOWNLOAD #######################################
  
  #########################################################################
  
  
  output$save <- downloadHandler(
    filename = function() {
      "result_manualfitting.Rda"
    },
    content = function(file) {
      
      showModal(modalDialog(div(tags$b("Loading...", style = "color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      
      res <- reslocal
      bc <- bclocal
      save(res, bc, file=file)
      #save(res, file=file)
      
    }
  )
  output$savetxt <- downloadHandler(
    filename = function(){
      "Profiles_txt.zip"
    },
    content = function(file){
      
      showModal(modalDialog(div(tags$b("Loading...", style = "color: steelblue;")), footer=NULL))
      on.exit(removeModal())
      
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      
      for ( i in 1:length(reslocal$allProfiles.refitted.manual)){
        
        
        
        write.table(reslocal$allProfiles.refitted.manual[[i]],
                    quote=F,
                    sep="\t",
                    col.names=T,
                    row.names=F,
                    file=paste0(temp_directory,"/",
                                paste0(names(reslocal$allTracks.processed)[i],"-manual_refit"),
                                ".ASCAT.scprofile.txt"))
      }
      
      zip::zip(
        zipfile = file,
        files = dir(temp_directory),
        root = temp_directory
      )
      
    },
    contentType = "application/zip"
    
  )
  
  
  #########################################################################
  
  ######################## SUNRISE ########################################
  
  #########################################################################
  
  
  observeEvent(input$sunrise2_click,{
   
    
    purity <- NULL
    ploidy <- NULL
    
    
    tryCatch(
      {
        
        d <- reslocal$manual$distance_matrix[[1]]
        #print(d)
        #print("####### Matrix before refit #######")
        #print(reslocal$manual$distance_matrix[[1]])
        
        purity <- coords$y
        ploidy <- coords$x
        
        ploidy_min<-as.numeric(rownames(d)[1])
        ploidy_max<-as.numeric(rownames(d)[nrow(d)])
        purity_min<-as.numeric(colnames(d)[1])
        purity_max<-as.numeric(colnames(d)[ncol(d)])
        
        ploidy <- as.numeric(ploidy*(ploidy_max - ploidy_min) + ploidy_min)
        purity <- as.numeric(purity*(purity_max - purity_min) + purity_min)
        
        ### automatic optima ###
        
        if(optValue()){
          
          best <- 1
          #print(localOpt)
          #print(localOpt$bao)
          
          dist <- crossdist(ploidy, purity, as.numeric(rownames(d)[localOpt[best, 1]]), as.numeric(colnames(d)[localOpt[best, 2]]))
          #print(dist)
          
          for (i in 1:nrow(localOpt)){
            
           
            dist2 <- crossdist(ploidy, purity, as.numeric(rownames(d)[localOpt[i, 1]]), as.numeric(colnames(d)[localOpt[i, 2]]))
           
            
            if (dist >= dist2){
              
              dist <- dist2
              best <- i
            }
            
          }
          
          
          ploidy <- as.numeric(rownames(d)[localOpt[best, 1]])
          purity <- as.numeric(colnames(d)[localOpt[best, 2]])
          
        }
       
        newprofile <- ascat.runAscat(bclocal, rho_manual=purity, psi_manual=ploidy)
        
        reslocal$manual$ploidy[[1]] <- ploidy
        reslocal$manual$purity[[1]] <- purity
     
        reslocal$manual$segments <<- newprofile$segments
      
        reslocal$manual$psi <<- newprofile$psi
        
        #print(newprofile$distance_matrix[[1]])
        #reslocal$manual$distance_matrix[[1]] <<- newprofile$distance_matrix[[1]]
        
        #print("####### Matrix after refit #######")
        #print(reslocal$manual$distance_matrix[[1]])
        #print(reslocal$manual$distance_matrix)
        
        output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(newprofile, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg())) })
        
          output$sunrise2 <- renderPlot({isolate( localOpt <<- ascat.plotSunrise(reslocal$manual$distance_matrix[[1]], psi_opt1 = reslocal$manual$ploidy[[1]], rho_opt1 = reslocal$manual$purity[[1]], optima=TRUE))})
        
      
      coords$y <<- NULL
      coords$x <<- NULL
      },
      error=function(e) {
        message('An Error Occurred')
        print(e)
        shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type = "error")
      },
      warning=function(w) {
        message('A Warning Occurred')
        print(w)
        shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
      })
  
    
  })
  
  #########################################################################
  
  ######################## MODIFY SEGMENTS ON GRAPH #######################
  
  #########################################################################
  

  
  observeEvent(input$profile2_click,{
    
    
    tryCatch(
      {
      
        #input$profile2_click
        points_x <- coords$x
        points_y <- coords$y
        #print(points_x)
        
       
        if(length(points_x) == 1){
          
          highlight$x <<- points_x
          highlight$y <<- points_y
          
          
          output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg(), highlight = highlight))})
            
          
        }
        else if (length(points_x) == 2){
          
          output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg(), highlight=highlight)) 
            
        
       
              
              x <- points_x[2]
              y <- points_y[2]
              
              points(x, y, col = "chartreuse", pch = 4, cex = 4)
              
          })
         
        }
        
        else if(length(points_x) == 3){
          
          
          highlight$x <<- c(highlight$x,points_x[3])
          highlight$y <<- c(highlight$y,points_y[3])
          
          output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg(), highlight = highlight))
          x <- points_x[2]
          y <- points_y[2]
          
          points(x, y, col = "chartreuse", pch = 4, cex = 4)
          
          })
        }
        
     
        if(length(points_x)==4){
          
    
          output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg(), highlight=highlight)) 
            

          x <- points_x[c(2,4)]
          y <- points_y[c(2,4)]
            
          points(x, y, col = "chartreuse", pch = 4, cex = 4)
            
       
          })
          
        }
        if(length(points_x)>4){
          coords$y <<- NULL
          coords$x <<- NULL
          highlight <<- NULL
        }
          })
    
})

  
  
  observeEvent(input$refit,{
    
      input$profile2_click
      chr1 <- NULL
      chr2 <- NULL
      #breaks <- c(0, cumsum(sapply(reslocal$allTracks.processed[[1]]$lSegs, function(x) max(x$output$loc.end))/1e+06))
      #coords$x
      tryCatch(
        {
          output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg())) })
          
         #if (REF=="hg38") {
          REF=data.frame(chrom=c(1:22, "X"),
                           start=rep(1, 23),
                           end=c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                                 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                                 58617616, 64444167, 46709983, 50818468, 156040895))
          
          REF$size=REF$end-REF$start+1
          # REF$middle=0
          # for (i in 1:nrow(REF)) {
          #   if (i==1) {
          #     REF$middle[i]=REF$size[i]/2
          #   } else {
          #     REF$middle[i]=sum(as.numeric(REF$size[1:(i-1)]))+REF$size[i]/2
          #   }
          # }; rm(i)
          #REF$cumul=cumsum(as.numeric(REF$size))
          
          REF$add=cumsum(as.numeric(c(0, REF$size[1:(nrow(REF)-1)])))
          
           
          chr <- length(REF[REF$add <= as.numeric(coords$x[1]), "chrom"])
          
          real_x_pos <- as.numeric(coords$x[1] - REF$add[which(REF$chrom==chr)])
          
          
          rescopy <- reslocal
          rescopy$original <- rescopy$manual
          
          segment_info <- reslocal$manual$segments[reslocal$manual$segments$startpos <= real_x_pos & reslocal$manual$segments$endpos >= real_x_pos & reslocal$manual$segments$chr == chr, ]
          #chr <- segment_info$chr
          #y2 <- round(coords$y[1], digits=0)
          y1 <- coords$y[1]
          y2 <- round(coords$y[2], digits=0)
          y3 <- coords$y[3]
          y4 <- round(coords$y[4], digits=0)
          
          refMajor <- NULL
          refMinor <- NULL
          
          if (segment_info$nMajor == segment_info$nMinor) {
            
            refMajor <- max(y2, y4)
            refMinor <- min(y2, y4)
            
          }
          else {
            
            if(abs(segment_info$nMajor -  y1) <=  abs(segment_info$nMinor - y1)) {
              refMajor <- y2
            }
            else{
              refMinor <- y2
            }
            
            if(abs(segment_info$nMajor -  y3) <=  abs(segment_info$nMinor - y3)) {
              refMinor <- y4
            }
            else{
              refMajor <- y4
            }
            
            
          }
          
          
          seg <- reslocal$manual$segments[reslocal$manual$segments$"chr"==chr, ]
          ind <- which(seg$startpos==segment_info$startpos)
           
          
          rescopy$manual <- refitProfileASCAT(rescopy$manual, refMinor=refMinor, refMajor=refMajor, chr= chr, gamma=0.55, wm=ind) 
          
          rescopy$original <- reslocal$original
          
          reslocal <<- rescopy
          
          output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg())) })
          
          output$sunrise2 <- renderPlot({isolate(localOpt <<- ascat.plotSunrise(reslocal$manual$distance_matrix[[1]], psi_opt1 = reslocal$manual$ploidy[[1]], rho_opt1 = reslocal$manual$purity[[1]]))})
          
          
          coords$y <<- NULL
          coords$x <<- NULL
        
          highlight <<- NULL
        },
        error=function(e) {
          message('An Error Occurred')
          print(e)
          shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type = "error")
        },
        warning=function(w) {
          message('A Warning Occurred')
          print(w)
          shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
        })
   
  
  })
  
  #########################################################################
  
  ######################## MODIFY SEGMENTS ################################
  
  #########################################################################
  
  observeEvent(input$modify,{
    
    vals <- as.numeric(chrs())
    
    
    if(! is.null(chrs())){
      #index <- getIndex(sampleName())
      tryCatch(
        {
         
        rescopy <- reslocal
        rescopy$original <- rescopy$manual
        
       
        rescopy$manual <- refitProfileASCAT(rescopy$manual, refMinor=vals[[3]], refMajor=vals[[2]], chr= vals[[1]], gamma=0.55) 
        
        rescopy$original <- reslocal$original
       
        reslocal <<- rescopy
       
        output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg())) })
        
        output$sunrise2 <- renderPlot({isolate(localOpt <<- ascat.plotSunrise(reslocal$manual$distance_matrix[[1]], psi_opt1 = reslocal$manual$ploidy[[1]], rho_opt1 = reslocal$manual$purity[[1]]))})
        
        
        },
        error=function(e) {
          message('An Error Occurred')
          print(e)
          shinyalert("Error", "Cannot fit profile: ploidy<0 or purity ∉ [0,1]. Please choose different values", type = "error")
        },
        warning=function(w) {
          message('A Warning Occurred')
          print(w)
          shinyalert("Warning", "New solution is ambiguous: reverted to old one", type = "error")
        })
    }
    else{
      return (NULL)
    }
    
  })


#########################################################################

######################## VIEW RAW DATA ################################

#########################################################################

observeEvent(input$logR,{
  vals(input$samples)
  #sample <- sampleName()
  #if(! is.null(sampleName())){
    plotlog(input$logR)
    #index <- getIndex(sampleName())
    if(input$logR) {
     
        output$profile <- renderPlot({newPlot <<- myascat.plotSegmentedData(bclocal)})

    }
    
    else {
     
      output$profile <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$original, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38")) })
      
    }
   
  #}
 
})
  

  observeEvent(input$BAF,{
    vals(input$samples)
    #sample <- sampleName()
    #if(! is.null(sampleName())){
      plotBAF(input$BAF)
      #index <- getIndex(sampleName())
      if(input$BAF) {
        

          output$profile <- renderPlot({newPlot <<- myascat.plotSegmentedData(bclocal, plotLogR=FALSE)})

      }
      
      else {
        
        output$profile <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$original, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38")) })
        
      }
      
    #}
   
  })

  #########################################################################
  
  ######################## PLOT RAW SEGMENTS ################################
  
  #########################################################################
  
  observeEvent(input$plotRaw,{
    vals(input$samples)
    #sample <- sampleName()
    #if(! is.null(sampleName())){
      plotRawSeg(input$plotRaw)
      #index <- getIndex(sampleName())
   
        output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg())) })
        
    
      
    #}
    
  })
  
  #########################################################################
  
  ######################## HIDE CN ################################
  
  #########################################################################
  observeEvent(input$hideSeg,{
    vals(input$samples)
    #sample <- sampleName()
    #if(! is.null(sampleName())){
      hideCN(input$hideSeg)
      #index <- getIndex(sampleName())
  
        
        
        output$profile2 <- renderPlot({isolate(plot2 <- ascat.plotAdjustedAscatProfile(reslocal$manual, SAMPLE=names(reslocal$original$purity)[[1]], REF="hg38", hideSeg=hideCN(), plot_unrounded = plotRawSeg())) })
        
      
    #}
   
  })
  
} 

shinyApp(ui = ui, server = server)
