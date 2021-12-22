#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

#devtools::install_github("kcha193/infoDecompuTE")

library(shiny)
library(infoDecompuTE)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  des <- reactive({

    if(input$phase1Des == "CRD")
      read.csv(paste0("OptimalDesigns/", input$phase1Des, 
             paste0("/design", input$phase1Des, input$trt, 
                    input$bioRep, 2, input$tag, ".csv")))
    else 
      read.csv(paste0("OptimalDesigns/", input$phase1Des, 
                      paste0("/design", input$phase1Des, input$trt, 
                             input$bioRep, input$cageRep, 2, input$tag, ".csv")))
    
  })
  
  output$desTablePhase1 <- renderTable({
    des <- des()
    
    if(input$phase1Des == "CRD")
      desTable <- data.frame(matrix(sort(unique(paste0(des$Ani, des$Trt))), ncol = 4))
    else 
      desTable <- data.frame(matrix(sort(unique(paste0(des$Cag, des$Ani, des$Trt))),
                                    ncol = 4))
    
    colnames(desTable) <- NULL
    
    desTable
  })
  
  
  output$desTable <- renderTable({
    des <- des()
   
    Run <- unique(des$Run)
    Tag <- unique(des$Tag)
    
    if(input$phase1Des == "CRD")
      desTable <- data.frame(matrix(paste0(des$Ani, des$Trt), 
                                    nrow = length(Run), ncol = length(Tag)))
    else 
      desTable <- data.frame(matrix(paste0(des$Cag, des$Ani, des$Trt),
                                    nrow = length(Run), ncol = length(Tag)))
    
    rownames(desTable) <- Run
    colnames(desTable) <- Tag
    
    desTable
  }, rownames = TRUE)
  
  output$anovaTable <- renderPrint({
    
    if(input$phase1Des == "CRD")
      summaryAovTwoPhase(des(), blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Trt", list.sep = FALSE)
    else 
      summaryAovTwoPhase(des(), blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Trt", list.sep = FALSE)
  })
  
  
  
})
