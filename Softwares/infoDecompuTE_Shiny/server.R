
#devtools::install_github("kcha193/infoDecompuTE")

library(shiny)
library(infoDecompuTE)
library(knitr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  des <- reactive({
    
    if(input$checkbox){
      des <- read.csv("design2.csv", row.names = NULL)
      
    } else {
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      des <- read.csv(inFile$datapath, row.names = NULL)
    }
    
    des
  })
  
  output$desDf <- renderTable({
    req(des())
  })
  
  output$expStr <- renderUI({
    if(input$phase == 1){
      inputPanel(  
        textInput("p1blk", label = h3("Block Structure"), value = "Ani/Sam"),
        textInput("trt", label = h3("Treatment structuret"), value = "Tag + Trt")
      )
    } else {
      inputPanel(  
        textInput("p2blk", label = h3("Phase 2 Block Structure"), value = "Run"),
        textInput("p1blk", label = h3("Phase 1 Block Structure"), value = "Ani/Sam"),
        textInput("trt", label = h3("Treatment structuret"), value = "Tag + Trt")
      )
    } 
  })
  
  outputPrint <- reactive({
    if(input$phase == 1){
      summaryAovOnePhase(req(des()), blk.str = input$p1blk, trt.str = input$trt, 
                                             list.sep = FALSE) 
      
    } else {
      summaryAovTwoPhase(req(des()), blk.str1 = input$p1blk, blk.str2 = input$p2blk, 
                                             trt.str = input$trt, list.sep = FALSE)   
    }
    
  })
  
  latexPrint <- reactive({
    if(input$phase == 1){
      summaryAovOnePhase(req(des()), blk.str = input$p1blk, trt.str = input$trt, 
                         list.sep = FALSE, latex = TRUE) 
      
    } else {
      summaryAovTwoPhase(req(des()), blk.str1 = input$p1blk, blk.str2 = input$p2blk, 
                         trt.str = input$trt, list.sep = FALSE, latex = TRUE)   
    }
    
    
  })
  
  output$anovaTable <- renderPrint({

      outputPrint()
  })
  
  output$latexUI <- renderUI({
    textAreaInput("latexOutput", "Latex code:", value = latexPrint(),
                  width = '1000px', height = '400px')
  })
  
  output$downloadOutput <- downloadHandler(
    filename = function() { 
      paste('output.txt', sep='') 
    },
    content = function(file) {
      write.table(capture.output(outputPrint()), file, row.names = FALSE, quote = FALSE)
    }
  )
  
  output$downloadLatexCode <- downloadHandler(
    filename = function() { 
      paste('latexCode.txt', sep='') 
    },
    content = function(file) {
      write.table(latexPrint(), file, row.names = FALSE, quote = FALSE)

    }
  )
  
  output$latexTable = downloadHandler(
    filename = 'table.pdf',
    
    content = function(file) {
      out = knit2pdf('input.Rnw', clean = TRUE)
      file.rename(out, file) # move pdf to file for downloading
    },
    
    contentType = 'application/pdf'
  )
  
})
