
#devtools:::install_github("kcha193/infoDecompuTE")

library(shiny)
library(infoDecompuTE)
library(knitr)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Application title
  titlePanel("infoDecompuTE Shiny App"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "file1",
        "Upload the design as a CSV File",
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv")
      ),
      checkboxInput("checkbox", label = "Use example design", value = FALSE),
      radioButtons(
        "phase",
        label = h3("Single-phase/Two-phase"),
        choices = list("Single-phase" = 1, "Two-phase" = 2),
        selected = 2
      ),
      uiOutput("expStr"),
      downloadButton('downloadOutput', 'Download Table'),
      downloadButton('downloadLatexCode', 'Download Latex Code'),
      downloadButton('latexTable', 'Download Latex Table')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(
      tabPanel("Theoretical ANOVA table", 
               verbatimTextOutput("anovaTable"), 
               uiOutput("latexUI")),
      tabPanel("Design Data Frame", tableOutput('desDf'))
    ))
  )
))
