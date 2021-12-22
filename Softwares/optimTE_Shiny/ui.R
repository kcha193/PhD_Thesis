#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(infoDecompuTE)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Optimal Phase 2 MudPIT-iTRAQ Design:"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      radioButtons("phase1Des", label = h3("Phase 1 Design"),
                   choices = list("CRD" = "CRD", "RCBD" = "RCBD", "BIBD" = "BIBD"), 
                   selected = "CRD"),
      selectInput("trt", label = h3("Treatment number:"), 
                  choices = 2:8, 
                  selected = 2),
      selectInput("bioRep", label = h3("Biological Replicate number:"), 
                  choices = 2:10, 
                  selected = 2),
      selectInput("cageRep", label = h3("Cage number:"), 
                  choices = 2:10, 
                  selected = 2),
      selectInput("tag", label = h3("Tag number:"), 
                  choices = c(4,8), 
                  selected = 4)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h2("Phase 1 Design:"),
      tableOutput("desTablePhase1"),
      h2("Optimal Phase 2 Design:"),
      tableOutput("desTable"),
      h2("Theatrical ANOVA Table:"),
      verbatimTextOutput("anovaTable")
    )
  )
))
