library(shiny)

shinyUI(fluidPage(
  titlePanel("2010-2018 WPD Arrest Data"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("year", label="Year", 
                  min=2010,max=2018,value=2010,sep="",animate=animationOptions(interval=500,loop=TRUE)),      
      radioButtons("data",label="Data:",c("Total Arrests","White Only Arrests","Black Only Arrests")),
      radioButtons("adj",label="Data Adjustment:",c("None","As a Percent of the Population","As a Percent of Total Arrests",
                                                    "Standardized Incidence Ratio","Poisson Regression")),
      br(),
      # div(img(src = "UNCWlogo.png", height = 70, width = 150), style="text-align: center;"),
      # br(),
      # div(p(a("Cape Fear Collective", 
      #         href = "https://capefearcollective.org/")), style="text-align: center;")
    ),
    mainPanel(
      textOutput("text"),
      plotOutput("map"))#,
      #tableOutput("table"))
  )))
