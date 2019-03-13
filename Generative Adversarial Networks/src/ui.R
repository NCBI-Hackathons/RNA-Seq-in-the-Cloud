fluidPage(
  titlePanel("Check The Lattice Points"),
  fluidRow(
    
    column(3, wellPanel(
      sliderInput("pos",
                   "Choose an Index to change",
                   min = 1,
                   max = 200,
                   value = 50)
      ),
      uiOutput("ui"),
      actionButton(inputId = "run_script",
                   label = "Refresh")
    ),
    

    column(7,
           plotOutput("distPlot"),
           verbatimTextOutput("r")
    )
  )
)