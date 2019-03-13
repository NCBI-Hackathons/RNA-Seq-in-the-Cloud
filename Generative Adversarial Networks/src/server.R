library(shiny)
library(ggplot2)
df <- data.frame(matrix(data = 0,ncol = 1, nrow = 200))
colnames(df) <- "value"

function(input, output, session) {
  
  output$ui <- renderUI({
    sliderInput("val", 
                       "Choose the input Value",
                       min = -3, 
                       max = 3, 
                       value = df[input$pos,1],
                      step = 0.01)
           
    
  })
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    df[input$pos,1] <<- input$val
    ggplot(df, aes(x=as.numeric(row.names(df)),y=value))+
      geom_bar(stat="identity")
  })
  
  text_value <- eventReactive(input$run_script, {
    runif(input$pos)
  })
  output$r <- renderPrint({text_value()})
  
}