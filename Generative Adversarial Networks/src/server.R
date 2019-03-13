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
    system("ls")
    myvec <- read.csv("/tmp/gan_visualization.csv", header=FALSE) # (n + 1) * 1
    case_control <- myvec[1,length(myvec)]
    normalized_counts <- myvec[1,1:length(myvec)-1]
    print(case_control)
    print(normalized_counts[,0:10])
    
    mymat <- read.csv("/tmp/gan_pca_mat.csv", header = FALSE) # n * 2
    coors <- as.matrix(normalized_counts) %*% t(as.matrix(mymat))
  })
  output$r <- renderPrint({text_value()})
  
}