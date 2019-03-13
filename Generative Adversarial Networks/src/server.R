library(shiny)
library(ggplot2)
library(data.table)
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
    latent_vec <- paste(df[, 1], collapse = ",")
    system(paste("~/gan/run.sh ", latent_vec))
    myvec <- fread("~/gan/gan_output.csv", sep=",", header = FALSE) # (n + 1) * 1
    case_control <- myvec[1, length(myvec)]
    normalized_counts <- myvec[1, 1:length(myvec)-1]
    mymat <- fread("~/gan/gan_pca_mat.csv", sep=",", header = FALSE) # n * 2
    coors <- as.matrix(mymat) %*% as.matrix(normalized_counts[1:15])
    write.csv(coors, "~/gan/coordinates.csv", row.names = FALSE)
    coors
  })
  
  coors <- read.csv("~/gan/coordinates.csv")
  output$r <- renderPlot({
    ggplot(data = coors, aes(x=coors[1,1],y=coors[2,1], size=4))+
      geom_point()+
      coord_cartesian(xlim = c(-1500, 1500), ylim = c(-1500, 1500))
  })

  output$r2 <- renderPrint({text_value()})
  
}