# source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")

ui <- fluidPage(
  rHandsontableOutput('hot')
)

server <- function(input, output, session) {
  
  # output$alignment_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
    data <- readAlignment("/Users/bust0037/Desktop/Kalanchoe_ITS_1.fa", type = "DNA")
    data <- pivot_wider(data, names_from = "position", values_from = "state")
    rownames(data) <- data$name
    data <- as.matrix(data[,-c(1)])
  
  output$hot <- renderRHandsontable({
    rhandsontable(data, 
                  # readOnly = TRUE
                  # width = 750, 
                  # height = 300
      ) %>%
      hot_cols(renderer = myrenderer)
  })
  
  myrenderer <- "function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.TextRenderer.apply(this, arguments);
                
                // Color the cells based on the value
                if (value == 'a') {
                    td.style.background = 'pink';
                } else if (value == 't') {
                    td.style.background = 'lightgreen';
                } else if (value == 'c') {
                    td.style.background = 'lightblue';
                } else if (value == 'g') {
                    td.style.background = 'lightyellow';
                }
  }"      
}

shinyApp(ui = ui, server = server)
