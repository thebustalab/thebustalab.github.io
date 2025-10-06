library(shiny)
library(ggplot2)
library(tibble)
library(dplyr)

empty_plot <- function(title) {
  ggplot() +
    theme_void(base_size = 14) +
    labs(title = title)
}

ui <- fluidPage(
  titlePanel("Island Water Wells"),
  fluidRow(
    column(
      width = 4,
      h3("North Aquifer"),
      numericInput(
        "north_mean",
        label = "Mean potassium (ppm)",
        value = 600,
        min = 0,
        max = 1000,
        step = 10
      ),
      numericInput(
        "north_sd",
        label = "SD (ppm)",
        value = 80,
        min = 1,
        max = 500,
        step = 1
      ),
      plotOutput("north_plot", height = "400px")
    ),
    column(
      width = 4,
      h3("South Aquifer"),
      numericInput(
        "south_mean",
        label = "Mean potassium (ppm)",
        value = 400,
        min = 0,
        max = 1000,
        step = 10
      ),
      numericInput(
        "south_sd",
        label = "SD (ppm)",
        value = 60,
        min = 1,
        max = 500,
        step = 1
      ),
      plotOutput("south_plot", height = "400px")
    ),
    column(
      width = 4,
      h3("Sample"),
      actionButton(
        "sample_button",
        label = "Sample!"
      ),
      numericInput(
        "asdf",
        label = "SD (ppm)",
        value = 80,
        min = 1,
        max = 500,
        step = 1
      ),
      plotOutput("sample_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  north_samples <- reactive({
    tibble(
      aquifer = "North Aquifer",
      value = rnorm(250, mean = input$north_mean, sd = input$north_sd)
    )
  })

  south_samples <- reactive({
    tibble(
      aquifer = "South Aquifer",
      value = rnorm(250, mean = input$south_mean, sd = input$south_sd)
    )
  })

  sampled_points <- reactiveVal(tibble(aquifer = character(), value = numeric()))

  observeEvent(input$sample_button, {
    sample_df <- bind_rows(
      north_samples() %>% slice_sample(n = 10, replace = FALSE),
      south_samples() %>% slice_sample(n = 10, replace = FALSE)
    )
    sampled_points(sample_df)
  }, ignoreNULL = FALSE)

  output$north_plot <- renderPlot({
    ggplot(north_samples(), aes(x = aquifer, y = value)) +
      geom_dotplot(
        binaxis = "y",
        stackdir = "center",
        stackratio = 1,
        dotsize = 0.8,
        fill = "maroon",
        color = "black"
      ) +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)"
      ) +
      theme_minimal(base_size = 14) +
      scale_y_continuous(limits = c(0,1000)) +
      theme(axis.text.x = element_blank())
  })

  output$south_plot <- renderPlot({
    ggplot(south_samples(), aes(x = aquifer, y = value)) +
      geom_dotplot(
        binaxis = "y",
        stackdir = "center",
        stackratio = 1,
        dotsize = 0.8,
        fill = "forestgreen",
        color = "black"
      ) +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)"
      ) +
      theme_minimal(base_size = 14) +
      scale_y_continuous(limits = c(0,1000)) +
      theme(axis.text.x = element_blank())
  })

  output$sample_plot <- renderPlot({
    sample_df <- sampled_points()
    req(nrow(sample_df) > 0)

    ggplot(sample_df, aes(x = aquifer, y = value, fill = aquifer)) +
      geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = 21, outlier.fill = "white") +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)",
        title = "Sampled wells (n = 10 per aquifer)"
      ) +
      scale_fill_manual(values = c("North Aquifer" = "#3182bd", "South Aquifer" = "#de2d26")) +
      theme_minimal(base_size = 14) +
      scale_y_continuous(limits = c(0,1000)) +
      theme(legend.position = "none")
  })
}

shinyApp(ui, server)
