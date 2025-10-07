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
        value = 60,
        min = -1000,
        max = 1000,
        step = 1
      ),
      numericInput(
        "north_sd",
        label = "SD (ppm)",
        value = 8,
        min = 0.1,
        max = 500,
        step = 0.1
      ),
      plotOutput("north_plot", height = "400px")
    ),
    column(
      width = 4,
      h3("South Aquifer"),
      numericInput(
        "south_mean",
        label = "Mean potassium (ppm)",
        value = 45,
        min = -1000,
        max = 1000,
        step = 1
      ),
      numericInput(
        "south_sd",
        label = "SD (ppm)",
        value = 6,
        min = 0.1,
        max = 500,
        step = 0.1
      ),
      plotOutput("south_plot", height = "400px")
    ),
    column(
      width = 4,
      h3("Sample"),
      actionButton("sample_button", "Sample!"),
      actionButton("rare_button", "Rare Sample"),
      plotOutput("sample_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  north_samples <- reactive({
    tibble(
      idx = seq_len(1000),
      aquifer = "North Aquifer",
      value = rnorm(1000, mean = input$north_mean, sd = input$north_sd)
    )
  })

  south_samples <- reactive({
    tibble(
      idx = seq_len(1000),
      aquifer = "South Aquifer",
      value = rnorm(1000, mean = input$south_mean, sd = input$south_sd)
    )
  })

  sampled_points <- reactiveVal(NULL)

  observeEvent(input$sample_button, {
    sample_df <- bind_rows(
      north_samples() %>% slice_sample(n = 10, replace = FALSE),
      south_samples() %>% slice_sample(n = 10, replace = FALSE)
    ) %>% mutate(mode = "Random")

    message("[Sample] Random draw created with ", nrow(sample_df), " wells")
    sampled_points(sample_df)
  }, ignoreNULL = FALSE)

  observeEvent(input$rare_button, {
    north_data <- north_samples()
    south_data <- south_samples()

    north_threshold <- quantile(north_data$value, probs = 0.8)
    south_threshold <- quantile(south_data$value, probs = 0.2)

    north_pool <- north_data %>% filter(value >= north_threshold)
    south_pool <- south_data %>% filter(value <= south_threshold)

    north_pick <- if (nrow(north_pool) > 0) {
      north_pool %>% slice_sample(n = min(10, nrow(north_pool)), replace = FALSE)
    } else {
      north_pool
    }

    south_pick <- if (nrow(south_pool) > 0) {
      south_pool %>% slice_sample(n = min(10, nrow(south_pool)), replace = FALSE)
    } else {
      south_pool
    }

    sample_df <- bind_rows(north_pick, south_pick)

    if (nrow(sample_df) > 0) {
      sample_df <- sample_df %>% mutate(mode = "Rare")
      message("[Sample] Rare draw created with ", nrow(sample_df), " wells")
      sampled_points(sample_df)
    } else {
      showNotification("Rare sample unavailable: adjust mean/SD to create extreme values.", type = "warning")
    }
  })

  output$north_plot <- renderPlot({
    sampled <- sampled_points()
    sample_idx <- if (is.null(sampled)) integer(0) else sampled$idx[sampled$aquifer == "North Aquifer"]

    data <- north_samples() %>%
      mutate(is_sample = idx %in% sample_idx)

    ggplot(data, aes(x = aquifer, y = value, fill = is_sample)) +
      geom_dotplot(
        binaxis = "y",
        stackdir = "center",
        stackratio = 0.6,
        dotsize = 0.8,
        binwidth = max(input$north_sd / 6, 0.5),
        color = NA
      ) +
      scale_fill_manual(values = c(`TRUE` = "#DAA520", `FALSE` = "#3182bd"), guide = "none") +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)"
      ) +
      coord_flip() +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_blank())
  })

  output$south_plot <- renderPlot({
    sampled <- sampled_points()
    sample_idx <- if (is.null(sampled)) integer(0) else sampled$idx[sampled$aquifer == "South Aquifer"]

    data <- south_samples() %>%
      mutate(is_sample = idx %in% sample_idx)

    ggplot(data, aes(x = aquifer, y = value, fill = is_sample)) +
      geom_dotplot(
        binaxis = "y",
        stackdir = "center",
        stackratio = 0.6,
        dotsize = 0.8,
        binwidth = max(input$south_sd / 6, 0.5),
        color = NA
      ) +
      scale_fill_manual(values = c(`TRUE` = "#DAA520", `FALSE` = "#de2d26"), guide = "none") +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)"
      ) +
      coord_flip() +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_blank())
  })

  output$sample_plot <- renderPlot({
    sample_df <- sampled_points()
    req(!is.null(sample_df))
    req(nrow(sample_df) > 0)

    counts_df <- as.data.frame(table(sample_df$aquifer), stringsAsFactors = FALSE)
    names(counts_df) <- c("aquifer", "count")

    title_suffix <- counts_df %>%
      mutate(label = sprintf("%s: %s", aquifer, count)) %>%
      pull(label) %>%
      paste(collapse = ", ")

    message("[Sample plot] Rendering with title suffix: ", title_suffix)

    ggplot(sample_df, aes(x = aquifer, y = value, fill = aquifer)) +
      geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = 21, outlier.fill = "white") +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)",
        title = paste0(sample_df$mode[1], " samples (", title_suffix, ")")
      ) +
      coord_flip() +
      scale_fill_manual(values = c("North Aquifer" = "#3182bd", "South Aquifer" = "#de2d26")) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })
}

shinyApp(ui, server)
