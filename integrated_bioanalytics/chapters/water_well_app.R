library(shiny)
library(tidyverse)

set.seed(123)

island <- tibble(
  theta = seq(0, 2 * pi, length.out = 301)
) %>% mutate(
  radius = 1 + 0.15 * sin(3 * theta) + 0.1 * sin(5 * theta),
  x = 2 * radius * cos(theta),
  y = 2 * radius * sin(theta)
)

make_wells <- function(aquifer_name, y_shift) {
  tibble(
    id = 1:100,
    theta = runif(100, 0, 2 * pi),
    radius = sqrt(runif(100)) * 0.65,
    z_score = rnorm(100)
  ) %>% mutate(
    x = 2 * radius * cos(theta),
    y = 2 * radius * sin(theta) + y_shift,
    aquifer = aquifer_name
  ) %>% select(aquifer, id, x, y, z_score)
}

well_locations <- bind_rows(
  make_wells("Northern aquifer", y_shift = 1.1),
  make_wells("Southern aquifer", y_shift = -1.1)
) %>% mutate(well_key = paste(aquifer, id, sep = "_"))

sample_default <- well_locations %>%
  group_by(aquifer) %>%
  slice_sample(n = 10) %>%
  ungroup() %>%
  pull(well_key)

ui <- fluidPage(
  titlePanel("Island Water Wells"),
  fluidRow(
    column(
      width = 4,
      h4("Aquifer distributions"),
      numericInput(
        "north_mean",
        label = "Northern aquifer mean (ppm)",
        value = 60,
        min = 0,
        max = 200,
        step = 1
      ),
      numericInput(
        "north_sd",
        label = "Northern aquifer sd (ppm)",
        value = 8,
        min = 0.1,
        max = 50,
        step = 0.1
      ),
      numericInput(
        "south_mean",
        label = "Southern aquifer mean (ppm)",
        value = 45,
        min = 0,
        max = 200,
        step = 1
      ),
      numericInput(
        "south_sd",
        label = "Southern aquifer sd (ppm)",
        value = 6,
        min = 0.1,
        max = 50,
        step = 0.1
      ),
      actionButton("sample", "Sample 10 From Each Aquifer")
    ),
    column(
      width = 4,
      plotOutput("abundance_plot", height = "520px")
    ),
    column(
      width = 4,
      plotOutput("island_plot", height = "520px")
    )
  )
)

server <- function(input, output, session) {
  sampled_keys <- reactiveVal(sample_default)

  wells_data <- reactive({
    north_mean <- input$north_mean
    north_sd <- input$north_sd
    south_mean <- input$south_mean
    south_sd <- input$south_sd

    well_locations %>%
      mutate(
        potassium_ppm = case_when(
          aquifer == "Northern aquifer" ~ north_mean + north_sd * z_score,
          aquifer == "Southern aquifer" ~ south_mean + south_sd * z_score
        )
      )
  })

  observeEvent(input$sample, {
    new_sample <- well_locations %>%
      group_by(aquifer) %>%
      slice_sample(n = 10) %>%
      ungroup() %>%
      pull(well_key)

    sampled_keys(new_sample)
  })

  sampled_data <- reactive({
    req(sampled_keys())

    wells_data() %>%
      filter(well_key %in% sampled_keys())
  })

  output$island_plot <- renderPlot({
    wells_df <- wells_data()
    sample_df <- sampled_data()

    ggplot() +
      geom_polygon(
        data = island,
        aes(x = x, y = y),
        fill = "#d9f0a3",
        color = "#66a61e",
        linewidth = 0.5
      ) +
      geom_point(
        data = wells_df,
        aes(x = x, y = y, color = aquifer),
        alpha = 0.5,
        size = 2
      ) +
      geom_point(
        data = sample_df,
        aes(x = x, y = y, color = aquifer),
        size = 3.6,
        stroke = 1.2
      ) +
      scale_color_manual(values = c("Northern aquifer" = "#3182bd", "Southern aquifer" = "#de2d26")) +
      coord_equal() +
      labs(
        x = NULL,
        y = NULL,
        subtitle = "Sampled wells are highlighted; adjust distributions on the left",
        color = "Aquifer"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        plot.subtitle = element_text(margin = margin(b = 10))
      )
  })

  output$abundance_plot <- renderPlot({
    wells_df <- wells_data()
    sample_df <- sampled_data()

    params <- tibble(
      aquifer = c("Northern aquifer", "Southern aquifer"),
      mean = c(input$north_mean, input$south_mean),
      sd = c(input$north_sd, input$south_sd)
    )

    hist_df <- wells_df %>% mutate(panel = "Distribution")

    density_df <- params %>%
      pmap_dfr(function(aquifer, mean, sd) {
        tibble(
          aquifer = aquifer,
          x = seq(mean - 4 * sd, mean + 4 * sd, length.out = 150),
          density = dnorm(x, mean = mean, sd = sd)
        )
      }) %>%
      mutate(panel = "Distribution")

    range_ppm <- range(wells_df$potassium_ppm, na.rm = TRUE)
    min_axis <- range_ppm[1]
    axis_span <- diff(range_ppm)
    if (!is.finite(axis_span) || axis_span == 0) {
      axis_span <- max(1, abs(min_axis) * 0.1 + 1)
    }
    label_pos <- min_axis - 0.05 * axis_span

    sample_vis <- sample_df %>%
      group_by(aquifer) %>%
      mutate(
        panel = "Sampled wells",
        sample_rank = as.integer(rank(potassium_ppm, ties.method = "first")),
        min_axis = min_axis,
        label_pos = label_pos,
        well_label = paste0("Well ", id)
      ) %>%
      ungroup() %>%
      arrange(aquifer, sample_rank)

    ggplot() +
      geom_histogram(
        data = hist_df,
        aes(x = potassium_ppm, fill = aquifer),
        bins = 20,
        color = "white",
        boundary = 0,
        alpha = 0.85
      ) +
      geom_point(
        data = density_df,
        aes(x = x, y = density, color = aquifer),
        size = 1.2,
        alpha = 0.9
      ) +
      geom_segment(
        data = sample_vis,
        aes(
          x = min_axis,
          xend = potassium_ppm,
          y = sample_rank,
          yend = sample_rank,
          color = aquifer
        ),
        linewidth = 0.6,
        alpha = 0.4
      ) +
      geom_point(
        data = sample_vis,
        aes(x = potassium_ppm, y = sample_rank, color = aquifer),
        size = 3
      ) +
      geom_text(
        data = sample_vis,
        aes(x = label_pos, y = sample_rank, label = well_label, color = aquifer),
        hjust = 1,
        size = 3.2,
        show.legend = FALSE
      ) +
      facet_grid(panel ~ aquifer, scales = "free_y", switch = "y") +
      scale_color_manual(values = c("Northern aquifer" = "#3182bd", "Southern aquifer" = "#de2d26")) +
      scale_fill_manual(values = c("Northern aquifer" = "#9ecae1", "Southern aquifer" = "#fcbba1")) +
      labs(
        x = "Potassium abundance (ppm)",
        y = NULL,
        title = "Aquifer potassium distributions and sampled wells"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        strip.text = element_text(face = "bold"),
        strip.placement = "outside",
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  })
}

shinyApp(ui, server)
