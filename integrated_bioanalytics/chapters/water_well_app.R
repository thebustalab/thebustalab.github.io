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
      plotOutput("sample_plot", height = "400px"),
      hr(),
      h3("Null Hypothesis View"),
      uiOutput("null_description"),
      plotOutput("null_distribution_plot", height = "300px"),
      plotOutput("null_extreme_samples_plot", height = "300px")
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

  observed_stats <- reactive({
    sample_df <- sampled_points()

    if (is.null(sample_df) || nrow(sample_df) == 0) {
      return(NULL)
    }

    counts <- table(sample_df$aquifer)
    required_groups <- c("North Aquifer", "South Aquifer")

    if (!all(required_groups %in% names(counts))) {
      return(NULL)
    }

    n_north <- counts[["North Aquifer"]]
    n_south <- counts[["South Aquifer"]]

    if (n_north == 0 || n_south == 0) {
      return(NULL)
    }

    mean_north <- mean(sample_df$value[sample_df$aquifer == "North Aquifer"])
    mean_south <- mean(sample_df$value[sample_df$aquifer == "South Aquifer"])

    list(
      mean_north = mean_north,
      mean_south = mean_south,
      diff = mean_north - mean_south,
      n_north = n_north,
      n_south = n_south
    )
  })

  null_results <- reactive({
    stats <- observed_stats()
    req(stats)

    pooled_mean <- mean(c(input$north_mean, input$south_mean))
    pooled_sd <- sqrt((input$north_sd^2 + input$south_sd^2) / 2)
    replicates <- 2000

    sim_list <- replicate(replicates, {
      north_draw <- if (stats$n_north > 0) rnorm(stats$n_north, mean = pooled_mean, sd = pooled_sd) else numeric(0)
      south_draw <- if (stats$n_south > 0) rnorm(stats$n_south, mean = pooled_mean, sd = pooled_sd) else numeric(0)

      list(
        diff = if (stats$n_north > 0 && stats$n_south > 0) mean(north_draw) - mean(south_draw) else NA_real_,
        north = north_draw,
        south = south_draw
      )
    }, simplify = FALSE)

    diff_values <- vapply(sim_list, function(x) x$diff, numeric(1))
    diffs_df <- tibble(
      replicate = seq_along(diff_values),
      diff = diff_values
    )

    extreme_mask <- abs(diff_values) >= abs(stats$diff)
    extreme_count <- sum(extreme_mask)
    p_value <- (extreme_count + 1) / (length(diff_values) + 1)

    extreme_ids <- which(extreme_mask)
    extreme_ids <- head(extreme_ids, 6)

    extreme_df <- if (length(extreme_ids) > 0) {
      bind_rows(lapply(extreme_ids, function(rep_id) {
        draw <- sim_list[[rep_id]]
        tibble(
          replicate = rep_id,
          aquifer = c(rep("North Aquifer", length(draw$north)), rep("South Aquifer", length(draw$south))),
          value = c(draw$north, draw$south),
          diff = draw$diff
        )
      }))
    } else {
      tibble()
    }

    list(
      diffs = diffs_df,
      p_value = p_value,
      observed_diff = stats$diff,
      extreme_samples = extreme_df,
      pooled_mean = pooled_mean,
      pooled_sd = pooled_sd,
      n_north = stats$n_north,
      n_south = stats$n_south
    )
  })

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

  output$null_description <- renderUI({
    stats <- observed_stats()
    results <- null_results()
    req(stats)
    req(results)

    mean_same <- sprintf("%.1f", results$pooled_mean)
    sd_same <- sprintf("%.1f", results$pooled_sd)
    diff_obs <- sprintf("%.2f", results$observed_diff)
    p_val <- sprintf("%.3f", results$p_value)

    HTML(paste0(
      "<p>Assuming both aquifers are really the same population (mean = ", mean_same,
      " ppm, sd = ", sd_same, " ppm), we repeatedly redraw <strong>",
      stats$n_north, "+", stats$n_south,
      " wells</strong>. The plots show how often those null samples create a difference in means as extreme as the observed Δ = ",
      diff_obs, " ppm. Only about ", p_val,
      " of those null worlds reach or exceed that gap — that’s your p-value story.</p>"
    ))
  })

  output$null_distribution_plot <- renderPlot({
    results <- null_results()
    diffs <- results$diffs

    if (is.null(diffs) || nrow(diffs) == 0) {
      return(empty_plot("Null sampling distribution unavailable"))
    }

    diffs <- diffs %>%
      mutate(extreme = abs(diff) >= abs(results$observed_diff))

    binwidth <- diff(range(diffs$diff)) / 30
    if (!is.finite(binwidth) || binwidth <= 0) {
      sd_est <- sd(diffs$diff, na.rm = TRUE)
      if (!is.finite(sd_est) || sd_est <= 0) {
        sd_est <- 0.5
      }
      binwidth <- sd_est / 5
    }
    if (!is.finite(binwidth) || binwidth <= 0) {
      binwidth <- 0.1
    }

    ggplot(diffs, aes(x = diff, fill = extreme)) +
      geom_histogram(color = "white", binwidth = binwidth, boundary = 0) +
      geom_vline(xintercept = 0, color = "#636363", linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = results$observed_diff, color = "#fd8d3c", linewidth = 0.9) +
      geom_vline(xintercept = -results$observed_diff, color = "#fd8d3c", linewidth = 0.9, linetype = "dotted") +
      scale_fill_manual(values = c(`TRUE` = "#feb24c", `FALSE` = "#c6dbef"), guide = "none") +
      labs(
        x = "Difference in sample means (North − South, ppm)",
        y = "Null samples",
        title = "If the aquifers are identical, how extreme is our sample?"
      ) +
      theme_minimal(base_size = 14)
  })

  output$null_extreme_samples_plot <- renderPlot({
    results <- null_results()
    samples <- results$extreme_samples

    if (is.null(samples) || nrow(samples) == 0) {
      return(empty_plot("No null draws reached the observed gap yet — try sampling again."))
    }

    samples <- samples %>%
      mutate(
        label = sprintf("Null draw %d\nΔ = %.2f ppm", replicate, diff)
      )

    ggplot(samples, aes(x = aquifer, y = value, color = aquifer)) +
      geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.9) +
      facet_wrap(~ label) +
      scale_color_manual(values = c("North Aquifer" = "#3182bd", "South Aquifer" = "#de2d26"), guide = "none") +
      labs(
        x = NULL,
        y = "Potassium abundance (ppm)",
        title = "What would null-world samples look like to beat our gap?"
      ) +
      theme_minimal(base_size = 14)
  })
}

shinyApp(ui, server)
