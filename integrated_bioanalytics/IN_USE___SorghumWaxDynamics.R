######## SHINY APP



#################################
###### Sorghum Wax Dynamics #####
## Mady Larson and Lucas Busta ##
#################################

#### Setup ####
  #Install these packages, they are needed for image analysis but don't come with phylochemistry (only need to do this once)
    # install.packages("exifr")
    # remotes::install_github("brianwdavis/quadrangle"), INSTALL_opts = "--no-multiarch")
    # devtools::install_github("hrbrmstr/qrencoder")

#### Create QR codes ####

  ## Read in the Google docs spreadsheet

    # plant_data <- read_sheet("https://docs.google.com/spreadsheets/d/1euesYJkds5e3idYqZoXEpoY-9chUUMv8mr8LMVtyzKU/edit?usp=sharing")

  ## Create QR codes for the plants

  # plot_list <- list()
  # for (i in 1:dim(plant_data)[1]) {
  #   metadata <- paste(plant_data$Plant_Number[i], plant_data$Experiment[i], plant_data$Planting_Date[i], plant_data$Treatment[i], sep = "_")
  #   
  #   png(paste0("/Users/bust0037/Desktop/mady/qrcodes/", metadata, ".png"))
  #   image(
  #     1 - qrencoder::qrencode(metadata), 
  #     asp = 1, xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1),
  #     col = c("black", "white")
  #   )
  #   dev.off()
  #   
  #   plot_list[[i]] <- cowplot::plot_grid(
  #     ggpubr::text_grob(metadata, size = 5),
  #     grid::rasterGrob(png::readPNG(paste0("/Users/bust0037/Desktop/mady/qrcodes/", metadata, ".png"))),
  #     ncol = 1, rel_heights = c(1,30)
  #   )
  #   
  # }

  ## Create print page for 20 images
  
  # do.call("grid.arrange", c(plot_list))
  
  ## then save as pdf (letter dimensions), print.

#### Analyze the photos ####

files <- googledrive::drive_ls(path = "https://drive.google.com/file/d/1S3zRFYFiCtf2Q5DaEVn98lgC9VpIXjDD/view?usp=sharing")$id
metadata_total <- list()
for (image_number in 1:length(files)) {
#### Import Image ####
  cat("Importing image: ")
  cat(image_number)
  cat(" of ")
  cat(length(files))
  cat("\n")
  
  #### Path to image

    googledrive::drive_download(
      # file = files[image_number],
      file = "https://drive.google.com/file/d/1uUfZSMFpjaYJc4ztvbWE67MXzStbo5pP/view?usp=sharing",
      path = "/Users/bust0037/Desktop/img1.JPEG",
      overwrite = TRUE
    )

    path_to_image <- "/Users/bust0037/Desktop/img1.JPEG"
  
  #### Import image to analyze
  
    image <- imager::load.image(path_to_image)
    image <- as_tibble(as.data.frame(as.cimg(image)))
    image$y <- -as.numeric(image$y)
    image$y <- image$y + -min(image$y)
    y <- image$y
    x <- image$x

#### Rotation analysis ####
    
    cat("Running rotation analysis...\n")
    rotation_analysis <- list()
    degrees_to_test <- seq(-4,4,0.5)
    for(i in 1:length(degrees_to_test)) {
      
      temp_image <- image
      
      out <- rotate_coord(
        x = temp_image[,c(1)][[1]],
        y = temp_image[,c(2)][[1]],
        center = c((max(temp_image[,c(1)])/2),(max(temp_image[,c(2)])/2)),
        angle = degrees_to_test[i]
      )
      
      temp_image$x <- round(out[,1], 0)
      temp_image$y <- round(out[,2], 0)
      
      center_x <- mean(temp_image$x)
      
      temp_image %>%
        group_by(x) %>%
        summarize(total_value = sum(value)/length(value)) %>%
        filter(x < center_x) -> output
      
      out1 <- min(output$total_value)
      left_boundary_x <- output$x[which.min(output$total_value)]
      
      temp_image %>%
        group_by(x) %>%
        dplyr::summarize(total_value = sum(value)/length(value)) %>%
        filter(x > center_x) -> output
      
      out2 <- min(output$total_value)
      right_boundary_x <- output$x[which.min(output$total_value)]
    
      rotation_analysis[[i]] <- data.frame(
        degree = degrees_to_test[i],
        out1 = out1,
        out2 = out2
      )

    }
    
    rotation_analysis <- do.call(rbind, rotation_analysis)
    
    ggplot() +
      geom_line(data = rotation_analysis, aes(x = degree, y = out1)) +
      geom_line(data = rotation_analysis, aes(x = degree, y = out2)) +
      geom_line(data = rotation_analysis, aes(x = degree, y = out1+out2))
  
    degrees_to_rotate <- rotation_analysis$degree[which.min(rotation_analysis$out1 + rotation_analysis$out2)]
    
  #### Apply rotation
  
    out <- rotate_coord(
      x = image[,c(1)][[1]],
      y = image[,c(2)][[1]],
      center = c((max(image[,c(1)])/2),(max(image[,c(2)])/2)),
      angle = degrees_to_rotate
    )
    
    image$x <- round(out[,1], 0)
    image$y <- round(out[,2], 0)
    
  # Cut off some bits to make the image square on top and bottom after rotation
    
    rotation_offset_1 <- min(filter(image, x == image$x[which.max(image$y)])$y)
    rotation_offset_2 <- max(filter(image, x == image$x[which.min(image$y)])$y)
    
    if (rotation_offset_1 > rotation_offset_2) {
      image <- filter(image, y < rotation_offset_1 & y > rotation_offset_2)
    } else {
      image <- filter(image, y < rotation_offset_2 & y > rotation_offset_1)
    }
      
    
    plot <- ggplot() +
      geom_tile(data = image, aes(x = x, y = y, fill = value)) +
      facet_grid(.~cc)

    png("/Users/bust0037/Desktop/img1_analyzed.png", units = "in", height = 3, width = 9, res = 300)
      plot
    dev.off()

  #### Reflect image if necessary

    # image$y <- abs(image$y-max(image$y))

#### Find the left and right boundaries of the template ####
  cat("Identifying image features...\n")
  center_x <- mean(image$x)
  
  image %>%
    group_by(x) %>%
    summarize(total_value = sum(value)/length(value)) %>%
    filter(x < center_x) -> output
  
  left_boundary_x <- output$x[which.min(output$total_value)]
  
  image %>%
    group_by(x) %>%
    dplyr::summarize(total_value = sum(value)/length(value)) %>%
    filter(x > center_x) -> output
  
  right_boundary_x <- output$x[which.min(output$total_value)]
  
  center_x <- round(mean(c(right_boundary_x, left_boundary_x)))

#### Find left and right of sample box ####

  image %>%
    group_by(x) %>%
    dplyr::summarize(total_value = sum(value)) %>%
    filter(x < center_x) %>%
    filter(x > left_boundary_x*1.5) -> output
  
  left_sample_boundary_x <- output$x[which.min(output$total_value)]
  
  image %>%
    group_by(x) %>%
    dplyr::summarize(total_value = sum(value)) %>%
    filter(x > center_x) %>%
    filter(x < right_boundary_x*0.9) -> output
  
  right_sample_boundary_x <- output$x[which.min(output$total_value)]

#### Find top and bottom of sample box ####  

  center_y <- mean(image$y)

  image %>%
    group_by(y) %>%
    dplyr::summarize(total_value = sum(value)) -> output
    
    y_sample_boundary_1 <- output$y[which.min(output$total_value)]
  
  image %>%
    group_by(y) %>%
    summarize(total_value = sum(value)) %>%
    filter(y < y_sample_boundary_1-100 | y > y_sample_boundary_1+100) -> output
  
    y_sample_boundary_2 <- output$y[which.min(output$total_value)]
    
  image %>%
      group_by(y) %>%
      summarize(total_value = sum(value)) %>%
      filter(y < y_sample_boundary_1-100 | y > y_sample_boundary_1+100) %>%
      filter(y < y_sample_boundary_2-100 | y > y_sample_boundary_2+100) -> output
    
  y_sample_boundary_3 <- output$y[which.min(output$total_value)]
  
  y_boundaries <- c(y_sample_boundary_1, y_sample_boundary_2, y_sample_boundary_3)
  top_most_boundary_y <- max(y_boundaries)
  bottom_sample_boundary_y <- min(y_boundaries)
  top_sample_boundary_y <- y_boundaries[!y_boundaries %in% c(top_most_boundary_y, bottom_sample_boundary_y)]
  
#### Find RGB boxes ####
    
    box_buffer <- (top_sample_boundary_y - bottom_sample_boundary_y)/20
    
    image %>% 
      filter(cc == 1) %>%
      filter(y < top_most_boundary_y) %>%
      filter(y > top_sample_boundary_y + box_buffer) %>%
      group_by(y) %>%
      summarize(total_value = sum(value)) -> output
    
    R_bottom_boundary_y <- round(min(output$y[output$total_value > 1000]) + box_buffer)
    R_top_boundary_y <- round(max(output$y[output$total_value > 1000]) - box_buffer)

    image %>% 
      filter(cc == 2) %>%
      filter(y < top_most_boundary_y) %>%
      filter(y > top_sample_boundary_y + box_buffer) %>%
      group_by(y) %>%
      summarize(total_value = sum(value)) -> output
    
    G_bottom_boundary_y <- round(min(output$y[output$total_value > 950]) + box_buffer)
    G_top_boundary_y <- round(max(output$y[output$total_value > 950]) - box_buffer)
    
    image %>% 
      filter(cc == 3) %>%
      filter(y < top_most_boundary_y) %>%
      filter(y > top_sample_boundary_y + box_buffer) %>%
      group_by(y) %>%
      summarize(total_value = sum(value)) -> output
    
    B_bottom_boundary_y <- round(min(output$y[output$total_value > 320]) + box_buffer)
    B_top_boundary_y <- round(max(output$y[output$total_value > 320]) - box_buffer)
    
    plot(output$total_value)
    
    
#### asf ####
    
    # image %>% 
    #   filter(x == left_sample_boundary_x) %>%
    #   filter(value < 0.1) -> output
    # 
    # y_boundary_top <- max(output$y)
    # 
    # image %>%
    #   group_by(y) %>%
    #   filter(y < y_boundary_top*0.99) %>%
    #   dplyr::summarize(total_value = sum(value)) -> output
    # 
    # y_boundary_bottom <- which.min(output$total_value)

#### Sample box buffers and standard scale ####

  horizontal_buffer <- round(0.04*(left_sample_boundary_x + (right_sample_boundary_x-left_sample_boundary_x)/2))
  vertical_buffer <- round(0.015*(bottom_sample_boundary_y + (top_sample_boundary_y-bottom_sample_boundary_y)/2))
  
  image %>%
    filter(x > left_boundary_x + horizontal_buffer) %>%
    filter(x < right_boundary_x - horizontal_buffer) %>%
    filter(y > top_sample_boundary_y + vertical_buffer) %>%
    filter(y < B_bottom_boundary_y - vertical_buffer) -> image_standard_scale
  
  image_standard_scale %>%
    group_by(x, cc) %>%
    summarize(mode = mean(value)) -> image_standard_scale

  # image %>%
  #   filter(x == standard_x_right) %>%
  #   filter(y > standard_y_bottom) %>%
  #   filter(y < standard_y_top) %>%
  #   filter(value > 0) %>%
  #   select(value) -> max_std_values
  # max_std_value <- mean(max_std_values$value)
  # 
  # image %>%
  #   filter(x == standard_x_left) %>%
  #   filter(y > standard_y_bottom) %>%
  #   filter(y < standard_y_top) %>%
  #   filter(value > 0) %>%
  #   select(value) -> min_std_values
  # min_std_value <- mean(min_std_values$value)
  # 
  # slope <- (max_std_value-min_std_value)/(standard_x_right-standard_x_left)
  # intercept <- min_std_value-(slope*standard_x_left)

## Find average value in sample box

  image %>%
    filter(x > left_sample_boundary_x + horizontal_buffer) %>%
    filter(x < right_sample_boundary_x - horizontal_buffer) %>%
    filter(y > bottom_sample_boundary_y + vertical_buffer) %>%
    filter(y < top_sample_boundary_y - vertical_buffer) -> image_sample
  
  image_sample %>%
    ggplot(aes(x = value)) + 
    geom_histogram() +
    facet_grid(.~cc)

  image_sample %>%
    group_by(cc) %>%
    dplyr::summarize(
      mode = mode(value, ignore_zero = TRUE),
      sd = sd(value)
    ) -> sample_output
    
  sample_output %>%
    filter(cc == 3) -> output

  x_coordinate_of_sample_mode <- filter(image_standard_scale, cc == 3)$x[which.min(abs(image_standard_scale$mode - output$mode))]
  
#### Summarize color data and metadata ####
  cat("Extracting metadata...\n")
  metadata <- data.frame(
    cc = image_standard_scale$cc,
    standard_scale_x = image_standard_scale$x,
    standard_scale_mode = image_standard_scale$mode
  )
  
  reference_R <- image %>%
    filter(y < R_top_boundary_y - vertical_buffer) %>%
    filter(y > R_bottom_boundary_y + vertical_buffer) %>%
    filter(x > left_sample_boundary_x + horizontal_buffer) %>%
    filter(x < right_sample_boundary_x - horizontal_buffer) %>%
    group_by(cc) %>%
    dplyr::summarize(
      sample_mode = mode(value, ignore_zero = TRUE)
    ) %>%
    filter(cc == 1) %>% select(sample_mode) %>% as.numeric()
  
  reference_G <- image %>%
    filter(y < G_top_boundary_y - vertical_buffer) %>%
    filter(y > G_bottom_boundary_y + vertical_buffer) %>%
    filter(x > left_sample_boundary_x + horizontal_buffer) %>%
    filter(x < right_sample_boundary_x - horizontal_buffer) %>%
    group_by(cc) %>%
    dplyr::summarize(
      sample_mode = mode(value, ignore_zero = TRUE)
    ) %>%
    filter(cc == 2) %>% select(sample_mode) %>% as.numeric()
  
  reference_B <- image %>%
    filter(y < B_top_boundary_y - vertical_buffer) %>%
    filter(y > B_bottom_boundary_y + vertical_buffer) %>%
    filter(x > left_sample_boundary_x + horizontal_buffer) %>%
    filter(x < right_sample_boundary_x - horizontal_buffer) %>%
    group_by(cc) %>%
    dplyr::summarize(
      sample_mode = mode(value, ignore_zero = TRUE)
    ) %>%
    filter(cc == 3) %>% select(sample_mode) %>% as.numeric()
  
  metadata$sample_mode <- rep(sample_output$mode, dim(metadata)[1] / dim(sample_output)[1])
  metadata$reference_mode <- rep(c(reference_R,reference_G,reference_B), dim(metadata)[1] / dim(sample_output)[1])
  
  
  #### Read QR code
  
  qr_location_data <- quadrangle::qr_scan_cpp(magick::image_read(path_to_image), lighten = TRUE, darken = TRUE)
  metadata$name <- quadrangle::qr_scan_js_from_corners(magick::image_read(path_to_image), qr_location_data$points)$data
  metadata$time_stamp <- exifr::read_exif(path_to_image)$SubSecCreateDate

  metadata_total[[image_number]] <- metadata
  
#### Plot ####
  cat("\nCreating output image...\n\n")
  plot <- ggplot() +
    geom_tile(data = image, aes(x = x, y = y, fill = value)) +
    geom_hline(yintercept = top_most_boundary_y, color = "yellow") +
    geom_hline(yintercept = bottom_sample_boundary_y, color = "yellow") +
    geom_hline(yintercept = top_sample_boundary_y, color = "yellow") +
    geom_hline(yintercept = R_top_boundary_y, color = "red") +
    geom_hline(yintercept = R_bottom_boundary_y, color = "red") +
    geom_hline(yintercept = G_top_boundary_y, color = "green") +
    geom_hline(yintercept = G_bottom_boundary_y, color = "green") +
    geom_hline(yintercept = B_top_boundary_y, color = "blue") +
    geom_hline(yintercept = B_bottom_boundary_y, color = "blue") +
    # # geom_hline(yintercept = standard_y_top, color = "yellow") +
    # # geom_hline(yintercept = standard_y_bottom, color = "yellow") +
    # # geom_hline(yintercept = y_boundary_1-(max(y)*0.02)) +
    # # geom_hline(yintercept = y_boundary_1+(max(y)*0.02)) +
    geom_vline(aes(xintercept = left_boundary_x), color = "red") +
    geom_vline(aes(xintercept = right_boundary_x), color = "red") +
    # # geom_vline(aes(xintercept = standard_x_left), color = "yellow") +
    # # geom_vline(aes(xintercept = standard_x_right), color = "yellow") +
    geom_vline(xintercept = left_sample_boundary_x, color = "yellow") +
    geom_vline(xintercept = right_sample_boundary_x, color = "yellow") +
    geom_vline(aes(xintercept = x_coordinate_of_sample_mode), color = "pink") +
    # facet_grid(.~cc)
    theme_bw()
    
  png(paste0("/Users/bust0037/Desktop/img1_analyzed", image_number, ".png"), units = "in", height = 3, width = 3, res = 100)
    print(plot)
  dev.off()
  
}

metadata_total <- do.call(rbind, metadata_total)

#### Analyze the metadata ####

  metadata_total %>%
  group_by(time_stamp, cc) %>%
  mutate(
    sample_x_coord = standard_scale_x[which.min(abs(sample_mode - standard_scale_mode))]
  ) -> metadata_total

  metadata_total %>%
    ggplot() +
    geom_line(
      aes(x = standard_scale_x, y = standard_scale_mode, color = time_stamp),
      alpha = 0.5) +  
    geom_hline(aes(yintercept = reference_mode, color = time_stamp)) +
    geom_vline(aes(xintercept = sample_x_coord, color = time_stamp)) +
    facet_grid(cc~.) + 
    # scale_color_brewer(palette = "Set1") +
    theme_bw()
    
#### Analyze image metadata from mady out.R ####

  out <- read_csv("/Users/bust0037/Desktop/out.csv")
  
  out %>%
    filter(QR == "17_Light_Manipulation_z8-12-21_Single_Shade_Cloth") %>%
    filter(type == "sample") %>%
    group_by(Date) %>%
    summarize(
      R = mode(R),
      G = mode(G),
      B = mode(B)
    ) %>%
    pivot_longer(cols = 2:4, names_to = "color", values_to = "value") %>%
    ggplot(aes(x = Date, y = value)) + geom_point() + facet_grid(color~., scales = "free") + geom_smooth(method = "lm")
    
    
  unique(out$QR)
  

#### Old code #### ######

# output <- output[-c(round(seq(
# 	y_boundary_1-(max(y)*0.068),
# 	y_boundary_1+(max(y)*0.068),
# 	1
# ))),]

# y_boundary_bottom <- output$y[which.min(output$total_value)]
# y_boundary_top
# y_boundary_bottom

# if (y_boundary_top < y_boundary_bottom) {
# 	temp <- y_boundary_top
# 	y_boundary_top <- y_boundary_bottom
# 	y_boundary_bottom <- temp
# }

image %>%
  group_by(x) %>%
  dplyr::summarize(total_value = sum(value)) %>%
  filter(x > center_x) -> output

right_boundary_x <- image$x[which.min(output$total_value)+center_x]

# ggplot()+
# 	geom_point(data = output, aes(x = x, y = total_value))

# center_y = round(mean(image$y))
#    image_along_center_y <- dplyr::filter(image, y == center_y)

#    ggplot() +
# 	geom_point(data = image_along_center_y, aes(x = x, y = value))

#    ## Left boundary
#     x_of_50_darkest_spots <- filter(image_along_center_y, x < center_x)[order(filter(image_along_center_y, x < center_x)$value),][1:50,]$x
#     left_boundary_x <- mean(x_of_50_darkest_spots)

## Right boundary
# x_of_50_darkest_spots <- filter(image_along_center_y, x > center_x)[order(filter(image_along_center_y, x > center_x)$value),][1:50,]$x
# right_boundary_x <- mean(x_of_50_darkest_spots)

## Plot it

# ggplot() +
# 	geom_point(data = image_along_center_y, aes(x = x, y = value)) +
# 	geom_vline(aes(xintercept = left_boundary_x), color = "blue") +
# 	geom_vline(aes(xintercept = right_boundary_x), color = "red")    


plot <- ggplot() +
  geom_tile(data = image, aes(x = x, y = y, fill = value)) +
  geom_vline(aes(xintercept = left_boundary_x), color = "blue") +
  geom_vline(aes(xintercept = right_boundary_x), color = "red")

png("/Users/bust0037/Desktop/mady/analysis.png", units = "in", height = 3, width = 3, res = 300)
plot
dev.off()

## Find center_x

center_x <- round(mean(c(right_boundary_x, left_boundary_x)))
image_along_center_x <- dplyr::filter(image, x == center_x)

ggplot() +
  geom_point(data = image_along_center_x, aes(x = value, y = y))

y_of_50_darkest_spots <- filter(image_along_center_x, y > center_y)[order(filter(image_along_center_x, y > center_y)$value),][1:50,]$y
experimental_box_y <- mean(y_of_50_darkest_spots)

## Plot it

plot <- ggplot() +
  geom_tile(data = image, aes(x = x, y = y, fill = value)) +
  geom_vline(aes(xintercept = left_boundary_x), color = "blue") +
  geom_vline(aes(xintercept = right_boundary_x), color = "red") +
  geom_vline(aes(xintercept = center_x), color = "red") +
  geom_hline(aes(yintercept = experimental_box_y), color = "red")

png("/Users/bust0037/Desktop/mady/analysis.png", units = "in", height = 3, width = 3, res = 300)
plot
dev.off()

## Find the middle of the color palette

summary <- image %>%
  group_by(y) %>%
  dplyr::summarize(total_value = sum(value))



center_x = round(mean(image$x))
image_along_center_x <- dplyr::filter(image, x == center_x)


image_along_center_y[which.min(filter(image_along_center_y, x < center_x)$value),]
image_along_center_y[which.min(filter(image_along_center_y, x > center_x)$value),]

ma_center_y <- as.numeric(movingAverage(center_y$value, n = 100))
ma_center_y <- dropNA(ma_center_y)
x = as.numeric(seq(1:length(ma_center_y)))
y = ma_center_y

dat <- data.frame(
  x = x,
  y = y,
  deriv1_y = as.vector(stats::predict(pspline::sm.spline(x, y), x, 1))
)

deriv1_extreme_indices <- c(cumsum(rle(diff(dat$deriv1_y) > 0.002)$lengths)+1) # This from https://stackoverflow.com/questions/46029090/how-to-find-changing-points-in-a-datasetderiv1_extreme_indices <- c(cumsum(rle(abs(diff(deriv1_y)) > 0.002)$lengths)+1) # This from https://stackoverflow.com/questions/46029090/how-to-find-changing-points-in-a-dataset
deriv1_extreme_indices <- deriv1_extreme_indices[diff(deriv1_extreme_indices) > 5]

center_indices <- c(
  deriv1_extreme_indices[1] + ((deriv1_extreme_indices[2] - deriv1_extreme_indices[1])/2),
  deriv1_extreme_indices[3] + ((deriv1_extreme_indices[4] - deriv1_extreme_indices[3])/2)
)
bar_ys <- y[center_indices]

ggplot() +
  geom_point(data = dat, aes(x = x, y = y)) +
  geom_vline(aes(xintercept = deriv1_extreme_indices)) +
  geom_vline(aes(xintercept = center_indices), color = "red")

bars <- image[image$x %in% bar_ys,]

plot <- ggplot() +
  geom_tile(data = image, aes(x = x, y = y, fill = value)) +
  geom_hline(aes(yintercept = deriv1_extreme_indices)) +
  geom_hline(aes(yintercept = center_indices), color = "red")

png("/Users/bust0037/Desktop/mady/analysis.png", units = "in", height = 3, width = 3, res = 300)
plot
dev.off()

## Find the middle of the image and thus the y coord of the test bar

middle_y <- bar_ys[2] - bar_ys[1]
middle_y <- dplyr::filter(image, y == middle_y)

ma_y <- as.numeric(movingAverage(middle_y$value, n = 100))
ma_y <- dropNA(ma_y)
x = as.numeric(seq(1:length(ma_y)))
y = ma_y

dat <- data.frame(
  x = as.numeric(seq(1:length(ma_y))),
  y = ma_y,
  deriv1_y = as.vector(stats::predict(pspline::sm.spline(x, y), x, 1))
)

deriv1_extreme_indices <- c(cumsum(rle(diff(dat$deriv1_y) > 0.002)$lengths)+1) # This from https://stackoverflow.com/questions/46029090/how-to-find-changing-points-in-a-datasetderiv1_extreme_indices <- c(cumsum(rle(abs(diff(deriv1_y)) > 0.002)$lengths)+1) # This from https://stackoverflow.com/questions/46029090/how-to-find-changing-points-in-a-dataset
deriv1_extreme_indices <- deriv1_extreme_indices[diff(deriv1_extreme_indices) > 5]

center_indices <- c(
  deriv1_extreme_indices[1] + ((deriv1_extreme_indices[2] - deriv1_extreme_indices[1])/2),
  deriv1_extreme_indices[3] + ((deriv1_extreme_indices[4] - deriv1_extreme_indices[3])/2)
)
bar_ys <- x[center_indices]

ggplot() +
  geom_point(data = dat, aes(x = x, y = y)) +
  geom_vline(aes(xintercept = deriv1_extreme_indices)) +
  geom_vline(aes(xintercept = center_indices), color = "red")

bars <- image[image$y %in% bar_ys,]


## Extract QR code information and time stamp

qr_location_data <- quadrangle::qr_scan_cpp(magick::image_read(path_to_image), lighten = TRUE, darken = TRUE)
data.frame(
  name = quadrangle::qr_scan_js_from_corners(magick::image_read(path_to_image), qr_location_data$points)$data,
  time_stamp = exifr::read_exif(path_to_image)$SubSecCreateDate
)

qr_location_data$points$y <- magick::image_info(magick::image_read(path_to_image))$height - qr_location_data$points$y
magick::image_ggplot(magick::image_read(path_to_image)) + theme_classic() + geom_point(data = qr_location_data$points, aes(x = x, y = y))


