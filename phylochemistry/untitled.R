source("")
library(viridis)

input_frame <- data.frame(
  Genus_species = "kalanchoe_species",
  GFF_in_path = "/Users/bust0037/Desktop/blocks_coords.gff",
  region_reach = 20000,
  concatenation_spacing = 10000,
  concatenate_by_species = FALSE,
  subset_mode = "none"
)

output <- readGFFs(GFFs, subset_mode = "none")

which_chr <- grep("Chr", output$seqnames)
length(unique(output[which_chr,]$seqnames))
which_scaff <- grep("Scaffold", output$seqnames)
length(unique(output[which_scaff,]$seqnames))
which_tig <- grep("tig", output$seqnames)
length(unique(output[which_tig,]$seqnames))
# which_other <- which(apply(cbind(
#   !seq(1,dim(output)[1],1) %in% which_chr,
#   !seq(1,dim(output)[1],1) %in% which_scaff
# ), 1, all))

# chr <- output[,]
# scaff <- output[,]

to_plot <- output[which_chr,]

to_plot <- to_plot[order(to_plot$start, decreasing = FALSE),]
to_plot %>%
  ggplot() +
    # geom_rect(aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = seqnames)) +
    geom_rect(aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = seqnames)) +
    scale_y_continuous(limits = c(-1,1)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    # facet_wrap(seqnames~.) +
    facet_grid(seqnames~.) +
    scale_fill_viridis(discrete=TRUE) +
    theme(
      axis.text.y = element_blank(),
      strip.text.y = element_text(angle = 0)
    )

pre <- readFasta("/Users/bust0037/Desktop/k_fed.contigs.fa")
post <- readFasta("/Users/bust0037/Desktop/k_fed.contigs.scaffolded.fa")

data <- rbind(
  data.frame(
    length = pre@ranges@width,
    file = "pre"
  ),
  data.frame(
    length = post@ranges@width,
    file = "post"
  )
)

library(treemapify)

ggplot(filter(data, file == "pre"), aes(area = length, fill = file)) +
  geom_treemap(color = "black")

ggplot(filter(data, file == "post"), aes(area = length, fill = file)) +
  geom_treemap(color = "black")

ggplot(data, aes(area = length, fill = file)) +
  geom_treemap(color = "black") +
  facet_grid(.~file)