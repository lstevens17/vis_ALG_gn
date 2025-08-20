#!/usr/bin/env Rscript

library(optparse)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(ggplot2)
library(gtools)



option_list = list(
  make_option(c("-b", "--busco"), type="character", default=NULL,
              help="busco full_table.tsv file", metavar="file.tsv"),
  make_option(c("-n", "--nigon"), type="character", default="gene2Nigon_busco20200927.tsv.gz",
              help="busco id assignment to Nigons [default=%default]", metavar="file.tsv"),
  make_option(c("-w", "--windowSize"), type="integer", default=5e5,
              help="window size to bin the busco genes [default=%default]. Sequences shorter than twice this integer will not be shown in the plot", metavar="integer"),
  make_option(c("-m", "--minimumGenesPerSequence"), type="integer", default=15,
              help="sequences (contigs/scaffolds) with less than this number of busco genes will not be shown in the plot [default=%default]", metavar="integer"),
  make_option(c("-o", "--outPlot"), type="character", default="Nigons.png",
              help="output image [default=%default]. Should include one of the following extensions: eps, ps, tex, pdf, jpeg, tiff, png, bmp or svg", metavar="file"),
  make_option(c("--height"), type="integer", default=6,
              help="height of plot. Increase this value according to the number of ploted sequences [default=%default]", metavar="integer"),
  make_option(c("--width"), type="integer", default=6,
              help="width of plot [default=%default]", metavar="integer"),
  make_option(c("-s", "--species"), type="character", default="",
              help="Title to be italicized in the plot [default=%default]", metavar="Genus_species")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load data
nigonDict <- read_tsv(opt$nigon,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(opt$busco,
                  col_names = c("Busco_id", "Status", "Sequence",
                                "start", "end", "strand", "Score", "Length",
                                "OrthoDB_url", "Description"),
                  col_types = c("ccciicdicc"),
                  comment = "#")) %>% mutate(Sequence = sub(":.*", "", Sequence))

windwSize <- opt$windowSize
minimumGenesPerSequence <- opt$minimumGenesPerSequence
spName <- opt$species
if(grepl(".", opt$species)) {
  spName <- paste0("*", sub("_", " ", opt$species), "*")
}

# Specify Nigon colors
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Filter data
fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")

consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence, mxGpos > windwSize * 2)

# Plot
plNigon <- group_by(consUsco, Sequence) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), windwSize),
                                            labels = seq(windwSize, max(stPos), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  count(ints, nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
  ggplot(aes(fill=nigon, y=n, x=(ints-windwSize)/1e6)) +
  facet_grid(scaffold_f ~ .) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  scale_x_continuous() +
  scale_fill_manual(values = cols) +
  xlab("Position (Mb)") +
  ylab("BUSCO count (n)") +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme(strip.text.y.right = element_text(angle = 0), text = element_text(size=9),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "#ededed")
  )

ggsave(opt$out, plNigon, width = opt$width, height = opt$height)
ggsave(opt$out, plNigon, width = opt$width, height = opt$height)
