
## Load packages
library(tidyverse)
library(stringr)
library(glue)
library(rlang)

arguments <- commandArgs(trailingOnly = TRUE)
Path <- arguments[1]

# load New Gene Name file
Gene_name <- "../Data/metadata/PF_GeneName.tsv"
Gene_name <- read_tsv(Gene_name)

## Let's combined all results in on single file
## Single country population

## I can also use map_dfr to combine all dataframes into one

## Load files
tibble_files <- list.files(path = Path, pattern = ".xlsx", full.names = TRUE)

## Arrange files by date of sample collection
tibble_files <- tibble(tibble_files) %>% 
    mutate(prefix = as.numeric(str_extract(tibble_files, "[0-9]{1,}"))) %>% 
    arrange(prefix) %>% 
    pull(tibble_files)

names(tibble_files) <- str_replace(string = tibble_files,
                                   pattern = "../results/Tables/(.*)_tajima.xlsx",
                                   replacement = "\\1")

Data <- map_dfr(.x = tibble_files, .f = read_tsv, .id = "Population") %>% 
        separate(col = `rownames(ff)`, 
                into = c("Start", "End"),
                sep = " - ") %>% 
        rename(Chromosome = `geneFile$V1[i]`)


Data <- Data %>% 
    # Convert Columns type
    type_convert() %>% 
    
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len = max(Start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(Data, ., by=c("Chromosome" = "Chromosome")) %>%
    
    # Convert Columns type
    type_convert() %>% 
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Start) %>%
    mutate( BPcum = Start + tot) %>% 
    
    # Add gene name
    inner_join(Gene_name) %>% 
    select(-c(tot, Chromosome, Gene_ID, NewGeneName, Chromosome))

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead.

axis_set <- Data %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

ylim <- abs(floor(max(Data$Tajima.D))) + 1
ylimits <- c(floor(min(Data$Tajima.D)), ylim)

# Ready to make the manhattan plot using ggplot2:
tajima.manhattan <- function(data, xcolumn, ycolumn, color, axis_set, ylimits)
{
    data %>% ggplot(aes(x = {{xcolumn}}, 
                        y = {{ycolumn}}, 
                        color = as_factor({{color}})
                        )) + 
        # Show all points
        geom_point(alpha = 0.75) +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
        
        # custom X axis:
        scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
        scale_y_continuous(expand = c(0,0), limits = {{ylimits}}) +    # remove space between plot area and x axis
        scale_color_manual(values = rep(c("black", "red") , 
                                        unique(length(axis_set$Chr))
                                        )) +
        scale_size_continuous(range = c(0.5,3)) +
        labs(x = "Chromosomes", 
             y = "Tajima's D") + 
        
        # Custom the theme:
        theme_bw() +
        # theme_minimal() +
        theme( 
            legend.position = "none",
            axis.text = element_text(face = "bold"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title = element_text(face = "bold.italic")) +
        facet_wrap(~Population, ncol = 5)
}

tajima.manhattan(Data, BPcum, Tajima.D, Chr, axis_set, ylimits)
ggsave("../results/Figures/tajima.png", width = 30, height = 20, units = "cm")
