
#===============
# Load libraries
#===============
library(ggridges)
library(tidyverse)
library(forcats)

Path <- "../results/Tables"

couleurs <- c("#FFFF00","#EEC591","#912CEE","#00FFFF","#87CEFF",
              "#006400","#7CFC00","#40E0D0","#828282","#000000")

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

#==================
##  RIDGELINE PLOT
#==================
Data %>%
    # Create Plot
    ggplot( aes(y = Population, x = Tajima.D,  fill = Population)) +
    geom_density_ridges_gradient(scale = 1) +
    scale_fill_manual(values = couleurs) + theme_minimal() +
    theme(
            legend.position="none",
            panel.spacing = unit(0.1, "lines"),
            axis.text = element_text(size = 10, color = 'black', face = "bold"),
            axis.title = element_text(size = 10, color = 'black', face = "bold")
        ) +
    xlab("Tajima's D") +
    ylab("")

#===========
# Save plot
#===========
ggsave("../results/Figures/Tajima_distribution.png")
