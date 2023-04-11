
library(tidyverse)

palette_color_factor <- c("Fresh" = "#66C2A5", "Fixed" ="#FC8D62", "Cryo" ="#8DA0CB")


### HEYN LAB FACS ###

FACSdata <- read.csv(here::here("05_experimental_results/CNAG_FACS/FACS_labelling_PBMC_cells.csv"), sep = ",", check.names = F)

FACSdata <- gather(FACSdata, key = "Celltypes", value = "Percentage", `B cells`:`Other cells`)
FACSdata$Condition <- factor(x = FACSdata$Condition,
                                 levels = c("Cryo", "Fixed"))

ggFACS <- ggplot(FACSdata, aes_string(x = "Condition", y = "Percentage", fill = "Condition", colour="Condition")) + 
  geom_boxplot(alpha = 0.5) + 
  geom_point(aes(group=Sample, shape=as.factor(Sample)), size =3, position = position_dodge(0.2)) +
  geom_line(aes(group = Sample), position = position_dodge(0.2), color = "darkgrey") +
  facet_wrap(~Celltypes, ncol = 6, strip.position = "bottom") + # scales = "free_y",
  scale_fill_manual(values = palette_color_factor) +
  scale_color_manual(values = palette_color_factor) +
  ggpubr::stat_compare_means(label = "p.format", paired = TRUE, label.y = 48, label.x.npc = "center", size=4.5) +
  labs(x="") +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16)) 
ggFACS

ggsave("FACS_PBMC_3replicates_plot.svg",
       plot = ggFACS,
       width = 15, height = 5,
       dpi = 300)



FACS_general <- read.csv(here::here("05_experimental_results/CNAG_FACS/FACS_labelling_PBMC_cells_general.csv"), sep = ",", check.names = F)

FACS_general <- gather(FACS_general, key = "Celltypes", value = "Percentage", `B cells`:`Other cells`)
FACS_general$Condition <- factor(x = FACS_general$Condition,
                             levels = c("Cryo", "Fixed"))

ggFACS_general <- ggplot(FACS_general, aes_string(x = "Condition", y = "Percentage", fill = "Condition", colour="Condition")) + 
  geom_boxplot(alpha = 0.5) + 
  geom_point(aes(group=Sample, shape=as.factor(Sample)), size =3, position = position_dodge(0.2)) +
  geom_line(aes(group = Sample), position = position_dodge(0.2), color = "darkgrey") +
  facet_wrap(~Celltypes, ncol = 6, strip.position = "bottom") + # scales = "free_y",
  scale_fill_manual(values = palette_color_factor) +
  scale_color_manual(values = palette_color_factor) +
  ggpubr::stat_compare_means(label = "p.format", paired = TRUE, label.y = 0, label.x.npc = 0.4, size=4.5) +
  labs(x="") +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16)) 
ggFACS_general

ggsave("FACS_PBMC_3replicates_general_plot.svg",
       plot = ggFACS_general,
       width = 8, height = 4,
       dpi = 300)



FACS_CD3 <- read.csv(here::here("05_experimental_results/CNAG_FACS/FACS_labelling_PBMC_cells_CD3.csv"), sep = ",", check.names = F)

FACS_CD3 <- gather(FACS_CD3, key = "Celltypes", value = "Percentage", `CD4+CD8+ T cells`:`CD4-CD8- T cells`)
FACS_CD3$Condition <- factor(x = FACS_CD3$Condition,
                                                   levels = c("Cryo", "Fixed"))

ggFACS_CD3 <- ggplot(FACS_CD3, aes_string(x = "Condition", y = "Percentage", fill = "Condition", colour="Condition")) + 
  geom_boxplot(alpha = 0.5) + 
  geom_point(aes(group=Sample, shape=as.factor(Sample)), size =3, position = position_dodge(0.2)) +
  geom_line(aes(group = Sample), position = position_dodge(0.2), color = "darkgrey") +
  facet_wrap(~Celltypes, ncol = 6, strip.position = "bottom") + #, scales = "free_y") +
  scale_fill_manual(values = palette_color_factor) +
  scale_color_manual(values = palette_color_factor) +
  ggpubr::stat_compare_means(label = "p.format", paired = TRUE, label.y = 70, label.x.npc = 0.4, size=4.5) +
  labs(x="") +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16)) 
ggFACS_CD3

ggsave("FACS_PBMC_3replicates_CD3_plot.svg",
       plot = ggFACS_CD3,
       width = 10, height = 4,
       dpi = 300)


### HEYN LAB MFI ###

MFIdata <- read.csv(here::here("05_experimental_results/CNAG_FACS/FACS_HumanPBMC_triplicate.csv"), sep = ",", check.names = F)
MFIdata <- gather(MFIdata, key = "Celltypes", value = "MFI", `CD3+`:`CD19+`)
MFIdata$Condition <- factor(x = MFIdata$Condition, levels = c("Cryo", "Fixed"))

ggMFI <- ggplot(MFIdata, aes_string(x = "Condition", y = "MFI", fill = "Condition", colour="Condition")) + 
  geom_bar(stat = "summary" , alpha = 0.5, aes(fill = Condition, colour=Condition)) + 
  geom_point(aes(group=Sample, shape=as.factor(Sample)), size =3, position = position_dodge(0.2)) +
  geom_line(aes(group = Sample), position = position_dodge(0.2), color = "darkgrey") +
  facet_wrap(~Celltypes, ncol = 6, strip.position = "bottom") + #, scales = "free_y") +
  scale_fill_manual(values = palette_color_factor) +
  scale_color_manual(values = palette_color_factor) +
  ggpubr::stat_compare_means(label = "p.format", paired = TRUE, label.y = 35000, label.x.npc = 0.4, size=4.5) +
  labs(x="") +
  theme_classic() +
  theme(legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16)) 
ggMFI

ggsave("MFI_PBMC_3replicates_plot.svg",
       plot = ggMFI,
       width = 10, height = 4,
       dpi = 300)



### PARISH LAB FACS ###

FACSdata <- read.csv(here::here("05_experimental_results/CNAG_FACS/FACS_labelling_PBMC_Parish_cells.csv"), sep = ",", check.names = F)

FACSdata <- gather(FACSdata, key = "Celltypes", value = "Percentage", `CD4+CD8+ T cells`:`Non T cells`)
FACSdata$Condition <- factor(x = FACSdata$Condition,
                             levels = c("Fresh", "Fixed"))
FACSdata$Celltypes <- factor(x = FACSdata$Celltypes,
                             levels = c("Non T cells", "CD4-CD8- T cells", "CD4-CD8+ T cells", "CD4+CD8- T cells", "CD4+CD8+ T cells" ))
palette_color_cells <- c("CD4+CD8+ T cells" = "#184e77",
                         "CD4+CD8- T cells" = "#1a759f",
                         "CD4-CD8+ T cells" = "#34a0a4",
                         "CD4-CD8- T cells" = "#76c893",
                         "Non T cells" = "#d9ed92")

ggFACS_p <- ggplot(FACSdata, 
                   aes(x = Condition, y = Percentage, fill = Celltypes)) + 
  geom_bar(stat="identity", alpha = 0.75) + 
  scale_fill_manual(values = palette_color_cells) +
  labs(x="", y = "Percentage of cells") +
  #geom_text(aes(label=paste0(sprintf("%.f", Percentage),"%")), position=position_stack(vjust=0.5)) +
  theme_classic() +
  coord_flip() +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size=16),
        axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5),
        strip.text = element_text(size=18)) 
ggFACS_p

ggsave("FACS_PBMC_Parish_plot.svg",
       plot = ggFACS_p,
       width = 8, height = 4,
       dpi = 300)

