###Data analysis for Wild Orca SRKW sequences
###Sofia Kaiaua and Amy Van Cise
###Fall 2024

### Set up environment ---------------------------------------------------------
library(tidyverse)
library(kableExtra)
library(phyloseq)
library(Biostrings)
library(RColorBrewer)
library(PNWColors)
library(patchwork)
library(ggbeeswarm)
library(viridis)

theme_set(theme_minimal())

setwd("G:/My Drive/00 UW/00.5 W.A.D.E. lab resources/Intern Projects/SRKW diet metabarcoding/04 Data analysis")

#Load data
load("G:/My Drive/00 UW/00.5 W.A.D.E. lab resources/Intern Projects/SRKW diet metabarcoding/04 Data analysis/SRKW_WO_16Sdiet_dada2out.Rdata")

### Merge taxa to species, remove O. orca seqs ---------------------------------
ps.sp <- tax_glom(ps, taxrank="Species") %>% 
  subset_taxa(Species != "Orcinus orca") %>% 
  prune_samples(sample_sums(.) > 0, .)
nsamples(ps.sp)

### Convert to proportional, create subset data -----------
ps.sp_proportional <- transform_sample_counts(ps.sp, function(x) x / sum(x))

ps.sp_subset <- subset_samples(ps.sp_proportional,
                               grepl("002-024|002-002|002-014",
                                     sample_names(ps.sp_proportional)))

ps.sp_proportional <- subset_samples(ps.sp_proportional, 
                                     !sample_names(ps.sp_proportional) %in% 
                                       c("WADE-002-024-nc", "WADE-002-002", "WADE-002-014-nc"))


### Filter out taxa with <1% proportional reads --------------------------------
#must be at least 1% of diet in 1 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.sp_proportional, f1, A=1)
ps.sp_major <- prune_taxa(lowcount.filt, ps.sp_proportional)

#Relative abundance Stacked Barplot for Individual Samples 
species_prop_bar <- plot_bar(ps.sp_major, x = "Sample", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Diet Composition of Samples",
       x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(limits = c("Atheresthes stomias","Oncorhynchus keta","Oncorhynchus kisutch","Oncorhynchus mykiss", 
                                 "Oncorhynchus tshawytscha", "Ophiodon elongatus", "Sebastes ensifer",
                                 "Carcharodon carcharias", "Gadus chalcogrammus",
                                 "Microstomus pacificus", "Hippoglossus stenolepis",
                                 "Oncorhynchus nerka"),
                      labels = c("Arrowtooth flounder", "Chum salmon", "Coho salmon", "Steelhed salmon", 
                                 "Chinook salmon", "Lingcod","Swordspine rockfish",
                                 "Great white shark", "Alaska polluck",
                                 "Pacific Dover sole", "Pacific halibut",
                                 "Sockeye salmon"))

### Compare Cleaned and Not Cleaned (nc) samples -------------------------------

#Final DNA concentration
dna_con_subset <- samdf %>% 
  rownames_to_column("SampleID") %>% 
  filter(!is.na(Index_concentration)) %>% 
  mutate(PCR_clean = case_when(grepl("nc", SampleID)~FALSE,
                               TRUE~TRUE)) %>% 
  separate(SampleID, sep = "-", into = c("LAb","Proj","Samp", NA)) %>% 
  unite(LAb:Samp, col = "SampleID", sep = "-") %>% 
  ggplot(aes(x=PCR_clean, y = Index_concentration, color = SampleID)) +
  geom_violin(color = "grey50") +
  geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) +
  theme(text = element_text(size = 14), legend.position = "none") +
  scale_color_viridis(discrete = TRUE) +
  ylab("Final DNA concentration (ng/ul)")

#Read count
reads_subset <- track %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  mutate(SampleID = case_when(SampleID == "WADE-002-014-rep"~"WADE-002-014-nc",
                              TRUE~SampleID)) %>% 
  mutate(PCR_clean = case_when(grepl("nc", SampleID)~FALSE,
                               TRUE~TRUE)) %>% 
  separate(SampleID, sep = "-", into = c("LAb","Proj","Samp", NA)) %>% 
  unite(LAb:Samp, col = "SampleID", sep = "-") %>% 
  select(SampleID, PCR_clean, nonchim) %>% 
  add_row(SampleID = c("WADE-002-063","WADE-002-054","WADE-002-056",
                       "WADE-002-028","WADE-002-019"),
          PCR_clean = TRUE, nonchim = 0) %>%
  group_by(SampleID) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  filter(SampleID != "WADE-002-002") %>% 
  ggplot(aes(x=PCR_clean, y = nonchim, color = SampleID)) +
  geom_violin(color = "grey50") +
  geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) +
  theme(text = element_text(size = 14), legend.position = "none") +
  scale_color_viridis(discrete = TRUE) +
  ylab("Read count")+
  ylim(0,80000)

#bar plot
species_prop_subset <- plot_bar(ps.sp_subset, x = "Sample", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~FieldID, scales = "free_x", 
             labeller = labeller(FieldID = c("F21SEP15.05A"="024",
                                            "F22JUL28.05A"="002",
                                            "F22JUN23.01A"="014"))) +
  labs(x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(limits = c("Atheresthes stomias","Oncorhynchus keta","Oncorhynchus kisutch","Oncorhynchus mykiss", 
                                 "Oncorhynchus tshawytscha", "Ophiodon elongatus", "Sebastes ensifer",
                                 "Carcharodon carcharias", "Gadus chalcogrammus",
                                 "Microstomus pacificus", "Hippoglossus stenolepis",
                                 "Oncorhynchus nerka"),
                      labels = c("Arrowtooth flounder", "Chum salmon", "Coho salmon", "Steelhed salmon", 
                                 "Chinook salmon", "Lingcod","Swordspine rockfish",
                                 "Great white shark", "Alaska polluck",
                                 "Pacific Dover sole", "Pacific halibut",
                                 "Sockeye salmon")) +
  guides(color=guide_legend(ncol=2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

  
###Compare Diet in Samples in Summer & Non-Summer Months -----------------------
#Create Field ID Column
ps_summer <- subset_samples(ps.sp_major, 
                            grepl("JUN|JUL|AUG|SEP|Aug|Sep|SEA", 
                                  sample_data(ps.sp_proportional)$FieldID))
ps_other_months <- subset_samples(ps.sp_major, 
                                  !grepl("JUN|JUL|AUG|SEP|Aug|Sep|SEA", 
                                         sample_data(ps.sp_proportional)$FieldID))
sample_data(ps_summer)$FieldID_group <- "Summer Months"
sample_data(ps_other_months)$FieldID_group <- "Non-Summer Months"
ps_combined <- merge_phyloseq(ps_summer, ps_other_months)

#Create Stacked Bar Chart to Compare Seasonal Diet Patterns
season_bar <- plot_bar(ps_combined, x = "Sample", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~ FieldID_group, scales = "free_x", space = "free_x") +
  scale_fill_discrete(limits = c("Atheresthes stomias","Oncorhynchus keta",
                                 "Oncorhynchus kisutch","Oncorhynchus mykiss", 
                                 "Oncorhynchus tshawytscha", "Ophiodon elongatus", 
                                 "Sebastes ensifer",
                                 "Carcharodon carcharias", "Gadus chalcogrammus",
                                 "Microstomus pacificus", "Hippoglossus stenolepis",
                                 "Oncorhynchus nerka"),
                      labels = c("Arrowtooth flounder", "Chum salmon", 
                                 "Coho salmon", "Steelhed salmon", 
                                 "Chinook salmon", "Lingcod",
                                 "Swordspine rockfish",
                                 "Great white shark", "Alaska polluck",
                                 "Pacific Dover sole", "Pacific halibut",
                                 "Sockeye salmon"))


#### Proportional abundance by month -------------------------------------------

### Convert phlyoseq object to dataframe

species_prop <- ps.sp_major@otu_table %>% 
  as.data.frame()
names(species_prop) <- ps.sp_major@tax_table[, 7]

sample_meta <- sample_data(ps.sp_major) %>%
  as("data.frame") %>%
  rownames_to_column("Sample") %>%
  mutate(
    Month = case_when(
      grepl("JAN", FieldID, ignore.case = TRUE) ~ 1,
      grepl("FEB", FieldID, ignore.case = TRUE) ~ 2,
      grepl("MAR", FieldID, ignore.case = TRUE) ~ 3,
      grepl("APR", FieldID, ignore.case = TRUE) ~ 4,
      grepl("MAY", FieldID, ignore.case = TRUE) ~ 5,
      grepl("JUN", FieldID, ignore.case = TRUE) ~ 6,
      grepl("JUL", FieldID, ignore.case = TRUE) ~ 7,
      grepl("AUG", FieldID, ignore.case = TRUE) ~ 8,
      grepl("SEP", FieldID, ignore.case = TRUE) ~ 9,
      grepl("OCT", FieldID, ignore.case = TRUE) ~ 10,
      grepl("NOV", FieldID, ignore.case = TRUE) ~ 11,
      grepl("DEC", FieldID, ignore.case = TRUE) ~ 12,
    )
  )
sample_meta <- sample_meta %>%
  mutate(
    Month = ifelse(Sample == "WADE-002-056-rep-nc", 9, Month)
  )
species_prop_meta <- species_prop %>%
  rownames_to_column("Sample") %>%
  left_join(sample_meta, by = "Sample") %>%
  relocate(14:length(.), .after = Sample)

species_prop_meta_long <- species_prop_meta %>%
  pivot_longer(
    cols = 21:length(.),
    names_to = "Species",
    values_to = "Proportion"
  ) %>% 
  filter(!(Species %in% c("Oncorhynchus mykiss","Gadus chalcogrammus", "Microstomus pacificus",
                          "Oncorhynchus nerka", "Carcharodon carcharias", "Sebastes ensifer")))


#plot proportional abundance
monthly_prop <- ggplot(species_prop_meta_long, aes(x = as.numeric(as.character(Month)), y = Proportion)) +
  geom_jitter(width = 0.1, size = 2, color = "slategray") +
  geom_smooth(method = "loess", span = 0.4, alpha = 0.4, color = "darkcyan", fill = "cyan4")+
  scale_x_continuous(breaks = 1:12, labels = 1:12) +
  xlab("Month") +
  ylab("Proportion of Diet") +
  coord_cartesian(ylim=c(0, 1)) +
  facet_wrap(~Species,
             ncol = 3,
             scales = "free",
             labeller = labeller(Species = c(
    "Atheresthes stomias" = "Arrowtooth flounder",
    "Oncorhynchus keta" = "Chum salmon",
    "Oncorhynchus kisutch" = "Coho salmon",
    "Oncorhynchus mykiss" = "Steelhead salmon",
    "Oncorhynchus tshawytscha" = "Chinook salmon",
    "Ophiodon elongatus" = "Lingcod",
    "Sebastes ensifer" = "Swordspine rockfish",
    "Carcharodon carcharias" = "Great white shark",
    "Gadus chalcogrammus" = "Alaska polluck",
    "Microstomus pacificus" = "Pacific Dover sole", 
    "Hippoglossus stenolepis" = "Pacific halibut",
    "Oncorhynchus nerka" = "Sockeye salmon"
  ))) +
  theme_minimal()
monthly_prop

monthly_prey_proportion <- species_prop_meta_long %>% 
  group_by(Month, Species) %>% 
  summarise(mean_monthly_prop = mean(Proportion))

### Proportional abundance and diet diversity by year --------------------------

samdf_year <- samdf %>% 
  mutate(
    Year = case_when(
      grepl("18", FieldID, ignore.case = TRUE) ~ 2018,
      grepl("19", FieldID, ignore.case = TRUE) ~ 2019,
      grepl("20", FieldID, ignore.case = TRUE) ~ 2020,
      grepl("21", FieldID, ignore.case = TRUE) ~ 2021,
      grepl("22", FieldID, ignore.case = TRUE) ~ 2022,
      grepl("23", FieldID, ignore.case = TRUE) ~ 2023,
      grepl("24", FieldID, ignore.case = TRUE) ~ 2024)
    )

sample_data(ps.sp_major) <- samdf_year
sample_data(ps.sp) <- samdf_year

#barplot
year_bar <- plot_bar(ps.sp_major, fill = "Species") +
  facet_wrap(~Year, scale = "free_x") +
  theme(axis.text.x = element_blank()) +
  scale_fill_discrete(limits = c("Atheresthes stomias","Oncorhynchus keta","Oncorhynchus kisutch","Oncorhynchus mykiss", 
                                 "Oncorhynchus tshawytscha", "Ophiodon elongatus", "Sebastes ensifer",
                                 "Carcharodon carcharias", "Gadus chalcogrammus",
                                 "Microstomus pacificus", "Hippoglossus stenolepis",
                                 "Oncorhynchus nerka"),
                      labels = c("Arrowtooth flounder", "Chum salmon", "Coho salmon", "Steelhed salmon", 
                                 "Chinook salmon", "Lingcod","Swordspine rockfish",
                                 "Great white shark", "Alaska polluck",
                                 "Pacific Dover sole", "Pacific halibut",
                                 "Sockeye salmon"))+
  theme(legend.position = "top")

#boxplot
year_box <- plot_richness(ps.sp, x = "Year", measures=c("Observed","Shannon","Simpson")) + 
  geom_boxplot(aes(group = as.factor(Year),
                   fill = as.factor(Year))) +
  theme(legend.position = "none") +
  scale_fill_manual(values=pnw_palette(n=8,name="Shuksan")) +
  scale_x_continuous(breaks=seq(2018,2024,1))


year_bar / year_box

### Save plots -----------------------------------------------------------------
save(year_box, year_bar, monthly_prey_proportion, monthly_prop, season_bar, 
     dna_con_subset, reads_subset, species_prop_bar, species_prop_subset, file = "SRKW_WO_16Sdiet_phyloseqout.Rdata")

