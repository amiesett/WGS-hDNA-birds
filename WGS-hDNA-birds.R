# R code fore analyses in manuscript: #
# "Library preparation method and DNA source influence endogenous DNA recovery #
# from 100-year-old avian museum specimens" #
# Author: AE Settlecowski #


library("tidyverse")
library("RColorBrewer")
library(ggpubr)
library(rstatix)
library(viridis)
library(ggforce)
library("lme4")

setwd("C:/Users/amies/Documents/Ethiopia-birds/Analyses/WGS-hDNA-birds")

###############################################################################
# Specify plot features #
###############################################################################
colors_method <- c("#009E73","#0072B2","#CC79A7")
colors_dig <- c("#E69F00","#CC79A7")
shapes_type <- c(25,22)

###############################################################################
# INPUT DATA AND PARSE #
###############################################################################
# SEQUENCE METRICS FROM BIOINFORMATICS ANALYSES # 
seq_metrics <- read_csv("Turdus_shallow_stats.csv", col_names = TRUE) %>%
  select(-c(`Avg_read_length`,`Total_depth`)) %>%
  mutate(Ref_dif/Unclipped_bases) %>%
  mutate(Cleaned_bases/Raw_bases) %>%
  mutate(Unclipped_bases/Cleaned_bases) %>%
  mutate(Unclipped_bases/Raw_bases)
## Create sample-only (no negatives) sequence metrics dataframe
seq_samples <- seq_metrics %>% 
  filter(str_detect(Sample, '_T_|_S_')) %>%
  extract(Sample, c("Sample","Type","Method"),
          regex = "(\\w{3}\\_\\w{3}\\_\\w+)_(\\w{1})_(\\w+)") %>%
  mutate(Type = str_replace(Type,"T", "toepad")) %>%
  mutate(Type = str_replace(Type,"S", "skin"))
seq_samples$Type <- factor(seq_samples$Type, levels=c("toepad","skin"))
seq_samples$Method <- factor(seq_samples$Method, levels=c("IDT","KAPA","SRSLY"))

### Plot genetic distance from reference genome to check for contamination 
p_contam <- ggplot(seq_samples %>% 
                     mutate(Sample = str_replace(Sample,"Tur_mig_162188A","Tur_mig_162188")) %>%
                     mutate(Sample = str_replace(Sample,"Tur_aby_83114A","Tur_aby_83114*")) %>%
                     mutate(Sample = str_replace(Sample,"Tur_aby_83114B","Tur_aby_83114")) %>%
                     mutate(Sample = str_replace(Sample,"Tur_mig_162188B","Tur_mig_162188*")),
                   mapping=aes(x=Sample, y=Ref_dif/Unclipped_bases,
                               color=Method, shape=Type)) +
  #geom_col(position="dodge") +
  geom_point() +
  scale_color_manual(values=colors_method) +
  #facet_grid(Type~.)+
  theme(plot.background=element_blank(), panel.background=element_blank(), 
        panel.border=element_rect(fill=NA), 
        legend.key=element_rect(fill=NA), legend.background=element_rect(fill=NA),
        text=element_text(size=12), axis.text.x=element_text(angle=30, vjust=1, hjust=1),
        plot.margin = margin(5,5,5,5, unit="pt"), legend.key.size = unit(7.5,"pt"),
        legend.box.spacing = unit(0,"pt"), legend.spacing = unit(-3,"pt")) +
  ylab(str_c("Reference differences"))
p_contam
ggsave("FigureS4.png", height=2.5, width=6.5, units="in", dpi=600)

## Create dataframe for predigest test replicates only
seq_predigest <- seq_samples %>%
  filter(Sample=="Tur_aby_83114A" | Sample=="Tur_mig_162188A" |
           Sample=="Tur_aby_83114B" | Sample=="Tur_mig_162188B") %>%
  mutate(`Predigest` = case_when(Sample=="Tur_aby_83114A" ~ "control",
                                  Sample=="Tur_aby_83114B" ~ "predigested",
                                  Sample=="Tur_mig_162188B" ~ "control",
                                  Sample=="Tur_mig_162188A" ~ "predigested")) %>%
  mutate(Sample = str_replace(Sample,"[A-B]", ""))
seq_predigest$Predigest <- factor(seq_predigest$Predigest, levels=c("control","predigested"))

### Remove control replicates from samples sequence dataframe
seq_samples <- seq_samples %>% 
  filter(Sample!="Tur_aby_83114A" & Sample!="Tur_mig_162188B") %>%
  mutate(Sample = str_replace(Sample,"[A-B]", ""))
## Create negatives-only sequence metrics dataframe
seq_negatives <- seq_metrics %>% filter(str_detect(Sample,"negative"))

# READ LENGTH METRICS #
rl_metrics <- read_csv("Turdus_shallow_rl_metrics.csv", col_names = TRUE)
# setwd("./read-lengths")
# rl_metrics=tibble(Sample=character(), Mean_rl=double(),
#                   Med_rl=double(), Std_rl=double())
# rl_comb_counts=tibble(Length=30:187)
# for (s in 1:length(seq_metrics$Sample)) {
#   rl_files <- list.files(path = ".", pattern=str_c(seq_metrics$Sample[s],"_lengths"),
#                          full.names=FALSE) %>% set_names()
#   rl_comb <- tibble(Count=double())
#   for (x in 1:length(rl_files)) {
#       rl_counts <- read_csv(rl_files[x], col_names = FALSE) %>%
#       rename(Count=X1, Length=X2) %>% uncount(Count)
#       rl_mean <- mean(rl_counts$Length) 
#       rl_med <- median(rl_counts$Length)
#       rl_std <- sd(rl_counts$Length)
#       rl_metrics <- rl_metrics %>% add_row(Sample=rl_files[x], Mean_rl=rl_mean,
#               Med_rl=rl_med, Std_rl=rl_std)
#       rl_comb <- rl_comb %>% rbind(rl_counts)
#   }
#   rl_mean <- mean(rl_comb$Length) 
#   rl_med <- median(rl_comb$Length)
#   rl_std <- sd(rl_comb$Length)
#   # Since printing rl_samples.csv have next line to concatenate "_lengthsC" to sample name
#   rl_metrics <- rl_metrics %>% add_row(Sample=str_c(seq_metrics$Sample[s], "_lengthsC")
#                                        Mean_rl=rl_mean, Med_rl=rl_med, Std_rl=rl_std)
#   rl_comb <- rl_comb %>% count(Length) %>% dplyr::rename(!!seq_metrics$Sample[s]:=n)
#   rl_comb_counts <- rl_comb_counts %>% left_join(rl_comb)
# }
# #### Can delete after confirming above edit to concatenate "_lengthsC" works
# rl_metrics <- rl_metrics %>% 
#   mutate(Sample = if_else(str_detect(Sample, regex(".*[^v]$")) == TRUE,
#                                      paste0(Sample, "_lengthsC"), Sample))
# ## Write dataframe to csv for future import
# write_csv(rl_metrics, "Turdus_shallow_rl_metrics.csv")

## Create sample-only read length metrics dataframe
rl_samples <- rl_metrics %>%
  filter(!str_detect(Sample,"neg")) %>%
  mutate(Sample=str_replace(Sample,".csv","")) %>%
  extract(Sample, c("Sample","Type","Method","Read"),
          regex = "(\\w{3}\\_\\w{3}\\_\\w+)_(\\w{1})_(\\w+)_lengths(\\w{1})") %>%
  mutate(Type = str_replace_all(Type, c("T"="toepad","S"="skin"))) %>%
  mutate(Read = str_replace_all(Read, c("1"="R1", "2"="R2",
                                        "U"="Merged_orphaned_R1",
                                        "C"="All_reads")))
rl_samples$Type <- factor(rl_samples$Type, levels=c("toepad","skin"))
rl_samples$Method <- factor(rl_samples$Method, levels=c("IDT","KAPA","SRSLY"))

## Create dataframe of predigest test replicates-only read length metrics
rl_predigest <- rl_samples %>%
  filter(Sample=="Tur_aby_83114A" | Sample=="Tur_mig_162188A" |
           Sample=="Tur_aby_83114B" | Sample=="Tur_mig_162188B") %>%
  mutate(`Predigest` = case_when(Sample=="Tur_aby_83114A" ~ "control",
                                  Sample=="Tur_aby_83114B" ~ "predigested",
                                  Sample=="Tur_mig_162188B" ~ "control",
                                  Sample=="Tur_mig_162188A" ~ "predigested")) %>%
  mutate(Sample = str_replace(Sample,"[A-B]", ""))
rl_predigest$Predigest <- factor(rl_predigest$Predigest, levels=c("control","predigested"))
### Remove control replicates from samples readlength dataframe
rl_samples <- rl_samples %>% 
  filter(Sample!="Tur_aby_83114A" & Sample!="Tur_mig_162188B") %>%
  mutate(Sample = str_replace(Sample,"[A-B]", ""))

## Negatives-only sequence metrics dataframe
rl_negatives <- rl_metrics %>% filter(str_detect(Sample,"negative"))

# LABWORK METRICS #
lab_metrics <- read_csv("WGS-hDNA-birds-lab-data.csv", col_names = TRUE) %>%
  select(c(Sample,Type,`Post-extract_yield`,`Post-repair_yield`)) %>%
  gather(`Post-extract_yield`,`Post-repair_yield`, key="Stage", value="Yield") %>%
  mutate(Stage = case_when(Stage=="Post-extract_yield" ~ "extract",
                           Stage=="Post-repair_yield" ~ "repaired")) %>% distinct()
## Create sample-only lab metrics dataframe
lab_samples <- lab_metrics %>% filter(!str_detect(Sample, "neg"))

## Create dataframe of only predigest replicates lab data
lab_predig <- lab_samples %>%
  filter(Sample=="Tur_aby_83114A" | Sample=="Tur_mig_162188A" |
           Sample=="Tur_aby_83114B" | Sample=="Tur_mig_162188B") %>%
  mutate(Predigest = case_when(Sample=="Tur_aby_83114A" ~ "control",
                                  Sample=="Tur_aby_83114B" ~ "predigested",
                                  Sample=="Tur_mig_162188B" ~ "control",
                                  Sample=="Tur_mig_162188A" ~ "predigested")) %>%
  mutate(Sample = str_replace(Sample,"[A-B]", ""))
### Remove control replicates from samples readlength dataframe
lab_samples <- lab_samples %>% 
  filter(Sample!="Tur_aby_83114A" & Sample!="Tur_mig_162188B") %>%
  mutate(Sample = str_replace(Sample,"[A-B]", ""))

# BIOANALYZER DATA #

#install.packages(install.packages("https://github.com/jwfoley/fraglyzeR/releases/download/v0.10.0/fraglyzeR_0.10.0.tar.gz", repos = NULL))
library("bioanalyzeR")
library(gridExtra)
library(grid)
library(ggpubr)

setwd("C:/Users/amies/Documents/Ethiopia-birds/Analyses/WGS-hDNA-birds/bioanalyzer-xml")

## Input traces
frag_f <- list.files(pattern=".xml")
frag_data <- read.electrophoresis(frag_f[1])
for (x in 2:length(frag_f)) {
  print(frag_f[x])
  frag_data <- frag_data %>% rbind(read.electrophoresis(frag_f[x]))
}

## Remove samples that were excluded prior to lib preps due to low conc.
frag_data <- frag_data %>% subset(!(str_detect(sample.name, '175416' ))) 

## Move negatives to their own dataframe and double check that all are negatives 
frag_negs <- frag_data %>% subset(str_detect(sample.name, 'negative')) 
frag_negs$samples$sample.name

## Create subset of extraction data
frag_dna <- frag_data %>% 
  subset(str_detect(sample.name, '_T_|_S_') 
         & str_detect(sample.name, 'extract|repaired'))
### Split `sample.name` and calculate mean molecule length >=100 and <100
frag_dna$samples <- frag_dna$samples %>%  
  extract(sample.name, c("sample.name","Type","Stage"),
          regex = "(\\w{3}\\_\\w{3}\\_\\w+)_(\\w{1})_(\\w+)") %>%
  mutate(Type = str_replace(Type,"T", "toepad")) %>%
  mutate(Type = str_replace(Type,"S", "skin")) %>%
  bind_cols(ratio_gt100=region.ratio(frag_dna, c(36,99), c(100,500)),
            summarize.custom(frag_dna, 36, 10000))
frag_dna$samples$Type <- factor(frag_dna$samples$Type, levels=c("toepad","skin"))
frag_dna$samples$Stage <- factor(frag_dna$samples$Stage, levels=c("extract","repaired"))

### Check for repeat traces
repeat_frag_dna <- frag_dna$samples %>% count(sample.name, Type, Stage)
#### 2 traces for 83114A&B, skin, repaired samples 
#### keep redo of each from batch: bioanalyzer-negatives-extract-redo-13-Apr-2022_001
frag_dna$samples <- frag_dna$samples %>%
  subset(!(batch=="bioanalyzer-repaired-extract1-extract2-12-13-3-Mar-2022_001"& 
             Type=="skin" & Stage=="repaired" & str_detect(sample.name, "83114")))

## Create subset of library data
frag_lib<- frag_data %>% 
  subset(str_detect(sample.name, '_T_|_S_') 
         & str_detect(sample.name, 'IDT|KAPA|SRSLY'))
### Split `sample.name` and calculate mean molecule length >=100
frag_lib$samples <- frag_lib$samples %>%  
  extract(sample.name, c("sample.name","Type","Method"),
          regex = "(\\w{3}\\_\\w{3}\\_\\w+)_(\\w{1})_(\\w+)") %>%
  mutate(Type = str_replace(Type,"T", "toepad")) %>%
  mutate(Type = str_replace(Type,"S", "skin")) %>%
  bind_cols(summarize.custom(frag_lib, 100, 10000))
frag_lib$samples$Type <- factor(frag_lib$samples$Type, levels=c("toepad","skin"))
frag_lib$samples$Method <- factor(frag_lib$samples$Method, levels=c("IDT","KAPA","SRSLY"))

### Check for repeat traces (should be 4 for each sample for combo of type & stage)
repeat_frag_lib <- frag_lib$samples %>% count(sample.name, Type, Method) %>%
  filter(n==2) %>% left_join(frag_lib$samples, by=c("sample.name","Type","Method"))
#### write_csv(repeat_frag_lib, "repeat-fragment-libraries.csv",col_names = TRUE)
#### keep 1 of each of the 7 repeated samples
frag_lib$samples <- frag_lib$samples %>% 
  subset(!(batch=="bioanalyzer-extract4-IDT-libs-01-03-4-Apr-2022_001" &
             Type=="toepad" & Method=="IDT"
           & str_detect(sample.name, "83107|83109"))) %>%
  subset(!(batch=="bioanalyzer-extract4-IDT-libs-01-03-4-Apr-2022_001" &
             Type=="skin" & Method=="IDT"
           & str_detect(sample.name, "83115"))) %>%
  subset(!(batch=="bioanalyzer-KAPA-libs31-40-neg02-22-Apr-2022" &
             Type=="skin" & Method=="KAPA"
           & str_detect(sample.name, "175413|162188A"))) %>%
  subset(!(batch=="bioanalyzer-IDT-libs04-14" &
             Type=="toepad" & Method=="IDT"
           & str_detect(sample.name, "83114"))) %>%
  subset(!(batch=="bioanalyzer-IDT-libs04-14" &
             Type=="skin" & Method=="IDT"
           & str_detect(sample.name, "162188A")))

## Create subset of predigestion replicates
### DNA
predig_frag_dna <- frag_dna %>%
  subset(str_detect(sample.name, "83114|162188"))
predig_frag_dna$samples <- predig_frag_dna$samples %>%
  mutate(`Predigest` = case_when(sample.name=="Tur_aby_83114A" ~ "control",
                                 sample.name=="Tur_aby_83114B" ~ "predigested",
                                 sample.name=="Tur_mig_162188B" ~ "control",
                                 sample.name=="Tur_mig_162188A" ~ "predigested")) %>%
  mutate(sample.name = str_replace(sample.name,"[A-B]", ""))

### Libraries
predig_frag_lib <- frag_lib %>%
  subset(str_detect(sample.name, "83114|162188"))
predig_frag_lib$samples <- predig_frag_lib$samples %>%
  mutate(`Predigest` = case_when(sample.name=="Tur_aby_83114A" ~ "control",
                                 sample.name=="Tur_aby_83114B" ~ "predigested",
                                 sample.name=="Tur_mig_162188B" ~ "control",
                                 sample.name=="Tur_mig_162188A" ~ "predigested")) %>%
  mutate(sample.name = str_replace(sample.name,"[A-B]", ""))

### Remove non-predigested replicates from frag_dna and frag_lib data frames
frag_dna <- frag_dna %>% 
  subset(sample.name!="Tur_aby_83114A" & sample.name!="Tur_mig_162188B")
frag_dna$samples <- frag_dna$samples %>% 
  mutate(sample.name = str_replace(sample.name,"[A-B]", ""))
#extracts <- subset(samples, str_detect(Stage,'extract'))
#repaired <- subset(samples, str_detect(Stage,'repaired'))
frag_lib <- frag_lib %>% 
  subset(sample.name!="Tur_aby_83114A" & sample.name!="Tur_mig_162188B")
frag_lib$samples <- frag_lib$samples %>% 
  mutate(sample.name = str_replace(sample.name,"[A-B]", ""))

# Plot bioanalyzer traces 
## Predigestion extract/repair traces
traces_predig_frag_dna <- qplot.electrophoresis(predig_frag_dna,
                                                linetype=Type,
                                                region.alpha=NA, 
                                                include.markers=TRUE,
                                                show.peaks="none",
                                                log="x",
                                                xlim=c(36,500), 
                                                scales="free_y") +
  geom_vline(data=predig_frag_dna$samples, aes(xintercept=`Mean in 36-10000`, 
                                               linetype=Type), color="grey") +
  facet_grid(Stage+Predigest~sample.name, scales="free") +
  theme(plot.background=element_blank(),
        panel.background=element_blank(), 
        panel.border=element_rect(fill=NA),
        legend.position = "none",
        legend.key=element_rect(fill=NA), legend.background=element_rect(fill=NA),
        aspect.ratio=1,
        plot.margin = margin(5,5,5,5, unit="pt"),
        text=element_text(size=10),
        axis.text.x=element_text(angle=45),
        legend.box.spacing = unit(0,"pt"), legend.spacing = unit(-3,"pt")) +
  ylab("molarity (pM) per length")
traces_predig_frag_dna
ggsave("FigureS2A.png", height=5.83, width=4, units="in",dpi=600)

## Predigestion library traces
traces_predig_frag_lib <- qplot.electrophoresis(predig_frag_lib,
                                                linetype=Type, color=Method,
                                                region.alpha=NA, 
                                                include.markers=TRUE,
                                                show.peaks="none",
                                                log="x",
                                                xlim=c(36,500), 
                                                scales="free_y") +
  geom_vline(data=predig_frag_lib$samples, aes(xintercept=`Mean in 100-10000`, 
                                               linetype=Type), color="grey") +
  facet_grid(Method+Predigest~sample.name, scales="free") +
  scale_color_manual(values=colors_method) +
  theme(plot.background=element_blank(),
        panel.background=element_blank(), 
        panel.border=element_rect(fill=NA),
        legend.position = "none",
        legend.key=element_rect(fill=NA), legend.background=element_rect(fill=NA),
        aspect.ratio=1,
        plot.margin = margin(5,5,5,5, unit="pt"),
        text=element_text(size=10),
        axis.text.x=element_text(angle=45),
        legend.box.spacing = unit(0,"pt"), legend.spacing = unit(-3,"pt")) +
  ylab("molarity (pM) per length")
traces_predig_frag_lib
ggsave("FigureS2B.png", height=8.5, width=4, units="in",dpi=600)

## Sample extract and library traces
### Create each sample facet plot with a loop
individuals <- frag_lib$samples %>% distinct(sample.name) %>% pull()
traces_frag_dna <- list()
traces_frag_lib <- list()
x <- 1
for (s in individuals){
  temp_frag_dna <- frag_dna %>% subset(sample.name==s)
  temp_frag_lib <- frag_lib %>% subset(sample.name==s)
  traces_frag_dna[[x]] <- qplot.electrophoresis(temp_frag_dna,
                                                linetype=Type,
                                                region.alpha=NA, 
                                                include.markers=TRUE,
                                                show.peaks="none",
                                                log="x",
                                                xlim=c(36,500), 
                                                scales="free_y") +
    geom_vline(data=temp_frag_dna$samples, aes(xintercept=`Mean in 36-10000`, 
                                               linetype=Type), color="grey") +
    facet_grid(Stage~sample.name, scales="free") +
    theme(plot.background=element_blank(),
          panel.background=element_blank(), 
          panel.border=element_rect(fill=NA),
          legend.position = "none",
          text=element_text(size=8),
          axis.text.x=element_text(angle=45),
          aspect.ratio=1,
          plot.margin = margin(2.5,2.5,-2.5,2.5, unit="pt"),
          axis.title=element_blank(),
          strip.text.y=element_blank())
  traces_frag_lib[[x]] <- qplot.electrophoresis(temp_frag_lib,
                                                linetype=Type,
                                                color=Method,
                                                region.alpha=NA, 
                                                include.markers=FALSE,
                                                show.peaks="none",
                                                log="x",
                                                xlim=c(36,500), 
                                                scales="free_y") +
    geom_vline(data=temp_frag_lib$samples, aes(xintercept=`Mean in 100-10000`, 
                                               linetype=Type), color="grey") +
    facet_grid(Method~sample.name, scales="free") +
    scale_color_manual(values=colors_method) +
    theme(plot.background=element_blank(),
          panel.background=element_blank(), 
          panel.border=element_rect(fill=NA),
          legend.position = "none",
          text=element_text(size=8),
          axis.text.x=element_text(angle=45),
          aspect.ratio=1,
          plot.margin = margin(2.5,2.5,-2.5,2.5, unit="pt"),
          axis.title=element_blank(),
          strip.text.y=element_blank())
  x<-x+1
}
### Gather all sample facet plots into single plot
y <- textGrob(expression(paste("molarity (pM) per length")), 
              rot = 90, gp = gpar(fontsize = 10)) 
x <- textGrob(expression(paste("length (bp)")), 
              gp = gpar(fontsize = 10))#, vjust=-5)
p_traces_frag_dna <- grid.arrange(grobs=traces_frag_dna,nrow=1,
                                  left=y, bottom=x)
ggsave("FigureS1A.png", p_traces_frag_dna, width=10, height=2.5, units="in", dpi=600)
p_traces_frag_lib <- grid.arrange(grobs=traces_frag_lib,nrow=1,
                                  left=y, bottom=x)
ggsave("FigureS1B.png", p_traces_frag_lib, width=10, height=3.5, units="in", dpi=600)

setwd("C:/Users/amies/Documents/Ethiopia-birds/Analyses/WGS-hDNA-birds")

## Rename sample.name column of bioanalyzer dataframes for consistency w/ others
frag_dna$samples <- frag_dna$samples %>% rename(Sample=sample.name)
predig_frag_dna$samples <- predig_frag_dna$samples %>% rename(Sample=sample.name)
frag_lib$samples <- frag_lib$samples %>% rename(Sample=sample.name)
predig_frag_lib$samples <- predig_frag_lib$samples %>% rename(Sample=sample.name)

## Add fragment size metrics to samples lab metrics dataframe
lab_samples <- lab_samples %>% 
  inner_join(frag_dna$samples, by=c("Sample","Type","Stage")) %>%
  select(Sample, Type, Stage, ratio_gt100, `Median in 36-10000`,
         `Mean in 36-10000`, `SD in 36-10000`, Yield)
lab_samples$Type <- factor(lab_samples$Type, levels=c("toepad","skin"))
## Add fragment size metrics to predigest lab metrics dataframe
lab_predig <- lab_predig %>%
  inner_join(predig_frag_dna$samples, by=c("Sample","Type","Stage", "Predigest")) %>%
  select(Sample, Type, Stage, Yield, Predigest, ratio_gt100, `Median in 36-10000`,
         `Mean in 36-10000`, `SD in 36-10000`)
lab_predig$Type <- factor(lab_predig$Type, levels=c("toepad","skin"))
lab_predig$Predigest <- factor(lab_predig$Predigest, levels=c("control","predigested"))

###############################################################################
# ANALYSIS OF DNA SOURCE EFFECT ON DNA (EXTRACTION/REPAIR) YIELD AND SIZE #
###############################################################################
# Check for outliers (83107 shows some contamination)
lab_samples %>% group_by(Stage, Type) %>% identify_outliers(Yield) -> outliers
lab_samples %>% group_by(Stage, Type) %>% identify_outliers(`Mean in 36-10000`) -> outliers
## Mean in 36-10000 outliers: 66823_S extract & repaired, 83114_T repaired, 83114_S repaired
lab_samples %>% group_by(Stage, Type) %>% identify_outliers(ratio_gt100) -> outliers
## Mean in ratio_gt100 outliers: 66823_S repaired, 83114_T repaired, 83114_S repaired

# Check normality
lab_samples %>% group_by(Stage, Type) %>% shapiro_test(Yield)
lab_samples %>% rename(Mean_length=`Mean in 36-10000`) %>% 
  group_by(Stage, Type) %>% shapiro_test(Mean_length)
## Mean in 36-10000 not normal: skin extract & repaired
lab_samples %>% group_by(Stage, Type) %>% shapiro_test(ratio_gt100)
## Mean in ratio_gt100 not normal: skin repaired
ggqqplot(lab_samples, "Yield", ggtheme = theme_bw()) +
  facet_grid(Stage ~ Type, labeller = "label_both")
ggqqplot(lab_samples, "Mean in 36-10000", ggtheme = theme_bw()) +
  facet_grid(Stage ~ Type, labeller = "label_both")
ggqqplot(lab_samples, "ratio_gt100", ggtheme = theme_bw()) +
  facet_grid(Stage ~ Type, labeller = "label_both")

# T-tests
ttest_extract_yield <-lab_samples %>% filter(Stage=="extract") %>%
  t_test(Yield ~ Type, paired = TRUE) %>% add_significance() %>%
  add_column(Stage="extract", .before=1)
ttest_repair_yield <- lab_samples %>% filter(Stage=="repaired") %>%
  t_test(Yield ~ Type, paired = TRUE) %>% add_significance() %>%
  add_column(Stage="repaired", .before=1)
ttest_extract_size <- lab_samples %>% filter(Stage=="extract") %>%
  t_test(`Mean in 36-10000` ~ Type, paired = TRUE) %>% add_significance() %>%
  add_column(Stage="extract", .before=1)
ttest_repair_size <- lab_samples %>% filter(Stage=="repaired") %>%
  t_test(`Mean in 36-10000` ~ Type, paired = TRUE) %>% add_significance() %>%
  add_column(Stage="repaired", .before=1)
ttest_extract_ratio <- lab_samples %>% filter(Stage=="extract") %>%
  t_test(ratio_gt100 ~ Type, paired = TRUE) %>% add_significance() %>%
  add_column(Stage="extract", .before=1)
ttest_repair_ratio <- lab_samples %>% filter(Stage=="repaired") %>%
  t_test(ratio_gt100 ~ Type, paired = TRUE) %>% add_significance() %>%
  add_column(Stage="repaired", .before=1)

## Combine results of ttest into 1 table and write to csv
source_ttest <- bind_rows(ttest_extract_size, ttest_extract_yield,
                          ttest_extract_ratio, ttest_repair_size,
                          ttest_repair_yield, ttest_repair_ratio)
write_csv(source_ttest, "source-ttests.csv")

# Output summary stats of lab metrics
lab_summ <- lab_samples %>% 
  select(-c("Median in 36-10000","SD in 36-10000")) %>% group_by(Type, Stage) %>%
  get_summary_stats(show=c("n","min","max","mean","sd"))
write_csv(lab_summ, "lab-summ-stats.csv")

# Plot DNA yield and size by source
## Yield - repeat for extract (Figs 2A, 2C) and repaired (Figs 2B, 2D)
temp_df <- lab_samples %>% select(Sample, Type, Stage, Yield) %>% spread(Type,Yield)
p_lab_bxp <- ggpaired(temp_df %>% filter(Stage=="extract"),
                      cond1="toepad", cond2="skin",
                      line.color="gray", fill=NA, width=0.75,
                      xlab="Source", ylab="DNA yield (ng/mg)") +
  stat_summary(fun=mean, geom="crossbar", width=0.75, size=0.4, linetype="dashed") +
  stat_summary(fun.data="mean_sd", geom="errorbar", width=0.2) +
  theme(plot.background=element_blank(), panel.background=element_blank(), 
        panel.border=element_rect(fill=NA), 
        text=element_text(size=10), plot.title = element_text(size=10),
        plot.margin = margin(5,5,5,5, unit="pt")) +
  labs(title="Extracted DNA")
p_lab_bxp
ggsave("Figure2A.png", height=3, width=2.25, units="in", dpi=600)

## Size
temp_df <- lab_samples %>% select(Sample, Type, Stage, `Mean in 36-10000`) %>% spread(Type,`Mean in 36-10000`)
p_lab_bxp <- ggpaired(temp_df %>% filter(Stage=="repaired") %>%
                        filter(Sample!="Tur_mig_66823"), cond1="toepad", cond2="skin",
                      line.color="gray", fill=NA, width=0.75,
                      xlab="Source", ylab="Mean DNA size (bp)\n") +
  stat_summary(fun=mean, geom="crossbar", width=0.75, size=0.4, linetype="dashed") +
  stat_summary(fun.data="mean_sd", geom="errorbar", width=0.2) +
  theme(plot.background=element_blank(), panel.background=element_blank(), 
        panel.border=element_rect(fill=NA), 
        text=element_text(size=10), plot.title = element_text(size=10),
        plot.margin = margin(5,5,5,5, unit="pt")) +
  labs(title="Repaired DNA")
p_lab_bxp
ggsave("Figure2D.png", height=3, width=2.25, units="in", dpi=600)

###############################################################################
# ANALYSIS OF METHOD AND DNA SOURCE EFFECTS ON SEQUENCING METRICS #
###############################################################################
# Add repaired fragment length, ratio_gt100, and avg_read length to seq_samples
seq_samples <- seq_samples %>%
  left_join(lab_samples %>% filter(Stage=="repaired") %>% 
    select(c(Sample,Type,ratio_gt100,`Mean in 36-10000`)), 
  by=c("Sample","Type")) %>% 
  left_join(rl_samples %>% filter(Read=="All_reads") %>%
              select(c(Sample,Type,Method,Mean_rl)),
            by=c("Sample","Type","Method"))

# Repeated, 2-way ANOVAs: Influence of method & source on sequencing efficiency 
# endogenous DNA content, mean read length
aov_map_raw <- anova_test(data=seq_samples, dv=Unclipped_bases/Raw_bases, 
                          wid=Sample, within=c(Type, Method)) # Method significant
aov_map_raw <- aov_map_raw[[1]] %>% 
  add_column(Variable="Unclipped_bases/Raw_bases", .before=1)
aov_clean_raw <- anova_test(data=seq_samples, dv=Cleaned_bases/Raw_bases, 
                          wid=Sample, within=c(Type, Method)) # Method significant
aov_clean_raw <- aov_clean_raw[[1]] %>%
  add_column(Variable="Cleaned_bases/Raw_bases", .before=1)
aov_map_clean <- anova_test(data=seq_samples, dv=Unclipped_bases/Cleaned_bases, 
                          wid=Sample, within=c(Type, Method)) # Method significant
aov_map_clean <- aov_map_clean[[1]] %>%
  add_column(Variable="Unclipped_bases/Cleaned_bases", .before=1)
aov_rl <- anova_test(data=seq_samples, dv=Mean_rl, 
                     wid=Sample, within=c(Type, Method)) # Interaction & method significant
aov_rl <- aov_rl[[1]] %>%
  add_column(Variable="Mean_rl", .before=1)
                              
## Significant interaction: 1-way Anova of Mean_rl
aov_rl_method <- seq_samples %>%
  group_by(Method) %>%
  anova_test(dv = Mean_rl, wid = Sample, within = Type) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH") %>% as_tibble()
aov_rl_method <- rename(aov_rl_method, Variable=Method)

## Combine all anova test results into 1 dataframe for output
method_source_aov <- bind_rows(aov_clean_raw, aov_map_clean, aov_map_raw,
                               aov_rl, aov_rl_method)
write_csv(method_source_aov, "method-source-aov.csv")

## Insignificant interaction w/ significant effect - Paired t-tests of method
ttest_map_raw <- seq_samples %>% 
  pairwise_t_test(Unclipped_bases/Raw_bases ~ Method, 
                  paired=TRUE, p.adjust.method = "BH") 
ttest_clean_raw <- seq_samples %>% 
  pairwise_t_test(Cleaned_bases/Raw_bases ~ Method, 
                  paired=TRUE, p.adjust.method = "BH")
ttest_map_clean <- seq_samples %>% 
  pairwise_t_test(Unclipped_bases/Cleaned_bases ~ Method, 
                  paired=TRUE, p.adjust.method = "BH")
ttest_rl <- seq_samples %>% 
  pairwise_t_test(Mean_rl ~ Method, 
                  paired=TRUE, p.adjust.method = "BH")
method_ttest <- bind_rows(ttest_map_raw, ttest_clean_raw, ttest_map_clean,
                                ttest_rl)
write_csv(method_ttest, "method_bases_ttest.csv")

# Plot boxplots of sequencing metrics by DNA source and method
seq_samples <- seq_samples %>%
  distinct(Sample,Method) %>% add_column(paired=1:24) %>%
  right_join(seq_samples)

## Repeat for each sequencing metric dependent variable
### Sequencing efficiency - Unclipped_bases/Raw_bases (3A)
### Endogenous DNA content - Unclipped_bases/Cleaned_bases (3B)
### Mean_rl (3C)
p_seq_bxp <- ggplot(seq_samples, mapping=aes(Type, Unclipped_bases/Cleaned_bases,
                                             color=Method, fill=Method)) +
  facet_grid(.~Method) +
  geom_boxplot(outlier.shape = NA, fill=NA, coef=0) +
  geom_point(aes(shape=Type, group=paired)) +
  geom_line(aes(group=paired),color="gray") +
  stat_summary(fun=mean, geom="crossbar", width=0.75, size=0.4, linetype="dashed") +
  stat_summary(fun.data="mean_sd", geom="errorbar", width=0.2) +
  scale_color_manual(values=colors_method) + 
  scale_fill_manual(values=colors_method) + 
  scale_shape_manual(values=shapes_type, ) +
  theme(plot.background=element_blank(), panel.background=element_blank(), 
        panel.border=element_rect(fill=NA), #aspect.ratio=1,
        text=element_text(size=12), #axis.text.x=element_text(angle=90),
        plot.margin = margin(5,5,5,5, unit="pt"),
        strip.background =element_rect(fill=NA, color="black"),
        legend.key=element_rect(fill=NA), legend.background=element_rect(fill=NA),
        legend.box.spacing = unit(0,"pt"), legend.spacing = unit(-3,"pt")) +
  xlab("Source") +
  ylab("Endogenous content") +
  guides(shape="none",color="none",fill="none",linetype="none")
p_seq_bxp
ggsave("Figure3B.png", height=2, width=4, units="in", dpi=600)

################################################################################
# ANALYSIS OF INPUT DNA SIZE ON SEQUENCING METRICS #
################################################################################
# Change variable name `Mean in 36-10000` bc lmer doesn't like it
seq_samples <- seq_samples %>% rename(Mean_frag=`Mean in 36-10000`)

# Run LMMs for each sequencing metric
## Sequencing efficiency - Mapped/raw bases
summary(lmer(Unclipped_bases/Raw_bases ~ Type + Method + (1|Sample),
             data=seq_samples))
summary(lmer(Unclipped_bases/Raw_bases ~ Type + Method 
             + Mean_frag + (1|Sample), data=seq_samples))
lmm_map_raw_null <- lmer(Unclipped_bases/Raw_bases ~ Type + Method + (1|Sample),
                         data=seq_samples, REML=FALSE)
lmm_map_raw_frag <- lmer(Unclipped_bases/Raw_bases ~ Type + Method 
                         + Mean_frag + (1|Sample), data=seq_samples, REML=FALSE)
anova(lmm_map_raw_null,lmm_map_raw_frag)
## Endogenous DNA content - Mapped/cleaned bases
summary(lmer(Unclipped_bases/Cleaned_bases ~ Type + Method + (1|Sample),
             data=seq_samples))
summary(lmer(Unclipped_bases/Cleaned_bases ~ Type + Method 
             + Mean_frag + (1|Sample), data=seq_samples))
lmm_map_clean_null <- lmer(Unclipped_bases/Cleaned_bases ~ Type + Method + (1|Sample),
                           data=seq_samples, REML=FALSE)
lmm_map_clean_frag <- lmer(Unclipped_bases/Cleaned_bases ~ Type + Method 
                           + Mean_frag + (1|Sample), data=seq_samples, REML=FALSE)
anova(lmm_map_clean_null,lmm_map_clean_frag)
## Mean read length
summary(lmer(Mean_rl ~ Type + Method + (1|Sample),
             data=seq_samples))
summary(lmer(Mean_rl ~ Type + Method 
             + Mean_frag + (1|Sample), data=seq_samples))
lmm_rl_null <- lmer(Mean_rl ~ Type + Method + (1|Sample),
                    data=seq_samples, REML=FALSE)
lmm_rl_frag <- lmer(Mean_rl ~ Type + Method 
                    + Mean_frag + (1|Sample), data=seq_samples, REML=FALSE)
anova(lmm_rl_null,lmm_rl_frag)

# Plot LMMs of sequencing metrics vs mean repaired fragment length
## Repeat for each sequencing metric dependent variable
### Sequencing efficiency - Unclipped_bases/Raw_bases (S3A)
### Endogenous DNA content - Unclipped_bases/Cleaned_bases (S3B)
### Mean_rl (S3C)
p_rl_v_frag <- ggplot(seq_samples %>% filter(Sample!="Tur_mig_66823"),
                         mapping=aes(x=`Mean in 36-10000`, y=Mean_rl,
                                     shape=Type, linetype=Type,
                                     color=Method, fill=Method)) +
  geom_point(alpha=0.75) +
  geom_smooth(method="lm",fill=NA) +
  scale_color_manual(values=colors_method) + 
  scale_fill_manual(values=colors_method) + 
  scale_shape_manual(values=shapes_type) +
  theme(plot.background=element_blank(), panel.background=element_blank(), 
        panel.border=element_rect(fill=NA), aspect.ratio=1,
        plot.margin = margin(5,5,5,5, unit="pt"), text=element_text(size=10), 
        legend.key=element_rect(fill=NA), legend.background=element_rect(fill=NA),
        legend.box.spacing = unit(0,"pt"), legend.spacing = unit(-3,"pt")) +
  xlab("Mean DNA size (bp)") +
  ylab("Mean read length (bp)") # Sequencing efficiency, Endogenous DNA content, Mean read length (bp)
p_rl_v_frag
ggarrange(p_seq_v_frag, p_end_v_frag, p_rl_v_frag, ncol=1,
          common.legend = TRUE, legend="right")
ggsave("FigureS3.png", height=6, width=3, units="in", dpi=600)

###############################################################################
# INFLUENCE OF PREDIGESTION ON DNA & SEQUENCING METRICS # 
###############################################################################
seq_predigest <- seq_predigest %>%
  left_join(lab_predig %>% filter(Stage=="repaired") %>% 
              select(c(Sample,Type,Predigest,Yield,ratio_gt100,`Mean in 36-10000`)), 
            by=c("Sample","Type","Predigest")) %>% 
  left_join(rl_predigest %>% filter(Read=="All_reads") %>%
              select(c(Sample,Type,Method,Mean_rl,Predigest)),
            by=c("Sample","Type","Method","Predigest"))

seq_predig_plot <- seq_predigest %>%
  select(c(Sample,Type,Method,Predigest,Yield,`Mean in 36-10000`,ratio_gt100,
           Mean_rl,`Unclipped_bases/Cleaned_bases`,`Unclipped_bases/Raw_bases`)) %>%
  pivot_wider(id_cols=c(Sample,Type,Method), names_from=Predigest,
              values_from=c(`Mean in 36-10000`,ratio_gt100, Yield, Mean_rl,
                            `Unclipped_bases/Cleaned_bases`,`Unclipped_bases/Raw_bases`))

# Plot difference between control and predigested replicates
## Yield, size, endogenous DNA, sequencing efficiency, mean read length
p_predig_yield <- ggplot(seq_predig_plot, mapping=aes(x=Sample,
                                                      y=`Yield_control`-`Yield_predigested`,
                                                      shape=Type)) +
  ggtitle("DNA yield (ng/mg)") +
  #ggtitle("DNA size (bp)") +
  #ggtitle("Sequencing efficiency") +
  #ggtitle("Endogenous DNA content") +
  #ggtitle("Mean read length (bp)") +
  geom_point() +
  #geom_point(mapping=aes(color=Method, fill=Method)) +
  #scale_color_manual(values=colors_method) + 
  scale_fill_manual(values=colors_method) + 
  scale_shape_manual(values=shapes_type) +
  theme(plot.background=element_blank(), panel.background=element_blank(), 
        panel.border=element_rect(fill=NA),plot.title=element_text(size=10),
        plot.margin = margin(5,5,5,5, unit="pt"), text=element_text(size=10), 
        axis.text.x=element_text(angle=30, vjust=1, hjust=1),
        legend.key=element_rect(fill=NA), legend.background=element_rect(fill=NA),
        legend.box.spacing = unit(0,"pt"), legend.spacing = unit(-3,"pt")) +
  ylab(NULL) +
  xlab(NULL)
FigS5 <- ggarrange(p_predig_yield, p_predig_size, p_predig_end, p_predig_eff, p_predig_rl,
          nrow=1, common.legend = TRUE, legend="top")
FigS5 <- annotate_figure(FigS4,
                         left=textGrob(expression(paste("control - predigested")), 
                                              rot = 90, gp = gpar(fontsize = 10)))
FigS5
ggsave("FigureS5.png", width=9, height=3, units="in", dpi=600)

###############################################################################
# OUTPUT SUMMARY STATISTICS #
###############################################################################
# All samples
## Seq metrics
write_csv(seq_predigest, "seq_metrics_predig.csv")
write_csv(lab_predig, "lab_metrics_predig.csv")
write_csv(seq_samples, "seq_metrics_samples.csv")
write_csv(rl_negatives %>% 
            filter(str_detect(Sample,"lengthsC")) %>%
            select(Sample,Mean_rl), "rl_neg.csv")
## lab metrics
write_csv(lab_samples %>%
            select(c(Sample,Type,Stage,`Mean in 36-10000`,Yield)),
          "lab_metrics_samples.csv")
write_csv(seq_negatives, "seq_metrics_neg.csv")


# By method and type
## seq metrics 
seq_summ_samples <- seq_samples %>% group_by(Type,Method) %>% 
  get_summary_stats(Raw_reads, Cleaned_reads, Mapped_reads, Mean_rl, 
                    `Unclipped_bases/Cleaned_bases`, `Unclipped_bases/Raw_bases`,
                    show = c("min","max","median","mean","sd","se")) %>%
  bind_rows(get_summary_stats(seq_samples, Raw_reads, Cleaned_reads, Mapped_reads, Mean_rl, 
                              `Unclipped_bases/Cleaned_bases`, 
                              `Unclipped_bases/Raw_bases`,
                              show = c("min","max","median","mean","sd","se")))
write_csv(seq_summ_samples, "summ_seq_samples.csv")

## lab metrics 
lab_summ_samples <- lab_samples %>% group_by(Type,Stage) %>% 
  get_summary_stats(Yield,`Mean in 36-10000`,
                    show = c("min","max","median","mean","sd","se"))
write_csv(lab_summ_samples, "summ_lab_samples.csv")

write_csv(seq_predig_plot, "summ_predig.csv")

seq_predigest <- seq_predigest %>%
  left_join(rl_predigest %>% filter(Read=="All_reads") %>%
              select(c(Sample,Type,Method,Mean_rl,Predigest)),
            by=c("Sample","Type","Method","Predigest")) %>% 
  filter(Predigest=="control")



