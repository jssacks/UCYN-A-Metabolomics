




library(tidyverse)
library(broom)


###Define inputs
vit.file <- "Raw_Data/UCYNA_B12_TQS_2.csv"

#meta_data
tqs.transitions.file <- "Meta_Data/TQS_Vit_Transitions.csv"




####Define QC Thresholds
blk.ratio <- 3
min.area <- 40000





#__________Load in data____________

####Load in SRM transition monitioring information (identifies which ion to quantify)
tqs.transitions <- read_csv(tqs.transitions.file)

#Load in vitamin dataset
vit.dat <- read_csv(vit.file) %>%
  #filter(Peptide %in% c("OH-Pseudocob", "CN-Pseudocob", "OH-B12", "CN-B12",)) #%>%
  rename("Compound" = "Precursor Ion Name",
         "SampID" = "Replicate Name",
         "Product_mz" = "Product Mz") %>%
  filter(Product_mz %in% tqs.transitions$Product_mz) %>%
  select(Compound, SampID, Area)



##________Define blanks and determine blank thresholds__________

#make list of sample names of blanks
blk.list <- c("241121_Smp_UCYNA_Z-2",
              "241121_Smp_UCYNA_Z-4",
              "241121_Smp_UCYNA_Z-7",
              "241121_Smp_UCYNA_Z-12",
              "241121_Smp_UCYNA_Z-13",
              "241121_Smp_UCYNA_Z-17")

#Pull out blank data
blk.dat <- vit.dat %>%
  filter(SampID %in% blk.list) %>%
  group_by(Compound) %>%
  mutate(mean.blk.area = mean(Area),
         blk.threshold = mean.blk.area*blk.ratio)


##_____Run QC for blk ratio threshold and minimum area
vit.dat.qc <- vit.dat %>%
  left_join(., blk.dat %>% select(Compound, blk.threshold) %>% unique()) %>%
  mutate(min.area.flag = case_when(Area > min.area ~ NA,
                                   TRUE ~ "Flag"),
         blk.ratio.flag = case_when(Area > blk.threshold ~ NA,
                                    TRUE ~ "Flag"))


# Quantify Vitamins using Calibration Curve ----------------------

# Separate calibration curve samples and pull out std concentration from sample name
cal.curve.dat <- vit.dat %>%
  filter(str_detect(SampID, "Std")) %>%
  filter(!str_detect(SampID, "Mix234")) %>%
  separate(SampID, into = c("Date", "Smp_Type", "Conc_nM", "UCYNA", "B12", "Rep"), remove = FALSE) %>% #separate SampID to get std concentration
  mutate(Conc_nM = as.numeric(str_remove(Conc_nM, "nM")))

# Visualize calibration curves:
ggplot(cal.curve.dat, aes(x = Conc_nM, y = Area)) +
  geom_smooth(method = lm) +
  geom_point(aes(color = Rep)) +
  facet_wrap(.~Compound, scales = "free")

# Run linear models to determine calibration curve slopes
lm.dat <- cal.curve.dat %>%
  group_by(Compound)

cal.curve.lms <- do(lm.dat,
                   tidy(
                     lm(Area ~ Conc_nM, data = .)))

cal.curve.outputs <- cal.curve.lms %>%
  filter(term == "Conc_nM") %>%
  select(Compound, estimate) %>%
  rename("slope" = estimate)



#___Combine calibration curve data with sample data to quantify vitamins in each sample
vit.quant.dat <- vit.dat.qc %>%
  left_join(., cal.curve.outputs) %>%
  mutate(nM_in_vial = Area/slope)



#_____ Quantify DMB separately using 1 point calibratin curve (response factor)
dmb.dat <- vit.dat.qc %>%
  filter(Compound == "DMB") %>%
  mutate(slope = Area[SampID == "241121_Std_50nM_RP-Mix234_InH2O"]/50,
         nM_in_vial = Area/slope)
        

#Combine quantified vitamin data
vit.quant.final.dat <- vit.quant.dat %>%
  filter(!Compound == "DMB") %>%
  rbind(dmb.dat)


#Export quantified vitamin data:
write_csv(vit.quant.final.dat, file = "Intermediates/Vit_Quant_QCed.csv")

