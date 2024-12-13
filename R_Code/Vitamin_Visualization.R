



library(tidyverse)



#define inputs
vit.quant.file <- "Intermediates/Vit_Quant_QCed.csv"
sample.info.file <- "Meta_Data/Sample_Information_12112024.csv"



#Read in data 
vit.dat <- read_csv(vit.quant.file) %>%
  filter(!str_detect(SampID, "Std"),
         !str_detect(SampID, "Poo")) %>%
  separate(SampID, into = c("Date", "Smp_Type", "UCYNA", "ID"), sep = "_", remove = FALSE) 



#Read in sample information
smp.info <- read_csv(sample.info.file)


##Combine sample information and data:
vit.dat.info <- vit.dat %>%
  left_join(., smp.info) %>%
  mutate(nmoles_in_vial = nM_in_vial*400/1e6,
         Smp_Biovolume_um3 = Cells_on_Filter*Cell_Volume_um3,
         Smp_Biovolume_L = (Cells_on_Filter*Cell_Volume_um3)/1E15,
         nmole_per_cell = nmoles_in_vial/Cells_on_Filter,
         nmole_per_biovolume_L = nmoles_in_vial/Smp_Biovolume_L)




#________Write data into nice neat form for export______________________
vit.dat.export <- vit.dat.info %>%
  filter(Compound %in% c("DMB", "Me-B12", "Ado-B12", "OH-B12", "CN-B12", 
                         "CN-Pseudocob", "OH-Pseudocob", "Me-Pseudocob", "Ado-Pseudocob")) %>%
  select(Compound, SampID, ID, type, replicate, Condition_B12, Condition_TOD, Condition_Culture, 
         Cells_on_Filter, Cell_Volume_um3, Smp_Biovolume_um3, Smp_Biovolume_L, Area, nM_in_vial, 
         nmoles_in_vial, nmole_per_cell, nmole_per_biovolume_L, min.area.flag, blk.ratio.flag) %>%
  rename("PeakArea" = Area)

write_csv(vit.dat.export, file = "Intermediates/UCYNA_Vitamin_Data_121324.csv")



#_________Visualize Data_________________

##Select data that was analyzed:
b12.dmb.dat <- vit.dat.info %>%
  filter(Compound %in% c("DMB", "Me-B12", "Ado-B12", "OH-B12", "CN-B12"))


#Plot all data:

#
ggplot(b12.dmb.dat, aes(x = SampID, y = nmoles_in_vial)) +
  geom_col(aes(color = Condition_B12, fill = Condition_TOD), size = 1) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  facet_grid(Compound~type, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

#
ggplot(b12.dmb.dat, aes(x = SampID, y = nmoles_in_vial)) +
  geom_col(aes(color = Condition_B12, fill = Condition_TOD), size = 1) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  facet_grid(Compound~type, scales="free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))








#Just culture data:
dat.cultures <- b12.dmb.dat %>%
  filter(Condition_Culture %in% c("Culture", "UCYN-A"))

ggplot(dat.cultures, aes(x = SampID, y = nmole_per_cell)) +
  geom_col(aes(color = Condition_B12, fill = Condition_TOD), size = 1) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  facet_grid(Compound~type, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))


ggplot(dat.cultures, aes(x = SampID, y = nmole_per_biovolume_L)) +
  geom_col(aes(color = Condition_B12, fill = Condition_TOD), size = 1) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  facet_grid(Compound~type, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))




####_____Pseudocobalamin data:
##Combine sample information and data:
pseudo.dat <- vit.dat %>%
  filter(Compound %in% c("CN-Pseudocob", "OH-Pseudocob", "Me-Pseudocob", "Ado-Pseudocob")) %>%
  left_join(., smp.info) #%>%
  # mutate(nmoles_in_vial = nM_in_vial*400/1e6,
  #        Smp_Biovolume_um3 = Cells_on_Filter*Cell_Volume_um3,
  #        Smp_Biovolume_L = (Cells_on_Filter*Cell_Volume_um3)/1E15,
  #        nmole_per_cell = nmoles_in_vial/Cells_on_Filter,
  #        nmole_per_biovolume_L = nmoles_in_vial/Smp_Biovolume_L)


pb12.fig <- ggplot(pseudo.dat, aes(x = SampID, y = Area)) +
  geom_col(aes(color = Condition_B12, fill = Condition_TOD), size = 1) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  facet_grid(Compound~type, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 6385) +
  geom_hline(yintercept = 40000, color = "blue") +
  ylab("Peak Area")
pb12.fig

ggsave(pb12.fig, file = "Figures/pB12_fig_12122024.pdf", height = 7, width = 6, scale = 1.2)


