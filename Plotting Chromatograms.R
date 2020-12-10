####
####    SCRIPT TO VISUALISE CHROMATOGRAMS
####    
####      This script reads chromatograms exported from ChromaTOF
####        (right click in the chromatogram, export as .csv)
####      We had some language encoding issues with the the csv-file
####      so we saved it as excel file. 
####      In theory the unmodified csv can be used as well. 
####    
####      This is an example for a particular setting, 
####      it needs to be adjusted for other situations. 
#### 
####      Matthias Pietzke - 09.12.2020
####

library(tidyverse)
library(readxl)
library(zoo)      # for rolling mean

## ACTION:  define the path
setwd("Z:/GC-TOF/Pegasus_RAW/Lizzy_RAW/2020/20201125_SB_TBDMS test derivatisations-V")


### ACTION: read in the file  
#  "Aspartate.xlsx"   / "Asparagine.xlsx"  /  "Glutamine.xlsx"
input = read_excel("Glutamine.xlsx") 


### ACTION: define the extracted mass used here and extracted in the file
#   Apartate =    418
#   Asparagine =  417
#   Glutamine =   431
extracted_mass = 431

## improve the input  # ACTION: Adapt to your experiment
input2 = input %>% 
  rename("intensity" = paste(extracted_mass )) %>% 
  # replace labels   
  mutate(Label = case_when(Sample == "e20331SB_12_Sample 1_2" ~ "CTRL - 0h",
                           Sample == "e20331SB_12_Sample_1_2)" ~ "CTRL - 0h",
                           Sample == "e20331SB_2_Sample 7_2" ~ "16h Gln-starvation",
                           Sample == "e20331SB_4_Sample 13_2" ~ "+24h NH3",
                           Sample == "e20331SB_6_Sample 19_2" ~ "+24h Asn",
                           Sample == "e20331SB_8_Sample 25_2" ~ "+72h NH3",
                           Sample == "e20331SB_10_Sample 31_2" ~ "+72h Asn"
                           )) 

## add smoothing when needed
input2 = input2 %>% 
  mutate(smoothed = rollmean(intensity, 7,    # 7 = number of points used for smoothing
                             fill = NA))

# double- check the input and the results
unique(input2$Sample)
unique(input2$Label)


##  define plotting order   ACTION: Define for your project
Order = c("CTRL - 0h", "16h Gln-starvation", "+24h NH3", "+72h NH3", "+24h Asn", "+72h Asn")

##  define colours   ACTION: define for ypur project
ColourList = c("CTRL - 0h" = "#7CB342",
               "16h Gln-starvation" = "red",
               "+24h NH3" = "goldenrod1",
               "+24h Asn" = "gainsboro",
               "+72h NH3" = "goldenrod3",
               "+72h Asn" = "dimgrey"
               )

## EXTRA: RI CALCULATION  
#  ACTION: Fill in RT - RI pairs for your experiment
RI_table = tribble(
   ~RT, ~RI,
   640,  1200,  
   1048, 1500,  
   1256, 1700,
   1436, 1900,
   1631, 2200,
   1883, 2800,
   2015, 3200,
   2159, 3600,
   1684, 2322,
   1607, 2163,
   1739, 2456  
   )

# show the plot
ggplot(RI_table, aes(x= RT, y=RI)) +
  geom_point() + 
  geom_smooth(se = FALSE) 

# calculate curve fit
fit = loess(RI ~ RT, RI_table)

###  calculate predicted RIs
###  this runs a loess fit, however the results may not be identical to the ChromaTOF processing
RI_table$RI = predict(fit, RI_table$RT)
input2$RI = predict(fit, input2$Time)


### generate final plot
ggplot(input2, aes(x = RI,         # ACTION:  use Time for time / RI when you calculated the RI 
                   y = smoothed,   # ACTION:  use intensity for raw data  / 
                                   #              smoothed for the smoothed data ! 
                  colour = factor(Label, levels = Order))) +
  geom_line(size = 1.2) +
  theme_bw(base_size = 16) + 
  scale_color_manual(NULL, values = ColourList) + 
  labs(title = "Glutamine TBDMS",        # ACTION: change to name of compound!
       x = "Retention Index",            # ACTION: change to retention time when using the rt.
       y = "Intensity of Quant-Mass") +
  #coord_cartesian(ylim = c(0, 40000)) +  # ACTION: comment or uncomment this line for cropping the plot
   theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = c(0.80, 0.80),
        legend.background = element_rect(colour = "black", size = 0.5) ,
        axis.text = element_text(colour = "black")
        )

# save the plot, remember to change the filename, too!
# fontsizes and lines are currently adjusted for 20 cm, you may adjust this when changing the expected size
ggsave("Glutamine_cropped_smoothed.png", width = 20, height = 20, units = "cm")

