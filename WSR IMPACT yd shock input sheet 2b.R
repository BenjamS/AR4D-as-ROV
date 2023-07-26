#===========================================================================
# Define country groupings to facilitate categorization into yd shock groups
# Note South Asia is left out b/c experts expressed high confidence that WSR
# not an issue there
# Important omissions to mention in paper:
#E. Eur/Balkans: "Romania", "Bulgaria", "Hungary", "Czechia", "Slovakia"
# Central & West Africa -- esp Sudan (Nigeria, Senegal, Mali produce little)
# To a lesser extent Finland, Norway, Georgia, Armenia, and Azerbaijan, NZ, Japan
# East Africa
eAfr_vec <- c("Ethiopia", "Kenya", "Uganda", "Tanzania")
eAfrName <- "East Africa"
# Southern Africa
sAfr_vec <- c("South Africa", "Zambia", "Zimbabwe")
sAfrName <- "Southern Africa"
# Other Balkans
otherBalks_vec <- c("Montenegro", "Bosnia and Herzegovina", "Croatia")
# Mediterranean Basin excl. France
medBasinExclFrance_vec <- c("Spain", "Italy", "Slovenia", "Albania", "Greece", "Turkey", otherBalks_vec,
                  "Egypt", "Israel", "Lebanon", "Libya",
                  "Morocco", "Syria", "Tunisia", "Algeria")
medBasinExclFranceName <- "Mediterranean Basin (excl. France)"
# Baltic States
baltStates_vec <- c("Latvia", "Lithuania", "Estonia")
# Western and Northern Europe, excl. Spain and Italy
wNeurExclSpainItaly_vec <- c("France", "Germany", "Belgium-Luxembourg",
           "UK", "Ireland", "Netherlands", "Austria",
           "Denmark", "Switzerland", "Poland", "Sweden", baltStates_vec)
wNeurExclSpainItalyName <- "Western & Northern Europe (excl. Spain)"
# WSR-vulnerable Central Asia
cAsiaVuln_vec <- c("Uzbekistan", "Tajikistan", "Kazakhstan")
cAsiaVulnName <- "WSR-vulnerable Central Asia"
# WSR-resilient Central Asia
cAsiaResi_vec <- c("Afghanistan", "Kyrgyzstan", "Turkmenistan", "Iran")
cAsiaResiName <- "WSR-resilient Central Asia"
# South & Central America
latAm_vec <- c("Argentina", "Brazil", "Uruguay", "Paraguay",
               "Chile", "Peru", "Colombia", "Ecuador", "Venezuela",
               "Bolivia", "Mexico", "Honduras", "El Salvador",
               "Nicaragua", "Panama", "Costa Rica")
latAmName <- "South & Central America"
#----------------------------------------------------------------------
# Now assign countries to yd shock groups as indicated by the expert consultation
groupYshkHi_vec <- c(eAfr_vec)
groupYshkMe_vec <- c(medBasinExclFrance_vec, cAsiaVuln_vec, latAm_vec,
                     "Iraq", "Canada", "Australia", "Russia")
groupYshkLo_vec <- c(wNeurExclSpainItaly_vec, cAsiaResi_vec, sAfr_vec,
                     "Ukraine", "USA", "China")
all_countries <- c(groupYshkHi_vec, groupYshkMe_vec, groupYshkLo_vec)
# And get sub region names ready for table display purposes
shkHiRegionName_vec <- eAfrName
shkMeRegionName_vec <- c(medBasinExclFranceName, cAsiaVulnName,
                         latAmName, "Iraq", "Canada", "Australia", "Russia")
shkLoRegionName_vec <- c(wNeurExclSpainItalyName, cAsiaResiName, sAfrName,
                         "Ukraine", "USA", "China")
# Assign the yd shock to each of these groups
# Direct yield loss (fractional) due to WSR outbreaks
# Best case
yShkHi_best <- 0.10
yShkMe_best <- 0.05
yShkLo_best <- 0
# Worst case
yShkHi_wrst <- 2 * yShkHi_best
yShkMe_wrst <- 2 * yShkMe_best
yShkLo_wrst <- 0.05
# Now assign outbreak intervals based on expert consultation and literature (in yrs)
# Note outbreak intervals assigned at a much more geographically aggregate level
# For IMPACT csv input file
freq_SSA <- 1
freq_WANA <- 3
freq_LAC <- 5 # Uruguay3-5, CWANA3
freq_RUSCA <- 5 
freq_NAM <- 10
freq_EUR <- 10
freq_EA <- 10
freq_AUS <- 10
freqVec <- c(freq_SSA, freq_WANA, freq_LAC, freq_RUSCA, freq_NAM,
             freq_EUR, freq_EA, freq_AUS)
freqRegNames <- c("Eastern & Southern Africa",
                  "West Asia & North Africa",
                  "South & Central America",
                  "Russia & Central Asia",
                  "North America",
                  "Europe",
                  "China",
                  "Australia")
#---
# For display table
freq_SSA <- 1
freq_WANA <- 3
freq_LAC_RUSCA <- 5
freq_NAM_EUR_EA_AUS <- 10
freqTabVec <- c(freq_SSA, freq_WANA, freq_LAC_RUSCA, freq_NAM_EUR_EA_AUS)
# Get region names ready for frequency table
freqRegTabNames <- c("Eastern & Southern Africa",
                    "West Asia & North Africa",
                    "South & Central America, Russia & Central Asia",
                    "Europe, North America, China, Australia")
#=======================================================================
# Get ydShk df together for IMPACT csv input file
dfYdShk <- data.frame(Country = all_countries, yShkBest = 0, yShkWrst = 0)
u <- dfYdShk$Country
dfYdShk$yShkBest[which(u %in% groupYshkHi_vec)] <- 1 - yShkHi_best
dfYdShk$yShkBest[which(u %in% groupYshkMe_vec)] <- 1 - yShkMe_best
dfYdShk$yShkBest[which(u %in% groupYshkLo_vec)] <- 1 - yShkLo_best
dfYdShk$yShkWrst[which(u %in% groupYshkHi_vec)] <- 1 - yShkHi_wrst
dfYdShk$yShkWrst[which(u %in% groupYshkMe_vec)] <- 1 - yShkMe_wrst
dfYdShk$yShkWrst[which(u %in% groupYshkLo_vec)] <- 1 - yShkLo_wrst
# Get outbreak frequency table together for IMPACT csv input file
dfFreq <- data.frame(Region = freqRegNames, freq = freqVec)
#------------------------------------------------------------------------
# Now, to merge the ydShk and frequency dfs, need a country-region key
# Can copy paste it from IMPACT ReportGen.xlsx file
# Also, IMPACT requires FPUs not country names
# FPUs can be got from an IMPACT chem fert price elasticity file downloaded from github
# The 3 letter country code that forms the second part of the FPU can be found in the ReportGen.xlsx file
# So the FPU file can be merged to the country-key file via this 3 letter code
# NOTE: The IMPACT 3 letter country code is same as the ISO3 code for many countries
# But it differs from the ISO3 code in many cases
# So don't use the ISO3 code as a guide!
#---
#this_folder <- "E:/BSCHIEK/2Blades WSR IMPACT/New run/"
this_folder <- "D:/OneDrive - CGIAR/Documents 1/CIAT 2/FnM Initiative/WSR 2Blades/"
this_file <- "FAO-IMPACT country-region key.csv"
this_filepath <- paste0(this_folder, this_file)
dfCtyRegKey <- read.csv(this_filepath, stringsAsFactors = F)
# ReportGen.xlsx includes Croatia, but in IMPACT documentation Croatia is part
# of "Other Balkans", so remove it.
# (It will soon be added back in for merge with dfYdShk, but we want the pristine
# IMPACT countries and regions here.)
dfCtyRegKey <- subset(dfCtyRegKey, IMPACTctyFull != "Croatia")
# Keep only relevant columns
dfCtyRegKey <- dfCtyRegKey[, setdiff(colnames(dfCtyRegKey), c("FAO1", "FAO2"))]
# Assign new names to regions
# This is just to merge ydShk and freq dfs
# These are not the impact regions to be analyzed in paper
dfCtyRegKey$Region <- NA
u <- dfCtyRegKey$IMPACTregion
dfCtyRegKey$Region[which(u == "FSU")] <- "Russia & Central Asia"
dfCtyRegKey$Region[which(u == "EUR")] <- "Europe"
dfCtyRegKey$Region[which(u == "NAM")] <- "North America"
dfCtyRegKey$Region[which(u == "SAS")] <- "South Asia"
dfCtyRegKey$Region[which(u == "SSA")] <- "Eastern & Southern Africa"
dfCtyRegKey$Region[which(u == "LAC")] <- "South & Central America"
dfCtyRegKey$Region[which(u == "EAP")] <- "China"
dfCtyRegKey$Region[which(u == "MEN")] <- "West Asia & North Africa"
dfCtyRegKey$Region[which(dfCtyRegKey$IMPACTctyFull == "Australia")] <- "Australia"
# Fine tune region groupings
# Put Ukraine and Belarus in Europe
dfCtyRegKey$Region[which(dfCtyRegKey$IMPACTctyFull == "Ukraine")] <- "Europe"
dfCtyRegKey$Region[which(dfCtyRegKey$IMPACTctyFull == "Belarus")] <- "Europe"
# Move Afghanistan from South Asia to Central Asia
dfCtyRegKey$Region[which(dfCtyRegKey$IMPACTctyFull == "Afghanistan")] <- "Russia & Central Asia"
# Move Iran from West Asia to Central Asia
dfCtyRegKey$Region[which(dfCtyRegKey$IMPACTctyFull == "Iran")] <- "Russia & Central Asia"
# # Put Mexico in North America
# dfCtyRegKey$IMPACT[which(dfCtyRegKey$Country == "Mexico")] <- "NAM"
dfCtyRegKey$IMPACTregion <- NULL
# Add the "Other Balkans" and "Baltic States" countries, which IMPACT does
# not have because it groups them as one
dfAdd <- data.frame(IMPACTcty3letter = NA,
                    IMPACTctyFull = c(otherBalks_vec, baltStates_vec),
                    Region = "Europe")
# The IMPACT 3 letter country codes for "Other Balkans" and "Baltic States"
# are OBN and BLT
dfAdd$IMPACTcty3letter[which(dfAdd$IMPACTctyFull %in% otherBalks_vec)] <- "OBN"
dfAdd$IMPACTcty3letter[which(dfAdd$IMPACTctyFull %in% baltStates_vec)] <- "BLT"
dfCtyRegKey <- as.data.frame(rbind(dfCtyRegKey, dfAdd))
# Rename country column as these are no longer all IMPACT country names
colnames(dfCtyRegKey)[2] <- "Country"
# Now merge with the chem fert price elasticity file by 3 letter code to get FPUs 
this_file <- "PF elast IMPACT.txt"
this_filepath <- paste0(this_folder, this_file)
dfFPU <- read.table(this_filepath, sep = ",", stringsAsFactors = F)
colnames(dfFPU) <- c("crop", "FPU", "sys", "type", "PF")
dfFPU <- dfFPU[-which(duplicated(dfFPU$FPU)), ]
dfFPU$IMPACTcty3letter <- gsub(".*_", "", dfFPU$FPU)
dfFPU <- dfFPU[, c("FPU", "IMPACTcty3letter")]
setdiff(dfCtyRegKey$IMPACTcty3letter, dfFPU$IMPACTcty3letter)
dfCtyRegKey <- merge(dfCtyRegKey, dfFPU)
#------------------------------------------------------------------------
# Now that the country-region-fpu key is sorted out, merge it with the ydShk df by country
# And then merge ydShk with frequency by region
setdiff(dfYdShk$Country, dfCtyRegKey$Country)
dfYdShk <- merge(dfYdShk, dfCtyRegKey)
dfYdShk <- merge(dfYdShk, dfFreq, by = "Region")
#setdiff(unique(dfYdShk$Country), all_countries)
this_file <- "outbreakSpatDist.csv"
this_filepath <- paste0(this_folder, this_file)
dfOutbreak <- read.csv(this_filepath, stringsAsFactors = F)
dfOutbreak$X <- NULL
setdiff(dfOutbreak$Region, dfYdShk$Region)
u <- dfOutbreak$Region
dfOutbreak$Region[grep("South America", u)] <- "South & Central America"
dfOutbreak$Region[grep("East Asia", u)] <- "China"
dfOutbreak$Region[grep("Africa South", u)] <- "Eastern & Southern Africa"
setdiff(dfOutbreak$Region, dfYdShk$Region)
colnames(dfOutbreak) <- gsub("X", "", colnames(dfOutbreak))
df <- merge(dfYdShk, dfOutbreak, by = "Region")
setdiff(df$Country, dfYdShk$Country)
setdiff(dfYdShk$Country, df$Country)
nrow(dfYdShk)
nrow(df)
yrCols <- colnames(df)[8:ncol(df)]
#---
keepCols <- c("FPU", "Country", "yShkBest", yrCols)
dfBest <- df[, keepCols]
dfBest[, yrCols] <- dfBest[, yrCols] * dfBest$yShkBest
dfBest$yShkBest <- NULL
dfBest <- dfBest[-which(duplicated(dfBest$FPU)), ] #Remove duplicate FPUs due to "Other Balkans" and "Baltic States"
dfBest$Country <- NULL
dfBest[dfBest == 0] <- 1
#---
keepCols <- c("FPU", "Country", "yShkWrst", yrCols)
dfWrst <- df[, keepCols]
dfWrst[, yrCols] <- dfWrst[, yrCols] * dfWrst$yShkWrst
dfWrst$yShkWrst <- NULL
dfWrst <- dfWrst[-which(duplicated(dfWrst$FPU)), ] #Remove duplicate FPUs due to "Other Balkans" and "Baltic States"
dfWrst$Country <- NULL
dfWrst[dfWrst == 0] <- 1
#---
setdiff(dfBest$Country, dfYdShk$Country)
setdiff(dfWrst$Country, df$Country)
setdiff(dfYdShk$Country, dfBest$Country)
#---------------------------------------------------------------------
this_file <- "WSR_best2.csv"
this_filepath <- paste0(this_folder, this_file)
write.csv(dfBest, this_filepath, row.names = F)
this_file <- "WSR_worst2.csv"
this_filepath <- paste0(this_folder, this_file)
write.csv(dfWrst, this_filepath, row.names = F)
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
#========================================================================
# Display tables for paper
# WSR Yield shock
dfTab <- data.frame(c("Region", shkHiRegionName_vec, paste(shkMeRegionName_vec, collapse = ", "), paste(shkLoRegionName_vec, collapse = ", ")),
                    c("Best case",
                      yShkHi_best * 100,
                      yShkMe_best * 100,
                      yShkLo_best * 100
                    ),
                    c("Worst case",
                      yShkHi_wrst * 100,
                      yShkMe_wrst * 100,
                      yShkLo_wrst * 100
                    )
)
colnames(dfTab) <- dfTab[1, ]
dfTab <- dfTab[-1, ]
library(flextable)
footText1 <- paste(c(medBasinExclFranceName, paste(medBasinExclFrance_vec, collapse = ", ")), collapse = " = ")
footText2 <- paste(c(sAfrName, paste(sAfr_vec, collapse = ", ")), collapse = " = ")
footText3 <- paste(c(cAsiaVulnName, paste(cAsiaVuln_vec, collapse = ", ")), collapse = " = ")
footText4 <- paste(c(cAsiaResiName, paste(cAsiaResi_vec, collapse = ", ")), collapse = " = ")
footText5 <- paste(c(wNeurExclSpainItalyName, paste(wNeurExclSpainItaly_vec, collapse = ", ")), collapse = " = ")
footText6 <- paste(c(latAmName, paste(latAm_vec, collapse = ", ")), collapse = " = ")
footText7 <- paste(c(eAfrName, paste(eAfr_vec, collapse = ", ")), collapse = " = ")
headerText <- "Table 1: Summary of expert assessment of WSR yield shock (%)"
dfTab %>% regulartable() %>%
  add_header_row(values = headerText, colwidths = c(3)) %>%
  add_footer_row(values = footText1, colwidths = c(3)) %>%
  add_footer_row(values = footText2, colwidths = c(3)) %>%
  add_footer_row(values = footText3, colwidths = c(3)) %>%
  add_footer_row(values = footText4, colwidths = c(3)) %>%
  add_footer_row(values = footText5, colwidths = c(3)) %>%
  add_footer_row(values = footText6, colwidths = c(3)) %>%
  add_footer_row(values = footText7, colwidths = c(3)) %>%
  fontsize(size = 8, part = "footer") %>%
  #line_spacing(space = 0.15, part = "footer") %>%
  autofit()
# Export to csv to then copy paste into Word
# Note have to copy paste the table footnotes manually
this_file <- "WSR impact expert summary table_ydShk.csv"
this_filepath <- paste0(this_folder, this_file)
write.csv(dfTab, this_filepath)
#------------------------------------------------------------------------
# WSR outbreak frequency display table
dfTab <- data.frame(c("Region", freqRegTabNames),
                    c("Outbreak frequency", freqTabVec))
# Export to csv to then copy paste into Word
#this_folder <- "D:/OneDrive - CGIAR/Documents 1/CIAT 2/FnM Initiative/WSR 2Blades/"
this_file <- "WSR impact expert summary table_outbreakInterval.csv"
this_filepath <- paste0(this_folder, this_file)
write.csv(dfTab, this_filepath)
#---
# Flextable table
colnames(dfTab) <- dfTab[1, ]
dfTab <- dfTab[-1, ]
balkanVec <- c("Slovenia", "Albania", "Greece", otherBalks_vec)
WANA_vec <- c(setdiff(medBasinExclFrance_vec, c(balkanVec, wNeurExclSpainItaly_vec, "Spain")), "Iraq")
NAM_vec <- c("Canada", "USA")
EUR_vec <- c(gsub(" \\(excl\\. Spain\\)", "", wNeurExclSpainItalyName), balkanVec)
footText1 <- paste(c(freqRegNames[2], paste(WANA_vec, collapse = ", ")), collapse = " = ")
footText2 <- paste(c("North America", paste(NAM_vec, collapse = ", ")), collapse = " = ")
footText3 <- paste(c("Europe", paste(EUR_vec, collapse = ", ")), collapse = " = ")
headerText <- "Table 2: Summary of expert assessment of WSR outbreak frequency (every x years)"
dfTab %>% regulartable() %>%
  add_header_row(values = headerText, colwidths = c(2)) %>%
  align(j = 2, align = "center", part = "body") %>%
  add_footer_row(values = footText1, colwidths = c(2)) %>%
  add_footer_row(values = footText2, colwidths = c(2)) %>%
  add_footer_row(values = footText3, colwidths = c(2)) %>%
  fontsize(size = 8, part = "footer") %>%
  #line_spacing(space = 0.15, part = "footer") %>%
  autofit()
#===========================================================================