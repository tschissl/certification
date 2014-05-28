####################
# Merging data for OKUM, GAS and MUH-1 certification project
# updated 20140519
# Thomas Meisel
####################

GOM <- data.frame(OKUM, MUH[6:60], GAS[6:60]) # colums 1:5 are identical for all RM and contain the factors
GOM <- merge(GOM, OKUM.methods) # merging all data and all methods
#defining factors 
GOM$Lab <- as.factor(GOM$Lab)
GOM$Packet <- as.factor(GOM$Packet)
GOM$Prep <- as.factor(GOM$Prep)
GOM$day <- as.factor(GOM$day)
GOM$names <- as.factor(GOM$names)
GOMorig <- GOM  #keeping the original file