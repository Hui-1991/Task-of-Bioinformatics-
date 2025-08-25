# Load required libraries
library(reshape2)
library(ggplot2)
library(tidyr)
library(extrafont)

# my theme 
my_theme = theme(
  strip.text = element_text(size = 20, face = "bold"),
  plot.title = element_text(size = 36, family = "Arial", face = 'bold'),
  axis.title.x = element_text(size = 36, family = "Arial", face = 'bold', margin = margin(t = 10)),
  axis.title.y = element_text(size = 36, family = "Arial", face = 'bold', margin = margin(r = 10)),
  axis.text.x = element_text(size = 18, family = "Arial", colour = 'black', face = 'bold'),
  axis.text.y = element_text(size = 18, family = "Arial", colour = 'black', face = 'bold'),
  legend.text = element_text(size = 20, family = "Arial", face = 'bold'),
  legend.key = element_rect(fill = 'white', color = 'white'),
  panel.background = element_rect(fill = 'white', colour = 'black')
)

##########################################################################
#                       Mexicana: Boxplots                               #
##########################################################################

# Load datasets
mex = read.table('Kennedy_mexicana.csv', header = TRUE, sep = ',',check.names = FALSE)

# Convert to long format (excluding sample and label columns)
mex_long = pivot_longer(
  mex,
  cols = -label,
  names_to = "Metabolite",
  values_to = "Value"
)

# Set label order: No_Drug appears first
mex_long$label = factor(mex_long$label, levels = c('No_Drug', 'DNDi5426'))

# Create boxplot for all metabolites
box_mex = ggplot(mex_long, aes(x = label, y = Value, fill = label)) +
  geom_boxplot() +
  facet_wrap(~ Metabolite, scales = "free_y") +
  scale_fill_manual(values = c("No_Drug" = "#525151ff",
                               "DNDi5426" = "#945ac7ff"))+
  labs(
    title = "Mexicana: No_Drug vs DNDi5426",
    x = "Group",
    y = "Intensity"
  ) +
  my_theme


# show and save 
box_mex
ggsave("box_mex.png", 
       plot = box_mex, width = 28,
       height = 15, dpi = 300)


##########################################################################
#                         Major: Boxplots                                #
##########################################################################
  
maj = read.csv("Kennedy_major.csv", header = TRUE, sep = ",", check.names = FALSE)


# Convert to long format (excluding sample and label columns)
maj_long = pivot_longer(
  maj,
  cols = 2:ncol(maj),
  names_to = "Metabolite",
  values_to = "Value"
)


# Set label order: No_Drug appears first
maj_long$label = factor(maj_long$label, levels = c('No_Drug', 'DNDi5426'))

# Create boxplot for all metabolites
box_maj = ggplot(maj_long, aes(x = label, y = Value, fill = label)) +
  geom_boxplot() +
  facet_wrap(~ Metabolite, scales = "free_y") +
  scale_fill_manual(values = c("No_Drug" = "#525151ff", 
                               "DNDi5426" = "#6495ED"))+
  labs(
    title = "Mexicana: No_Drug vs DNDi5426",
    x = "Group",
    y = "Intensity"
  ) +
  my_theme

# show and save 
box_maj
ggsave("box_maj.png", 
       plot = box_maj, width = 28,
       height = 15, dpi = 300)



##########################################################################
#                Flod-change & combined Barplots                         #
##########################################################################

# -------------------------- Mexicana subset --------------------------- #

# subset for mexicana 
sub_mex = mex[, 1:6]
data_mex = sub_mex[,2:ncol(sub_mex)]

#Seperate groups
nodrug_mex = data_mex[sub_mex$label == "No_Drug", ]
drug_mex   = data_mex[sub_mex$label == "DNDi5426", ]

# Calculate mean value for each metabolite of each group
mean_nodrug_mex = colMeans(nodrug)
mean_drug_mex = colMeans(drug)

#Calculate fold change = DNDi5426/No_drug
ratio = mean_drug/mean_nodrug

# Create a new dataframe for foldchange
mex_df = data.frame(
  Metabolite = names(ratio),
  DNDi5426_L.mexicana = as.numeric(ratio)
)


# -------------------------- Major subset --------------------------- #

# Subset for major 
sub_maj = maj[, 1:7]
data_maj = sub_maj[,2:ncol(sub_maj)]

# Seperate groups
nodrug_maj = data_maj[sub_maj$label=="No_Drug",]
drug_maj = data_maj[sub_maj$label=="DNDi5426",]

# Calculate mean values
mean_nodrug_maj = colMeans(nodrug_maj)
mean_drug_maj = colMeans(drug_maj)

# Fold change
ratio_maj = mean_drug_maj/mean_nodrug_maj

# Create a new dataframe for foldchange
maj_df = data.frame(
  Metabolite = names(ratio_maj),
  DNDi5426_L.major = as.numeric(ratio_maj)
)

# -------------Combine Mexicana and Major dataframe ------------------- #

# Merge two dataframe
combind = merge(mex_df,maj_df, by = "Metabolite", all=TRUE)

# Handle missing value for CDP-ehthanolamine 
combind[is.na(combind$DNDi5426_L.mexicana), "DNDi5426_L.mexicana"] = 2

# Convert to long format for plotting
long_combind = pivot_longer(
  combind, 
  cols = c("DNDi5426_L.mexicana","DNDi5426_L.major"),
  names_to = "Type",
  values_to = "Foldchange"
  )


# set value for No_drug 
baseline = data.frame(
  Metabolite = unique(long_combind$Metabolite),
  Type = "No_Drug",
  FoldChange = 1
)

# combine all dataframe
all = rbind(long_combind,baseline)
  
# set order 
all$Type = factor(
  all$Type, 
  levels = c("No_Drug","DNDi5426_L.mexicana","DNDi5426_L.major"))


# Create bar plots 
combine_bar = ggplot(all, aes(x=Type,y=Foldchange,fill = Type))+
  geom_bar(stat = "identity")+
   facet_wrap(~Metabolite, nrow =1, scales = "free_y")+
  labs(title = "Fold Change per Metabolite (with No_Drug baseline)",
       x = "Group",
       y = "Fold Change") +
  scale_fill_manual(values = c(
    "No_Drug" = "#525151ff",
    "DNDi5426_L.mexicana" = "#945ac7ff",
    "DNDi5426_L.major" = "#6495ED"
  )) +
  my_theme +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

# display
combine_bar

# save plots
ggsave("combine_bar2.png", 
       plot = combine_bar, width = 25,
       height = 10, dpi = 300) 
  

