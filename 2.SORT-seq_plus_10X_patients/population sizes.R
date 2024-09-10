#=========================================================================================================================
## Population size comparisons between patients
##========================================================================================================================
## Define cell subsets for ease of use
# Retrieve tissue specific cells
pbmc.cells   <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "PBMC")
plaque.cells <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "plaque")

# Retrieve patient specific cells
P1.cells <- WhichCells(pbmc_plaque.seurat, expression = Patient == "P1")
P2.cells <- WhichCells(pbmc_plaque.seurat, expression = Patient == "P2")
P3.cells <- WhichCells(pbmc_plaque.seurat, expression = Patient == "P3")

# Retrieve patient specific cells per tissue
P1.pbmc.cells   <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "PBMC"   & Patient == "P1")
P2.pbmc.cells   <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "PBMC"   & Patient == "P2")
P3.pbmc.cells   <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "PBMC"   & Patient == "P3")

P1.plaque.cells <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "plaque" & Patient == "P1")
P2.plaque.cells <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "plaque" & Patient == "P2")
P3.plaque.cells <- WhichCells(pbmc_plaque.seurat, expression = Tissue == "plaque" & Patient == "P3")


#=========================================================================================================================
## Population size comparisons between patients
##========================================================================================================================
# Compare cluster size per patient
freq_table <- prop.table(x = table(Idents(pbmc_plaque.seurat), pbmc_plaque.seurat@meta.data[, "Patient"]), margin = 2)
m <- melt(freq_table)
m$Var2 <- as.factor(m$Var2)
ggplot(m, aes(x = Var2, fill = factor(Var1), y = value)) +
  geom_col(position = position_fill(reverse = T)) +
  scale_fill_manual(values = cluster_colours) +
  coord_flip() +
  theme_light() + 
  ylab("% of cells") +
  xlab("Patient") +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.50,0.75,1.00), labels = c(0,25,50,75,100)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(), aspect.ratio = 1:1
  )
ggsave("results/freq patient.pdf",width = 10, height = 10)


#Plot number of cells per category
#Get the number of cells per category
cell.numbers <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i))
}

# SPlit per patient
cell.numbers.p1 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.p1[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P1"))
}

cell.numbers.p2 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.p2[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P2"))
}

cell.numbers.p3 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.p3[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P3"))
}

# Split per tissue
cell.numbers.plaque <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.plaque[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Tissue == "plaque"))
}

cell.numbers.pbmc <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.pbmc[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Tissue == "PBMC"))
}

# Split per patient and tissue
cell.numbers.pbmc.p1 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.pbmc.p1[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P1" & Tissue == "PBMC"))
}

cell.numbers.pbmc.p2 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.pbmc.p2[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P2" & Tissue == "PBMC"))
}

cell.numbers.pbmc.p3 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.pbmc.p3[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P3" & Tissue == "PBMC"))
}

cell.numbers.plaque.p1 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.plaque.p1[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P1" & Tissue == "plaque"))
}

cell.numbers.plaque.p2 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.plaque.p2[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P2" & Tissue == "plaque"))
}

cell.numbers.plaque.p3 <- list()
for(i in unique(pbmc_plaque.seurat$predicted.celltype.l1)){
  cell.numbers.plaque.p3[i] <- length(WhichCells(pbmc_plaque.seurat, expression = predicted.celltype.l1 == i & Patient == "P3" & Tissue == "plaque"))
}

#Get the total number of cells
total.cells <- ncol(pbmc_plaque.seurat)

#Sanity check
total.cells == Reduce("+", cell.numbers)

#Plot bars for number or percentage of cells
d <- melt(cell.numbers)
colnames(d) <- c("Number","Category")
d$Percentage <- (d$Number / total.cells) * 100

#Sanity check
Reduce("+", d$Percentage)

#Plot relative values with numbers
ggplot(d, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "Dodge") +
  geom_text(aes(label=round(Percentage, 1)), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_x_discrete(limits = d$Category[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), 
                     name = "% of cells", 
                     breaks = c(0,10,20,30,40,50,60), 
                     limits = c(00,60),
                     sec.axis = sec_axis(trans = ~ (. / 100) * total.cells, name = "# of cells", breaks = c(0,500,1000,2500,5000,10000,15000,20000))) +
  theme_pubr(legend = "right") +
  ggtitle("Percentage and number of cells per Category") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_L1_category_relative_show_values.pdf")


# Plot bars for number or percentage of cells
# Split per patient
d.1 <- melt(cell.numbers.p1)
d.2 <- melt(cell.numbers.p2)
d.3 <- melt(cell.numbers.p3)

d.1$Patient <- "P1"
d.2$Patient <- "P2"
d.3$Patient <- "P3"

colnames(d.1) <- c("Number","Category", "Patient")
colnames(d.2) <- c("Number","Category", "Patient")
colnames(d.3) <- c("Number","Category", "Patient")

d.1$Percentage <- (d.1$Number / length(P1.cells)) * 100
d.2$Percentage <- (d.2$Number / length(P2.cells)) * 100
d.3$Percentage <- (d.3$Number / length(P3.cells)) * 100

d <- rbind(rbind(d.1, d.2), d.3)

#Sanity check
Reduce("+", d$Percentage)

#Plot relative values with numbers
ggplot(d, aes(x = Category, y = Percentage, fill = Patient, group = Patient)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_x_discrete(limits = d$Category[order(d.1$Number, decreasing = T)]) +
  geom_text(aes(label=round(Percentage, 1)), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_y_continuous(expand = c(0,0), 
                     name = "% of cells", 
                     breaks = c(0,10,20,30,40,50,60,70), 
                     limits = c(00,70),
                     sec.axis = sec_axis(trans = ~ (. / 100) * total.cells, name = "# of cells", breaks = c(0,500,1000,2500,5000,10000,15000,20000))) +
  theme_pubr(legend = "right") +
  ggtitle("Percentage and number of cells per Category") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_patient_per_L1_category_relative_show_values.pdf")

# Split per tissue
d.pbmc   <- melt(cell.numbers.pbmc)
d.plaque <- melt(cell.numbers.plaque)

d.pbmc$Tissue   <- "PBMC"
d.plaque$Tissue <- "Plaque"

colnames(d.pbmc)   <- c("Number","Category", "Tissue")
colnames(d.plaque) <- c("Number","Category", "Tissue")

d.pbmc$Percentage   <- (d.pbmc$Number / length(pbmc.cells)) * 100
d.plaque$Percentage <- (d.plaque$Number / length(plaque.cells)) * 100

d <- rbind(d.pbmc, d.plaque)

#Sanity check
Reduce("+", d$Percentage)

#Plot relative values with numbers
ggplot(d, aes(x = Category, y = Percentage, fill = Tissue, group = Tissue)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_x_discrete(limits = d$Category[order(d.pbmc$Number, decreasing = T)]) +
  geom_text(aes(label=round(Percentage, 1)), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_y_continuous(expand = c(0,0), 
                     name = "% of cells", 
                     breaks = c(0,10,20,30,40,50,60,60), 
                     limits = c(00,60),
                     sec.axis = sec_axis(trans = ~ (. / 100) * total.cells, name = "# of cells", breaks = c(0,500,1000,2500,5000,10000,15000,20000))) +
  theme_pubr(legend = "right") +
  ggtitle("Percentage and number of cells per Category") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_tissue_per_L1_category_relative_show_values.pdf")


# Split per patient and tissue
# PBMC
d.1 <- melt(cell.numbers.pbmc.p1)
d.2 <- melt(cell.numbers.pbmc.p2)
d.3 <- melt(cell.numbers.pbmc.p3)

d.1$Patient <- "P1"
d.2$Patient <- "P2"
d.3$Patient <- "P3"

colnames(d.1) <- c("Number","Category", "Patient")
colnames(d.2) <- c("Number","Category", "Patient")
colnames(d.3) <- c("Number","Category", "Patient")

d.1$Percentage <- (d.1$Number / length(P1.pbmc.cells)) * 100
d.2$Percentage <- (d.2$Number / length(P2.pbmc.cells)) * 100
d.3$Percentage <- (d.3$Number / length(P3.pbmc.cells)) * 100

d <- rbind(rbind(d.1, d.2), d.3)

#Sanity check
Reduce("+", d$Percentage)

#Plot relative values with numbers
ggplot(d, aes(x = Category, y = Percentage, fill = Patient, group = Patient)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_x_discrete(limits = d$Category[order(d.1$Number, decreasing = T)]) +
  geom_text(aes(label=round(Percentage, 1)), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_y_continuous(expand = c(0,0), 
                     name = "% of cells", 
                     breaks = c(0,10,20,30,40,50,60,70), 
                     limits = c(00,70),
                     sec.axis = sec_axis(trans = ~ (. / 100) * total.cells, name = "# of cells", breaks = c(0,500,1000,2500,5000,10000,15000,20000))) +
  theme_pubr(legend = "right") +
  ggtitle("Percentage and number of cells per Category") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_patient_in_PBMCs_per_L1_category_relative_show_values.pdf")


# Plaque
d.1 <- melt(cell.numbers.plaque.p1)
d.2 <- melt(cell.numbers.plaque.p2)
d.3 <- melt(cell.numbers.plaque.p3)

d.1$Patient <- "P1"
d.2$Patient <- "P2"
d.3$Patient <- "P3"

colnames(d.1) <- c("Number","Category", "Patient")
colnames(d.2) <- c("Number","Category", "Patient")
colnames(d.3) <- c("Number","Category", "Patient")

d.1$Percentage <- (d.1$Number / length(P1.plaque.cells)) * 100
d.2$Percentage <- (d.2$Number / length(P2.plaque.cells)) * 100
d.3$Percentage <- (d.3$Number / length(P3.plaque.cells)) * 100

d <- rbind(rbind(d.1, d.2), d.3)

#Sanity check
Reduce("+", d$Percentage)

#Plot relative values with numbers
ggplot(d, aes(x = Category, y = Percentage, fill = Patient, group = Patient)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_x_discrete(limits = d$Category[order(d.1$Number, decreasing = T)]) +
  geom_text(aes(label=round(Percentage, 1)), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_y_continuous(expand = c(0,0), 
                     name = "% of cells", 
                     breaks = c(0,10,20,30,40,50,60,70), 
                     limits = c(00,70),
                     sec.axis = sec_axis(trans = ~ (. / 100) * total.cells, name = "# of cells", breaks = c(0,500,1000,2500,5000,10000,15000,20000))) +
  theme_pubr(legend = "right") +
  ggtitle("Percentage and number of cells per Category") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_patient_in_plaques_per_L1_category_relative_show_values.pdf")
