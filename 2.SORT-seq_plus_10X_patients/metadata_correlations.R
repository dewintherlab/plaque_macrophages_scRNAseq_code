#=========================================================================================================================
## Analyze correlations with metadata
##========================================================================================================================
dir.create("metadata_correlations/Individual trait UMAPs", showWarnings = F, recursive = T)
dir.create("metadata_correlations/Individual trait distribution", showWarnings = F, recursive = T)
dir.create("metadata_correlations/Individual trait UMAPs binned", showWarnings = F, recursive = T)
dir.create("metadata_correlations/Individual trait distribution binned", showWarnings = F, recursive = T)


## Process metadata
# Read the file
meta.data <- read.xlsx(file = "raw_data/2021-09-16 AtheroExpress Database Marie.xlsx", sheetIndex = 1, header = T, as.data.frame = T)
dim(meta.data)

# Fix NA values
# Work-around weird 'charToDate()' error when checking if a string is NA while there is a date like structure in the string by forcing those columns to 'factor'
for(theCol in colnames(meta.data)[grep("date", colnames(meta.data), ignore.case = T)]){
  meta.data[, theCol] <- factor(meta.data[, theCol])
}
meta.data[meta.data == "NA"] <- NA

# Filter on quantity of NAs (keep only columns with less than 10% NA)
max.na    <- floor(nrow(meta.data) * 0.10)
meta.data <- meta.data[,colSums(apply(meta.data, 1:2, is.na)) <= max.na]
dim(meta.data)

# Filter out columns with the same value for all rows
meta.data <- meta.data[,apply(meta.data, 2, function(x)length(unique(x))) > 1]
dim(meta.data)

# Filter out columns with the same value or NA for all rows
meta.data <- meta.data[,!(apply(meta.data, 2, function(x)length(unique(x))) == 2 & colSums(apply(meta.data, 1:2, is.na)) > 0)]
dim(meta.data)

## Add metadata to 43p seurat object (as we don't have htis info for the 10X guys)
# Fetch patient information per cell
md.df <- data.frame(Patient=mye.43p.seurat$Patient)

# Expand metadata to all cells
md.df <- merge(md.df, meta.data, by.x = 1, by.y = 2, all.x = T)
md.df$Patient <- NULL

# And add to the Seurat object
for (i in names(md.df)){
  mye.43p.seurat <- AddMetaData(mye.43p.seurat, md.df[,i], col.name = i)
}

##=========================================================================================================================
## Visualise
# Plot individual UMAPs per trait
i     <- 1
i.max <- length(colnames(meta.data[-2])) #-2 because we merged on that column
for(theTrait in colnames(meta.data)[-2]){
  cat(paste("Plotting ", theTrait, " (", i, " of ", i.max, ")\n", sep = ""))
  i <- i + 1
  customUMAP(object    = mye.43p.seurat, 
             group.by  = theTrait, 
             pt.size   = 4, 
             label     = F, 
             title     = theTrait,
             shuffle   = T, 
             file.name = paste("metadata_correlations/Individual trait UMAPs/UMAP ", theTrait, ".pdf"))
}



## Plot individual relative distributions per trait
# Get distributions per trait, per population
trait.dist <- list()
i          <- 1
i.max      <- length(colnames(meta.data[-2])) #-2 because we merged on that column
for(theTrait in colnames(meta.data[-2])){
  cat(paste("Working on trait: ", theTrait, " (", i, " of ", i.max, ")\n", sep =""))
  
  # make a df in the list and populate with the distribution for this trait over all cells
  trait.dist[[theTrait]]            <- data.frame(table(mye.43p.seurat@meta.data[,theTrait]))
  row.names(trait.dist[[theTrait]]) <-  trait.dist[[theTrait]]$Var1
  trait.dist[[theTrait]]$Var1       <- NULL
  colnames(trait.dist[[theTrait]])  <- "Total"
  
  # Then add the distribution per population
  for(thePop in unique(mye.43p.seurat$method.int.idents)){
    cat(paste("\tAdding distribution for pop: ", thePop, "...\n", sep =""))
    
    # Subset out the desired population and get the distribution
    tmp.pop.seurat        <- subset(mye.43p.seurat, method.int.idents == thePop)
    tmp.pop.df            <- data.frame(table(tmp.pop.seurat@meta.data[,theTrait]))
    row.names(tmp.pop.df) <- tmp.pop.df$Var1
    tmp.pop.df$Var1       <- NULL
    colnames(tmp.pop.df)  <- thePop
    
    ## Add this pop's distribution to the overall df in the list
    # We need to merge instead of cbind just in case a pop doesn't have a certain category of the trait.
    trait.dist[[theTrait]]            <- merge(trait.dist[[theTrait]], tmp.pop.df, by.x = 0, by.y = 0, all.x = T)
    row.names(trait.dist[[theTrait]]) <- trait.dist[[theTrait]]$Row.names
    trait.dist[[theTrait]]$Row.names  <- NULL
  }
  
  # Get rid of NAs introduced by trait cats not in certain populations
  trait.dist[[theTrait]][is.na(trait.dist[[theTrait]])] <- 0
  
  i <- i + 1
}
head(trait.dist)

# Clean up: remove traits with one observation (shouldn't be there to begin with)
for(theTrait in names(trait.dist)){
  if(dim(trait.dist[[theTrait]])[1] < 2){
    cat(paste("Removed trait: ", theTrait, "\n", sep = ""))
    trait.dist[[theTrait]] <- NULL
  }
}

## Compute statistical significance
# Run fisher test per trait. Compare each population to the total.
sig.trait.dist <- data.frame()
for(theTrait in names(trait.dist)){
  cat(paste("Working on", theTrait, "...\n", sep = " "))
  
  # First populate a new df with the first comparison (Total vs. column 2)
  tmp.trait.df           <- data.frame(fisher.test(x = trait.dist[[theTrait]][,c(1,2)], simulate.p.value = T)$p.value, row.names = theTrait)
  colnames(tmp.trait.df) <- colnames(trait.dist[[theTrait]])[2]
  
  # Iterate over the rest of the comparisons (Total vs. column n)
  for(thePop in colnames(trait.dist[[theTrait]])[c(-1, -2)]){
    # Make a fresh df so we can more easily manipulate the names
    tmp.pop.df           <- data.frame(fisher.test(x = cbind(trait.dist[[theTrait]][,1], trait.dist[[theTrait]][,thePop]), simulate.p.value = T)$p.value, row.names = theTrait)
    colnames(tmp.pop.df) <- colnames(trait.dist[[theTrait]][, thePop, drop = F])
    
    # Add back to the tmp trait df
    tmp.trait.df <- cbind(tmp.trait.df, tmp.pop.df)
  }

  # Add back to the main sig trait df
  if(nrow(sig.trait.dist) > 0 ){
    sig.trait.dist <- rbind(sig.trait.dist, tmp.trait.df)
  } else{
    sig.trait.dist <- tmp.trait.df
  }
}
tail(sig.trait.dist)

# Adjust p values for the amount of tests we ran we ran within each trait (number of populations)
padj.sig.trait.dist <- apply(sig.trait.dist, c(1,2), function(x){p.adjust(p = x, method = "BH", n = length(levels(mye.43p.seurat$method.int.idents)))})

# anything significant?
sum(padj.sig.trait.dist <= 0.05)

# Which are significant?
padj.sig.trait.dist[rowSums(padj.sig.trait.dist <= 0.05) > 0,]

## Now compare traits as a whole, not on the population level)
# Run fisher test per trait.
per.trait.sig.trait.dist <- data.frame()
for(theTrait in names(trait.dist)){
  cat(paste("Working on", theTrait, "...\n", sep = " "))
  tmp.trait.df <- data.frame(fisher.test(x = trait.dist[[theTrait]][,-1], simulate.p.value = T)$p.value, row.names = theTrait)
  colnames(tmp.trait.df) <- "pvalue"
  
  # Add back to the main per trait sig trait df
  if(nrow(per.trait.sig.trait.dist) > 0 ){
    per.trait.sig.trait.dist <- rbind(per.trait.sig.trait.dist, tmp.trait.df)
  } else{
    per.trait.sig.trait.dist <- tmp.trait.df
  }
}
tail(per.trait.sig.trait.dist)

# Adjust p values for the amount of tests we ran ( traits )
padj.per.trait.sig.trait.dist <- p.adjust(p = per.trait.sig.trait.dist$pvalue, method = "BH", n = nrow(per.trait.sig.trait.dist))
names(padj.per.trait.sig.trait.dist) <- row.names(per.trait.sig.trait.dist)

# anything significant?
sum(padj.per.trait.sig.trait.dist <= 0.05)

# Which are significant?
names(padj.per.trait.sig.trait.dist[padj.per.trait.sig.trait.dist < 0.05])


## Make a bar plot per trait dist
for(theTrait in names(trait.dist)){
  m            <- trait.dist[[theTrait]]
  m            <- as.data.frame(t(apply(m, 1, function(x){x/colSums(m)})))
  m[,theTrait] <- row.names(m)
  m            <- melt(m)
  colnames(m)  <- c(theTrait, "Population", "Frequency")
  m$Population <- factor(m$Population, levels = c("Total",names(M.int_refined.pop.colors[c(-1,-2,-3)])))
  levels(m$Population)
  label_width <- 1 + max(strwidth(m$Population, units = "inches"))
  
  if(theTrait %in% row.names(padj.sig.trait.dist[rowSums(padj.sig.trait.dist <= 0.05) > 0,])){ # If there are sig diff populations, plot significance
    stat.test            <- melt(padj.sig.trait.dist[theTrait,])
    stat.test$group1     <- "Total"
    stat.test$group2     <- row.names(stat.test)
    stat.test            <- stat.test[,c(2,3,1)]
    colnames(stat.test)  <- c("group1", "group2", "p.adj")
    stat.test            <- stat.test[stat.test$p <= 0.05,]
    row.names(stat.test) <- NULL
    stat.test            <- stat.test %>% add_significance()
    #stat.test           <- stat.test %>% add_y_position(fun = "max", formula = Frequency ~ Population, data = m) 
    # y pos is always the same, since we have a stacked bar
    stat.test$y.position <- 1.25

  ggplot(m, aes(x = Population, y = Frequency)) + 
      geom_bar(position = position_stack(), stat = "identity", mapping = aes_string(fill = theTrait)) +
      coord_flip() +
      scale_fill_startrek() + 
      theme_pubr(base_size = 14) +
      scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0", "25%", "50%", "75%", "100%"),oob = scales::oob_keep) +
      annotate(x=0, xend=0, y=0, yend=1, colour="black", lwd=1, geom="segment") +
      stat_pvalue_manual(coord.flip = T, remove.bracket = T, data = stat.test, x = "group2", label = "p.adj.signif", size = 6, vjust = 0.75) +
      ggtitle(label = theTrait, subtitle = paste("Fisher exact FDR:", round(padj.per.trait.sig.trait.dist[theTrait], digits = 3), sep = " ")) +
      theme(axis.text.y  = element_text(angle = 0, size = 14, colour = M.int_refined.pop.colors[c(-1,-2)]),
            axis.text.x  = element_text(size = 10),
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank(),
            axis.line    = element_blank(),
            legend.title = element_blank(),
            aspect.ratio = 1.2)
    ggsave(paste("metadata_correlations/Individual trait distribution/", theTrait, ".pdf", sep = ""), width = label_width + 4, height = 4)
    
  }else{ # No sig diff populations
  ggplot(m, aes(x = Population, y = Frequency)) + 
    geom_bar(position = position_stack(), stat = "identity", mapping = aes_string(fill = theTrait)) + 
    coord_flip() +
    scale_fill_startrek() + 
    theme_pubr(base_size = 14) + 
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0", "25%", "50%", "75%", "100%"),oob = scales::oob_keep) +
    annotate(x=0, xend=0, y=0, yend=1, colour="black", lwd=1, geom="segment") +
    ggtitle(label = theTrait, subtitle = paste("Fisher exact FDR:", round(padj.per.trait.sig.trait.dist[theTrait], digits = 3), sep = " ")) +
    theme(axis.text.y = element_text(angle = 0, size = 14, colour = M.int_refined.pop.colors[c(-1,-2)]),
          axis.text.x  = element_text(size = 10),
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.line    = element_blank(), 
          legend.title = element_blank(),
          aspect.ratio = 1.2)
  ggsave(paste("metadata_correlations/Individual trait distribution/", theTrait, ".pdf", sep = ""), width = label_width + 4, height = 4)
  }
}


##=========================================================================================================================
## Let's bin the traits with more than 10 options
binned.meta.data <- meta.data[,1:2]
for(theTrait in colnames(meta.data[,c(-1,-2)])){
  # Keep only traits with over 10 levels
  if(length(levels(factor(meta.data[,theTrait]))) > 10){
    # Keep only numeric traits
    if(check.numeric(meta.data[1,theTrait])){
      cat(paste("Binning:", theTrait,"\n"))
     
      # Define the quartiles of this trait
      q <- quantile(as.numeric(as.character(meta.data[, theTrait])), na.rm = T)
      q.levels <-             paste("n <", round(q[2], digits = 0))
      q.levels <- c(q.levels, paste(round(q[2], digits = 0), "< n <", round(q[3], digits = 0)))
      q.levels <- c(q.levels, paste(round(q[3], digits = 0), "< n <", round(q[4], digits = 0)))
      q.levels <- c(q.levels, paste("n >", round(q[4], digits = 0)))
      
      # Keep a working copy of the trait were we replace NA values
      tmp.df                <- meta.data[,theTrait, drop = F]
      tmp.df[is.na(tmp.df)] <- 0
      
      # Build a results data frame, case by case
      results.df            <- tmp.df
      results.df[,theTrait] <- as.character(results.df[,theTrait])
        
      results.df[as.numeric(as.character(tmp.df[,theTrait]))  < q[2]                                                     , theTrait] <- q.levels[1]
      results.df[as.numeric(as.character(tmp.df[,theTrait])) >= q[2] & as.numeric(as.character(tmp.df[,theTrait])) < q[3], theTrait] <- q.levels[2]
      results.df[as.numeric(as.character(tmp.df[,theTrait])) >= q[3] & as.numeric(as.character(tmp.df[,theTrait])) < q[4], theTrait] <- q.levels[3]
      results.df[as.numeric(as.character(tmp.df[,theTrait])) >= q[4]                                                     , theTrait] <- q.levels[4]
      
      results.df[,theTrait] <- factor(results.df[,theTrait], levels = q.levels)
      
      # Add results back to the main df
      binned.meta.data <- cbind(binned.meta.data, results.df)
    }
  }
}

## Add binned metadata to the seurat object
colnames(binned.meta.data) <- paste0(colnames(binned.meta.data), ".binned", sep = "")
head(binned.meta.data)

# Fetch patient information per cell
md.df <- data.frame(Patient=mye.43p.seurat$Patient)

# Expand metadata to all cells
md.df <- merge(md.df, binned.meta.data, by.x = 1, by.y = 2, all.x = T)
md.df$Patient <- NULL

# And add to the Seurat object
for (i in names(md.df)){
  mye.43p.seurat <- AddMetaData(mye.43p.seurat, md.df[,i], col.name = i)
}


## Visualise
# Plot individual UMAPs per trait
i     <- 1
i.max <- length(colnames(binned.meta.data[-2])) #-2 because we merged on that column
for(theTrait in colnames(binned.meta.data)[-2]){
  cat(paste("Plotting ", theTrait, " (", i, " of ", i.max, ")\n", sep = ""))
  i <- i + 1
  customUMAP(object    = mye.43p.seurat, 
             group.by  = theTrait, 
             pt.size   = 4, 
             label     = F, 
             title     = theTrait,
             shuffle   = T, 
             file.name = paste("metadata_correlations/Individual trait UMAPs binned/UMAP ", theTrait, ".pdf"))
}


## Plot individual relative distributions per trait
# Get distributions per trait, per population
trait.dist <- list()
i          <- 1
i.max      <- length(colnames(binned.meta.data[-2])) #-2 because we merged on that column
for(theTrait in colnames(binned.meta.data[-2])){
  cat(paste("Working on trait: ", theTrait, " (", i, " of ", i.max, ")\n", sep =""))
  
  # make a df in the list and populate with the distribution for this trait over all cells
  trait.dist[[theTrait]]            <- data.frame(table(mye.43p.seurat@meta.data[,theTrait]))
  row.names(trait.dist[[theTrait]]) <-  trait.dist[[theTrait]]$Var1
  trait.dist[[theTrait]]$Var1       <- NULL
  colnames(trait.dist[[theTrait]])  <- "Total"
  
  # Then add the distribution per population
  for(thePop in unique(mye.43p.seurat$method.int.idents)){
    cat(paste("\tAdding distribution for pop: ", thePop, "...\n", sep =""))
    
    # Subset out the desired population and get the distribution
    tmp.pop.seurat        <- subset(mye.43p.seurat, method.int.idents == thePop)
    tmp.pop.df            <- data.frame(table(tmp.pop.seurat@meta.data[,theTrait]))
    row.names(tmp.pop.df) <- tmp.pop.df$Var1
    tmp.pop.df$Var1       <- NULL
    colnames(tmp.pop.df)  <- thePop
    
    ## Add this pop's distribution to the overall df in the list
    # We need to merge instead of cbind just in case a pop doesn't have a certain category of the trait.
    trait.dist[[theTrait]]            <- merge(trait.dist[[theTrait]], tmp.pop.df, by.x = 0, by.y = 0, all.x = T)
    row.names(trait.dist[[theTrait]]) <- trait.dist[[theTrait]]$Row.names
    trait.dist[[theTrait]]$Row.names  <- NULL
  }
  
  # Get rid of NAs introduced by trait cats not in certain populations
  trait.dist[[theTrait]][is.na(trait.dist[[theTrait]])] <- 0
  
  i <- i + 1
}
head(trait.dist)


## Compute statistical significance
# Run fisher test per trait. Compare each population to the total.
sig.trait.dist <- data.frame()
for(theTrait in names(trait.dist)){
  cat(paste("Working on", theTrait, "...\n", sep = " "))
  
  # First populate a new df with the first comparison (Total vs. column 2)
  tmp.trait.df           <- data.frame(fisher.test(x = trait.dist[[theTrait]][,c(1,2)], simulate.p.value = T)$p.value, row.names = theTrait)
  colnames(tmp.trait.df) <- colnames(trait.dist[[theTrait]])[2]
  
  # Iterate over the rest of the comparisons (Total vs. column n)
  for(thePop in colnames(trait.dist[[theTrait]])[c(-1, -2)]){
    # Make a fresh df so we can more easily manipulate the names
    tmp.pop.df           <- data.frame(fisher.test(x = cbind(trait.dist[[theTrait]][,1], trait.dist[[theTrait]][,thePop]), simulate.p.value = T)$p.value, row.names = theTrait)
    colnames(tmp.pop.df) <- colnames(trait.dist[[theTrait]][, thePop, drop = F])
    
    # Add back to the tmp trait df
    tmp.trait.df <- cbind(tmp.trait.df, tmp.pop.df)
  }
  
  # Add back to the main sig trait df
  if(nrow(sig.trait.dist) > 0 ){
    sig.trait.dist <- rbind(sig.trait.dist, tmp.trait.df)
  } else{
    sig.trait.dist <- tmp.trait.df
  }
}
tail(sig.trait.dist)

# Adjust p values for the amount of tests we ran we ran within each trait (number of populations)
padj.sig.trait.dist <- apply(sig.trait.dist, c(1,2), function(x){p.adjust(p = x, method = "BH", n = length(levels(mye.43p.seurat$method.int.idents)))})

# anything significant?
sum(padj.sig.trait.dist <= 0.05)

# Which are significant?
padj.sig.trait.dist[rowSums(padj.sig.trait.dist <= 0.05) > 0,]

## Now compare traits as a whole, not on the population level)
# Run fisher test per trait.
per.trait.sig.trait.dist <- data.frame()
for(theTrait in names(trait.dist)){
  cat(paste("Working on", theTrait, "...\n", sep = " "))
  tmp.trait.df <- data.frame(fisher.test(x = trait.dist[[theTrait]][,-1], simulate.p.value = T)$p.value, row.names = theTrait)
  colnames(tmp.trait.df) <- "pvalue"
  
  # Add back to the main per trait sig trait df
  if(nrow(per.trait.sig.trait.dist) > 0 ){
    per.trait.sig.trait.dist <- rbind(per.trait.sig.trait.dist, tmp.trait.df)
  } else{
    per.trait.sig.trait.dist <- tmp.trait.df
  }
}
tail(per.trait.sig.trait.dist)

# Adjust p values for the amount of tests we ran ( traits )
padj.per.trait.sig.trait.dist <- p.adjust(p = per.trait.sig.trait.dist$pvalue, method = "BH", n = nrow(per.trait.sig.trait.dist))
names(padj.per.trait.sig.trait.dist) <- row.names(per.trait.sig.trait.dist)

# anything significant?
sum(padj.per.trait.sig.trait.dist <= 0.05)

# Which are significant?
names(padj.per.trait.sig.trait.dist[padj.per.trait.sig.trait.dist < 0.05])


## Make a bar plot per trait dist
for(theTrait in names(trait.dist)){
  m            <- trait.dist[[theTrait]]
  m            <- as.data.frame(t(apply(m, 1, function(x){x/colSums(m)})))
  m[,theTrait] <- row.names(m)
  m            <- melt(m)
  colnames(m)  <- c(theTrait, "Population", "Frequency")
  m$Population <- factor(m$Population, levels = c("Total",names(M.int_refined.pop.colors[c(-1,-2,-3)])))
  levels(m$Population)
  label_width <- 1 + max(strwidth(m$Population, units = "inches"))
  
  if(theTrait %in% row.names(padj.sig.trait.dist[rowSums(padj.sig.trait.dist <= 0.05) > 0,])){ # If there are sig diff populations, plot significance
    stat.test            <- melt(padj.sig.trait.dist[theTrait,])
    stat.test$group1     <- "Total"
    stat.test$group2     <- row.names(stat.test)
    stat.test            <- stat.test[,c(2,3,1)]
    colnames(stat.test)  <- c("group1", "group2", "p.adj")
    stat.test            <- stat.test[stat.test$p <= 0.05,]
    row.names(stat.test) <- NULL
    stat.test            <- stat.test %>% add_significance()
    #stat.test           <- stat.test %>% add_y_position(fun = "max", formula = Frequency ~ Population, data = m) 
    # y pos is always the same, since we have a stacked bar
    stat.test$y.position <- 1.25
    
    ggplot(m, aes(x = Population, y = Frequency)) + 
      geom_bar(position = position_stack(), stat = "identity", mapping = aes_string(fill = theTrait)) +
      coord_flip() +
      scale_fill_startrek() + 
      theme_pubr(base_size = 14) +
      scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0", "25%", "50%", "75%", "100%"),oob = scales::oob_keep) +
      annotate(x=0, xend=0, y=0, yend=1, colour="black", lwd=1, geom="segment") +
      stat_pvalue_manual(coord.flip = T, remove.bracket = T, data = stat.test, x = "group2", label = "p.adj.signif", size = 6, vjust = 0.75) +
      ggtitle(label = theTrait, subtitle = paste("Fisher exact FDR:", round(padj.per.trait.sig.trait.dist[theTrait], digits = 3), sep = " ")) +
      theme(axis.text.y  = element_text(angle = 0, size = 14, colour = M.int_refined.pop.colors[c(-1,-2)]),
            axis.text.x  = element_text(size = 10),
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank(),
            axis.line    = element_blank(),
            legend.title = element_blank(),
            aspect.ratio = 1.2)
    ggsave(paste("metadata_correlations/Individual trait distribution binned/", theTrait, ".pdf", sep = ""), width = label_width + 4, height = 4)
    
  }else{ # No sig diff populations
    ggplot(m, aes(x = Population, y = Frequency)) + 
      geom_bar(position = position_stack(), stat = "identity", mapping = aes_string(fill = theTrait)) + 
      coord_flip() +
      scale_fill_startrek() + 
      theme_pubr(base_size = 14) + 
      scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0", "25%", "50%", "75%", "100%"),oob = scales::oob_keep) +
      annotate(x=0, xend=0, y=0, yend=1, colour="black", lwd=1, geom="segment") +
      ggtitle(label = theTrait, subtitle = paste("Fisher exact FDR:", round(padj.per.trait.sig.trait.dist[theTrait], digits = 3), sep = " ")) +
      theme(axis.text.y = element_text(angle = 0, size = 14, colour = M.int_refined.pop.colors[c(-1,-2)]),
            axis.text.x  = element_text(size = 10),
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank(),
            axis.line    = element_blank(), 
            legend.title = element_blank(),
            aspect.ratio = 1.2)
    ggsave(paste("metadata_correlations/Individual trait distribution binned/", theTrait, ".pdf", sep = ""), width = label_width + 4, height = 4)
  }
}
