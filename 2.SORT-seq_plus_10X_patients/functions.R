#========================================================================================================================
## Correlation plot
##========================================================================================================================
plot.cor <- function(object = final.pop.call.integrated.full.seurat, assay = "RNA", features = "GAPDH", cor.feature = "CLEC4A", group.by = "ident", ncol = 2, file.name = "cor.plot.pdf", width = 10, height = 10){
  # Get average expression for the requested features per group
  avg.expression.gene.set <- as.data.frame(AverageExpression(object   = object, 
                                                             assays   = assay, 
                                                             features = c(features, cor.feature), 
                                                             group.by = group.by)[[assay]])

  # Scale from 0 to 1
  avg.expression.gene.set <- as.data.frame(apply(avg.expression.gene.set, 1, rescale))
  
  
  # Iterate over all comparisons with CLEC4A
  ggscatter(data           = avg.expression.gene.set,
            x              = cor.feature,
            y              = colnames(avg.expression.gene.set)[-length(colnames(avg.expression.gene.set))],
            combine        = T,
            ncol           = ncol,
            add            = "reg.line",
            add.params     = list(color = "red", fill = "grey"),
            conf.int       = TRUE, 
            cor.coef       = TRUE, 
            cor.method     = "spearman",
            xlab           = paste(cor.feature,"(AU)", sep = " "),
            ylab           = "(AU)", 
            cor.coef.coord = c(0,1.1), 
            ylim           = c(0,1.1),
            xlim           = c(0,1),
  ) +
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
    theme_pubr(border = T)+
    theme(axis.title = element_text(face = "italic"),
          aspect.ratio = 1)
  ggsave(file.name, width = width, height = height)
}


##========================================================================================================================
## DGIdb drugability plots helper functions
##========================================================================================================================
shape.DGI.gene.pathway.overlap <- function(x = intersected.pathways.foamy.ont.unique.Mye.type.markers){
  # First reduce the list of dataframes to a list of unqiue genes per pathway (we are now only interested in the gene/pathway count, not the drug count)
  tmp.list <- list()
  for(thePathway in names(x)){
    tmp.list[[thePathway]] <- unique(x[[thePathway]]$Gene)
  }
  
  # Get the maxium number of genes per pathway so we can pad out the df
  n <- max(unlist(lapply(tmp.list, length)))
  df <- data.frame(row.names = 1:n)
  
  # Add each list element vector of unique genes as a column of the new df, padding with NA as needed
  for(thePathway in names(tmp.list)){
    if(length(tmp.list[[thePathway]] < n)){
      df[,thePathway] <- c(tmp.list[[thePathway]], rep(NA, (n - length(tmp.list[[thePathway]]))))
    }else{
      df[,thePathway] <- tmp.list[[thePathway]]
    }
  }
  return(df)
}

concat.DGI.score_and_abundance <- function(x = intersected.pathways.foamy.ont.unique.Mye.type.markers){
  # Prep the data
  x <- do.call("rbind", x)[,1:3]
  x$Pathway <- gsub("\\.[0-9]+", "", row.names(x))
  row.names(x) <- NULL
  x <- x[order(x$Gene, x$Pathway, x$Score, decreasing = T),]
  
  # Get high score per gene
  d <- duplicated(x$Gene)
  x.score <- x[!d,c(1,3)]
  x.score$Score <- as.numeric(x.score$Score)
  row.names(x.score) <- NULL
  
  # Get abundance of drugs per gene
  x <- x[order(x$Gene, x$Drug, decreasing = T),]
  x.abundance <- data.frame()
  for(theGene in unique(x$Gene)){
    s  <- subset(x, Gene == theGene)
    df <- data.frame("Gene" = theGene, "Abundance" = length(unique(s$Drug)))
    x.abundance <- rbind(x.abundance, df)
  }
  x.abundance$Abundance <- as.numeric(x.abundance$Abundance)
  x <- merge(x.score, x.abundance, by =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "Gene")
  return(x)
}

##========================================================================================================================
## Make DGIdb drugability plots
##========================================================================================================================
plot.DGI.results <- function(dgi.list = Mye.markers.dgi.filtered, n = 15, save.dir = "drugability_results/markers_plus_LRI_pairs", name = "potential drugs"){
  for( i in names(dgi.list)){
    m              <- dgi.list[[i]]
    cat(paste(i,"\n"))
    if(nrow(m != 0)){
      m$Drug       <- factor(m$Drug, levels = rev(unique(m$Drug)))
      m$Score      <- as.numeric(m$Score)
      m            <- m[order(m$Score, decreasing = T),]
      row.names(m) <- NULL # Reset row numbering
      
      # Check if we have n distinct genes
      if(length(unique(m$Gene)) > n){
        # Take top n unique genes from m. (Select on the first occurrence of n + 1, because if we select the last occurrence of n it breaks if gene[n]'s score is fragmented (not sorted together))
        m <- m[1:row.names(head(m[m$Gene == unique(m$Gene)[(n + 1)],,drop = F], n = 1)),]
        m <- m[1:(nrow(m) - 1),]  
      }else if (length(unique(m$Gene)) <= n){
        # Set n to the max number of genes we have, no need to alter m
        n <- length(unique(m$Gene))
      }
      
      m$GeneDrug <- paste(m$Gene, m$Drug, sep = " : ")
      m$Gene     <- factor(m$Gene, levels = unique(m$Gene))
      
      # Plot top hits per cluster in a barplot
      ggplot(m, aes(x = Drug, y = Score, fill = Score)) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() +
        scale_fill_gradient(low = "dodgerblue", high = "coral") +
        scale_x_discrete(breaks = m$Drug, labels = m$GeneDrug) +
        xlab("Gene : Drug") +
        ylab("Score") +
        ggtitle(paste(i, "potential drugs")) +
        theme(panel.background = element_blank(),
              legend.position = "none",
              axis.ticks.y     = element_blank(),
              plot.margin      = unit(c(1,1,1,1), "cm"),
              axis.text        = element_text(size = 14, face = "plain"), 
              aspect.ratio     = 2/1
        )
      ggsave(paste(save.dir, "/Cluster ", i, " target ", name, " top ", n, " bars.pdf", sep = ""), width = 10)
    
      # Plot top hits per cluster in a dotplot
      ggplot(m, aes(x = Gene, y = Drug, col = Score)) +
        geom_point(size = 3) +
        scale_color_gradient(low = "dodgerblue", high = "coral") +
        xlab("Gene") +
        ylab("Drug") +
        ggtitle(paste(i, "potential drugs")) +
        theme(panel.background = element_blank(),
              axis.ticks.y     = element_blank(),
              plot.margin      = unit(c(1,1,1,1), "cm"),
              axis.text        = element_text(size = 14, face = "plain"),
              axis.text.x      = element_text(angle = 45, hjust = 1),
              aspect.ratio     = 1.5/1
        )
      ggsave(paste(save.dir, "/Cluster ", i, " target ", name, " top ", n, " dots.pdf", sep = ""), width = 10)
    }
  }
}

plot.DGI.summary <- function(x = intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance, name = "Foamy populations", save.dir = "drugability_results/markers_plus_LRI_pairs"){
  ggplot(x, aes(x = Count, y = Score)) +
    geom_point(aes( size = Abundance, color = Aggregate)) +
    theme_pubr(border = T, legend = "right") +
    scale_color_continuous(name = "Aggregate score") +
    scale_size_continuous(name = "# of drugs per target") +
    ylab("DGIdb Score") +
    xlab("# of enriched pathways per gene") +
    theme(aspect.ratio = 1) +
    ggtitle(paste(name, "drugable targets")) +
    geom_label_repel(aes(label = Gene))
  ggsave(paste(save.dir, "/Cluster ", name, " hits summary.pdf", sep = ""), width = 12, height = 8)
}
##========================================================================================================================
## Make cellphonedb receptor ligand interaction plots
##========================================================================================================================
custom_LRI_plots <- function(object = exp.m.s.cpdb.relevant, save.dir = "cellphonedb_results", name = "my_lri", sub.re = NULL, subclus.colors = subclus.colors, top.score = 1, verbose = F){
  ## First round check the ligands
  if(verbose) cat("Plotting ligands...\n")
  sub.object <- object[grep(sub.re, object$ligand.cells),]

  # Plot all ligands and scores colored by receptor cells.
  if(verbose) cat("  All ligands...\n")
  ggplot(object, aes(x = ligand, y= interaction.score, col = receptor.cells)) +
    geom_point(size = 2) +
    geom_hline(yintercept = top.score, size = 0.5, col = "red") +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(name, "cell ligands", sep = " ")) +
    ylab("Interaction score") +
    theme_pubr(border = T) +
    theme(axis.line = element_blank(),
          panel.border = element_rect(size = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1/5
    )
  ggsave(paste(save.dir, "/", name, " cell plots - ALL ligands colored by receptor cell.pdf", sep = ""), width = 40, height = 10)
  
  # Filter: Keep only one instance of a ligand interacting with any 'sub.re' cells
  if(verbose) cat("  Filtering...\n")
  sub.object$receptor_int <- paste(sub.object$receptor, sub.object$receptor.cells, sep = ".")
  sub.object <- sub.object[order(sub.object$receptor_int),]
  d <- duplicated(sub.object$receptor_int)
  sub.object <- sub.object[!d,]
  
  # Sort on score
  sub.object <- sub.object[order(sub.object$interaction.score, decreasing = T),]

  # Get the frequencies of interactions per population
  sub.object.table <- sort(table(sub.object$receptor.cell), decreasing = T)

  # Plot the frequency of interactions per population
  if(verbose) cat("  Frequency of populations...\n")
  m <- melt(sub.object.table)
  colnames(m) <- c("Population", "Frequency")
  ggplot(m, aes(x = Population, y = Frequency, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    scale_fill_manual(values = subclus.colors) + 
    theme_pubr() +
    ylab("# of receptors") +
    xlab("Receptor population") +
    ggtitle(paste("Receptors interacting with ligands from", name, "cells", sep = " ")) +
    scale_y_continuous(limits = c(0,roundCeiling(max(as.vector(m$Frequency)))),expand = c(0,0)) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1)
  ggsave(paste(save.dir, "/", name, " cell plots - Frequency of receptors interacting with ligands from ", name,  " cells per cell type.pdf", sep = ""))
  
  # Get the frequencies of receptors
  sub.object.receptor.table <- sort(table(as.vector(sub.object$receptor)), decreasing = T)
  
  # Plot the frequency of receptors
  if(verbose) cat("  Frequencies...\n")
  m <- melt(sub.object.receptor.table)
  colnames(m) <- c("Gene", "Frequency")
  m <- m[!is.na(m$Gene),]
  ggplot(m, aes(x = Gene, y = Frequency)) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    theme_pubr() +
    ylab("Frequency of receptor") +
    xlab("Receptor") +
    ggtitle(paste("Receptors interacting with ligands from", name, "cells", sep = " ")) +
    scale_y_continuous(limits = c(0,roundCeiling(max(as.vector(m$Frequency)))), expand = c(0,0)) +
    theme(legend.position = "right",
          aspect.ratio = 1/7.5, 
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste(save.dir, "/", name, " cell plots - Frequency of receptors interacting with ligands from ", name,  " cells accumulative.pdf", sep = ""), width = 20, height = 8)
  
  # Plot all interactions
  if(verbose) cat("  All interaction matrix...\n")
  m <- melt(sub.object)
  colnames(m)
  m$variable  <- NULL
  colnames(m) <- c(colnames(m)[1:6], "interaction.score")
  ggplot(m, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score)) +
    geom_point() +
    geom_jitter() +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(name, "cell ligands", sep = " ")) +
    theme_light() +
    ylab("Receptors (in other cells)") +
    xlab(paste("Ligands (in", name, "cells)", sep = " ")) +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 16),
          aspect.ratio = 1)
  ggsave(paste(save.dir, "/", name, " cell plots - All ligands and receptors matrix ligand viewpoint.pdf", sep = ""), width = 25, height = 25)
  
  # Plot top scores only
  if(verbose) cat("  Top score matrix...\n")
  m <- melt(sub.object)
  colnames(m)
  m$variable  <- NULL
  colnames(m) <- c(colnames(m)[1:6], "interaction.score")
  m <- m[m$interaction.score > top.score,]
  ggplot(m, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score)) +
    geom_point() +
    geom_jitter() +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(name, "cell ligands", sep = " "), subtitle = paste("interaction scores >", top.score, sep = "")) +
    theme_light() +
    ylab("Receptors (in other cells)") +
    xlab(paste("Ligands (in", name, "cells)", sep = " ")) +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 16),
          aspect.ratio = 1)
  ggsave(paste(save.dir, "/", name, " cell plots - Interaction score over ", top.score, " - ligands and receptors matrix ligand viewpoint.pdf", sep = ""), width = 25, height = 25)
  
  
  #=============================================================================================================
  ## Second round check the receptors
  if(verbose) cat("\nPlotting receptors...\n")
  sub.object <- object[grep(sub.re, object$receptor.cells),]

  # Plot all receptors and scores colored by ligand cells.
  if(verbose) cat("  All receptors...\n")
  ggplot(object, aes(x = receptor, y= interaction.score, col = ligand.cells)) +
    geom_point(size = 2) +
    geom_hline(yintercept = top.score, size = 0.5, col = "red") +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(name, "cell receptors", sep = " ")) +
    ylab("Interaction score") +
    theme_pubr(border = T) +
    theme(axis.line = element_blank(),
          panel.border = element_rect(size = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1/5
    )
  ggsave(paste(save.dir, "/", name, " cell plots - ALL receptors colored by ligand cell.pdf", sep = ""), width = 40, height = 10)
  
  # Filter: Keep only one instance of a ligand interacting with any 'sub.re' cells
  if(verbose) cat("  Filtering...\n")
  sub.object$ligand_int <- paste(sub.object$ligand, sub.object$ligand.cells, sep = ".")
  sub.object <- sub.object[order(sub.object$ligand_int),]
  d <- duplicated(sub.object$ligand_int)
  sub.object <- sub.object[!d,]
  
  # Sort on score
  sub.object <- sub.object[order(sub.object$interaction.score, decreasing = T),]
  
  # Get the frequencies of interactions per population
  sub.object.table <- sort(table(sub.object$ligand.cell), decreasing = T)
  
  # Plot the frequency of interactions per population
  if(verbose) cat("  Frequencies per population...\n")
  m <- melt(sub.object.table)
  colnames(m) <- c("Population", "Frequency")
  ggplot(m, aes(x = Population, y = Frequency, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    scale_fill_manual(values = subclus.colors) + 
    theme_pubr() +
    ylab("# of ligands") +
    xlab("Ligand population") +
    ggtitle(paste("Ligands interacting with receptors on", name, "cells", sep = " ")) +
    scale_y_continuous(limits = c(0,roundCeiling(max(as.vector(m$Frequency)))),expand = c(0,0)) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1)
  ggsave(paste(save.dir, "/", name, " cell plots - Frequency of ligands interacting with receptors from ", name,  " cells per cell type.pdf", sep = ""))
  
  # Get the frequencies of ligands
  sub.object.ligand.table <- sort(table(as.vector(sub.object$ligand)), decreasing = T)
  
  # Plot the frequency of receptors
  if(verbose) cat("  Frequencies...\n")
  m <- melt(sub.object.ligand.table)
  colnames(m) <- c("Gene", "Frequency")
  m <- m[!is.na(m$Gene),]
  ggplot(m, aes(x = Gene, y = Frequency)) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    theme_pubr() +
    ylab("Frequency of ligand") +
    xlab("Ligand") +
    ggtitle(paste("Ligands interacting with receptors from", name, "cells", sep = " ")) +
    scale_y_continuous(limits = c(0,roundCeiling(max(as.vector(m$Frequency)))), expand = c(0,0)) +
    theme(legend.position = "right",
          aspect.ratio = 1/7.5, 
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste(save.dir, "/", name, " cell plots - Frequency of ligands interacting with receptors from ", name,  " cells accumulative.pdf", sep = ""), width = 20, height = 8)
  
  # Plot all interactions
  if(verbose) cat("  All interaction matrix...\n")
  m <- melt(sub.object)
  colnames(m)
  m$variable  <- NULL
  colnames(m) <- c(colnames(m)[1:6], "interaction.score")
  ggplot(m, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score)) +
    geom_point() +
    geom_jitter() +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(name, "cell receptors", sep = " ")) +
    theme_light() +
    ylab(paste("Receptors (in", name, "cells)", sep = " ")) +
    xlab("Ligands (in other cells)") +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 16),
          aspect.ratio = 1)
  ggsave(paste(save.dir, "/", name, " cell plots - All ligands and receptors matrix receptor viewpoint.pdf", sep = ""), width = 25, height = 25)
  
  # Plot top scores only
  if(verbose) cat("  Top score matrix...\n\n")
  m <- melt(sub.object)
  colnames(m)
  m$variable  <- NULL
  colnames(m) <- c(colnames(m)[1:6], "interaction.score")
  m <- m[m$interaction.score > top.score,]
  ggplot(m, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score)) +
    geom_point() +
    geom_jitter() +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(name, "cell receptors", sep = " "), subtitle = paste("interaction scores >", top.score, sep = "")) +
    theme_light() +
    ylab(paste("Receptors (in", name, "cells)", sep = " ")) +
    xlab("Ligands (in other cells)") +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 16),
          aspect.ratio = 1)
  ggsave(paste(save.dir, "/", name, " cell plots - Interaction score over ", top.score, " - ligands and receptors matrix receptor viewpoint.pdf", sep = ""), width = 25, height = 25)
}


#=========================================================================================================================
## Cehck if a charcter string can be interpreted as a date
##========================================================================================================================
check.date <- function(x){
  return(!is.na(as.Date(as.character(x), optional = T)))
}

#=========================================================================================================================
## Plot a cseries of ustomised plots for stratified expression
##========================================================================================================================
stratPlots <- function(object = samples.seurat, assay = "RNA", group.by = NULL, features = "GAPDH", cols = NULL, x.label = "Stratified", name = "plot_basename", ncol = 1){
   p1 <- VlnPlot(object  = object,
                 assay   = "RNA",  
                features = features, 
                group.by = group.by, 
                cols     = cols, 
                ncol     = ncol, combine = T) & theme(axis.title.x = element_blank()) 
  p1 <- p1 + plot_annotation(caption = x.label, theme = theme(plot.caption = element_text(hjust = 0.5)))
  
  ggsave(paste(name, " - violin plot.pdf", sep = ""), height = 5*(length(features)/ncol), plot = p1)

  DotPlot(object = object, assay = "RNA", features = features, group.by = group.by, cols = "Spectral", dot.scale = 15) + scale_colour_gsea()
  ggsave(paste(name, " - dot plot.pdf", sep = ""), width = 2*length(features))
  
  FeaturePlot(object = object, assay = "RNA", features = features, pt.size = 2, split.by = group.by, cols = c("grey", "firebrick"), order = T, ncol = ncol)
  ggsave(paste(name, " - feature plot.pdf", sep = ""), height = 5*(length(features)/ncol), width = 5*(length(features)/ncol))
}


#==================================================================================================================================
## Stratify the seurat object based on expresion quntile of a given gene, and Plot a series of custom plots for stratified data.
##=================================================================================================================================
stratifyByExpression <- function(object = all.seur.combined, assay = "RNA", strat.by = "GAPDH", gene.groups = list("Foamy genes" = c("ABCA1", "ABCG1", "OLR1"), "Resident genes" = c("LYVE1", "MRC1", "FOLR2"), "LAM genes" = c("TREM2", "CD9", "GPNMB"), "Inflammatory genes" = c("IL1B", "TNF", "NLRP3", "CASP1")), file.name = "myplot", return.object = T, do.plot = T, verbose = T, onlyUMAP = F, lower.quantile = 0.25, upper.quantile = 0.75){
  # Sanity check gene.groups
  if(!onlyUMAP & type(gene.groups) != "list"){
    cat("ERROR: gene.groups argument has to be a named list of vectors with gene names!\n")
    if(return.object){
      return(object)
    }else{
      return(invisible(NULL))
    }
  }
  
  # Get marker expression levels
  if(!strat.by %in% row.names(object)){
    cat("ERROR:", strat.by, "not present in this Seurat object's RNA assay!\n")
    if(return.object){
      return(object)
    }else{
      return(invisible(NULL))
    }
  }
  
  theMarker <- GetAssayData(object = object, assay = "RNA")[strat.by,]
  
  # Check some stats
  if(verbose){
    cat("Cells not expressing ",          strat.by, " ", sum(theMarker == 0), "\n", sep = "")
    cat("Cells expressing ",              strat.by, " ", sum(theMarker != 0), "\n", sep = "")
    cat("Perentage of cells expressing ", strat.by, " ", round((sum(theMarker != 0) / length(theMarker)) * 100), "%\n", sep = "")
  }
  
  # First subset the cells with expression of theMarker
  theMarker.expressed <- theMarker[theMarker != 0]
  
  # Then stratify the expressing cells
  if(verbose){
    cat("Quantiles of", strat.by, "expressing cells:\n", names(summary(theMarker.expressed)),"\n", summary(theMarker.expressed), "\n")
  }
  
  theMarker.low  <- names(theMarker.expressed)[theMarker.expressed <= quantile(theMarker.expressed, lower.quantile)]
  theMarker.med  <- names(theMarker.expressed)[theMarker.expressed >  quantile(theMarker.expressed, lower.quantile) & theMarker.expressed <= quantile(theMarker.expressed, upper.quantile)]
  theMarker.high <- names(theMarker.expressed)[theMarker.expressed >  quantile(theMarker.expressed, upper.quantile)]

  # Add as metadata to the seurat object
  marker.col <- paste(strat.by, "expr", sep = "_")
  object     <- AddMetaData(object, metadata = rep("Zero", ncol(object)), col.name = marker.col)

  object@meta.data[row.names(object@meta.data) %in% theMarker.low,  marker.col] <- "Low"
  object@meta.data[row.names(object@meta.data) %in% theMarker.med,  marker.col] <- "Medium"
  object@meta.data[row.names(object@meta.data) %in% theMarker.high, marker.col] <- "High"

    object@meta.data[,marker.col] <- factor( object@meta.data[,marker.col], levels = c("Zero", "Low", "Medium", "High"))
  
  # Plot the stratification
  if(do.plot){
    VlnPlot(object   = object, 
            assay    = assay, 
            features = strat.by, 
            group.by = marker.col, 
            cols     = c("grey", "bisque", "coral", "firebrick")) + xlab(paste("Stratified by", strat.by ,"expression", sep = " "))
    ggsave(filename  = paste(file.name, " - violin.pdf"))
    
    customUMAP(object    = object, 
               group.by  = marker.col, 
               pt.size   = 4, 
               cols      = c("grey", "bisque", "coral", "firebrick"), 
               title     = paste(strat.by, "expression", sep = " "),
               file.name = paste(file.name, " - UMAP.pdf"))
    
    # And some genes, if requested
    if(!onlyUMAP){
      for(theGroup in names(gene.groups)){
        stratPlots(object   = object, 
                   group.by = marker.col, 
                   features = gene.groups[[theGroup]], 
                   cols     = c("grey", "bisque", "coral", "firebrick"),
                   x.label  = paste("Stratified by", strat.by ,"expression", sep = " "), 
                   name     = paste(file.name, " - ", theGroup, sep = ""),
                   ncol     = 1)
      }
    }
  }
    
  if(return.object){
    return(object)
  }
}

#=========================================================================================================================
## Plot a customised UMAP (dimplot), keeping standard values as most often used in this project.
##========================================================================================================================
# Wrapper to exttract polt limits so we can dynamically alter the position of the 'axes arrows' on the umap plot
get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

customUMAP <- function(object = object, group.by = NULL, pt.size = 4, label = F,label.size = NULL, cols = NULL, title = NULL, font.size = 14, reduction = "umap", shuffle = T, legend.pos = "top", seed = 666, file.name = "myplot.pdf", plot.height = 10, plot.width = 10, sizes.highlight = 1, cells.highlight = NULL, cells = NULL){
  p1 <- DimPlot(object = object, group.by = group.by, pt.size = pt.size, label = label, label.size = label.size, cols = cols, reduction = reduction, shuffle = shuffle, seed = seed, sizes.highlight = sizes.highlight, cells.highlight = cells.highlight, cells = cells) +
        ggtitle(title) +
        xlab("UMAP 2") +
        ylab("UMAP 1") + 
        theme_pubr(base_size = font.size, legend = legend.pos) +
        theme(axis.line        = element_blank(), 
              axis.text        = element_blank(), 
              axis.ticks       = element_blank(),
              axis.title       = element_text(hjust = 0.025),
              panel.background = element_blank(),
              title            = element_text(size = (font.size + 2), face = "bold"))

  # Add dynamic 'axes arrows'
    p1 <- p1 + coord_cartesian(xlim   = c((floor(get_plot_limits(p1)$xmin) - 0.2), ceiling(get_plot_limits(p1)$xmax)),
                               ylim   = c((floor(get_plot_limits(p1)$ymin) - 0.2), ceiling(get_plot_limits(p1)$ymax)),
                               expand = F, 
                               clip   = "off") +
          annotate(geom = "segment", 
                   x    = floor(get_plot_limits(p1)$xmin), 
                   xend = (floor(get_plot_limits(p1)$xmin) + 1.5), 
                   y    = floor(get_plot_limits(p1)$ymin), 
                   yend = floor(get_plot_limits(p1)$ymin), 
                   lineend = "round", lwd = 2, arrow = grid::arrow(length = unit(15, "pt"), angle = 20)) +
          annotate(geom = "segment", 
                   x    = floor(get_plot_limits(p1)$xmin), 
                   xend = floor(get_plot_limits(p1)$xmin),
                   y    = floor(get_plot_limits(p1)$ymin), 
                   yend = (floor(get_plot_limits(p1)$ymin) + 1.5), 
                   lineend = "round", lwd = 2, arrow = grid::arrow(length = unit(15, "pt"), angle = 20))
ggsave(filename = file.name, height = plot.height, width = plot.width, plot = p1)
}

#=========================================================================================================================
## Plot a customised violinPlot, keeping standard values as most often used in this project.
##========================================================================================================================
customVln <- function(object = samples.seurat, group.by = NULL, idents = NULL, features = "GAPDH", assay = "RNA", draw.names = T, name = "plot_name.pdf", splitPlot = F, ncol = NULL, stack = F, pt.size = 0, width = 15, height = 15, split.by = NULL, cols = NULL){
  if(stack == T){
    # If we stack we also wanna flip
    VlnPlot(pt.size    = pt.size, 
            group.by   = group.by,
            idents     = idents,
            object     = object, 
            features   = features,
            assay      = assay,
            split.by   = split.by, 
            split.plot = splitPlot,
            stack      = T,
            flip       = T) &
      theme_pubr(base_size = 14, x.text.angle = 45) &
      theme(panel.background = element_blank(),
            plot.margin      = unit(c(1,1,1,5), units = "cm"),
             title           = element_text(size = 16, face = "bold")
      ) 
    ggsave(filename  = name, 
           width     = width, 
           height    = height,
           limitsize = F)
  }else{
    if(draw.names == T){
    VlnPlot(pt.size    = pt.size, 
            group.by   = group.by, 
            idents     = idents,
            object     = object, 
            features   = features,
            assay      = assay,
            ncol       = ncol,
            split.by   = split.by, 
            split.plot = splitPlot,
            cols       = cols) & 
      theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") &
      theme(panel.background = element_blank(),
            plot.margin      = unit(c(1,1,1,5), units = "cm"),
            title            = element_text(size = 16, face = "bold")
      )
    ggsave(filename  = name, 
           width     = width, 
           height    = height,
           limitsize = F)
    }else{
      VlnPlot(pt.size    = pt.size, 
              group.by   = group.by, 
              idents     = idents,
              object     = object, 
              features   = features,
              assay      = assay,
              ncol       = ncol,
              split.by   = split.by, 
              split.plot = splitPlot,
              cols       = cols) & 
        theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") &
        theme(panel.background = element_blank(),
              plot.margin      = unit(c(1,1,1,5), units = "cm"),
              title            = element_text(size = 16, face = "bold"),
              axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()
        )
      ggsave(filename  = name, 
             width     = width, 
             height    = height,
             limitsize = F)
    }
  }
}


#=========================================================================================================================
## Plot a customised DotPlot, keeping standard values as most often used in this project.
##========================================================================================================================
customDot <- function(object = samples.seurat, group.by = NULL, idents = NULL, features = "GAPDH", assay = "RNA", cluster.idents = T, name = "plot_name.pdf", split.by = NULL, dot.scale = 7, width = "auto", height = "auto"){
  if(width == "auto"){
    # Base width 5. Add 0.2 for every extra feature, and add the maximum string width of the cluster names.
    pop.width <- max(strwidth(levels(Idents(object)), units = "inches")) * 2.54 # Convert to cm
    if(length(features) <= 5){
      width <- 5 + pop.width
      } else{
        width <- 5 + ((length(features) - 5) * 0.5) + pop.width
      }
  }

  if(height == "auto"){
    # Determine number of categories on the y axis:
    # Number of idents...
    if(is.null(group.by)){
      num.cats <- length(unique(Idents(object)))
    # Or number of whatever we're grouping by...
    }else{
      num.cats <- length(unique(object@meta.data[,group.by]))
    }
    # multiplied by the number of categories we are splitting by
    if(!is.null(split.by)){
      num.cats <- num.cats * length(unique(object@meta.data[,split.by]))
    }
    
    # Base height 5. Add 0.25 for every extra category
    if(num.cats <= 5){
      height <- 5}
    else{
        height <- 5 + ((num.cats - 5) * 0.25)
      }
    }

    DotPlot(object       = object,
          group.by       = group.by, 
          idents         = idents,
          cluster.idents = cluster.idents,
          dot.scale      = dot.scale, 
          split.by       = split.by, 
          cols           = "Spectral", 
          assay          = assay,
          features       = features) + 
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")) +
    scale_color_gsea()
    ggsave(filename  = name, 
           width     = width, 
           height    = height,
           limitsize = F)
}


#=========================================================================================================================
## Plot a customised FeaturePlot, keeping standard values as most often used in this project.
##========================================================================================================================
customFeature <- function(object = samples.seurat, cols = c("grey", "blue"), features = "GAPDH", name = "plot_name.pdf", reduction = "umap", pt.size = 1, order = T, width = 15, height = 15, ncol = NULL){
  FeaturePlot(object    = object,
              features  = features,
              cols      = cols,
              reduction = reduction, 
              pt.size   = pt.size, 
              order     = order, 
              ncol      = ncol) & 
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") &
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold"), aspect.ratio = 1)
  ggsave(filename  = name, 
         width     = width, 
         height    = height,
         limitsize = F)
}

#=========================================================================================================================
## Plot a bunch of custom plots at once, using most commonly ussed variations.
## Don't add .pdf to the name here as we are making compund names based on the settings used!
## Implemented: feature, violins and DotPlots
##========================================================================================================================
bunchOfCustomPlots <- function(object = samples.seurat, idents = NULL, group.by = NULL, features = "GAPDH", assay = "RNA", reduction = "umap", Vln.draw.names = T, name = "plot_name", feature.pt.size = 3, Vln.pt.size = 0, dot.scale = 7, Vln.width = 15, Vln.height = 15, Vln.stack = FALSE, Vln.color = NULL, Dot.width = "auto", Dot.height = "auto", ncol = NULL){
  # Feature plot
  customFeature(object = object, features = features, name = paste(name, " - feature plot.pdf", sep = ""), ncol = ncol, pt.size = feature.pt.size, reduction = reduction)

  # Violin plots
  customVln(object = object, idents = idents, group.by = group.by, features = features, assay = assay, draw.names = Vln.draw.names, name = paste(name, " - violin plot.pdf", sep = ""), width = Vln.width, height = Vln.height, pt.size = Vln.pt.size, ncol = ncol, stack = Vln.stack, cols = Vln.color)
  
  # Dot plots
  customDot(object = object, idents = idents, group.by = group.by, features = features, assay = assay, name = paste(name, " - dot plot.pdf", sep = ""),    width = Dot.width, height = Dot.height, dot.scale = dot.scale)
}


#-----------------------------------------------------------------------------------------
# Define a function to plot the GO graphs
plot_GO_from_table <- function(go_data, name="GO_terms.pdf", sig_cut=0.05){
  #Prep our vars
  go_data[,1]       <-  gsub("_", " ", go_data[,1])
  go_data[,4]       <- -log10(go_data[,2])
  colnames(go_data) <-  c("term","padj", "direction","log10_FDR_qvalue")
  sig_cut           <- -log10(sig_cut)
  height            <-  (nrow(go_data) / 2) + 1
  cat(height, "\n")
  #Calculate how wide the labels will be
  label_width <- 1 + (max(strwidth(go_data$term, units = "inches")) * 2.54) #Converting to cm
  #print(label_width) #old debug line
  
  #Make the plot!
  ggplot(go_data,aes(x = reorder(term, log10_FDR_qvalue)  , y = log10_FDR_qvalue, fill = direction)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c(Down = "dodgerblue", Up = "coral1")) +
    coord_flip() +
    ylab("-log10(FDR q-value)") +
    xlab("") +
    ylim(0,roundCeiling(x = max(go_data$log10_FDR_qvalue))) +
    geom_hline(yintercept = sig_cut, col = "purple", linetype = "dashed", size = 1) +
    
    theme(legend.position  = "bottom",
          panel.background = element_blank(),
          axis.ticks.y     = element_blank(),
          plot.margin      = unit(c(1,1,1,1), "cm"),
          axis.text        = element_text(size = 16, face = "bold", colour = "black"),
    )
  
  #And save it
  ggsave(name,width = (label_width / 2.54) + 7.36, height = height)
} #plot_GO_from_table()

#-----------------------------------------------------------------------------------------
# Define a function to plot GO volcanoes
plot_GO_volcano_from_table <- function(go_data, name="GO_terms.pdf", sig_cut=0.05){
  #Prep our vars
  go_data[,1]       <-  gsub("_", " ", go_data[,1])
  go_data[,5]       <- -log10(go_data[,2])
  colnames(go_data) <-  c("term","padj", "direction", "avg.log2FC", "log10_FDR_qvalue")
  
  go_data$fill.type  <- "Down.sig"
  go_data[which(go_data$avg.log2FC < 0 & go_data$padj >  sig_cut), "fill.type"] <- "Down.ns"
  go_data[which(go_data$avg.log2FC > 0 & go_data$padj >  sig_cut), "fill.type"] <- "Up.ns"
  go_data[which(go_data$avg.log2FC > 0 & go_data$padj <= sig_cut), "fill.type"] <- "Up.sig"
  
  #Make the plot!
  ggplot(go_data,aes(x = avg.log2FC  , y = log10_FDR_qvalue, color = fill.type)) +
    geom_point(size = 2) + 
    scale_color_manual(values = c(Down.ns = "#C6E3FF", Up.ns = "#FFDBD4", Down.sig = "dodgerblue", Up.sig = "coral1")) +
    ylab("-log10(FDR q-value)") +
    xlab("log2(FC)") +
    ylim(0,roundCeiling(x = max(go_data$log10_FDR_qvalue))) +
    geom_hline(yintercept = -log10(sig_cut), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(legend.position  = "none",
          panel.background = element_blank(),
          axis.line        = element_line(),
          plot.margin      = unit(c(1,1,1,1), "cm"),
          axis.text        = element_text(size = 16, face = "bold", colour = "black"),
    )
  
  #And save it
  ggsave(name)
} #plot_GO_volcano_from_table()

#-----------------------------------------------------------------------------------------
#Define helper function for rounding the axis of the plots
roundCeiling <- function(x) {
  if(length(x) != 1) stop("'x' must be of length 1")
  if(x < 5){
    return(5)
  }
  else if(x < 10){
    return(10)
  }
  else if(x < 15){
    return(15)
  }
  else if(x < 100){
    round(x+5,-1)
  }
  else if(x < 1000){
    round(x+50,-2)
  }
  else{
    round(x+500,-3)
  }
} #roundCeiling()

#wrapper for write.table so that it makes a header without blank tab at the start when using row names
write.better.table <- function(data, file = "myfile.txt"){
  write.table(transform(data,Symbol = rownames(data)
  )[,c(length(colnames(data))+1,
       1:length(colnames(data))
  )
  ],
  file      = file, 
  quote     = F, 
  row.names = F,
  sep       = "\t"
  )
}


#-----------------------------------------------------------------------------------------
# Color spacer
colorSpacer <- function(startcolor = "red", middlecolors = NULL, endcolor = "red4", colors = NULL, steps = 4, return.colors = TRUE, per.line = 10, show.labels = T, show.grid = T, grid.color = "gray"){
  # If a premade list of colors is given, use that.
  if(!is.null(colors)){
    steps <- length(colors)
  }else{
    # If start, end point, and (optionally) mid point colors are given, use those instead
    if(is.null(middlecolors)){
      pal <- colorRampPalette(c(startcolor, endcolor))
    }else{
      pal <- colorRampPalette(c(startcolor, middlecolors, endcolor))
    }
    colors <- pal(steps)
  }
  if(return.colors){
    return(colors)
  }else{
    # Make a plot
    # Define number of lines to plot and how many colors will go on the last line
    if(steps <= per.line){
      lines      <- 1
      last.steps <- 0
      per.line  <- steps
    }else{
      lines      <- ceiling(steps/per.line)
      last.steps <- per.line - steps %% per.line
      if(last.steps == per.line){
        last.steps <- 0
      }
    }
    # Start building the plot
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "" )
    
    # Add colored rectangles
    if(show.grid){
      rect(     rep((0:(per.line - 1) / per.line), lines)[1:(per.line * lines - last.steps)],  
                sort(rep((0:(   lines - 1) / lines),    per.line), decreasing=T)[1:(per.line * lines - last.steps)],   
                rep((1:per.line       / per.line), lines)[1:(per.line * lines - last.steps)],
                sort(rep((1:lines          / lines),    per.line), decreasing=T)[1:(per.line * lines - last.steps)],  
                border = grid.color,
                col    = colors)
    }else{
      rect(     rep((0:(per.line - 1) / per.line), lines)[1:(per.line * lines - last.steps)],  
                sort(rep((0:(   lines - 1) / lines),    per.line), decreasing=T)[1:(per.line * lines - last.steps)],   
                rep((1:per.line       / per.line), lines)[1:(per.line * lines - last.steps)],
                sort(rep((1:lines          / lines),    per.line), decreasing=T)[1:(per.line * lines - last.steps)],  
                border = NA,
                col    = colors)
    }
    
    # Add axes (if requested)
    if(show.labels){
      axis(2, 
           at     = (0:(lines - 1) / lines) + (1 / (lines * 2)), 
           labels = rev(seq(0,(lines * per.line),per.line)[1:lines]), 
           tick   = F, lty = 6, pos = 0.01)
      axis(3,
           at     = 0:(per.line - 1) / per.line + (1 / (per.line * 2)),
           labels = seq(1,per.line),
           tick   = F, lty = 6, pos = 0.99)
      
      # Add labels
      mtext("Little step #", 3, line = 1)
      mtext("Big step #", 2, line = 1)
      mtext("Step #: Little + Big Step #", 1)
    }
  }
}
#-----------------------------------------------------------------------------------------

get_ontology <- function(res, name = "cluster", outdir = ".", return.data = F, full_GSEA = T, plot.top.n = 20, universe = all.seur.combined, volcano.plot = T){
  #Set up our data
  #Get background of all expressed genes in our dataset
  universe <- row.names(universe)
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")$ENTREZID
  
  #Extract genes from res
  genes <- data.frame(Symbol = res[,"gene"])

    #Return emptyhanded if empty cluster
  if (nrow(genes)==0){
    return("NA")
  }
  
  #Retrieve entrez IDs
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = as.vector(genes$Symbol), columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  genes <- genes[!duplicated(genes$SYMBOL),]
  colnames(genes) <- c("symbol", "entrezID")
  genes <- genes[,c(2,1)]
  
  #Add logFC column
  genes$logFC <- res[,"avg_log2FC"]
  
  #-------------------------------------------------------
  #First run the custom collection and GO to make neat bar graphs
  #CP & H
  #Create annotation table
  gs.annots <- buildCustomIdx(genes$entrezID,
                              species = "human",
                              gsets = ont_sets,
                              label = "CP_and_H",
                              name = "CP_and_H"
  )
  
  #Get pathway enrichments
  gsa.pathways <- egsea.ora(geneIDs = genes$entrezID,
                            title = "x", 
                            universe = universe,
                            gs.annots = gs.annots,
                            symbolsMap = genes,
                            logFC = genes$logFC,
                            sort.by = "p.adj",
                            report.dir = ".",
                            num.threads = 11,
                            report = F
  )
  
  #GO
  #Create annotation table
  gs.annots <- buildMSigDBIdx(entrezIDs = genes$entrezID, 
                              species = "human", 
                              geneSets = "c5"
  )
  
  #Get GO enrichments
  gsa.go <- egsea.ora(geneIDs = genes$entrezID,
                      gs.annots = gs.annots,
                      title = "x",
                      universe = universe,
                      symbolsMap = genes,
                      logFC = genes$logFC,
                      sort.by = "p.adj",
                      report.dir = ".",
                      num.threads = 11,
                      report = F
  )

  
  #Extract top n enriched sets
  if(length(gsa.pathways@results) > 0){
    top.pathways <- data.frame(name      = row.names(topSets(gsa.pathways, number = plot.top.n, names.only = F)),
                               padj      = topSets(gsa.pathways, number = plot.top.n, names.only = F)$p.adj,
                               direction = topSets(gsa.pathways, number = plot.top.n, names.only = F)$direction)
  }
  
  if(length(gsa.go@results) > 0){
    top.GO       <- data.frame(name      = row.names(topSets(gsa.go, number = plot.top.n, names.only = F)), 
                               padj      = topSets(gsa.go, number = plot.top.n, names.only = F)$p.adj,
                               direction = topSets(gsa.go, number = plot.top.n, names.only = F)$direction)
  }
  
  #Check if any sets are significant, otherwise return emptyhanded
  if(exists("top.pathways") && exists("top.GO")){
    if (!sum(unlist(lapply(c(top.pathways$padj, top.GO$padj), function(x) any(x<0.1))))){
      cat("No significant sets found! Exiting...\n")
      return("NA")
    }
  }else if(exists("top.GO")){
    if (!sum(unlist(lapply(c(top.GO$padj), function(x) any(x<0.1))))){
      cat("No significant sets found! Exiting...\n")
      return("NA")
    }
  }else if(exists("top.pathways")){
    if (!sum(unlist(lapply(c(top.pathways$padj), function(x) any(x<0.1))))){
      cat("No significant sets found! Exiting...\n")
      return("NA")
    }
  }
  
  #Plot!
  if(exists("top.pathways")){
   plot_GO_from_table(top.pathways, name = paste(outdir, "/", name, ".pathways.pdf", sep = ""))
  }
  if(exists("top.GO")){
   plot_GO_from_table(top.GO,       name = paste(outdir, "/", name, ".GO_terms.pdf", sep = ""))
  }
  
  #-------------------------------------------------------
  #Now run the full GSEA suite and make a report dir, if wanted
  if(full_GSEA){
    #Build annotation index
    gs.annots <- buildIdx(entrezIDs = genes$entrezID,
                          species = "human",
                          msigdb.gsets = "all",
                          gsdb.gsets = "all",
                          go.part = T,
                          kegg.updated = T
    )
    
    #Run EGSEA
    egsea.ora(geneIDs = genes$entrezID,
              title = name,
              universe = universe,
              gs.annots = gs.annots,
              logFC = genes$logFC,
              symbolsMap = genes,
              display.top = 20,
              sort.by = "p.adj",
              report.dir = paste(outdir, "/", name, "_EGSEA", sep = ""),
              kegg.dir = paste(outdir, "/", name, "_EGSEA/kegg-dir", sep = ""),
              num.threads = 4,
              interactive = F,
              report = T,
              verbose = F
    )
  }
  
  # Plot a volcano if wanted
  if(volcano.plot){
    if(length(gsa.pathways@results) > 0){
      pathways <- data.frame(name          = row.names(topSets(gsa.pathways, number = 500, names.only = F)),
                             padj          = topSets(gsa.pathways, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.pathways, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.pathways, number = 500, names.only = F)$avg.logfc.dir
      )
    }
    
    if(length(gsa.go@results) > 0){
      goterms  <- data.frame(name          = row.names(topSets(gsa.go, number = 500, names.only = F)), 
                             padj          = topSets(gsa.go, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.go, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.go, number = 500, names.only = F)$avg.logfc.dir)
    }
    
    #Plot!
    if(exists("pathways")){
      plot_GO_volcano_from_table(pathways, name = paste(outdir, "/", name, ".pathway volcano.pdf", sep = ""))
    }
    if(exists("goterms")){
      plot_GO_volcano_from_table(goterms,  name = paste(outdir, "/", name, ".GO_term volcano.pdf", sep = ""))
    }
  }
  
  
  #Return the tables if wanted
  if(return.data){
    cat("Returning Data...\n")
    if(length(gsa.pathways@results) > 0){
      pathways <- data.frame(name          = row.names(topSets(gsa.pathways, number = 500, names.only = F)),
                             padj          = topSets(gsa.pathways, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.pathways, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.pathways, number = 500, names.only = F)$avg.logfc.dir)
    }
    if(length(gsa.go@results) > 0){
      goterms  <- data.frame(name          = row.names(topSets(gsa.go, number = 500, names.only = F)), 
                             padj          = topSets(gsa.go, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.go, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.go, number = 500, names.only = F)$avg.logfc.dir)
    }
    
    if(exists("pathways") && exists("goterms")){
      returnList        <- list(pathways,goterms)
      names(returnList) <- c("Pathways", "GO_terms")
    }else if(exists("pathways")){
      returnList        <- list(pathways)
      names(returnList) <- c("Pathways")
    }else if(exists("goterms")){
      returnList        <- list(goterms)
      names(returnList) <- c("GO_terms")
    }else{
      returnList <- "NA"
    }
    return(returnList)
  }
} #get_ontology(res)

get_ontology_DE <- function(res, name = "cluster", outdir = ".", return.data = F, full_GSEA = T, plot.top.n = 20, universe = all.seur.combined, volcano.plot = T){
  #Set up our data
  #Get background of all expressed genes in our dataset
  universe <- row.names(universe)
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")$ENTREZID
  
  #Extract genes from res
  genes <- data.frame("Symbol" = row.names(res))

  #Return emptyhanded if empty or too small cluster
  if (nrow(genes)<5){
    return("NA")
  }

  #Retrieve entrez IDs
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = as.vector(genes$Symbol), columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  genes <- genes[!duplicated(genes$SYMBOL),]
  colnames(genes) <- c("symbol", "entrezID")
  genes <- genes[,c(2,1)]

    #Add logFC column
  genes$logFC <- res[,"log2FoldChange"]
  
  #-------------------------------------------------------
  #First run the custom collection and GO to make neat bar graphs
  #CP & H
  #Create annotation table
  gs.annots <- buildCustomIdx(genes$entrezID,
                              species = "human",
                              gsets = ont_sets,
                              label = "CP_and_H",
                              name = "CP_and_H"
  )
  
  #Get pathway enrichments
  gsa.pathways <- egsea.ora(geneIDs = genes$entrezID,
                            title = "x", 
                            universe = universe,
                            gs.annots = gs.annots,
                            symbolsMap = genes,
                            logFC = genes$logFC,
                            sort.by = "p.adj",
                            report.dir = ".",
                            num.threads = 11,
                            report = F
  )
  
  #GO
  #Create annotation table
  gs.annots <- buildMSigDBIdx(entrezIDs = genes$entrezID, 
                              species = "human", 
                              geneSets = "c5"
  )
  
  #Get GO enrichments
  gsa.go <- egsea.ora(geneIDs = genes$entrezID,
                      gs.annots = gs.annots,
                      title = "x",
                      universe = universe,
                      symbolsMap = genes,
                      logFC = genes$logFC,
                      sort.by = "p.adj",
                      report.dir = ".",
                      num.threads = 11,
                      report = F
  )
  
  
  #Extract top n enriched sets
  if(length(gsa.pathways@results) > 0){
    top.pathways <- data.frame(name      = row.names(topSets(gsa.pathways, number = plot.top.n, names.only = F)),
                               padj      = topSets(gsa.pathways, number = plot.top.n, names.only = F)$p.adj,
                               direction = topSets(gsa.pathways, number = plot.top.n, names.only = F)$direction)
  }
  
  if(length(gsa.go@results) > 0){
    top.GO       <- data.frame(name      = row.names(topSets(gsa.go, number = plot.top.n, names.only = F)), 
                               padj      = topSets(gsa.go, number = plot.top.n, names.only = F)$p.adj,
                               direction = topSets(gsa.go, number = plot.top.n, names.only = F)$direction)
  }
  
  #Check if any sets are significant, otherwise return emptyhanded
  if(exists("top.pathways") && exists("top.GO")){
    if (!sum(unlist(lapply(c(top.pathways$padj, top.GO$padj), function(x) any(x<0.1))))){
      cat("No significant sets found! Exiting...\n")
      return("NA")
    }
  }else if(exists("top.GO")){
    if (!sum(unlist(lapply(c(top.GO$padj), function(x) any(x<0.1))))){
      cat("No significant sets found! Exiting...\n")
      return("NA")
    }
  }else if(exists("top.pathways")){
    if (!sum(unlist(lapply(c(top.pathways$padj), function(x) any(x<0.1))))){
      cat("No significant sets found! Exiting...\n")
      return("NA")
    }
  }
  
  #Plot!
  if(exists("top.pathways")){
    plot_GO_from_table(top.pathways, name = paste(outdir, "/", name, ".pathways.pdf", sep = ""))
  }
  if(exists("top.GO")){
    plot_GO_from_table(top.GO,       name = paste(outdir, "/", name, ".GO_terms.pdf", sep = ""))
  }
  
  #-------------------------------------------------------
  #Now run the full GSEA suite and make a report dir, if wanted
  if(full_GSEA){
    #Build annotation index
    gs.annots <- buildIdx(entrezIDs = genes$entrezID,
                          species = "human",
                          msigdb.gsets = "all",
                          gsdb.gsets = "all",
                          go.part = T,
                          kegg.updated = T
    )
    
    #Run EGSEA
    egsea.ora(geneIDs = genes$entrezID,
              title = name,
              universe = universe,
              gs.annots = gs.annots,
              logFC = genes$logFC,
              symbolsMap = genes,
              display.top = 20,
              sort.by = "p.adj",
              report.dir = paste(outdir, "/", name, "_EGSEA", sep = ""),
              kegg.dir = paste(outdir, "/", name, "_EGSEA/kegg-dir", sep = ""),
              num.threads = 4,
              interactive = F,
              report = T,
              verbose = F
    )
  }
  
  # Plot a volcano if wanted
  if(volcano.plot){
    if(length(gsa.pathways@results) > 0){
      pathways <- data.frame(name          = row.names(topSets(gsa.pathways, number = 500, names.only = F)),
                             padj          = topSets(gsa.pathways, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.pathways, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.pathways, number = 500, names.only = F)$avg.logfc.dir
      )
    }
    
    if(length(gsa.go@results) > 0){
      goterms  <- data.frame(name          = row.names(topSets(gsa.go, number = 500, names.only = F)), 
                             padj          = topSets(gsa.go, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.go, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.go, number = 500, names.only = F)$avg.logfc.dir)
    }
    
    #Plot!
    if(exists("pathways")){
      plot_GO_volcano_from_table(pathways, name = paste(outdir, "/", name, ".pathway volcano.pdf", sep = ""))
    }
    if(exists("goterms")){
      plot_GO_volcano_from_table(goterms,  name = paste(outdir, "/", name, ".GO_term volcano.pdf", sep = ""))
    }
  }
  
  
  #Return the tables if wanted
  if(return.data){
    cat("Returning Data...\n")
    if(length(gsa.pathways@results) > 0){
      pathways <- data.frame(name          = row.names(topSets(gsa.pathways, number = 500, names.only = F)),
                             padj          = topSets(gsa.pathways, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.pathways, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.pathways, number = 500, names.only = F)$avg.logfc.dir)
    }
    if(length(gsa.go@results) > 0){
      goterms  <- data.frame(name          = row.names(topSets(gsa.go, number = 500, names.only = F)), 
                             padj          = topSets(gsa.go, number = 500, names.only = F)$p.adj,
                             direction     = topSets(gsa.go, number = 500, names.only = F)$direction,
                             avg.logfc.dir = topSets(gsa.go, number = 500, names.only = F)$avg.logfc.dir)
    }
    
    if(exists("pathways") && exists("goterms")){
      returnList        <- list(pathways,goterms)
      names(returnList) <- c("Pathways", "GO_terms")
    }else if(exists("pathways")){
      returnList        <- list(pathways)
      names(returnList) <- c("Pathways")
    }else if(exists("goterms")){
      returnList        <- list(goterms)
      names(returnList) <- c("GO_terms")
    }else{
      returnList <- "NA"
    }
    return(returnList)
  }
} #get_ontology(res)
