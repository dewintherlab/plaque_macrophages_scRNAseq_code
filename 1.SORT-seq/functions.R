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
bunchOfCustomPlots <- function(object = samples.seurat, idents = NULL, group.by = NULL, features = "GAPDH", assay = "RNA", Vln.draw.names = T, name = "plot_name", feature.pt.size = 3, Vln.pt.size = 0, dot.scale = 7, Vln.width = 15, Vln.height = 15, Vln.stack = FALSE, Vln.color = NULL, Dot.width = "auto", Dot.height = "auto", ncol = NULL){
  # Feature plot
  customFeature(object = object, features = features, name = paste(name, " - feature plot.pdf", sep = ""), ncol = ncol, pt.size = feature.pt.size)
  
  # Violin plots
  customVln(object = object, idents = idents, group.by = group.by, features = features, assay = assay, draw.names = Vln.draw.names, name = paste(name, " - violin plot.pdf", sep = ""), width = Vln.width, height = Vln.height, pt.size = Vln.pt.size, ncol = ncol, stack = Vln.stack, cols = Vln.color)
  
  # Dot plots
  customDot(object = object, idents = idents, group.by = group.by, features = features, assay = assay, name = paste(name, " - dot plot.pdf", sep = ""),    width = Dot.width, height = Dot.height, dot.scale = dot.scale)
}




#-----------------------------------------------------------------------------------------

merge.on.rownames <- function(x, y){
  m <- merge(x, y, by = 0)
  row.names(m) <- m$Row.names
  m$Row.names <- NULL
  return(m)
}

#-----------------------------------------------------------------------------------------

update.symbols <- function(my.counts.table = my.counts.table, n = n){
  #Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
  #Retrieve the mapping

  gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(my.counts.table), columns = c("SYMBOL","ENTREZID"), keytype = "ALIAS", multiVals = "first")
  
  # Mitochondrial genes work around: org.Hs.eg.db gets its mappings from NCBI, which weirdly does not support proper MT- gene names but instead converts it to some aliases of there own.
  # So here we just take the MT- symbls and apply them directly to the 'SYMBOL' column.
  gene.info[grep("^MT-", gene.info[,"ALIAS"]),"SYMBOL"] <- gene.info[grep("^MT-", gene.info[,"ALIAS"]),"ALIAS"]
 
  # In fact, turns out there is a couple 100 lncRNAs without a proper NCBI symbol so let's add those aliases back in as well
  gene.info[is.na(gene.info$SYMBOL),"SYMBOL"] <-  gene.info[is.na(gene.info$SYMBOL),"ALIAS"]
  
  #Keep only one mapping per alias
  d <- duplicated(gene.info$ALIAS)
  gene.info <- gene.info[!d,]
  
  #Add mappings to the table
  my.counts.table$SYMBOL <- gene.info$SYMBOL
  
  #Remove non-mappings (old LOCs and stuff that are not availble anymore)
  na <- is.na(my.counts.table$SYMBOL)
  my.counts.table <- my.counts.table[!na,]
  
  #Keep only highest expressed SYMBOL if multiple old aliases are now merged
  #Get maximum number of dups per gene
  num_dups <- max(table(my.counts.table$SYMBOL))
  
  #Loop until only one entry per gene left
  while(num_dups > 1){
    #First retrieve the list of duplicated symbols
    o <- order(my.counts.table$SYMBOL)
    my.counts.table <- my.counts.table[o,]
    d <- duplicated(my.counts.table$SYMBOL)
    
    #Then keep only the highest expressed one (this part will fail if we happen upon 2 rows with equal expression, but for now that hasn't happened yet so meh ;)
    h <- rowSums(my.counts.table[d,1:n]) > rowSums(my.counts.table[which(d==T)-1,1:n])
    h <- c(h, rowSums(my.counts.table[which(d==T)-1,1:n]) > rowSums(my.counts.table[d,1:n]))
    h <- h[which(!h)]
    my.counts.table <- my.counts.table[!row.names(my.counts.table) %in% names(h),]
    
    num_dups <- num_dups - 1 #One dup removed
  }
  
  #Overwrite old symbols in row names and clean up
  row.names(my.counts.table) <- my.counts.table$SYMBOL
  my.counts.table$SYMBOL <- NULL
  
  return(my.counts.table)
}
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

#Retrieve all cell names (not) expressing a certain gene
GetCellNames <- function(myObject = all.seur.combined, theGene, treshold = 0, invert = F){
  m        <- as.matrix(myObject@data)
  theRows  <- grep(theGene, 
                   row.names(m)
                   )
  if(invert){
    theCells <- colnames(m[,
                           m[theRows,] <= treshold
                          ]
                      )
  } else {
    theCells <- colnames(m[,
                           m[theRows,] > treshold
                          ]
    )
  }
  return(theCells)
}
