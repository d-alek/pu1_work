# data = wtkoaml_all (start from there) "Integrated logFC logRpkm wtkoaml ALL coded_01.txt" 
# only changed AML & CML sets to ..._KEGG to avoid file overwrite issue with aml
# sets = CaG (done as a test) + colnames(wtkoaml_all)[c(15,17,21:30,32:34,39,40,44,46:48,50,52:54,56,58,61,62,67,72,76:96)]

plot_HP_HT_perSetFunct <- function(data, sets, ...) {
  
  # here you can further limit logFC scope of wtkoaml_all, based on logRPKM difference ("wtkoaml_cut"):
  wtkoaml_cut <- dplyr::filter(data, 
                               abs((wt1+wt2)-(ko1+ko2))>=2.5 | 
                               abs((ko1+ko2)-(aml1+aml2))>2.5 | 
                               abs((wt1+wt2)-(aml1+aml2))>3.5)
  # filter on genes that are members of at least one set ("sbst"):
  wtkoaml_cut <- dplyr::mutate(wtkoaml_cut, set.no = rowSums(select(wtkoaml_cut, c(15:17, 21:96))))
  sbst <- dplyr::filter(wtkoaml_cut, set.no >= 1)
  # improve column situation and order logically ("sbst"):
  sbst <- dplyr::mutate(sbst, Pu1_motif = ifelse((RGAGGAARY == 1 | WGAGGAAG == 1), 1, 0))
  sbst <- dplyr::select(sbst, -wnt, -RGAGGAARY, -WGAGGAAG)
  sbst <- sbst[ , c(1:14, 95, 18, 15:17, 19:27, 37:38, 28:29, 32:36, 39:42, 30:31, 43:51, 53:54, 52, 55:58, 60, 59, 61:80, 84:87, 81, 83, 88:90, 82, 92, 91, 93:94)]
  sbst <- dplyr::rename(sbst, Pu1_ChIP = ChIP)
  
  
  for (i in seq_along(sets)) {
    
    set <- sets[i]
    
    ##1 KO dn part ##
    temp_sbst <- filter_( sbst, paste(set, "== 1 & ((wt1+wt2)-(ko1+ko2)) > 2.85 & adjPko <= 0.01") )
    
    if (nrow(temp_sbst) != 0) {
      
      #1 heatmap:
      all.de <- list(temp_sbst$Symbol) # can change it to $GeneName if required
      every.g <- temp_sbst[,3:8]
      rownames(every.g) <- temp_sbst$Symbol
      heat.ab <- every.g - rowMeans(every.g)
      
      rowseps <- cumsum(sapply(all.de, FUN=length))
      lwid <- c(4, 0.3*ncol(heat.ab) + 2)
      lhei <- c(2, 0.2*nrow(heat.ab) + 2)
      last <- 0L
      re.heat <- list()
      for (x in all.de) {
        current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
        if (length(x)!=1L) {
          ddr <- as.dendrogram(hclust(as.dist(1-cor(t(current))), method = "average")) # tried: average, single, complete, ward.D, ward.D2, mcquitty
          #ddr <- as.dendrogram(hclust(dist(current), method = "average"))
          #ddr <- reorder(ddr, -rowMeans(current))
          ddr <- reorder(ddr, -rowMeans(current[,1:2,drop=FALSE])) # see which makes more sense
          current <- current[order.dendrogram(ddr),,drop=FALSE] 
        }
        re.heat[[length(re.heat)+1L]] <- current
        last <- last + length(x) 
      }
      heat.ab <- do.call(rbind, re.heat)
      heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
      
      if (nrow(heat.abb) == 1) {
        heat.abb <- rbind(heat.abb, heat.abb)
      }
      
      cur.fname <- paste(set, "_ko_dn_HP.pdf", sep = "")
      pdf(cur.fname, height=sum(lhei), width=sum(lwid))
      heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
                col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
                rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10),
                main = paste(set, "(KO dn)", sep = " "))
      dev.off()
      
      #1 heattable:
      sbst_long <- gather(temp_sbst, key = set, value = true, range = c(15:94))
      set.order <- colnames(temp_sbst)[c(15:94)]
      gen.order <- rev(rownames(heat.abb))
      pdf(paste(set, "_ko_dn_HT.pdf", sep = ""), width = ((0.2*80)+2), height = ((0.2*nrow(temp_sbst))+2), useDingbats = FALSE)
      p <- ggplot(sbst_long, aes(set, Symbol)) + 
        geom_tile(aes(fill = true), colour = "white") + 
        scale_fill_gradient(low = "gray91", high = "steelblue") + 
        theme_classic(base_size = 16) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
        labs(x = NULL, y = NULL) + 
        scale_x_discrete(limits = set.order, expand = c(0,0)) + 
        scale_y_discrete(limits = gen.order, expand = c(0,0)) + 
        theme(legend.position = "none")
      print(p)
      dev.off()
      
      names1 <- select(temp_sbst, GeneID, Symbol, GeneName)
      names1 <- add_column(names1, Change = rep(paste(set, "_ko_dn", sep = ""), nrow(temp_sbst)))
      names1 <- names1[match(rownames(heat.abb), names1$Symbol), ]
    }
    
    
    ##2 KO up part ##
    temp_sbst <- filter_( sbst, paste(set, "== 1 & ((wt1+wt2)-(ko1+ko2)) < -2.9 & adjPko <= 0.01") )
    
    if (nrow(temp_sbst) != 0) {
      
      #2 heatmap:
      all.de <- list(temp_sbst$Symbol) # can change it to $GeneName if required
      every.g <- temp_sbst[,3:8]
      rownames(every.g) <- temp_sbst$Symbol
      heat.ab <- every.g - rowMeans(every.g)
      
      rowseps <- cumsum(sapply(all.de, FUN=length))
      lwid <- c(4, 0.3*ncol(heat.ab) + 2)
      lhei <- c(2, 0.2*nrow(heat.ab) + 2)
      last <- 0L
      re.heat <- list()
      for (x in all.de) {
        current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
        if (length(x)!=1L) {
          ddr <- as.dendrogram(hclust(as.dist(1-cor(t(current))), method = "average")) # tried: average, single, complete, ward.D, ward.D2, mcquitty
          #ddr <- as.dendrogram(hclust(dist(current), method = "average"))
          #ddr <- reorder(ddr, -rowMeans(current))
          ddr <- reorder(ddr, -rowMeans(current[,3:4,drop=FALSE])) # see which makes more sense
          current <- current[order.dendrogram(ddr),,drop=FALSE] 
        }
        re.heat[[length(re.heat)+1L]] <- current
        last <- last + length(x) 
      }
      heat.ab <- do.call(rbind, re.heat)
      heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
      
      if (nrow(heat.abb) == 1) {
        heat.abb <- rbind(heat.abb, heat.abb)
      }
      
      cur.fname <- paste(set, "_ko_up_HP.pdf", sep = "")
      pdf(cur.fname, height=sum(lhei), width=sum(lwid))
      heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
                col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
                rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10),
                main = paste(set, "(KO up)", sep = " "))
      dev.off()
      
      #2 heattable:
      sbst_long <- gather(temp_sbst, key = set, value = true, range = c(15:94))
      set.order <- colnames(temp_sbst)[c(15:94)]
      gen.order <- rev(rownames(heat.abb))
      pdf(paste(set, "_ko_up_HT.pdf", sep = ""), width = ((0.2*80)+2), height = ((0.2*nrow(temp_sbst))+2), useDingbats = FALSE)
      p <- ggplot(sbst_long, aes(set, Symbol)) + 
        geom_tile(aes(fill = true), colour = "white") + 
        scale_fill_gradient(low = "gray91", high = "steelblue") + 
        theme_classic(base_size = 16) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
        labs(x = NULL, y = NULL) + 
        scale_x_discrete(limits = set.order, expand = c(0,0)) + 
        scale_y_discrete(limits = gen.order, expand = c(0,0)) + 
        theme(legend.position = "none")
      print(p)
      dev.off()
      
      names2 <- select(temp_sbst, GeneID, Symbol, GeneName)
      names2 <- add_column(names2, Change = rep(paste(set, "_ko_up", sep = ""), nrow(temp_sbst)))
      names2 <- names2[match(rownames(heat.abb), names2$Symbol), ]
    }
    
    ##3 KO part together ##
    temp_sbst <- filter_( sbst, paste(set, "== 1 & (((wt1+wt2)-(ko1+ko2)) > 2.85 | 
                                      ((wt1+wt2)-(ko1+ko2)) < -2.9 ) & adjPko <= 0.01") )
    
    if (nrow(temp_sbst) != 0) {
      
      #3 heatmap:
      all.de <- list(temp_sbst$Symbol) # can change it to $GeneName if required
      every.g <- temp_sbst[,3:8]
      rownames(every.g) <- temp_sbst$Symbol
      heat.ab <- every.g - rowMeans(every.g)
      
      rowseps <- cumsum(sapply(all.de, FUN=length))
      lwid <- c(4, 0.3*ncol(heat.ab) + 2)
      lhei <- c(2, 0.2*nrow(heat.ab) + 2)
      last <- 0L
      re.heat <- list()
      for (x in all.de) {
        current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
        if (length(x)!=1L) {
          ddr <- as.dendrogram(hclust(as.dist(1-cor(t(current))), method = "average")) # tried: average, single, complete, ward.D, ward.D2, mcquitty
          #ddr <- as.dendrogram(hclust(dist(current), method = "average"))
          #ddr <- reorder(ddr, -rowMeans(current))
          ddr <- reorder(ddr, -rowMeans(current[,1:2,drop=FALSE])) # see which makes more sense
          current <- current[order.dendrogram(ddr),,drop=FALSE] 
        }
        re.heat[[length(re.heat)+1L]] <- current
        last <- last + length(x) 
      }
      heat.ab <- do.call(rbind, re.heat)
      heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
      
      if (nrow(heat.abb) == 1) {
        heat.abb <- rbind(heat.abb, heat.abb)
      }
      
      cur.fname <- paste(set, "_ko_HP.pdf", sep = "")
      pdf(cur.fname, height=sum(lhei), width=sum(lwid))
      heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
                col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
                rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10),
                main = paste(set, "(KO)", sep = " "))
      dev.off()
      
      #3 heattable:
      sbst_long <- gather(temp_sbst, key = set, value = true, range = c(15:94))
      set.order <- colnames(temp_sbst)[c(15:94)]
      gen.order <- rev(rownames(heat.abb))
      pdf(paste(set, "_ko_HT.pdf", sep = ""), width = ((0.2*80)+2), height = ((0.2*nrow(temp_sbst))+2), useDingbats = FALSE)
      p <- ggplot(sbst_long, aes(set, Symbol)) + 
        geom_tile(aes(fill = true), colour = "white") + 
        scale_fill_gradient(low = "gray91", high = "steelblue") + 
        theme_classic(base_size = 16) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
        labs(x = NULL, y = NULL) + 
        scale_x_discrete(limits = set.order, expand = c(0,0)) + 
        scale_y_discrete(limits = gen.order, expand = c(0,0)) + 
        theme(legend.position = "none")
      print(p)
      dev.off()
      
      names3 <- select(temp_sbst, GeneID, Symbol, GeneName)
      names3 <- add_column(names3, Change = rep(paste(set, "_ko", sep = ""), nrow(temp_sbst)))
      names3 <- names3[match(rownames(heat.abb), names3$Symbol), ]
    }
    
    
    ##4 AML dn part ##
    temp_sbst <- filter_( sbst, paste(set, "== 1 & 
                           ( (((ko1+ko2)-(aml1+aml2)) > 2.9 & adjPaml <= 0.01) | 
                               (((wt1+wt2)-(aml1+aml2)) > 4 & adjPaml <= 0.02) ) & 
                           !(((wt1+wt2)-(ko1+ko2)) > 2.85) & 
                           !(((wt1+wt2)-(ko1+ko2)) < -2.9)") )
    
    if (nrow(temp_sbst) != 0) {
      
      #4 heatmap:
      all.de <- list(temp_sbst$Symbol) # can change it to $GeneName if required
      every.g <- temp_sbst[,3:8]
      rownames(every.g) <- temp_sbst$Symbol
      heat.ab <- every.g - rowMeans(every.g)
      
      rowseps <- cumsum(sapply(all.de, FUN=length))
      lwid <- c(4, 0.3*ncol(heat.ab) + 2)
      lhei <- c(2, 0.2*nrow(heat.ab) + 2)
      last <- 0L
      re.heat <- list()
      for (x in all.de) {
        current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
        if (length(x)!=1L) {
          ddr <- as.dendrogram(hclust(as.dist(1-cor(t(current))), method = "average")) # tried: average, single, complete, ward.D, ward.D2, mcquitty
          #ddr <- as.dendrogram(hclust(dist(current), method = "average"))
          #ddr <- reorder(ddr, -rowMeans(current))
          ddr <- reorder(ddr, -rowMeans(current[,3:4,drop=FALSE])) # see which makes more sense
          current <- current[order.dendrogram(ddr),,drop=FALSE] 
        }
        re.heat[[length(re.heat)+1L]] <- current
        last <- last + length(x) 
      }
      heat.ab <- do.call(rbind, re.heat)
      heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
      
      if (nrow(heat.abb) == 1) {
        heat.abb <- rbind(heat.abb, heat.abb)
      }
      
      cur.fname <- paste(set, "_aml_dn_HP.pdf", sep = "")
      pdf(cur.fname, height=sum(lhei), width=sum(lwid))
      heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
                col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
                rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10),
                main = paste(set, "(AML dn)", sep = " "))
      dev.off()
      
      #4 heattable:
      sbst_long <- gather(temp_sbst, key = set, value = true, range = c(15:94))
      set.order <- colnames(temp_sbst)[c(15:94)]
      gen.order <- rev(rownames(heat.abb))
      pdf(paste(set, "_aml_dn_HT.pdf", sep = ""), width = ((0.2*80)+2), height = ((0.2*nrow(temp_sbst))+2), useDingbats = FALSE)
      p <- ggplot(sbst_long, aes(set, Symbol)) + 
        geom_tile(aes(fill = true), colour = "white") + 
        scale_fill_gradient(low = "gray91", high = "steelblue") + 
        theme_classic(base_size = 16) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
        labs(x = NULL, y = NULL) + 
        scale_x_discrete(limits = set.order, expand = c(0,0)) + 
        scale_y_discrete(limits = gen.order, expand = c(0,0)) + 
        theme(legend.position = "none")
      print(p)
      dev.off()
      
      names4 <- select(temp_sbst, GeneID, Symbol, GeneName)
      names4 <- add_column(names4, Change = rep(paste(set, "_aml_dn", sep = ""), nrow(temp_sbst)))
      names4 <- names4[match(rownames(heat.abb), names4$Symbol), ]
    }
    
    
    ##5 AML up part ##
    temp_sbst <- filter_( sbst, paste(set, "== 1 & 
                           ( (((ko1+ko2)-(aml1+aml2)) < -2.9 & adjPaml <= 0.01) | 
                               (((wt1+wt2)-(aml1+aml2)) < -4 & adjPaml <= 0.02) ) & 
                           !(((wt1+wt2)-(ko1+ko2)) > 2.85) & 
                           !(((wt1+wt2)-(ko1+ko2)) < -2.9)") )
    
    if (nrow(temp_sbst) != 0) {
      
      #5 heatmap:
      all.de <- list(temp_sbst$Symbol) # can change it to $GeneName if required
      every.g <- temp_sbst[,3:8]
      rownames(every.g) <- temp_sbst$Symbol
      heat.ab <- every.g - rowMeans(every.g)
      
      rowseps <- cumsum(sapply(all.de, FUN=length))
      lwid <- c(4, 0.3*ncol(heat.ab) + 2)
      lhei <- c(2, 0.2*nrow(heat.ab) + 2)
      last <- 0L
      re.heat <- list()
      for (x in all.de) {
        current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
        if (length(x)!=1L) {
          ddr <- as.dendrogram(hclust(as.dist(1-cor(t(current))), method = "average")) # tried: average, single, complete, ward.D, ward.D2, mcquitty
          #ddr <- as.dendrogram(hclust(dist(current), method = "average"))
          #ddr <- reorder(ddr, -rowMeans(current))
          ddr <- reorder(ddr, -rowMeans(current[,5:6,drop=FALSE])) # see which makes more sense
          current <- current[order.dendrogram(ddr),,drop=FALSE] 
        }
        re.heat[[length(re.heat)+1L]] <- current
        last <- last + length(x) 
      }
      heat.ab <- do.call(rbind, re.heat)
      heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
      
      if (nrow(heat.abb) == 1) {
        heat.abb <- rbind(heat.abb, heat.abb)
      }
      
      cur.fname <- paste(set, "_aml_up_HP.pdf", sep = "")
      pdf(cur.fname, height=sum(lhei), width=sum(lwid))
      heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
                col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
                rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10),
                main = paste(set, "(AML up)", sep = " "))
      dev.off()
      
      #5 heattable:
      sbst_long <- gather(temp_sbst, key = set, value = true, range = c(15:94))
      set.order <- colnames(temp_sbst)[c(15:94)]
      gen.order <- rev(rownames(heat.abb))
      pdf(paste(set, "_aml_up_HT.pdf", sep = ""), width = ((0.2*80)+2), height = ((0.2*nrow(temp_sbst))+2), useDingbats = FALSE)
      p <- ggplot(sbst_long, aes(set, Symbol)) + 
        geom_tile(aes(fill = true), colour = "white") + 
        scale_fill_gradient(low = "gray91", high = "steelblue") + 
        theme_classic(base_size = 16) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
        labs(x = NULL, y = NULL) + 
        scale_x_discrete(limits = set.order, expand = c(0,0)) + 
        scale_y_discrete(limits = gen.order, expand = c(0,0)) + 
        theme(legend.position = "none")
      print(p)
      dev.off()
      
      names5 <- select(temp_sbst, GeneID, Symbol, GeneName)
      names5 <- add_column(names5, Change = rep(paste(set, "_aml_up", sep = ""), nrow(temp_sbst)))
      names5 <- names5[match(rownames(heat.abb), names5$Symbol), ]
    }
    
    
    ##6 AML part together ##
    temp_sbst <- filter_( sbst, paste(set, "== 1 & 
                                      ( (abs((ko1+ko2)-(aml1+aml2)) > 2.9 & adjPaml <= 0.01) | 
                                      (abs((wt1+wt2)-(aml1+aml2)) > 4 & adjPaml <= 0.02) ) & 
                                      !(((wt1+wt2)-(ko1+ko2)) > 2.85) & 
                                      !(((wt1+wt2)-(ko1+ko2)) < -2.9)") )
    
    if (nrow(temp_sbst) != 0) {
      
      #6 heatmap:
      all.de <- list(temp_sbst$Symbol) # can change it to $GeneName if required
      every.g <- temp_sbst[,3:8]
      rownames(every.g) <- temp_sbst$Symbol
      heat.ab <- every.g - rowMeans(every.g)
      
      rowseps <- cumsum(sapply(all.de, FUN=length))
      lwid <- c(4, 0.3*ncol(heat.ab) + 2)
      lhei <- c(2, 0.2*nrow(heat.ab) + 2)
      last <- 0L
      re.heat <- list()
      for (x in all.de) {
        current <- heat.ab[last+1:length(x),,drop=FALSE] # drop=FALSE to preserve data frame
        if (length(x)!=1L) {
          ddr <- as.dendrogram(hclust(as.dist(1-cor(t(current))), method = "average")) # tried: average, single, complete, ward.D, ward.D2, mcquitty
          #ddr <- as.dendrogram(hclust(dist(current), method = "average"))
          #ddr <- reorder(ddr, -rowMeans(current))
          ddr <- reorder(ddr, -rowMeans(current[,3:4,drop=FALSE])) # see which makes more sense
          current <- current[order.dendrogram(ddr),,drop=FALSE] 
        }
        re.heat[[length(re.heat)+1L]] <- current
        last <- last + length(x) 
      }
      heat.ab <- do.call(rbind, re.heat)
      heat.abb <- as.matrix(heat.ab, rownames.force=TRUE)
      
      if (nrow(heat.abb) == 1) {
        heat.abb <- rbind(heat.abb, heat.abb)
      }
      
      cur.fname <- paste(set, "_aml_HP.pdf", sep = "")
      pdf(cur.fname, height=sum(lhei), width=sum(lwid))
      heatmap.2(heat.abb, dendrogram="none", trace="none", density.info="none",
                col=bluered, cexCol=1.2, symbreaks=TRUE, breaks=100, lwid=lwid, lhei=lhei, cexRow=1.2,
                rowsep=rowseps, Colv=FALSE, Rowv=FALSE, sepwidth=c(0.2, 0.2), margin=c(10,10),
                main = paste(set, "(AML)", sep = " "))
      dev.off()
      
      #6 heattable:
      sbst_long <- gather(temp_sbst, key = set, value = true, range = c(15:94))
      set.order <- colnames(temp_sbst)[c(15:94)]
      gen.order <- rev(rownames(heat.abb))
      pdf(paste(set, "_aml_HT.pdf", sep = ""), width = ((0.2*80)+2), height = ((0.2*nrow(temp_sbst))+2), useDingbats = FALSE)
      p <- ggplot(sbst_long, aes(set, Symbol)) + 
        geom_tile(aes(fill = true), colour = "white") + 
        scale_fill_gradient(low = "gray91", high = "steelblue") + 
        theme_classic(base_size = 16) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
        labs(x = NULL, y = NULL) + 
        scale_x_discrete(limits = set.order, expand = c(0,0)) + 
        scale_y_discrete(limits = gen.order, expand = c(0,0)) + 
        theme(legend.position = "none")
      print(p)
      dev.off()
      
      names6 <- select(temp_sbst, GeneID, Symbol, GeneName)
      names6 <- add_column(names6, Change = rep(paste(set, "_aml", sep = ""), nrow(temp_sbst)))
      names6 <- names6[match(rownames(heat.abb), names6$Symbol), ]
    }
    
    
    # writing table of full gene names for the set, ordered as in heatplots, heattables:
    CaG_order_genenames <- bind_rows(names1, names2, names3, names4, names5, names6)
    write.table(CaG_order_genenames, file=paste(set, "_order_gnames.txt", sep = ""), row.names=F, sep='\t', quote=FALSE)
    
  }
}






