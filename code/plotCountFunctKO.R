# data = wtkoaml_cut (start from there, see how is it "cut" in previous file)
# sets = sets_names
# colurs used: gray100, gray10, firebrick (alt. maroon), dodgerblue4, darkslategrey 
# notes: wnt column is empty; think about ~ comon limits for gradient across files 
# all_de (2.975,4) CaG..(2.25,4.75) KEGG&H (2.05,4.915) M (2.2,4.78); can size n() be normalised?
# notes: have to repeat: Apo(32), P53(34), PI3K(40), Notch(54) (png coerces to lowercase, pdf & svg coerce to caps) 
# (because overwriten with same-name diff-case file)

plotCountFunct <- function(data, sets, ...) {
  
  sets <- sets[59:79 , ]
  
  mi <- numeric()
  mx <- numeric()
  st <- character()
  ct <- integer()
  
  for (i in seq_along(sets[ , 1])) {
    
    set <- sets[i,1]
    set_name <- sets[i,2]
    
    selection <- filter(data, ko_change == "down" | ko_change == "up")
    selection <- filter_(selection, paste(set, "== 1")) # notice use of filter_ ; needed to fix "t sup" to "t_sup for this
    selection <- group_by(selection, ko_change, pu_motif_bind_2cat)
    selection <- mutate(selection, av_logPko = median(logPko)) # mutate does group statistics here!!!
    selection <- mutate(selection, cnt = n())
    
    counts = summarise(selection, cnt = mean(cnt))
    minc = min(counts$cnt)
    maxc = max(counts$cnt)
    if ((maxc - minc) > 3) {
      midc = round(mean(c(minc, maxc), 0))
    } else {
      midc = 0
    }
    
    png(paste(set,"_ko.png", sep = ""), width = 10, height = 5, units = "in", res = 300) # useDingbats = FALSE (for .pdf); units = "in", res = 300 (for .tiff, .png)
    g <- ggplot(selection, aes(pu_motif_bind_2cat, ko_change)) + 
      geom_point(aes(size = cnt, fill = av_logPko), shape = 21, colour = "gray70", alpha = 0.8) + 
      scale_size_area(max_size = 50, name = "Count", breaks = c(minc, midc, maxc)) + 
      labs(x ="Pu.1 motif or binding", y = "Gene DE", title = paste(set_name,"(KO)", sep = " ")) + 
      theme_bw() + theme(panel.grid.major = element_blank(), axis.ticks = element_blank(), legend.key = element_rect(colour = "white"), legend.box ="horizontal", legend.box.just = "right", legend.title.align = 0.5 , legend.text.align = 0, legend.margin = unit(0, "mm")) + 
      scale_x_discrete(limits = c("0", "1", "2"), labels = c("0" = "none", "1" = "motif", "2" = "binding"), expand = c(0.3, 0.6)) + 
      scale_y_discrete(limits = c("down", "up")) + 
      scale_fill_gradient(low = "gray100", high = "darkslategrey", name = "Median -log10 p-value", limits = c(2.2, 4.78))
    print(g)
    dev.off()
    
    mi[i] = min(selection$av_logPko)
    mx[i] = max(selection$av_logPko)
    st[i] = set
    ct[i] = nrow(selection)
    
  }
  
  minlogp = min(mi)
  maxlogp = max(mx)
  logp_range <- c(minlogp, maxlogp)
  
  stct <- data.frame(st, ct)
  
  extralist <- list("sets" = stct, "logp_range" = logp_range) # single argument returns only allowed, so a list is handy container
  return(extralist)
}

