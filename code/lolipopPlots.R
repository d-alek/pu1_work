# make lolipop type plots with hallmark, kegg, and motif set test logP values:
### Hallmarks
hallmarks <- read.table("lolipop_halmark28.txt", header=T, sep="\t", quote="", as.is=TRUE)
require(tidyverse)
hallmarks <- mutate(hallmarks, l.rom.mix.ko = -log10(rom.mix.ko), l.rom.dn.ko = log10(rom.dn.ko), l.rom.up.ko = -log10(rom.up.ko), l.cam.dn.ko = log10(cam.dn.ko), l.cam.up.ko = -log10(cam.up.ko), l.rom.mix.aml = -log10(rom.mix.aml), l.rom.dn.aml = log10(rom.dn.aml), l.rom.up.aml = -log10(rom.up.aml),  l.cam.dn.aml = log10(cam.dn.aml), l.cam.up.aml = -log10(cam.up.aml))
hallmarks <- arrange(hallmarks, desc(Id))
hallmarks <- mutate(hallmarks, set = factor(Set, levels = hallmarks$Set), class = factor(Class, levels = rev(unique(hallmarks$Class))))
# !!! this last bit factorises set and class variables according to my original Id order - to override alphabeical ordering

# ko H mix:
pdf("Hallmark_ko_mix.pdf", width = 7, height = 5, useDingbats = FALSE) # useDingbats = FALSE (for .pdf); units = "in", res = 300 (for .tiff, .png)
g <- ggplot(hallmarks, aes(l.rom.mix.ko, set)) + 
  geom_segment(aes(x = 0, y = set, xend = l.rom.mix.ko, yend = set, colour = class), size = 0.3) + 
  geom_point(aes(colour = class)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key = element_rect(colour = "white"), legend.title = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "WT-KO contrast \nHallmark sets with differentially expressed genes") + 
  geom_vline(xintercept = 1.3, linetype = "dotted", size = 0.25) + 
  scale_x_continuous(limits = c(0, 3.75), breaks = c(0, 3), labels = c("0", "3"))
print(g)
dev.off()

# ko H directional:
pdf("Hallmark_ko_directional.pdf", width = 8, height = 5, useDingbats = FALSE)
g <- ggplot(hallmarks, aes(y = set)) + 
  geom_segment(aes(x = 0, y = set, xend = l.rom.dn.ko, yend = set), size = 0.325, colour = "blue") + 
  geom_point(aes(x = l.rom.dn.ko), colour = "blue") + 
  geom_segment(aes(x = 0, y = set, xend = l.rom.up.ko, yend = set), size = 0.325, colour = "red") + 
  geom_point(aes(x = l.rom.up.ko), colour = "red") + 
  geom_segment(aes(x = 0, y = set, xend = l.cam.dn.ko, yend = set), size = 0.15, colour = "blue", linetype = "dashed") + 
  geom_point(aes(x = l.cam.dn.ko), colour = "blue", shape = 1) + 
  geom_segment(aes(x = 0, y = set, xend = l.cam.up.ko, yend = set), size = 0.15, colour = "red", linetype = "dashed") +
  geom_point(aes(x = l.cam.up.ko), colour = "red", shape = 1) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "WT-KO contrast \nHallmark sets with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  #geom_hline(yintercept = c(3.5, 12.5, 14.5, 19.5, 22.5, 27.5), size = 0.2) +
  scale_x_continuous(limits = c(-5.4025, 5.25), breaks = c(-3, 0, 3), labels = c("3", "0", "3")) # width ~ 10.5
print(g)
dev.off()

# aml H mix:
pdf("Hallmark_aml_mix.pdf", width = 7, height = 5, useDingbats = FALSE)
g <- ggplot(hallmarks, aes(l.rom.mix.aml, set)) + 
  geom_segment(aes(x = 0, y = set, xend = l.rom.mix.aml, yend = set, colour = class), size = 0.3) + 
  geom_point(aes(colour = class)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key = element_rect(colour = "white"), legend.title = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "KO-AML contrast \nHallmark sets with differentially expressed genes") + 
  geom_vline(xintercept = 1.3, linetype = "dotted", size = 0.25) + 
  scale_x_continuous(limits = c(0, 3.75), breaks = c(0, 3), labels = c("0", "3"))
print(g)
dev.off()

# aml H directional:
pdf("Hallmark_aml_directional.pdf", width = 11.5, height = 5, useDingbats = FALSE)
g <- ggplot(hallmarks, aes(y = set)) + 
  geom_segment(aes(x = 0, y = set, xend = l.rom.dn.aml, yend = set), size = 0.325, colour = "blue") + 
  geom_point(aes(x = l.rom.dn.aml), colour = "blue") + 
  geom_segment(aes(x = 0, y = set, xend = l.rom.up.aml, yend = set), size = 0.325, colour = "red") + 
  geom_point(aes(x = l.rom.up.aml), colour = "red") + 
  geom_segment(aes(x = 0, y = set, xend = l.cam.dn.aml, yend = set), size = 0.15, colour = "blue", linetype = "dashed") + 
  geom_point(aes(x = l.cam.dn.aml), colour = "blue", shape = 1) + 
  geom_segment(aes(x = 0, y = set, xend = l.cam.up.aml, yend = set), size = 0.15, colour = "red", linetype = "dashed") +
  geom_point(aes(x = l.cam.up.aml), colour = "red", shape = 1) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "KO-AML contrast \nHallmark sets with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  #geom_hline(yintercept = c(3.5, 12.5, 14.5, 19.5, 22.5, 27.5), size = 0.2) +
  scale_x_continuous(limits = c(-5, 11.706), breaks = c(-3, 0, 3, 6, 9), labels = c("3", "0", "3", "6", "9")) # width ~ 16.5
print(g)
dev.off()

#### KEGGs
keggs <- read.table("lolipop_KEGG37.txt", header=T, sep="\t", quote="", as.is=TRUE)
keggs <- mutate(keggs, l.up.ko = -log10(up.ko), l.dn.ko = log10(dn.ko), l.up.aml = -log10(up.aml), l.dn.aml = log10(dn.aml))
keggs <- arrange(keggs, desc(Id))
keggs <- mutate(keggs, path = factor(pathway, levels = keggs$pathway), class = factor(Class, levels = rev(unique(keggs$Class))))

# ko KEGG directional:
pdf("KEGG_ko_directional.pdf", width = 8, height = 6.5, useDingbats = FALSE)
g <- ggplot(keggs, aes(y = path)) + 
  geom_segment(aes(x = 0, y = path, xend = l.dn.ko, yend = path), size = 0.3, colour = "blue") + 
  geom_point(aes(x = l.dn.ko), colour = "blue") + 
  geom_segment(aes(x = 0, y = path, xend = l.up.ko, yend = path), size = 0.3, colour = "red") + 
  geom_point(aes(x = l.up.ko), colour = "red") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "WT-KO contrast \nKEGG pathways with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  scale_x_continuous(limits = c(-28, 28.5), breaks = c(-20, -10, 0, 10, 20), labels = c("20", "10", "0", "10", "20")) 
print(g)
dev.off()

# aml KEGG directional:
pdf("KEGG_aml_directional.pdf", width = 9.2, height = 6.5, useDingbats = FALSE)
g <- ggplot(keggs, aes(y = path)) + 
  geom_segment(aes(x = 0, y = path, xend = l.dn.aml, yend = path), size = 0.3, colour = "blue") + 
  geom_point(aes(x = l.dn.aml), colour = "blue") + 
  geom_segment(aes(x = 0, y = path, xend = l.up.aml, yend = path), size = 0.3, colour = "red") + 
  geom_point(aes(x = l.up.aml), colour = "red") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "KO-AML contrast \nKEGG pathways with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  scale_x_continuous(limits = c(-39.21, 28.52), breaks = c(-30, -20, -10, 0, 10, 20), labels = c("30", "20", "10", "0", "10", "20")) # width ~ 10.5
print(g)
dev.off()

# plus create one KEGG for 8 classes colour scheme:
pdf("KEGG_classes_sup.pdf", width = 9.2, height = 6.5, useDingbats = FALSE)
g <- ggplot(keggs, aes(y = path)) + 
  geom_segment(aes(x = 0, y = path, xend = l.dn.aml, yend = path, colour = class), size = 0.3) + 
  geom_point(aes(x = l.dn.aml, colour = class)) + 
  geom_segment(aes(x = 0, y = path, xend = l.up.aml, yend = path, colour = class), size = 0.3) + 
  geom_point(aes(x = l.up.aml, colour = class)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11), legend.key = element_rect(colour = "white"), legend.title = element_blank()) + 
  labs(x = "-log10 p-value", y = NULL, title = "KO-AML contrast \nKEGG pathways with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  scale_x_continuous(limits = c(-39.21, 28.52), breaks = c(-30, -20, -10, 0, 10, 20), labels = c("30", "20", "10", "0", "10", "20")) # width ~ 10.5
print(g)
dev.off()

#### Motifs
motifs <- read.table("lolipop_motif20.txt", header=T, sep="\t", quote="", as.is=TRUE)
motifs <- mutate(motifs, l.rom.up.ko = -log10(rom.up.ko), l.rom.dn.ko = log10(rom.dn.ko), l.rom.up.aml = -log10(rom.up.aml), l.rom.dn.aml = log10(rom.dn.aml))
motifs <- arrange(motifs, desc(Id))
motifs <- mutate(motifs, motif = factor(Motif, levels = motifs$Motif))

# ko M directional:
pdf("Motif_ko_directional.pdf", width = 4.5, height = 3.55, useDingbats = FALSE)
g <- ggplot(motifs, aes(y = motif)) + 
  geom_segment(aes(x = 0, y = motif, xend = l.rom.dn.ko, yend = motif), size = 0.3, colour = "blue") + 
  geom_point(aes(x = l.rom.dn.ko), colour = "blue") + 
  geom_segment(aes(x = 0, y = motif, xend = l.rom.up.ko, yend = motif), size = 0.3, colour = "red") + 
  geom_point(aes(x = l.rom.up.ko), colour = "red") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "WT-KO contrast \nMotif sets with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  scale_x_continuous(limits = c(-3.5, 3.4), breaks = c(-3, 0, 3), labels = c("3", "0", "3"))
print(g)
dev.off()

# aml M directional:
pdf("Motif_aml_directional.pdf", width = 4.5, height = 3.55, useDingbats = FALSE)
g <- ggplot(motifs, aes(y = motif)) + 
  geom_segment(aes(x = 0, y = motif, xend = l.rom.dn.aml, yend = motif), size = 0.3, colour = "blue") + 
  geom_point(aes(x = l.rom.dn.aml), colour = "blue") + 
  geom_segment(aes(x = 0, y = motif, xend = l.rom.up.aml, yend = motif), size = 0.3, colour = "red") + 
  geom_point(aes(x = l.rom.up.aml), colour = "red") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) + 
  labs(x = "-log10 p-value", y = NULL, title = "KO-AML contrast \nMotif sets with down- and up-regulated genes") + 
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dotted", size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  scale_x_continuous(limits = c(-3.5, 3.4), breaks = c(-3, 0, 3), labels = c("3", "0", "3"))
print(g)
dev.off()

# plot TF levels.. simple line plot
