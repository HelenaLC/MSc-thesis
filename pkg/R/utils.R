# utilities
# ==============================================================================

# aesthetics for plotting
methods <- c(
    "diffcyt-limma", "diffcyt-LMM", 
    "edgeR-rawCounts", "edgeR-normCounts", "edgeR-scaledCPM")
method_colors <- c("darkblue", "cornflowerblue", "red3", "coral", "gold")
method_colors <- setNames(method_colors, methods)

prettify <- theme(
    aspect.ratio=1,
    plot.margin=unit(c(0,.1,0,0), "cm"),
    panel.grid=element_blank(),
    axis.title=element_text(size=8),
    axis.text=element_text(color="black", size=6),
    strip.background=element_blank(),
    strip.text=element_text(face="bold", hjust=0, size=8),
    plot.title=element_text(face="bold", hjust=0, size=8),
    legend.margin=margin(0,0,0,0, "cm"),
    legend.title=element_text(size=8),
    legend.text=element_text(size=6),
    legend.key.size=unit(0.5, "cm"))