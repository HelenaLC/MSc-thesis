# DD paterns
# ==============================================================================
# Generates a schematic of the DD patterns proposed in
# Korthauer et al. "A statistical approach for identifying 
# differential distributions in single-cell rna-seq experiments." 
# Genome Biology, 17:222, 2016.
# ------------------------------------------------------------------------------

library(ggplot2)
library(reshape2)

x <- seq(0, 6, .01)
df <- data.frame(x=x,
  DE=c(dnorm(x,2,.5), dnorm(x,4,.5)),
  DP=c(.6*dnorm(x,2,.5)+.4*dnorm(x,4,.5), .4*dnorm(x,2,.5)+.6*dnorm(x,4,.5)),
  DM=c(.8*dnorm(x,2,.5), .4*dnorm(x,2,.5)+.6*dnorm(x,4,.5)),
  DB=c(.5*dnorm(x,2,.5)+.5*dnorm(x,4,.5), .7*dnorm(x,3,1)),
  EP=c(.6*dnorm(x,2,.5)+.4*dnorm(x,4,.5), .6*dnorm(x,2,.4)+.4*dnorm(x,4,.6)),
  EE=c(dnorm(x,3,.45), dnorm(x,3,.55)),
  condition=rep(c("A", "B"), each=length(x)))
df[, 2:7] <- t(t(df[, 2:7]) / colMaxs(as.matrix(df[, 2:7])))
df <- melt(df, id.vars=c("x", "condition"), variable.name="category")

p <- ggplot(df, aes(x=x, y=value, col=condition)) + 
  geom_line(size=.8) + guides(color=FALSE) +
  facet_wrap(~ category, scales="free", nrow=1) + 
  scale_x_continuous(limits=c(0, 6), expand=c(.04,0)) +
  scale_y_continuous(limits=c(0, 1), expand=c(.06,0)) +
  theme_classic() + theme(aspect.ratio=2/3,
    panel.background=element_blank(),
    plot.margin=unit(c(0,.1,0,0),"cm"),
    strip.background=element_blank(),
    strip.text=element_text(size=8, face="bold", hjust=0),
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    axis.title=element_blank())
ggsave("figures/dd_patters.pdf", 
    plot=p, width=14, height=2.1, units="cm")
