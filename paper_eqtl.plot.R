source("paper_eqt.main.R")

## plot for YRI vs. pigmentation
YRI.pigmentation.plot = ggplot(unique(all_gene.expression_norm.resid.df[,c("Sample","YRI","pigment.lm.residuals")]),aes(x=log(YRI),y=pigment.lm.residuals)) + geom_point() + geom_smooth(method=lm,se=F,fullrange=T) + theme_bw() + scale_y_continuous(expand=c(0.1,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(margin = margin(r=20)), axis.title.x = element_text(margin = margin(t=20))) + ggtitle("YRI ancestry vs. pigmentation") + ylab("Norm. expression")
