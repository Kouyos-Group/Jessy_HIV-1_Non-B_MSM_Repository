#### Combine all the trees:
library(grid)
library(ggpubr)

ggarrange(allA, allC, allF,
          allG,allAE, allAG,
          ncol = 2, nrow = 3,
          labels = c("A", "B","C","D","E","F"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))

ggarrange(Aplot, Cplot, Fplot,
          Gplot,AEplot, AGplot,
          ncol = 3, nrow = 2,
          labels = c("A", "B","C","D","E","F"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))


alltrees<-ggarrange(ATREE, CTREE, FTREE,
          GTREE,AETREE, AGTREE,
          ncol = 2, nrow = 3,
          #labels = c("A", "B","C","D","E","F"), 
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 18, color = "black"))

bottom_row <- plot_grid(leg3, leg1,leg2, ncol=3) 

plot_grid(alltrees, 
          bottom_row, ncol=1, rel_heights = c(2,1)) 


## Save the plot:
# setwd("~/OneDrive/Master Thesis/Figures/Trees")
# ggsave("Tree_combined_legendbottom.tiff", width = 50, height = 50, units = "cm", limitsize = FALSE)


