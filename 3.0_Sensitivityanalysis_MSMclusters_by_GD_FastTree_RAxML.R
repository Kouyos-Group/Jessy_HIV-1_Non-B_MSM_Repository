
#########
### Figure: Dependence of MSM cluster extraction on genetic distance threshold

library(grid)
library(ggpubr)
####FASTTREE

# Or specify the factor levels in the order you want
table$subtype <- factor(table$subtype, levels = c("CRF01_AE","CRF02_AG","A","C","F","G"))

plot_freq <-ggplot(data=table, aes(x=distance, y=countMSM, color= subtype)) +
  geom_point(size=5) +
  geom_path(aes(group = subtype)) + 
  scale_x_discrete(name ="Distance [%]",
                   breaks= paste0("d",seq(10,50,5)),
                   labels= seq(1.0,5.0,0.5)) +
  scale_y_continuous(breaks = c(1,seq(5,35,5)), limits = c(1,30))+
  labs(x="Distance",  y="Frequency")+
  scale_color_manual(values=c("darkorange3","darkgoldenrod2","violetred2","cyan4","olivedrab","mediumpurple3")) +
  guides(color=guide_legend(title = "Subtype", reverse = F, size= 14))+ 
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")
plot_freq

table$percMSM

plot_fraction <-ggplot(data=table, aes(x=distance, y=percMSM, color= subtype)) +
  geom_point(size=5) +
  geom_path(aes(group = subtype)) + 
  scale_x_discrete(name ="Distance [%]",
                   breaks= paste0("d",seq(10,50,5)),
                   labels= seq(1.0,5.0,0.5)) +
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0, 1))+
  labs(x="Distance",  y="Fraction")+
  scale_color_manual(values=c("darkorange3","darkgoldenrod2","violetred2","cyan4","olivedrab","mediumpurple3")) +
  guides(color=guide_legend(title = "Subtype", reverse = F, size= 14))+ 
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
plot_fraction

ggarrange(plot_freq, plot_fraction,ncol = 2, nrow = 1,
          labels = c("A", "B"), 
          common.legend = TRUE, legend = "right",
          font.label = list(size = 18, color = "black"))



####RAxML

# Or specify the factor levels in the order you want
tableMSM$subtype <- factor(tableMSM$subtype, levels = c("CRF01_AE","CRF02_AG","A","C","F","G"))

plot_freq_raxml <-ggplot(data=tableMSM, aes(x=distance, y=countMSM, color= subtype)) +
  geom_point(size=5) +
  geom_path(aes(group = subtype)) + 
  scale_x_discrete(name ="Distance [%]",
                   breaks= paste0("d",seq(10,50,5)),
                   labels= seq(1.0,5.0,0.5)) +
  scale_y_continuous(breaks = c(1,seq(5,35,5)), limits = c(1,20))+
  labs(x="Distance",  y="Frequency")+
  scale_color_manual(values=c("darkorange3","darkgoldenrod2","violetred2","cyan4","olivedrab","mediumpurple3")) +
  guides(color=guide_legend(title = "Subtype", reverse = F, size= 14))+ 
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")
plot_freq_raxml



plot_fraction_raxml <-ggplot(data=tableMSM, aes(x=distance, y=percMSM, color= subtype)) +
  geom_point(size=5) +
  geom_path(aes(group = subtype)) + 
  scale_x_discrete(name ="Distance [%]",
                   breaks= paste0("d",seq(10,50,5)),
                   labels= seq(1.0,5.0,0.5)) +
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0, 1))+
  labs(x="Distance",  y="Fraction")+
  scale_color_manual(values=c("darkorange3","darkgoldenrod2","violetred2","cyan4","olivedrab","mediumpurple3")) +
  guides(color=guide_legend(title = "Subtype", reverse = F, size= 14))+ 
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    #legend.position = c(0.9, 0.8),
    legend.text = element_text(colour="black", size=14, face="plain"),
    legend.title = element_text(colour="black", size=14, face="bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
plot_fraction_raxml

ggarrange(plot_freq_raxml, plot_fraction_raxml,ncol = 2, nrow = 1,
          labels = c("A", "B"), 
          common.legend = TRUE, legend = "right",
          font.label = list(size = 18, color = "black"))



### FastTree (A,B) and RAxML (C,D)
ggarrange(plot_freq, plot_fraction,
          plot_freq_raxml, plot_fraction_raxml,
          ncol = 2, nrow = 2,
          labels = c("A", "B","C","D"), 
          common.legend = TRUE, legend = "right",
          font.label = list(size = 18, color = "black"))



  
  