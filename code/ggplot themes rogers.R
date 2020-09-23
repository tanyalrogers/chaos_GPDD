# personal ggplot themes

classic=theme_classic() + theme(panel.background = element_rect(color="black", fill=NULL), axis.text = element_text(color="black"), axis.ticks = element_line(color="black"))
xlabvert=theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
xlabnull=theme(axis.text.x=element_blank()) 
legalpha=guides(colour = guide_legend(override.aes = list(alpha = 1)))

removefacets=theme(strip.background = element_blank(),strip.text = element_blank())
removefacetbackground=theme(strip.background = element_blank())

