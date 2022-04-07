plot <- ggplot() + geom_point(aes(x = N_present, y = N_detectable)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Leopards Present") +
  ylab("Leopards Detectable") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black'))