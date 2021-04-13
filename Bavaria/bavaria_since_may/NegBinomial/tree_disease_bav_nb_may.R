library(partykit)
library(ggplot2)
data = dat_bp_may_nb
dat = dat_bp_may_nb %>% mutate(t=1:n())
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  MASS::glm.nb(y ~ 0 + x,  start = start, ...)
}
#poisson with overdispersion:

dat$t2 = dat$t
segmented_1st1<- mob(backpro~t | t2,
                    data = dat, fit = logit,control = mob_control(minsize = 7))
plot(segmented_1st1)


#extracting breakpoints:
#not sure how to extract breakpoints more easily:
bp_data = segmented_1st1$fitted
bp_data$t = 1:nrow(bp_data)
colnames(bp_data)
colnames(bp_data) = c("fitted", "t")
colnames(bp_data)
breakpoints_mob = (bp_data %>%  group_by(fitted) %>% summarise(max_t = max(t)))$max_t
breakpoints_mob = breakpoints_mob[1:(length(breakpoints_mob)-1)]
breakpoints_mob #44 86 100
predictions = predict(segmented_1st1, type = "response")
#a view, where the change points are
dat[c(breakpoints_mob),] 

tree_nb_disease = ggplot() + 
  #ylim(0, 8000) +
  geom_col(data = dat, aes(x=t, y=backpro),
           col = "grey", fill = "lightgrey")+
  geom_line(aes(x=seq(1,max(dat$t)),y=predictions, 
                color = predictions),
            col = "red",  lwd = 1.0) +
  geom_vline(xintercept=breakpoints_mob, lty=2, col ="red", size=0.8)   +
  theme_bw() +
  scale_x_continuous(breaks = seq(max(dat$t), min(dat$t), by = -7), 
                     labels = strftime(seq(max(dat$date), min(dat$date), by = "-1 weeks"), 
                                       format = "%d.%m.")) + 
  scale_y_continuous(expand=expansion(mult = c(0, .1))) +
  # geom_text(mapping=aes(x=break_points$bp,  y = 20, label = strftime(break_points$bp_date, format = "%d.%m.")), 
  #           size = 5, angle = 90, vjust =-.5, hjust = -0.1) +
  labs(x = "Date", y = "Number of disease onsets with model-based tree(NegBinomial")
tree_nb_disease
png("D:/Bavaria/bavaria_since_may/NegBinomial/tree-disease-onset_NegBinomial-since-may.png", width = 1200, height = 480)
tree_nb_disease
dev.off()


