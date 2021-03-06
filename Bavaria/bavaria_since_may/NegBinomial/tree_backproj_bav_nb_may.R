library(partykit)
library(ggplot2)
data = dat_bp_bav_back_may_nb
dat = data %>% mutate(t=1:n())
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  MASS::glm.nb(y ~ 0 + x,  start = start, ...)
}
#poisson with overdispersion:

dat$t2 = dat$t

segmented_min<-mob(backpro~t | t2,
                   data = dat, fit = logit,control = mob_control(minsize = 7))

plot(segmented_min)


#extracting breakpoints:
#not sure how to extract breakpoints more easily:
bp_data = segmented_min$fitted
bp_data$t = 1:nrow(bp_data)
colnames(bp_data)
colnames(bp_data) = c("fitted", "t")
colnames(bp_data)
breakpoints_mob = (bp_data %>%  group_by(fitted) %>% summarise(max_t = max(t)))$max_t
breakpoints_mob = breakpoints_mob[1:(length(breakpoints_mob)-1)]
breakpoints_mob #16,23,30,80,94,121 

#a view, where the change points are
dat[c(breakpoints_mob),] 
predictions = predict(segmented_min, type = "response")



tree_nb_backproj = ggplot() + 
  #ylim(0, 8000) +
  geom_col(data = dat, aes(x=t, y=backpro),
           col = "grey", fill = "lightgrey")+
  geom_line(aes(x=seq(1,max(dat$t)),y=predictions, 
                color = predictions),
            col = "red",  lwd = 1.1) +
  
  geom_vline(xintercept=breakpoints_mob, lty=2, col ="red", size=0.8)   +
  theme_bw() +
  scale_x_continuous(breaks = seq(max(dat$t), min(dat$t), by = -7), 
                     labels = strftime(seq(max(dat$date), min(dat$date), by = "-1 weeks"), 
                                       format = "%d.%m.")) + 
  ylim(0,400)+
  # geom_text(mapping=aes(x=break_points$bp,  y = 20, label = strftime(break_points$bp_date, format = "%d.%m.")), 
  #           size = 5, angle = 90, vjust =-.5, hjust = -0.1) +
  labs(x = "Date", y = "Number of estimated infections")
tree_nb_backproj
png("D:/Bavaria/bavaria_since_may/NegBinomial/tree-backproj_NegBinomial-since-may.png", width = 1200, height = 480)
tree_nb_backproj
dev.off()