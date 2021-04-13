library(partykit)
library(ggplot2)
data = dat_bp
dat = data %>% mutate(t=1:n())
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = quasipoisson, start = start, ...)
}
#poisson with overdispersion:

dat$t2 = dat$t
segmented <- glmtree(backpro~t | t2,
                     data = dat, family = quasipoisson(link = "log"),minsize = 5)



#extracting breakpoints:
#not sure how to extract breakpoints more easily:
bp_data = segmented$fitted
bp_data$t = 1:nrow(bp_data)
colnames(bp_data)
colnames(bp_data) = c("fitted", "t")
colnames(bp_data)
breakpoints_mob = (bp_data %>%  group_by(fitted) %>% summarise(max_t = max(t)))$max_t
breakpoints_mob = breakpoints_mob[1:(length(breakpoints_mob)-1)]
breakpoints_mob
predictions = predict(segmented, newdata = data.frame(t=seq(1,max(dat$t), 0.1),
                                                      t2 = seq(1,max(dat$t), 0.1)), type = "response")

dat[c(breakpoints_mob),]
mobtree = ggplot() + 
  #ylim(0, 8000) +
  geom_col(data = dat, aes(x=t, y=backpro),
           col = "grey", fill = "lightgrey")+
  geom_line(aes(x=seq(1,max(dat$t), 0.1),y=predictions, 
                color = predictions),
            col = "red",  lwd = 1.1) +
  
  geom_vline(xintercept=breakpoints_mob, lty=2, col ="red", size=0.8)   +
  theme_bw() +
  scale_x_continuous(breaks = seq(max(dat$t), min(dat$t), by = -7), 
                     labels = strftime(seq(max(dat$date), min(dat$date), by = "-1 weeks"), 
                                       format = "%d.%m.")) + 
  scale_y_continuous(expand=expansion(mult = c(0, .1))) +
  # geom_text(mapping=aes(x=break_points$bp,  y = 20, label = strftime(break_points$bp_date, format = "%d.%m.")), 
  #           size = 5, angle = 90, vjust =-.5, hjust = -0.1) +
  labs(x = "Date", y = "Number of disease onsets with model-based tree")
mobtree
png("D:/Bavaria/before_may (quasi-poisson)/results/disease_onset_mob.png", width = 1000, height = 460)
mobtree
dev.off()