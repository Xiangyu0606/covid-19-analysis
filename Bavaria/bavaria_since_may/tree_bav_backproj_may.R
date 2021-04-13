library(partykit)
library(ggplot2)
data_back_may = dat_bp_back_may
dat_back_may = data_back_may %>% mutate(t=1:n())
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = quasipoisson, start = start, ...)
}
#poisson with overdispersion:

dat_back_may$t2 = dat_back_may$t
segmented_backproj_may <- glmtree(backpro~t | t2,
                         data = dat_back_may, family = quasipoisson(link = "log"),minsize = 7
)
plot(segmented_backproj_may)

#extracting breakpoints:
#not sure how to extract breakpoints more easily:
bp_data_backproj_may = segmented_backproj_may$fitted
bp_data_backproj_may$t = 1:nrow(bp_data_backproj_may)
colnames(bp_data_backproj_may)
colnames(bp_data_backproj_may) = c("fitted", "t")
colnames(bp_data_backproj_may)
breakpoints_mob_may_backproj = (bp_data_backproj_may %>%  group_by(fitted) %>% summarise(max_t = max(t)))$max_t
breakpoints_mob_may_backproj = breakpoints_mob_may_backproj[1:(length(breakpoints_mob_may_backproj)-1)]
breakpoints_mob_may_backproj
predictions = predict(segmented_backproj_may, newdata = data.frame(t=seq(1,max(dat_back_may$t), 0.1),
                                                          t2 = seq(1,max(dat_back_may$t), 0.1)), type = "response"
)

#view for change points
dat_back_may[c(breakpoints_mob_may_backproj),]


mob4bp_backproj = ggplot() + 
  #ylim(0, 8000) +
  geom_col(data = dat_back_may, aes(x=t, y=backpro),
           col = "grey", fill = "lightgrey")+
  geom_line(aes(x=seq(1,max(dat_back_may$t), 0.1),y=predictions, 
                color = predictions),
            col = "red",  lwd = 1.1) +
  
  geom_vline(xintercept=breakpoints_mob_may_backproj, lty=2, col ="red", size=0.8)   +
  theme_bw()+
  scale_x_continuous(breaks = seq(max(dat_back_may$t), min(dat_back_may$t), by = -7), 
                     labels = strftime(seq(max(dat_back_may$date), min(dat_back_may$date), by = "-1 weeks"), 
                                       format = "%d.%m.")) + 
  ylim(0,400) +
  # geom_text(mapping=aes(x=break_points$bp,  y = 20, label = strftime(break_points$bp_date, format = "%d.%m.")), 
  #           size = 5, angle = 90, vjust =-.5, hjust = -0.1) +
  labs(x = "Date", y = "Number of estimated infections")
mob4bp_backproj
png("D:/Bavaria/bavaria_since_may/tree-backproj-since-may.png", width = 1200, height = 480)
mob4bp_backproj
dev.off()
