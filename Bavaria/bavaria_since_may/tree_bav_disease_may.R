library(partykit)
library(ggplot2)
data_may = dat_bp_may
dat_may = data_may %>% mutate(t=1:n())
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = quasipoisson, start = start, ...)
}
#poisson with overdispersion:

dat_may$t2 = dat_may$t
segmented_may <- glmtree(backpro~t | t2,
                     data = dat_may, family = quasipoisson(link = "log"),minsize= 6)
plot(segmented_may)


#extracting breakpoints:
#not sure how to extract breakpoints more easily:
bp_data_disease_may = segmented_may$fitted
bp_data_disease_may$t = 1:nrow(bp_data_disease_may)
colnames(bp_data_disease_may)
colnames(bp_data_disease_may) = c("fitted", "t")
colnames(bp_data_disease_may)
breakpoints_mob_may = (bp_data_disease_may %>%  group_by(fitted) %>% summarise(max_t = max(t)))$max_t
breakpoints_mob_may = breakpoints_mob_may[1:(length(breakpoints_mob_may)-1)]
breakpoints_mob_may
predictions = predict(segmented_may, newdata = data.frame(t=seq(1,max(dat_may$t), 0.1) , 
                                                          t2 = seq(1,max(dat_may$t), 0.1)), type = "response"
)

#view for change points   
dat_may[c(breakpoints_mob_may),]


mob4p = ggplot() + 
  #ylim(0, 8000) +
  geom_col(data = dat_may, aes(x=t, y=backpro),
           col = "grey", fill = "lightgrey")+
  geom_line(aes(x=seq(1,max(dat_may$t), 0.1),y=predictions, 
                color = predictions),
            col = "red",  lwd = 1.1) +
  
  geom_vline(xintercept=breakpoints_mob_may, lty=2, col ="red", size=0.8)   +
  theme_bw() +
  ylim(0,400)+
  scale_x_continuous(breaks = seq(max(dat_may$t), min(dat_may$t), by = -7), 
                     labels = strftime(seq(max(dat_may$date), min(dat_may$date), by = "-1 weeks"), 
                                       format = "%d.%m.")) + 
  labs(x = "Date", y = "Number of disease onsets with model-based tree")
mob4p
png("D:/Bavaria/bavaria_since_may/tree-disease-onset-since-may.png", width = 1200, height = 480)
mob4p
dev.off()

