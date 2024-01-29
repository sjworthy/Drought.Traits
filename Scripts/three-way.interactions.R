gbm.call <- BRT.final.sub.3.4$gbm.call
n.trees <- gbm.call$best.trees
depth <- gbm.call$interaction.depth
gbm.x <- gbm.call$gbm.x
n.preds <- length(gbm.x)
pred.names <- gbm.call$interaction.depthe$gbm.call$predictor.names
data <- gbm.call$dataframe[, gbm.x]

x.var.pred.1 = expand.grid(seq(min(data$log.leafN, na.rm = T), max(data$log.leafN, na.rm = T), length = 20),
                               seq(min(data$log.height, na.rm = T), max(data$log.height, na.rm = T), length = 20))

x.var.pred.2 = expand.grid(seq(min(data$log.rootN, na.rm = T), max(data$log.rootN, na.rm = T), length = 20),
                           seq(min(data$log.SLA, na.rm = T), max(data$log.SLA, na.rm = T), length = 20))
x.var.pred.3 = expand.grid(seq(min(data$log.root.depth, na.rm = T), max(data$log.root.depth, na.rm = T), length = 400))

x.var.pred.3 = expand.grid(seq(min(data$log.root.depth, na.rm = T), max(data$log.root.depth, na.rm = T), length = 20),
                           seq(min(data$log.RTD, na.rm = T), max(data$log.RTD, na.rm = T), length = 20))
x.var.pred.4 = expand.grid(seq(min(data$log.SRL, na.rm = T), max(data$log.SRL, na.rm = T), length = 20),
                           seq(min(data$log.root.diam, na.rm = T), max(data$log.root.diam, na.rm = T), length = 20))
                               
x.var.pred.frame = cbind(x.var.pred.1,x.var.pred.2,x.var.pred.3,x.var.pred.4) 
x.var.pred.frame = cbind(x.var.pred.2, x.var.pred.3)
colnames(x.var.pred.frame)=c("log.rootN","log.SLA","log.root.depth")

predict.set.1 = x.var.pred.frame[,c(1:3)]
for(i in 4:8){
predict.set.1[,i] <- mean(data[,i], na.rm = T)
}
colnames(predict.set.1)=pred.names

prediction.set.1 <- gbm::predict.gbm(BRT.final.all.variable, predict.set.1, 
                               n.trees = n.trees, type = "link")

prediction.set.1 <- gbm::predict.gbm(BRT.final.sub.3.4, x.var.pred.frame, 
                                     n.trees = n.trees, type = "link")

interaction.test.model.set.1 <- lm(prediction.set.1 ~ as.factor(x.var.pred.frame[, 1]) + 
                                     as.factor(x.var.pred.frame[, 2]) + as.factor(x.var.pred.frame[, 3]))

interaction.flag.set.1 <- round(mean(resid(interaction.test.model.set.1)^2) * 
                            1000, 2)

# repeat code above for each set of 3-way interactions
# create an index of the values in descending order

search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")

n.important <- max(2,round(0.1 * ((n.preds^2)/2),0))



