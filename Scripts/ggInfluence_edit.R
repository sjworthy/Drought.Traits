ggInfluence_test <-
function (gbm.object, col.bar = "skyblue3", show.signif = TRUE, 
          col.signif = "#8B3A3A", main = gbm.call$response.name, x.label = "Relative influence (%)", 
          signif = FALSE, plot = T, ...) {
  gbm.call <- gbm.object$gbm.call
  if (signif == F) {
    dfContr <- data.frame(gbm.object$contributions)
    Influence <- dfContr[2]
    ggInfluence <- ggplot(dfContr, aes(x = reorder(var, rel.inf), 
                                       y = rel.inf)) + geom_bar(fill = col.bar, stat = "identity") + 
      coord_flip() + theme(axis.title = element_text(size = 15), 
                           axis.text = element_text(face = "plain",size = 12), plot.title = element_text(size = 16, 
                                                                                               face = "bold"), axis.line.x = element_line(size = 0.3), 
                           panel.background = element_rect(fill = "white"), 
                           panel.grid.major.y = element_line(size = 0.2, linetype = "dotted", 
                                                             color = "grey"), panel.grid.major.x = element_line(linetype = "blank"), 
                           panel.grid.minor = element_line(linetype = "blank")) + 
      ylab(x.label) + xlab("") + ggtitle(main)
    if (show.signif == TRUE) {
      ggInfluence <- ggInfluence + geom_hline(yintercept = 100/length(dfContr$var), 
                                              colour = col.signif, linetype = "dashed")
    }
    if (plot == T) {
      print(ggInfluence)
    }
    return(Influence)
  }
  else {
    dfContr <- data.frame(gbm.object$contributions)
    dfContr <- subset(dfContr, dfContr[[2]] > (100/length(dfContr[[1]])))
    dfContr[[2]] <- dfContr[[2]] * 100/sum(dfContr[[2]])
    Influence <- dfContr[2]
    ggInfluence <- ggplot(dfContr, aes(x = reorder(var, rel.inf), 
                                       y = rel.inf)) + geom_bar(fill = col.bar, stat = "identity") + 
      coord_flip() + theme(axis.title = element_text(size = 13), 
                           axis.text = element_text(face = "plain"), plot.title = element_text(size = 14, 
                                                                                               face = "bold"), axis.line.x = element_line(size = 0.3), 
                           panel.background = element_rect(fill = "white"), 
                           panel.grid.major.y = element_line(size = 0.2, linetype = "dotted", 
                                                             color = "grey"), panel.grid.major.x = element_line(linetype = "blank"), 
                           panel.grid.minor = element_line(linetype = "blank")) + 
      ylab(x.label) + xlab("") + ggtitle(main)
    if (plot == T) {
      print(ggInfluence)
    }
    return(Influence)
  }
}
