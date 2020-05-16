find.col.match <- function(x, tp = 1){
  library(dplyr)
  #Filtering timepoint
  x <- x %>% filter(TP == tp)
  #removing missing values
  ind.n <- unlist(lapply(x, function(x){length(x[!is.na(x)])}))
  #creates list
  matchable <- list()
  #cluster comparable variables
  for(i in 1:length(unique(ind.n))){
    matchable[i] <- list(ind.n[ind.n == unique(ind.n)[i]])
  }
  return(matchable)
}


covbox <- function(x = covid, y = c("IL27"), cat = c("Test_date2", "gender"),
                   fun = c(mean = mean, sd = sd), facet = NULL, rm.na = T){

  library(tidyverse)
  library(lubridate)
  library(tidytext)
  library(hrbrthemes)
  library(viridis)
  library(kableExtra)
  library(dplyr)
  library(scales)
  library(knitr)
  library(plotly)
  library(plyr)


  #Remove missing?
  if(rm.na == T){
    x <- x[!is.na(x$Test_date2),]
  }

  #filter patient
  patients <- split(x, f = x$ID2)
  #patients <- patients[lapply(patients, function(x){nrow(x)}) > 1]

  patients <- lapply(patients,
                     function(x){
                       level_key <- 1:length(x$Test_date)
                       names(level_key) <-   x$Test_date
                       x$Test_date2 <- level_key
                       return(x)
                     })


  # my color palette
  # my.phase.palette <- hue_pal()(length(patients))

  bind <- bind_rows(patients, .id="df")
  bind$Test_date2 <- as.character(bind$Test_date2)


  # Plot
  p <- bind %>%
    ggplot(aes(x=Test_date2, y=get(y), fill=Test_date2)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("Var Kinetics") +
    xlab("Timepoint") +
    ylab(y)
  if(length(facet) == !0){
    p <- p + facet_wrap(~as.character(get(facet)), scale="free", )
  }

  return(list(plot = p, dt = bind))
}



make.bin<- function(x = x , t = "DaysfromSxOnset", bin = 2){
  library(dplyr)
  x <- x[!is.na(x[[t]]),]
  time <- x[[t]]

  bins <- cut(x = x[[t]], include.lowest = T,
              breaks = seq(min(time), max(time)+2, bin))

  x[["Timepoints"]] <- cut(x = x[[t]], include.lowest = T,
                           breaks = seq(min(time), max(time)+2, bin), labels = 1:length(levels(bins)))
  x[["Timepoints"]] <- as.numeric(x[["Timepoints"]])

  return(x)
}


na.perc <- function(x){

  p <- table(is.na(x))/length(x);
  p <- p*100

  if(length(p) == 1){
    if(names(p) == FALSE){
      p["TRUE"] <- 0}
    if(names(p) == TRUE){
      p["FALSE"] <- 0}
  }
  return(p[["TRUE"]])

}

filt.sh <- function(x = x, vars = y, cond1 = z){
  library(dplyr)
  x %>% dplyr::filter(.data[[vars[[1]]]] %in% cond1)}


test.anova <- function(y = "IL6", x = "Timepoints", group = "ICU", data = data){

  if(group %in% c("ICU", "gender")){
    mod1 <- lm(get(y) ~ get(x),  data = data)
    mod2 <- lm(get(y) ~ get(x) + get(group), data = data)
    return(anova(mod1, mod2))}

  if(group %in% c("genderICU")){
    mod1 <- lm(get(y) ~ get(x), data = data)
    mod2 <- lm(get(y) ~ get(x) + ICU, data = data)
    mod3 <- lm(get(y) ~ get(x) + ICU + gender, data = data)

    return(anova(mod1, mod2, mod3))}
}


r2.lm <- function(x = covid.trim.groups, c = "IL6"){
  library(dplyr)

  mod <- list()
  s <- list()
  adj.r.squared <- list()

  for(i in names(x)){
    mod[[i]] <- lm(get(c) ~ Timepoints, data = x[[i]]) }

  for(i in names(mod)){
    s[[i]] <- summary(mod[[i]]) }

  for(i in names(mod)){
    adj.r.squared[[i]] <- s[[i]]$adj.r.squared }

  val <- min(unlist(adj.r.squared))-max(unlist(adj.r.squared))

  return(val)
}

calc_corr <- function(x = covid.trim.groups,
                      names4.r2vals = names4.r2vals,
                      groups = groups,
                      cv.filt=cv.filt){
  r2vals <- list()
  aov <- list()

  for(i in names4.r2vals){
    r2vals[[i]] <- covid::r2.lm(x = x, c = i)
  }

  r2vals <- unlist(r2vals)


  if(groups == "genderICU"){ r2vals <- r2vals[-276]}

  r2vals[r2vals < 0] <- r2vals[r2vals < 0] * -1
  r2vals <- sort(r2vals, decreasing = T)

  for(i in names4.r2vals){

    aov.test <- test.anova(y = i, x = "Timepoints", group = groups, data = cv.filt)
    aov[[i]] <- aov.test$`Pr(>F)`[2]
  }

  aov <- unlist(aov)
  aov <- as.data.frame(aov)

  r2vals <- as.data.frame(r2vals)

  r2vals$meas <- rownames(r2vals)
  idx <- match(rownames(aov), rownames(r2vals))
  r2vals <- cbind(aov, r2vals[idx,])

  colnames(r2vals) <- c("aov_pvalue", "R2 Diff", "Observations")

  return(r2vals)

}









#  Rvals$r2vals <- as.data.frame(readRDS("r2vals.rds"))
