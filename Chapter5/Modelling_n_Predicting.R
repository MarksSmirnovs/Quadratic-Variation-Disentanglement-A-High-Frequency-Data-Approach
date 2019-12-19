# Source ----------------------------------------------------------------
source("functions.R")

#   -----------------------------------------------------------------------

data <- readRDS("Application/predictors.rds")



plot(data$return_1)

train_set <- data %>% filter(Date < min(tail(data$Date,round(nrow(data)*0.25))))
test_set <- data %>% filter(Date >= min(tail(data$Date,round(nrow(data)*0.25))))

md_fit <- lm(formula = QV ~ IV_1 + JV_1 + return_1 + 
     VIX_1 + IV_5 + JV_5 + return_5 +
     VIX_5 + IV_22 + I(IV_22^2) + JV_22 + return_22 +VIX_22+ I(VIX_22^2), train_set)


#Functions needed for making 10^ notation in the tables.
mantissa <- function(x) {
  if (x == 0) 0
  else {
    log <- log10(abs(x))
    10^(log - floor(log))
  }
}

exponent <- function(x) {
  if (x == 0) 0
  else floor(log10(abs(x)))
}



mantissas=apply( summary(md_fit)$coefficients, c(1,2), mantissa)%>%round(4)
exponents=apply( summary(md_fit)$coefficients, c(1,2), exponent)


table <- summary(md_fit)$coefficients%>%round(4)

index=which(abs(table)<=0.001)

table[index]=paste0("$",mantissas[index],"\\cdot10^{",exponents[index],"}","$")


rownames(table)[2:nrow(table)] <- paste0("$", rownames(table)[2:nrow(table)], "$")
kableExtra::kable(table, "latex", booktabs = T , escape = F) %>% kable_styling()

# Model validation for full model--------------------------------------------------------



resids <- md_fit$residuals %>% data.frame
qqfig <- ggplot(resids, aes(sample = .))
qqfig <- qqfig + stat_qq() + stat_qq_line(col = "red") +
  xlab("Theoretical Quantiles") + ylab("Residuals Quantiles")  + ggtitle("Q-Q Plot")

resid_fit <- ggplot() + geom_point(aes(y = md_fit$residuals, x = md_fit$fitted.values)) +
  ylab("Residuals") + xlab("Fitted Values") + ggtitle("Residuals and Fitted Values")

gridExtra::grid.arrange(
  qqfig, resid_fit, layout_matrix = rbind(c(1,2))
)


# Backwards elimination Based on AIC. --------------------------------------------------


bkwrd_select <- step(md_fit, direction = "backward",k=2)#AIC

mantissas=apply( summary(bkwrd_select)$coefficients, c(1,2), mantissa)%>%round(4)

exponents=apply( summary(bkwrd_select)$coefficients, c(1,2), exponent)


table <- summary(bkwrd_select)$coefficients%>%round(4)

index=which(abs(table)<=0.001)

table[index]=paste0("$",mantissas[index],"\\cdot10^{",exponents[index],"}","$")


rownames(table)[2:nrow(table)] <- paste0("$", rownames(table)[2:nrow(table)], "$")
kableExtra::kable(table, "latex", booktabs = T , escape = F) %>% kable_styling()

#Find change in AIC for reducued model
AICs = c(-3378.4,
         -3372.2,
         -3369.8,
         -3368.6,
         -3365.4,
         -3364.0,
         -3363.6)


AICs - -3378.54

# Model validation --------------------------------------------------------


resids <- bkwrd_select$residuals %>% data.frame
qqfig <- ggplot(resids, aes(sample = .))
qqfig <- qqfig + stat_qq() + stat_qq_line(col = "red") +
  xlab("Theoretical Quantiles") + ylab("Residuals Quantiles")  + ggtitle("Q-Q Plot")

resid_fit <- ggplot() + geom_point(aes(y = bkwrd_select$residuals, x = bkwrd_select$fitted.values)) +
  ylab("Residuals") + xlab("Fitted Values") + ggtitle("Residuals and Fitted Values")

gridExtra::grid.arrange(
  qqfig, resid_fit, layout_matrix = rbind(c(1,2))
)


# Backwards elimination Based on BIC. --------------------------------------------------

bkwrd_select <- step(md_fit, direction = "backward",k=log(nrow(train_set)))#BIC

mantissas=apply( summary(bkwrd_select)$coefficients, c(1,2), mantissa)%>%round(4)

exponents=apply( summary(bkwrd_select)$coefficients, c(1,2), exponent)


table <- summary(bkwrd_select)$coefficients%>%round(4)

index=which(abs(table)<=0.001)

table[index]=paste0("$",mantissas[index],"\\cdot10^{",exponents[index],"}","$")


rownames(table)[2:nrow(table)] <- paste0("$", rownames(table)[2:nrow(table)], "$")
kableExtra::kable(table, "latex", booktabs = T , escape = F) %>% kable_styling()

# Model prediction --------------------------------------------------------

QV_preds <- lapply(1:nrow(test_set), function(idx){
  
  if(idx == 0){
    browser()
    single_pred <- predict(bkwrd_select, test_set[idx,], interval = "predict")
    return(data.frame(Date = test_set[idx,]$Date, QV_pred = single_pred[1,1],
                      upper_pred = single_pred[1,3], lower_pred = single_pred[1,2]))
  }else{
    new_train <- rbind(train_set, test_set[1:(idx-1),])
    new_md <- lm(formula = QV ~ VIX_1 + return_5 +
                   IV_22 + I(IV_22^2) + JV_22 + return_22 + VIX_22, new_train)
    single_pred <- predict(new_md, test_set[idx,], interval = "predict")
    data.frame(Date = test_set[idx,]$Date, QV_pred = single_pred[1,1],
               upper_pred = single_pred[1,3], lower_pred = single_pred[1,2])
  }
}) %>% do.call(what = rbind)

plot_df <- merge(QV_preds, test_set, by = "Date")

ggplot(plot_df) + geom_line(aes(x = Date, y = QV_pred, col = "QV Prediction")) +
  geom_line(aes(x = Date, y = QV, col = "QV")) + 
  geom_ribbon(aes(x = Date, ymin = lower_pred, ymax = upper_pred, fill = "Prediction Interval"),
              alpha = 0.4) +
  scale_fill_manual(labels = c("Prediction Interval"),
                    breaks = "Prediction Interval",
                    values = "darkgrey",
                    name = "") + 
  scale_color_manual(labels = c("QV Prediction", "QV"),
                     breaks = c("QV Prediction", "QV"),
                     values = c("Black", "Red"),
                     name = "") + ylab("")

plot_df %>% summarise(
  mape = mean(abs((QV - QV_pred)/ QV)), 
  coverage = mean(QV >= lower_pred & QV <= upper_pred)
) %>% round(4)




