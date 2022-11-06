#' @name linreg
#' @param formula formula
#' @param data data
#' @author Mahnaz , Bita
#' @description Implementing Linear Regression
#' @title linreg
#' @export linreg
#' @import ggplot2
#' @import cowplot
#' @references \url{http://staff.www.ltu.se/~jove/courses/c0002m/least_squares.pdf}


linreg <- function(formula,data){
  col_names <- c(names(data))
  
  for(i in 1:length(all.vars(formula))){
    stopifnot(all.vars(formula)[i] %in% col_names)
  }
  
  x <- model.matrix(formula, data)
  m <- data[all.vars(formula)[1]]
  y <- as.matrix(m)
  Bhat <- as.vector(solve((t(x) %*% x))%*% (t(x) %*% y))
  yhat <- as.vector(c( x %*% Bhat))
  ehat <- as.vector(y-yhat)
  degrees_of_freedom <- nrow(x)-ncol(x)
  residual_variance <- c((t(ehat)%*%ehat ) / degrees_of_freedom)
  variance_regression_coefficients <- residual_variance * solve((t(x) %*% x))
  t_values <- Bhat / sqrt(diag(variance_regression_coefficients))
  betahat <- t_values*sqrt(diag(variance_regression_coefficients))
  
  
  #QR
  QR <- qr(x)
  QR_pivot <- QR$pivot[1:QR$rank]
  betahat_QR <- solve.qr(QR, y)
  yhat_QR <- as.vector(c(x[,QR_pivot] %*% betahat))
  resi_QR <- as.vector(y - yhat_QR)
  degrees_of_freedom_QR <- nrow(x) - QR$rank
  
  residual_variance_QR <- c((t(resi_QR)%*%resi_QR ) / degrees_of_freedom)
  variance_regression_coefficients_QR = sqrt(residual_variance_QR)
  
  return_obj <- list(coefficients = betahat, fitted_values = yhat,
                     residuals = ehat, degree_freedom = degrees_of_freedom,
                     rv = residual_variance, varreg_coef = variance_regression_coefficients, 
                     t_value = t_values
                     ,call=match.call(),coefficients_QR=betahat_QR,
                     varr_QR=variance_regression_coefficients_QR,data=data,formula=formula)
  
  
  
  class(return_obj) <- 'linreg'
  return(return_obj)
  
  
}

#' @name print
#' @param obj data
#' @author Mahnaz , Bita
#' @description Implementing print
#' @title print
#' @export
# print 
print <- function (obj) {
  UseMethod("print")
}
print.linreg <- function(obj){
  cat("Call:\n")
  print(obj$call)
  cat("\nCoefficients:\n")
  print(obj$coefficients)
  
}

#' @name resid
#' @param obj data
#' @author Mahnaz , Bita
#' @description Implementing resid
#' @title resid
#' @export
# resid
resid <- function (obj) {
  UseMethod("resid")
}
resid.linreg <- function(obj){
  return(obj$residuals)
  
}


#' @name pred
#' @param obj data
#' @author Mahnaz , Bita
#' @description Implementing pred
#' @title pred
#' @export

# pred 
pred <- function (obj) {
  UseMethod("pred")
}
pred.linreg <- function(obj){
  return(obj$fitted_values)
  
}


#' @name coef
#' @param obj data
#' @author Mahnaz , Bita
#' @description Implementing coef
#' @title coef
#' @export
# coef 
coef <- function (obj) {
  UseMethod("coef")
}
coef.linreg <- function(obj){
  return(obj$coefficients)
}



#' @name summary
#' @param obj data
#' @author Mahnaz , Bita
#' @description Implementing summary
#' @title summary
#' @export
# summary 
summary <- function (obj) {
  UseMethod("summary")
}
summary.linreg <- function(obj)
{
  p_val <- function(p_value) {
    stars <- c()
    for(i in 1:length(p_value)){
      if (p_value[i] > 0.1) stars[i] = (" ")
      if (p_value[i] > 0.05) stars[i] = (".")
      if (p_value[i] > 0.01) stars[i] = ("*")
      if (p_value[i] > 0.001) stars[i] = ("**")
      else stars[i] = ("***")
    }
    return(stars)
  }
  Coefficient_df <- data.frame(Estimate = obj$coefficients,
                               StdErr = sqrt(diag(obj$varreg_coef)),
                               t_values = obj$t_value,
                               p_value = 2*pt(-abs(obj$t_value),df=obj$degree_freedom),
                               stars = p_val(2*pt(-abs(obj$t_value),df=obj$degree_freedom))
  )
  colnames(Coefficient_df) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)','')
  cat("Call:\n")
  print(obj$call)
  cat("\nCoefficients: \n")
  print(Coefficient_df)
  cat("\nResidual standard error:" , sqrt(obj$rv) ,"on", obj$degree_freedom ,"degrees of freedom")
}

#' @name plot
#' @param obj data
#' @author Mahnaz , Bita
#' @description Implementing plot
#' @title plot
#' @export
plot <- function (obj) {
  UseMethod("plot")
}
plot.linreg <- function(object) {
  
  std <- sqrt(abs(object$residuals))
  
  Residuals_Fitted <- ggplot2::ggplot(data=object$data,aes(x = object$fitted_values, y = object$residuals,label=Species))+ 
    ggplot2::geom_point(shape=1,size=3) + ggplot2::geom_smooth(method="lm", colour="red",se = FALSE)+ 
    ggplot2::labs(title="Residuals vs Fitted",
                  x="Fitted valueslinreg(Petal.Length~Species, data = iris)", 
                  y="Residuals")+ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                                panel.background = element_rect(fill = "white", colour = "black"))+
    ggplot2::geom_text(aes(label=ifelse(object$residuals>1&object$residuals<(-1),
                                        as.character(Species),'')),hjust=0,vjust=0)
  
  Scale_Location <- ggplot2::ggplot(data=object$data,aes(x = object$fitted_values, 
                                                         y = std,label=Species))+ 
    ggplot2::geom_point(shape=1,size=3) + ggplot2::geom_smooth(method="lm", colour="red",
                                                               se = FALSE)+ 
    ggplot2::labs(title="ScaleLocation",
                  x="Fitted valueslinreg(Petal.Length~Species, data = iris)", 
                  y=expression(sqrt(abs("Standard resduals"))))+
    ggplot2::theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_rect(fill = "white", 
                                                   colour = "black"))
  
  
  plot_grid(Residuals_Fitted,Scale_Location,ncol = 1, nrow = 2)
  
}
