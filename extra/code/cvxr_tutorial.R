## ---- echo = FALSE, fig = TRUE, fig.align = "center", fig.cap="Source: https://www.solver.com"----
knitr::include_graphics("figures/02/convexchord.gif")


## ---- eval = FALSE-------------------------------------------------------
## required_packages  <- c("tidyr",
##                         "ggplot2",
##                         "nnls",
##                         "ISLR",
##                         "glmnet",
##                         "isotone",
##                         "profvis",
##                         "dplyr",
##                         "survey",
##                         "expm",
##                         "RColorBrewer",
##                         "kableExtra")
## install.packages(required_packages)


## ---- eval = FALSE-------------------------------------------------------
## optimx(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,
##             method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
##             control=list(),
##              ...)


## ---- echo = FALSE-------------------------------------------------------
library(nnls)
library(kableExtra)

#' Print a matrix in a stylized way using row and column names if specified
#' @param the matrix to be printed
#' @param row_names optional row names to use can be math
#' @param col_names optional col names to use can be math
print_matrix <- function(m, row_names = NULL, col_names = NULL) {
  if (!is.null(row_names)) rownames(m) <- row_names
  if (!is.null(col_names)) colnames(m) <- col_names  
  knitr::kable(m, format = "html") %>%
    kable_styling("striped") %>%
    column_spec(1:2, background = "#ececec")
}


## ------------------------------------------------------------------------
set.seed(123)
n <- 50; p <- 10;
beta <- -4:5    # beta is just -4 through 5.
X <- matrix(rnorm(n * p), nrow=n)
colnames(X) <- paste0("beta_", beta)
Y <- X %*% beta + rnorm(n)


## ------------------------------------------------------------------------
ls.model <- lm(Y ~ 0 + X)   # There is no intercept in our model above
m <- matrix(coef(ls.model), ncol = 1)


## ---- echo = FALSE-------------------------------------------------------
print_matrix(m, row_names = paste0("$\\beta_{", 1:p, "}$"))


## ---- message = FALSE----------------------------------------------------
library(CVXR, warn.conflicts=FALSE)


## ------------------------------------------------------------------------
beta <- Variable(p)


## ------------------------------------------------------------------------
objective <- Minimize(sum((Y - X %*% beta)^2))


## ------------------------------------------------------------------------
problem <- Problem(objective)


## ------------------------------------------------------------------------
result <- solve(problem)


## ---- echo = FALSE-------------------------------------------------------
solution_status <- result$status
objective_value <- result$value
solution <- result$getValue(beta)
cat(sprintf("OLS Solution Status: %s, OLS Objective value: %f\n", solution_status, objective_value))


## ------------------------------------------------------------------------
m <- cbind(result$getValue(beta), coef(ls.model))


## ---- echo = FALSE-------------------------------------------------------
print_matrix(m, row_names = paste0("$\\beta_{", 1:p, "}$"), col_names = c("CVXR est.", "lm est."))


## ------------------------------------------------------------------------
objective2 <- Minimize(sum(abs(Y - X %*% beta)))
problem2 <- Problem(objective2)
result2 <- solve(problem2)
cat(sprintf("LAD Solution Status: %s, LAD Objective value: %f\n", result2$status, result2$value))
m2 <- cbind(result2$getValue(beta), coef(ls.model))


## ---- echo = FALSE-------------------------------------------------------
print_matrix(m2, row_names = paste0("$\\beta_{", 1:p, "}$"), col_names = c("CVXR LAD est.", "lm est."))


## ------------------------------------------------------------------------
cat("LAD objective value: %f, OLS objective value: %f\n",
    result2$value, result$value)


## ------------------------------------------------------------------------
m3 <- cbind(result$getValue(beta), result2$getValue(beta))


## ---- echo = FALSE-------------------------------------------------------
print_matrix(m3, row_names = paste0("$\\beta_{", 1:p, "}$"), col_names = c("Problem 1 est.", "Problem 2 est."))


## ------------------------------------------------------------------------
problem <- Problem(objective, constraints = list(beta >= 0))
result <- solve(problem)
betaEstimate <- result$getValue(beta)


## ---- echo = FALSE-------------------------------------------------------
m <- matrix(betaEstimate, ncol = 1)
print_matrix(m, row_names = paste0("$\\beta_{", 1:p, "}$"))


## ------------------------------------------------------------------------
nnls.fit <- nnls::nnls(X, Y)$x
m <- cbind(betaEstimate, nnls.fit)


## ---- echo = FALSE-------------------------------------------------------
print_matrix(m, row_names = paste0("$\\beta_{", 1:p, "}$"), col_names = c("CVXR NNLS est.", "nnls est."))


## ---- eval = FALSE-------------------------------------------------------
## constraint1 <- beta[1] + beta[2] + beta[3] + beta[4] <= 0


## ------------------------------------------------------------------------
A <- matrix(c(rep(1, 4), rep(0, 6)), nrow = 1)


## ---- echo = FALSE-------------------------------------------------------
print_matrix(A, col_names =  paste0("$\\beta_{", 1:p, "}$"))


## ------------------------------------------------------------------------
constraint1 <- A %*% beta <= 0


## ------------------------------------------------------------------------
problem <- Problem(objective, constraints = list(constraint1))
ex1 <- solve(problem)


## ------------------------------------------------------------------------
betaEstimate <- ex1$getValue(beta)


## ---- echo = FALSE-------------------------------------------------------
m <- matrix(betaEstimate, ncol = 1)
print_matrix(m, row_names = paste0("$\\beta_{", 1:p, "}$"))


## ------------------------------------------------------------------------
B <- diag(c(rep(0, 4), rep(1, 6)))


## ---- echo = FALSE-------------------------------------------------------
print_matrix(B, row_names = paste0("$\\beta_{", 1:p, "}$"), col_names = paste0("$\\beta_{", 1:p, "}$"))


## ------------------------------------------------------------------------
constraint2 <- B %*% beta <= 4
problem2 <- Problem(objective, constraints = list(constraint1, constraint2))
ex2 <- solve(problem2)
betaEstimate <- ex2$getValue(beta)


## ---- echo = FALSE-------------------------------------------------------
m <- matrix(betaEstimate, ncol = 1)
print_matrix(m, row_names =  paste0("$\\beta_{", 1:p, "}$"))


## ------------------------------------------------------------------------
problem3 <- Problem(objective,
                   constraints = list(beta >= 0, diff(beta) >= 0))
ex3 <- solve(problem3)
betaEstimate <- ex3$getValue(beta)


## ---- echo = FALSE-------------------------------------------------------
m <- matrix(betaEstimate, ncol = 1)
print_matrix(m, row_names =  paste0("$\\beta_{", 1:p, "}$"))


## ------------------------------------------------------------------------
D1 <- D2 <- diag(0, p)
D1[1:5, 1:5] <- 1; D2[5:p, 5:p] <- 1;
D3 <- D4 <- diag(0, p - 1)
D3[1:4, 1:4] <- 1; D4[5:(p-1), 5:(p-1)] <- 1;
constraints = list(D3 %*% diff(D1 %*% beta) >= 0, D4 %*% diff(D2 %*% beta) <= 0)
problem4 <- Problem(objective, constraints)
ex4 <- solve(problem4)
betaEstimate <- ex4$getValue(beta)


## ---- echo = FALSE-------------------------------------------------------
m <- matrix(betaEstimate, ncol = 1)
print_matrix(m, row_names =  paste0("$\\beta_{", 1:p, "}$"))


## ---- eval = FALSE-------------------------------------------------------
## beta <- Variable(p)
## objective <- Minimize(sum((Y - X %*% beta)^2))
## constraints <- list(beta >= 0)
## problem <- Problem(objective, constraints)
## result <- solve(problem)
## solution_status <- result$status
## objective_value <- result$value
## beta_hat <- result$getValue(beta)


## ---- eval = FALSE-------------------------------------------------------
## objective <- Minimize(sum((Y - X %*% beta)^2))


## ---- eval = FALSE-------------------------------------------------------
## objective <- Minimize(sum_squares(Y - X %*% beta))


## ------------------------------------------------------------------------
result <- solve(problem, verbose = TRUE)


## ------------------------------------------------------------------------
library(ISLR)
data(Smarket)
glmfit <- stats::glm(formula = Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume,
                     data = Smarket, family = binomial)	       


## ------------------------------------------------------------------------
y <- as.integer(Smarket$Direction) - 1L
X <- cbind(1, as.matrix(Smarket[, c("Lag1", "Lag2", "Lag3", "Lag4", "Lag5", "Volume")]))
p <- ncol(X)
beta <- Variable(p)
objective <- -sum(X[y <= 0, ] %*% beta) - sum(logistic(-X %*% beta))
problem <- Problem(Maximize(objective))
result <- solve(problem)
beta_hat <- result$getValue(beta)


## ---- echo = FALSE-------------------------------------------------------
print_matrix(cbind(beta_hat, coef(glmfit)), row_names = paste0("$\\beta_{", seq.int(0, p-1L), "}$"), col_names = c("CVXR", "GLM"))


## ------------------------------------------------------------------------
fitted_values <- 1 / (1 + exp(result$getValue(-X %*% beta)))


## ------------------------------------------------------------------------
sum((fitted_values - glmfit$fitted.values))^2


## ------------------------------------------------------------------------
log_odds <- result$getValue(X %*% beta)


## ---- echo = FALSE, message = FALSE--------------------------------------
library(isotone)


## ------------------------------------------------------------------------
data("pituitary")
str(pituitary)


## ------------------------------------------------------------------------
res_p <- with(pituitary, gpava(age, size))


## ------------------------------------------------------------------------
x_p <- with(pituitary, {
      n <- length(size)
      x <- Variable(n)
      objective <- Minimize(cvxr_norm(size - x, 2))
      constraint <- list(diff(x) >= 0)
      problem <- Problem(objective, constraint)
      result <- solve(problem)
      result$getValue(x)
})


## ------------------------------------------------------------------------
print_matrix(cbind(res_p$x, x_p), col_names = c("isotone", "CVXR"))


## ------------------------------------------------------------------------
res_s <- with(pituitary, gpava(age, size, ties = "secondary"))
res_t <- with(pituitary, gpava(age, size, ties = "tertiary"))


## ------------------------------------------------------------------------
x_s <- with(pituitary, {
    n <- length(size)
    x <- Variable(n)
    objective <- Minimize(p_norm(size - x, 2))
    secondary_constraints <- lapply(base::split(x = seq_len(n),
                                                f = age),
                                    function(i) diff(x[i]) == 0)
    constraint <- c(diff(x) >= 0,
                    secondary_constraints)
    problem <- Problem(objective, constraint)
    solve(problem)$getValue(x)
})


## ---- echo = FALSE-------------------------------------------------------
m <- cbind(res_s$x, x_s)
print_matrix(m, col_names = c("Isotone (S)", "CVXR (S)"))


## ------------------------------------------------------------------------
x_t <- with(pituitary, {
    n <- length(size)
    x <- Variable(n)
    objective <- Minimize(p_norm(size - x, 2))
    blocks <- base::split(x = seq_len(n),
                          f = pituitary$age)
    block_means <- lapply(blocks, function(i) {
        v <- numeric(n)
        v[i] <- 1.0 / length(i)
        matrix(v, nrow = 1) %*% x
    })
    block_mean_vector <- do.call(vstack, block_means)
    constraint <- list(diff(block_mean_vector) >= 0)
    problem <- Problem(objective, constraint)
    solve(problem)$getValue(x)
})


## ---- echo = FALSE-------------------------------------------------------
m <- cbind(res_t$x, x_t)
print_matrix(m, col_names = c("Isotone (T)", "CVXR (T)"))


## ---- echo = FALSE, message = FALSE--------------------------------------
library(glmnet)


## ------------------------------------------------------------------------
#' Define the elastic penalty
#' @param beta the arg min variable
#' @param lambda the penalization parameter
#' @param alpha the elastic net parameter, 0 = ridge, 1 = lasso
elastic_penalty <- function(beta, lambda = 0, alpha = 0) {
    ridge <- (1 - alpha) / 2 * sum_squares(beta)
    lasso <- alpha * cvxr_norm(beta, 1)
    lambda * (lasso + ridge)
}


## ------------------------------------------------------------------------
## Problem data
set.seed(4321)
p <- 10
n <- 500
DENSITY <- 0.25    # Fraction of non-zero beta
beta_true <- matrix(rnorm(p), ncol = 1)
idxs <- sample.int(p, size = floor((1 - DENSITY) * p), replace = FALSE)
beta_true[idxs] <- 0
sigma <- 45
X <- matrix(rnorm(n * p, sd = 5), nrow = n, ncol = p)
eps <- matrix(rnorm(n, sd = sigma), ncol = 1)
Y <- X %*% beta_true + eps


## ------------------------------------------------------------------------
TRIALS <- 10
beta_vals <- matrix(0, nrow = p, ncol = TRIALS)
lambda_vals <- 10^seq(-2, log10(50), length.out = TRIALS)


## ------------------------------------------------------------------------
beta <- Variable(p)  
loss <- sum_squares(Y - X %*% beta) / (2 * n)
## Elastic-net regression LASSO
alpha <- 1
beta_vals <- sapply(lambda_vals,
                    function (lambda) {
                        obj <- loss + elastic_penalty(beta, lambda, alpha)
                        prob <- Problem(Minimize(obj))
                        result <- solve(prob)
                        result$getValue(beta)
                    })


## ---- echo = FALSE-------------------------------------------------------
print_matrix(round(beta_vals, 3),
             row_names = sprintf("$\\beta_{%d}$", seq_len(p)),
             col_names = sprintf("$\\lambda = %.3f$", lambda_vals))


## ------------------------------------------------------------------------
plot(0, 0, type = "n", main = "CVXR Regularization Path for Lasso Regression",
     xlab = "Log Lambda", ylab = "Coefficients",
     ylim = c(-1, 2), xlim = c(-4, 4))
matlines(log(lambda_vals), t(beta_vals))


## ------------------------------------------------------------------------
model_net <- glmnet(X, Y, family = "gaussian", alpha = alpha,
                    lambda = lambda_vals,
                    standardize = FALSE,
                    intercept = FALSE,
                    thresh = 1e-8)
## Reverse order to match beta_vals
coef_net <- as.data.frame(as.matrix(coef(model_net)[-1, seq(TRIALS, 1, by = -1)]))


## ---- echo = FALSE-------------------------------------------------------
print_matrix(round(coef_net, 3),
             row_names = sprintf("$\\beta_{%d}$", seq_len(p)),
             col_names = sprintf("$\\lambda = %.3f$", lambda_vals))


## ------------------------------------------------------------------------
library(glmnet)
data(QuickStartExample)
catn <- function(...) cat(..., "\n")
objective_value <- function(y, x, coefs, lambda, alpha) {
    n <- nrow(x)
    ridge <- sum(coefs^2) ; l1 <- sum(abs(coefs))
    sum((y - (x %*% coefs))^2) / (2 * n) + lambda * ((1 - alpha) / 2 * ridge + alpha * l1)
}
alpha  <- 0; lambda  <- 1;
fit <- glmnet(x, y, intercept=F, standardize=F, lambda=1, alpha=0)


## ------------------------------------------------------------------------
objective_value(y, x, coef(fit)[-1, ], lambda, alpha)


## ------------------------------------------------------------------------
beta <- Variable(20)
elastic_reg <- function(beta, lambda = 0, alpha = 0) {
   ridge <- (1 - alpha) * sum(beta^2) * .5
   lasso <- alpha * p_norm(beta, 1)
   lambda * (lasso + ridge)
}
loss <- sum((y - x %*% beta)^2)/(2*length(y))
obj <- loss + elastic_reg(beta, lambda = 1, 0)
prob <- Problem(Minimize(obj))
result <- solve(prob)
print(result$value)


## ------------------------------------------------------------------------
catn <- function(...) cat(..., "\n")
## Standardize the y
y_s <- local({
    n <- length(y)
    m <- mean(y); s <- as.numeric(sqrt(var(y) * (n - 1) / n));
    result <- (y - m) / s ## scale using 1/n
    attr(result, "scaled:center") <- m
    attr(result, "scaled:scale") <- s
    result
})


## ------------------------------------------------------------------------
## STANDARDIZED COMPARISON
fit_s <- glmnet(x, y_s, intercept=F, standardize=F, lambda = lambda, alpha=alpha)
catn("Glmnet objective (scaled y)",
     objective_value(y_s, x, coef(fit_s)[-1], lambda, alpha))


## ------------------------------------------------------------------------
elastic_reg <- function(beta, lambda = 0, alpha = 0) {
    ridge <- (1 - alpha) / 2 * sum_squares(beta)
    lasso <- alpha * p_norm(beta, 1)
    lambda * (lasso + ridge)
}
loss <- sum_squares(y_s - x %*% beta) / (2 * nrow(x))
obj <- loss + elastic_reg(beta, lambda = lambda, alpha)
prob <- Problem(Minimize(obj))
beta_est <- solve(prob)$getValue(beta)
catn("CVXR objective (scaled y):", objective_value(y_s, x, beta_est, lambda, alpha))


## ------------------------------------------------------------------------
## NONSTANDARDIZED COMPARISON
fit <- glmnet(x, y, intercept=F, standardize=F, lambda = lambda, alpha=alpha)
catn("Glmnet objective (unscaled y)", objective_value(y, x, coef(fit)[-1], lambda, alpha))


## ------------------------------------------------------------------------
loss <- sum_squares(y - x %*% beta) / (2 * nrow(x))
obj <- loss + elastic_reg(beta, lambda = lambda / attr(y_s, "scaled:scale"), alpha)
prob <- Problem(Minimize(obj))
beta_est <- solve(prob)$getValue(beta)
catn("CVXR objective (unscaled y)", objective_value(y, x, beta_est, lambda, alpha))


## ------------------------------------------------------------------------
print_matrix(round(cbind(beta_est, coef(fit)[-1]), 3), 
             row_names = sprintf("$\\beta_{%d}$", seq_len(20)),
             col_names = c("CVXR", "GLMNET")) 


## ------------------------------------------------------------------------
beta <- Variable(p)  
loss  <- sum(huber(Y - X %*% beta, M = 0.5))
## Elastic-net regression LASSO
alpha <- 1
beta_vals <- sapply(lambda_vals,
                    function (lambda) {
                        obj <- loss + elastic_penalty(beta, lambda, alpha)
                        prob <- Problem(Minimize(obj))
                        result <- solve(prob)
                        result$getValue(beta)
                    })


## ---- echo = FALSE-------------------------------------------------------
print_matrix(round(beta_vals, 3),
             row_names = sprintf("$\\beta_{%d}$", seq_len(p)),
             col_names = sprintf("$\\lambda = %.3f$", lambda_vals))


## ---- message = FALSE, echo = FALSE--------------------------------------
library(ggplot2)
library(boot)


## ------------------------------------------------------------------------
data(cdiac)
str(cdiac)


## ------------------------------------------------------------------------
neariso_fit <- function(y, lambda) {
    m <- length(y)
    beta <- Variable(m)
    obj <- 0.5 * sum_squares(y - beta) + lambda * sum(pos(diff(beta)))
    prob <- Problem(Minimize(obj))
    solve(prob)$getValue(beta)
}


## ------------------------------------------------------------------------
neariso_fit_stat <- function(data, index, lambda) {
    sample <- data[index,]                  # Bootstrap sample of rows
    sample <- sample[order(sample$year),]   # Order ascending by year
    neariso_fit(sample$annual, lambda)
}


## ---- eval = FALSE-------------------------------------------------------
## set.seed(123)
## boot.neariso <- boot(data = cdiac,
##                      statistic = neariso_fit_stat,
##                      R = 10, lambda = 0.44)
## ci.neariso <- t(sapply(seq_len(nrow(cdiac)),
##                        function(i) boot.ci(boot.out = boot.neariso, conf = 0.95,
##                                            type = "norm", index = i)$normal[-1]))
## data.neariso <- data.frame(year = cdiac$year,
##                            annual = cdiac$annual,
##                            est = boot.neariso$t0,
##                            lower = ci.neariso[, 1],
##                            upper = ci.neariso[, 2])


## ---- eval = FALSE-------------------------------------------------------
## (plot.neariso <- ggplot(data = data.neariso) +
##      geom_point(mapping = aes(year, annual), color = "red") +
##      geom_line(mapping = aes(year, est), color = "blue") +
##      geom_ribbon(mapping = aes(x = year, ymin = lower,ymax = upper),alpha=0.3) +
##      labs(x = "Year", y = "Temperature Anomalies")
## )


## ---- eval = FALSE-------------------------------------------------------
## nearconvex_fit <- function(y, lambda) {
##     m <- length(y)
##     beta <- Variable(m)
##     obj <- 0.5 * sum_squares(y - beta) + lambda * sum(pos(diff(beta, differences = 2)))
##     prob <- Problem(Minimize(obj))
##     solve(prob)$getValue(beta)
## }
## 
## nearconvex_fit_stat <- function(data, index, lambda) {
##     sample <- data[index,]                  # Bootstrap sample of rows
##     sample <- sample[order(sample$year),]   # Order ascending by year
##     nearconvex_fit(sample$annual, lambda)
## }
## 
## set.seed(987)
## boot.nearconvex <- boot(data = cdiac,
##                         statistic = nearconvex_fit_stat,
##                         R = 5,
##                         lambda = 0.44)
## 
## ci.nearconvex <- t(sapply(seq_len(nrow(cdiac)),
##                           function(i) boot.ci(boot.out = boot.nearconvex, conf = 0.95,
##                                               type = "norm", index = i)$normal[-1]))
## data.nearconvex <- data.frame(year = cdiac$year,
##                               annual = cdiac$annual,
##                               est = boot.nearconvex$t0,
##                               lower = ci.nearconvex[, 1],
##                               upper = ci.nearconvex[, 2])
## 


## ---- eval = FALSE-------------------------------------------------------
## (plot.nearconvex <- ggplot(data = data.nearconvex) +
##      geom_point(mapping = aes(year, annual), color = "red") +
##      geom_line(mapping = aes(year, est), color = "blue") +
##      geom_ribbon(mapping = aes(x = year, ymin = lower,ymax = upper),alpha=0.3) +
##      labs(x = "Year", y = "Temperature Anomalies")
## )


## ----prereqs, message = FALSE, echo = FALSE------------------------------
library(profvis)


## ------------------------------------------------------------------------
data(cdiac)
y <- cdiac$annual
m <- length(y)
lambda <- 0.44
beta <- Variable(m)
obj <- 0.5 * sum((y - beta)^2) + lambda * sum(pos(diff(beta)))
prob <- Problem(Minimize(obj))
soln <- solve(prob)
betaHat <- soln$getValue(beta)


## ------------------------------------------------------------------------
data(cdiac)
y <- cdiac$annual
profvis({
    beta <- Variable(m)
    obj <- Minimize(0.5 * sum((y - beta)^2) + lambda * sum(pos(diff(beta))))
    prob <- Problem(obj)
    soln <- solve(prob)
    betaHat <- soln$getValue(beta)
})


## ------------------------------------------------------------------------
prob_data <- get_problem_data(prob, solver = "ECOS")


## ---- eval = FALSE-------------------------------------------------------
## soln <- solve(prob, verbose = TRUE)


## ------------------------------------------------------------------------
solver_output <- ECOSolveR::ECOS_csolve(c = prob_data[["c"]],
                                        G = prob_data[["G"]],
                                        h = prob_data[["h"]],
                                        dims = prob_data[["dims"]],
                                        A = prob_data[["A"]],
                                        b = prob_data[["b"]])


## ------------------------------------------------------------------------
direct_soln <- unpack_results(prob, "ECOS", solver_output)


## ------------------------------------------------------------------------
profvis({
    beta <- Variable(m)
    obj <- Minimize(0.5 * sum((y - beta)^2) + lambda * sum(pos(diff(beta))))
    prob <- Problem(obj)
    prob_data <- get_problem_data(prob, solver = "ECOS")
    solver_output <- ECOSolveR::ECOS_csolve(c = prob_data[["c"]],
                                            G = prob_data[["G"]],
                                            h = prob_data[["h"]],
                                            dims = prob_data[["dims"]],
                                            A = prob_data[["A"]],
                                            b = prob_data[["b"]])
    direct_soln <- unpack_results(prob, "ECOS", solver_output)
})


## ------------------------------------------------------------------------
identical(betaHat, direct_soln$getValue(beta))


## ------------------------------------------------------------------------
## Simulation data.
set.seed(123)
N <- 100
K <- 4
p <- 50
X <- matrix(rnorm(n = N * p, mean = 0, sd = 1), nrow = N, ncol = p)
Z <- matrix(rbinom(n = N * K, size = 1, prob = 0.5), nrow = N, ncol = K)

## Response model.
beta <- rep(x = 0, times = p)
beta[1:4] <- c(2, -2, 2, 2)
coeffs <- cbind(beta[1], beta[2], beta[3] + 2 * Z[, 1], beta[4] * (1 - 2 * Z[, 2]))
mu <- diag(X[, 1:4] %*% t(coeffs))
y <- mu + 0.5 * rnorm(N, mean = 0, sd = 1)


## ---- eval = FALSE-------------------------------------------------------
## plasso_fit <- function(y, X, Z, lambda, alpha = 0.5, intercept = TRUE,
##                        ZERO_THRESHOLD= 1e-6, verbose = FALSE) {
##     N <- length(y)
##     p <- ncol(X)
##     K <- ncol(Z)
## 
##     beta0 <- 0
##     if (intercept) {
##         beta0 <- Variable(1) * matrix(1, nrow = N, ncol = 1)
##     }
##     ## Define_Parameters
##     ## Build_Penalty_Terms
##     ## Compute_Fitted_Value
##     ## Build_Objective
##     ## Define_and_Solve_Problem
##     ## Return_Values
## }
## 
## ## Fit pliable lasso using CVXR.
## # pliable <- pliable_lasso(y, X, Z, alpha = 0.5, lambda = lambda)


## ----Define_Parameters, eval = FALSE-------------------------------------
## beta <- Variable(p)
## theta0 <- Variable(K)
## theta <- Variable(p, K); theta_transpose <- t(theta)


## ----Build_Penalty_Terms, eval = FALSE-----------------------------------
## penalty_term1 <- sum(cvxr_norm(hstack(beta, theta), 2, axis = 1))
## penalty_term2 <- sum(cvxr_norm(theta, 2, axis = 1))
## penalty_term3 <- sum(cvxr_norm(theta, 1))


## ----Compute_Fitted_Value, eval = FALSE----------------------------------
## xz_theta <- lapply(seq_len(p),
##                    function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% theta_transpose[, j])
## XZ_term <- Reduce(f = '+', x = xz_theta)
## y_hat <- beta0 + X %*% beta + Z %*% theta0 + XZ_term


## ----Build_Objective, eval = FALSE---------------------------------------
## objective <- sum_squares(y - y_hat) / (2 * N) +
##     (1 - alpha) * lambda * (penalty_term1 + penalty_term2) +
##     alpha * lambda * penalty_term3


## ----Define_and_Solve_Problem, eval = FALSE------------------------------
## prob <- Problem(Minimize(objective))
## result <- solve(prob, verbose = verbose)
## beta_hat <- result$getValue(beta)


## ----Return_Results, eval = FALSE----------------------------------------
## theta0_hat <- result$getValue(theta0)
## theta_hat <- result$getValue(theta)
## 
## ## Zero out stuff before returning
## beta_hat[abs(beta_hat) < ZERO_THRESHOLD] <- 0.0
## theta0_hat[abs(theta0_hat) < ZERO_THRESHOLD] <- 0.0
## theta_hat[abs(theta_hat) < ZERO_THRESHOLD] <- 0.0
## list(beta0_hat = if (intercept) result$getValue(beta0)[1] else 0.0,
##      beta_hat = beta_hat,
##      theta0_hat = theta0_hat,
##      theta_hat = theta_hat,
##      criterion = result$value)


## ------------------------------------------------------------------------
plasso_fit <- function(y, X, Z, lambda, alpha = 0.5, intercept = TRUE,
                          ZERO_THRESHOLD= 1e-6, verbose = FALSE) {
    N <- length(y)
    p <- ncol(X)
    K <- ncol(Z)

    beta0 <- 0
    if (intercept) {
        beta0 <- Variable(1) * matrix(1, nrow = N, ncol = 1)
    }
    beta <- Variable(p)
    theta0 <- Variable(K)
    theta <- Variable(p, K); theta_transpose <- t(theta)
    penalty_term1 <- sum(cvxr_norm(hstack(beta, theta), 2, axis = 1))
    penalty_term2 <- sum(cvxr_norm(theta, 2, axis = 1))
    penalty_term3 <- sum(cvxr_norm(theta, 1))
    xz_theta <- lapply(seq_len(p),
                       function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% theta_transpose[, j])
    XZ_term <- Reduce(f = '+', x = xz_theta)
    y_hat <- beta0 + X %*% beta + Z %*% theta0 + XZ_term
    objective <- sum_squares(y - y_hat) / (2 * N) +
        (1 - alpha) * lambda * (penalty_term1 + penalty_term2) +
        alpha * lambda * penalty_term3
    prob <- Problem(Minimize(objective))
    result <- solve(prob, verbose = verbose)
    beta_hat <- result$getValue(beta)
    theta0_hat <- result$getValue(theta0)
    theta_hat <- result$getValue(theta)
    
    ## Zero out stuff before returning
    beta_hat[abs(beta_hat) < ZERO_THRESHOLD] <- 0.0
    theta0_hat[abs(theta0_hat) < ZERO_THRESHOLD] <- 0.0
    theta_hat[abs(theta_hat) < ZERO_THRESHOLD] <- 0.0
    list(beta0_hat = if (intercept) result$getValue(beta0)[1] else 0.0,
         beta_hat = beta_hat,
         theta0_hat = theta0_hat,
         theta_hat = theta_hat,
         criterion = result$value)
}


## ------------------------------------------------------------------------
result <- plasso_fit(y, X, Z, lambda = 0.6, alpha = 0.5, intercept = FALSE)


## ------------------------------------------------------------------------
cat(sprintf("Objective value: %f\n", result$criterion))


## ------------------------------------------------------------------------
index <- which(result$beta_hat != 0)
est.table <- data.frame(matrix(result$beta_hat[index], nrow = 1))
names(est.table) <- paste0("$\\beta_{", index, "}$")
knitr::kable(est.table, format = "html", digits = 3) %>%
    kable_styling("striped")


## ------------------------------------------------------------------------
est.table <- data.frame(matrix(result$theta0_hat, nrow = 1))
names(est.table) <- paste0("$\\theta_{0,", 1:K, "}$")
knitr::kable(est.table, format = "html", digits = 3) %>%
    kable_styling("striped")


## ------------------------------------------------------------------------
est.table <- data.frame(result$theta_hat[1:5, ])
names(est.table) <- paste0("$\\theta_{,", 1:K, "}$")
knitr::kable(est.table, format = "html", digits = 3) %>%
    kable_styling("striped")


## ---- message = FALSE, echo = FALSE--------------------------------------
library(dplyr)
library(survey)
build_df <- function(api, method, wt_value) {
    d <- data.frame(apisrs$stype, apisrs$sch.wide, wt_value )
    names(d) <- c("stype", "sch.wide", "weight")
    rownames(d) <- NULL
    d %<>%
        group_by(stype, sch.wide) %>%
        summarize(value = first(weight), frequency = n())
    names(d) <- c("stype", "sch.wide", paste(method, "wts."), "Frequency")
    d
}
build_table <- function(d1, d2, title) {
    d <- inner_join(d1, d2, by = c("stype", "sch.wide"))
    names(d) <- gsub("Frequency.x|Frequency.y", "Frequency", names(d))
    d %>%
    knitr::kable(format = "html", digits = 3, caption = title) %>%
    kable_styling("striped") %>%
    column_spec(1:6, background = "#ececec")
}


## ------------------------------------------------------------------------
data(api)
design_api <- svydesign(id = ~dnum, weights = ~pw, data = apisrs)
formula <- ~stype + sch.wide
T <- apply(model.matrix(object = formula, data = apipop),
           2,
           sum)

cal_api <- calibrate(design_api, formula, population = T, calfun = cal.raking)
w_survey <- weights(cal_api)


## ------------------------------------------------------------------------
di <- apisrs$pw
X <- model.matrix(object = formula, data = apisrs)
A <- di * X
n <- nrow(apisrs)
g <- Variable(n)
constraints <- list(t(A) %*% g == T)

## Raking
Phi_R <- Minimize(sum(di * (-entr(g) - g + 1)))
p <- Problem(Phi_R, constraints)
res <- solve(p)
w_cvxr <- di * res$getValue(g)


## ---- echo = FALSE-------------------------------------------------------
## Using functions in the *un echoed* preamble of this document...
build_table(d1 = build_df(apisrs, "Survey", w_survey),
            d2 = build_df(apisrs, "CVXR", w_cvxr),
            title = "Calibration weights from Raking")


## ------------------------------------------------------------------------
w_survey_q <- weights(calibrate(design_api, formula, population = T, calfun = cal.linear))


## ------------------------------------------------------------------------
## Quadratic
Phi_Q  <- Minimize(sum_squares(g - 1) / 2)
p <- Problem(Phi_Q, constraints)
res <- solve(p, solver = "SCS")
w_cvxr_q <- di * res$getValue(g)


## ---- echo = FALSE-------------------------------------------------------
build_table(d1 = build_df(apisrs, "Survey", w_survey_q),
            d2 = build_df(apisrs, "CVXR", w_cvxr_q),
            title = "Calibration weights from Quadratic metric")


## ------------------------------------------------------------------------
u <- 1.10; l <- 0.90
w_survey_l <- weights(calibrate(design_api, formula, population = T, calfun = cal.linear,
                                bounds = c(l, u)))


## ------------------------------------------------------------------------
Phi_L <- Minimize(sum(-entr((g - l) / (u - l))  -
                      entr((u - g) / (u - l)))) 
p <- Problem(Phi_L, c(constraints, list(l <= g, g <= u)))
res <- solve(p)
w_cvxr_l <- di * res$getValue(g)


## ---- echo = FALSE-------------------------------------------------------
build_table(d1 = build_df(apisrs, "Survey", w_survey_l),
            d2 = build_df(apisrs, "CVXR", w_cvxr_l),
            title = "Calibration weights from Logit metric")


## ------------------------------------------------------------------------
hellinger <- make.calfun(Fm1 = function(u, bounds)  ((1 - u / 2)^-2) - 1,
                         dF= function(u, bounds) (1 -u / 2)^-3 ,
                         name = "Hellinger distance")
w_survey_h <- weights(calibrate(design_api, formula, population = T, calfun = hellinger))


## ------------------------------------------------------------------------
Phi_h <- Minimize(sum((1 - g / 2)^(-2)))
p <- Problem(Phi_h, constraints)
res <- solve(p)
w_cvxr_h <- di * res$getValue(g)


## ---- echo = FALSE-------------------------------------------------------
build_table(d1 = build_df(apisrs, "Survey", w_survey_h),
            d2 = build_df(apisrs, "CVXR", w_cvxr_h),
            title = "Calibration weights from Hellinger distance metric")


## ------------------------------------------------------------------------
w_survey_s <- weights(calibrate(design_api, formula, population = T, calfun = cal.sinh,
                                bounds = c(l, u)))


## ------------------------------------------------------------------------
Phi_s <- Minimize(sum( 0.5 * (exp(g) + exp(-g))))
p <- Problem(Phi_s, c(constraints, list(l <= g, g <= u)))
res <- solve(p)
w_cvxr_s <- di * res$getValue(g)


## ---- echo = FALSE-------------------------------------------------------
build_table(d1 = build_df(apisrs, "Survey", w_survey_s),
            d2 = build_df(apisrs, "CVXR", w_cvxr_s),
            title = "Calibration weights from derivative of sinh metric")


## ---- message = FALSE, echo = FALSE--------------------------------------
library(CVXR)
library(ggplot2)
library(grid)
library(Matrix)
library(expm)
## 
## Reference: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

theme_bare <- theme(
    axis.line = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(), 
    legend.position = "none", 
    panel.background = element_rect(fill = "white"), 
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    )


## ------------------------------------------------------------------------
set.seed(1)
n <- 10      ## Dimension of matrix
m <- 1000    ## Number of samples

## Create sparse, symmetric PSD matrix S
A <- rsparsematrix(n, n, 0.15, rand.x = stats::rnorm)
Strue <- A %*% t(A) + 0.05 * diag(rep(1, n))    ## Force matrix to be strictly positive definite


## ------------------------------------------------------------------------
R <- base::solve(Strue)


## ------------------------------------------------------------------------
x_sample <- matrix(stats::rnorm(n * m), nrow = m, ncol = n) %*% t(expm::sqrtm(R))
Q <- cov(x_sample)    ## Sample covariance matrix


## ------------------------------------------------------------------------
alphas <- c(10, 8, 6, 4, 1)
S <- Semidef(n)    ## Variable constrained to positive semidefinite cone
obj <- Maximize(log_det(S) - matrix_trace(S %*% Q))

S.est <- lapply(alphas,
                function(alpha) {
                    constraints <- list(sum(abs(S)) <= alpha)
                    ## Form and solve optimization problem
                    prob <- Problem(obj, constraints)
                    result <- solve(prob)
                    
                    ## Create covariance matrix
                    R_hat <- base::solve(result$getValue(S))
                    Sres <- result$getValue(S)
                    Sres[abs(Sres) <= 1e-4] <- 0
                    Sres
                })


## ---- echo = FALSE, fig.width=6, fig.height=6----------------------------
## Plotting function
plotSpMat <- function(S, alpha) {
    n <- nrow(S)
    df <- expand.grid(j = seq_len(n), i = seq_len(n))
    df$z = as.character(as.numeric(S) != 0)
    p <- ggplot(data = df, mapping = aes(x = i, y = j, fill = z)) +
        geom_tile(color = "black") +
        scale_fill_brewer(type = "qual", palette = "Paired") +
        scale_y_reverse()
    if (missing(alpha)) {
        p <- p + xlab("Truth")
    } else {
        p <- p + xlab(parse(text=(paste0("alpha == ", alpha))))
    }
    p + theme_bare
}


## ---- fig.width=6, fig.height=4------------------------------------------
do.call(multiplot, args = c(list(plotSpMat(Strue)),
                            mapply(plotSpMat, S.est, alphas, SIMPLIFY = FALSE),
                            list(layout = matrix(1:6, nrow = 2, byrow = TRUE))))


## ---- message = FALSE, echo = FALSE--------------------------------------
library(CVXR)
library(ggplot2)
library(RColorBrewer)
library(tidyr)


## ------------------------------------------------------------------------
## Problem data
set.seed(10)
n <- 10
SAMPLES <- 100
mu <- matrix(abs(rnorm(n)), nrow = n)
Sigma <- matrix(rnorm(n^2), nrow = n, ncol = n)
Sigma <- t(Sigma) %*% Sigma

## Form problem
w <- Variable(n)
ret <- t(mu) %*% w
risk <- quad_form(w, Sigma)
constraints <- list(w >= 0, sum(w) == 1)

## Risk aversion parameters
gammas <- 10^seq(-2, 3, length.out = SAMPLES)
ret_data <- rep(0, SAMPLES)
risk_data <- rep(0, SAMPLES)
w_data <- matrix(0, nrow = SAMPLES, ncol = n)

## Compute trade-off curve
for(i in seq_along(gammas)) {
    gamma <- gammas[i]
    objective <- ret - gamma * risk
    prob <- Problem(Maximize(objective), constraints)
    result <- solve(prob)
    
    ## Evaluate risk/return for current solution
    risk_data[i] <- result$getValue(sqrt(risk))
    ret_data[i] <- result$getValue(ret)
    w_data[i,] <- result$getValue(w)
}


## ---- eval = FALSE-------------------------------------------------------
## result$getValue(risk)
## result$getValue(ret)


## ------------------------------------------------------------------------
cbPalette <- brewer.pal(n = 10, name = "Paired")
p1 <- ggplot() +
    geom_line(mapping = aes(x = risk_data, y = ret_data), color = "blue") +
    geom_point(mapping = aes(x = sqrt(diag(Sigma)), y = mu), color = "red")

markers_on <- c(10, 20, 30, 40)
nstr <- sprintf("gamma == %.2f", gammas[markers_on])
df <- data.frame(markers =  markers_on, x = risk_data[markers_on],
                 y = ret_data[markers_on], labels = nstr)

p1 + geom_point(data = df, mapping = aes(x = x, y = y), color = "black") +
    annotate("text", x = df$x + 0.2, y = df$y - 0.05, label = df$labels, parse = TRUE) +
    labs(x = "Risk (Standard Deviation)", y = "Return")


## ------------------------------------------------------------------------
w_df <- data.frame(paste0("grp", seq_len(ncol(w_data))),
                   t(w_data[markers_on,]))
names(w_df) <- c("grp", sprintf("gamma == %.2f", gammas[markers_on]))
tidyW <- gather(w_df, key = "gamma", value = "fraction", names(w_df)[-1], factor_key = TRUE)
ggplot(data = tidyW, mapping = aes(x = gamma, y = fraction)) +
    geom_bar(mapping = aes(fill = grp), stat = "identity") +
    scale_x_discrete(labels = parse(text = levels(tidyW$gamma))) +
    scale_fill_manual(values = cbPalette) +
    guides(fill = FALSE) +
    labs(x = "Risk Aversion", y = "Fraction of Budget")


## ---- eval = FALSE-------------------------------------------------------
## constr <- list(p_norm(w,1) <= Lmax, sum(w) == 1)


## ------------------------------------------------------------------------
lapply(list(LPSOLVE = "lpSolveAPI",
            GLPK = "Rglpk"),
       function(x) x %in% installed.packages()[, 1])


## conda install -c mosek mosek


## pip install -f https://download.mosek.com/stable/wheel/index.html Mosek


## ------------------------------------------------------------------------
installed_solvers()


## ------------------------------------------------------------------------
## Problem data
m <- 501
L <- 2
h <- L / (m - 1)

## Form objective
x <- Variable(m)
y <- Variable(m)
objective <- Minimize(sum(y))

## Form constraints
constraints <- list(x[1] == 0, y[1] == 1,
                    x[m] == 1, y[m] == 1,
                    diff(x)^2 + diff(y)^2 <= h^2)

## Solve the catenary problem
prob <- Problem(objective, constraints)
result <- solve(prob)


## ------------------------------------------------------------------------
cat("Solution status is", result$status)


## ------------------------------------------------------------------------
if ("MOSEK" %in% installed_solvers()) {
    result <- solve(prob, solver = "MOSEK")
    cat("Solution status is", result$status)
} else {
    cat("Solution not available as MOSEK is not installed!")
}


## ---- eval=FALSE---------------------------------------------------------
## setClass("QuadOverLin", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")


## ---- eval=FALSE---------------------------------------------------------
## quad_over_lin <- function(x, y) { new("QuadOverLin", x = x, y = y) }


## ---- eval=FALSE---------------------------------------------------------
## setMethod("initialize", "QuadOverLin", function(.Object, ..., x, y) {
##   .Object@x <- x
##   .Object@y <- y
##   callNextMethod(.Object, ..., args = list(.Object@x, .Object@y))
## })


## ---- eval=FALSE---------------------------------------------------------
## setMethod("initialize", "Atom", function(.Object, ..., args = list(), .size = NA_real_) {
##   .Object@args <- lapply(args, as.Constant)
##   validate_args(.Object)
##   .Object@.size <- size_from_args(.Object)
##   callNextMethod(.Object, ...)
## })


## ---- eval=FALSE---------------------------------------------------------
## setMethod("validate_args", "QuadOverLin", function(object) {
##   if(!is_scalar(object@args[[2]]))
##     stop("[QuadOverLin: validation] y must be a scalar")
## })


## ---- eval=FALSE---------------------------------------------------------
## setMethod(".domain", "QuadOverLin", function(object) { list(object@args[[2]] >= 0) })


## ---- eval=FALSE---------------------------------------------------------
## setMethod("to_numeric", "QuadOverLin", function(object, values) { sum(values[[1]]^2) / values[[2]] })


## ---- eval=FALSE---------------------------------------------------------
## setMethod("size_from_args", "QuadOverLin", function(object) { c(1,1) })
## setMethod("sign_from_args",  "QuadOverLin", function(object) { c(TRUE, FALSE) })
## setMethod("is_atom_convex", "QuadOverLin", function(object) { TRUE })
## setMethod("is_atom_concave", "QuadOverLin", function(object) { FALSE })


## ---- eval=FALSE---------------------------------------------------------
## setMethod("is_incr", "QuadOverLin", function(object, idx) { (idx == 1) && is_positive(object@args[[idx]]) })
## setMethod("is_decr", "QuadOverLin", function(object, idx) { ((idx == 1) && is_negative(object@args[[idx]])) || (idx == 2) })


## ---- eval=FALSE---------------------------------------------------------
## setMethod(".grad", "QuadOverLin", function(object, values) {
##   X <- values[[1]]
##   y <- as.numeric(values[[2]])
##   if(y <= 0)
##     return(list(NA_real_, NA_real_))
##   else {
##     # DX = 2X/y, Dy = -||X||^2_2/y^2
##     Dy <- -sum(X^2)/y^2
##     Dy <- Matrix(Dy, sparse = TRUE)
##     DX <- 2.0*X/y
##     DX <- Matrix(as.numeric(t(DX)), sparse = TRUE)
##     return(list(DX, Dy))
##   }
## })


## ---- eval=FALSE---------------------------------------------------------
## QuadOverLin.graph_implementation <- function(arg_objs, size, data = NA_real_) {
##   x <- arg_objs[[1]]
##   y <- arg_objs[[2]]   # Known to be a scalar.
##   t <- create_var(c(1,1))
##   two <- create_const(2, c(1,1))
##   constraints <- list(SOC(lo.sum_expr(list(y, t)),
##                           list(lo.sub_expr(y, t),
##                                lo.mul_expr(two, x, x$size))),
##                       create_geq(y))
##   list(t, constraints)
## }


## ---- eval=FALSE---------------------------------------------------------
## x <- Variable(n)
## y <- Variable()
## obj <- quad_over_lin(x,y)
## constr <- c(domain(obj), A*x == b)
## prob <- Problem(Minimize(obj), constr)
## solve(prob)

