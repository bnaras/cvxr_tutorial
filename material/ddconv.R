library(CVXR)

set.seed(1782)
n <- 5000
X1 <- rexp(n = n, rate = 0.447)
Z1 <- rnorm(n = n, mean = 0, sd = sqrt(3.2))
Y1 <- X1 + Z1

X2 <- rgamma(n = n, shape = 5, rate = 1)
Z2 <- rnorm(n = n, mean = 0, sd = sqrt(3.2))
Y2 <- X2 + Z2

K <- floor(min(200, 3 * sqrt(n)))

make_data <- function(Y, N, K, z.mean = 0, z.sd = 1,
                         penalty = c("d2", "g")) {
    penalty <- match.arg(penalty)
    f_Y <- hist(x = Y, breaks = K, plot = FALSE)$density
    x <- seq(from = min(Y), to = max(Y), length.out = K)
    arg_mat <- sapply(x, function(w) x - w)
    delta <- (x[K] - x[1]) / (K - 1)
    C_mat <- delta * dnorm(arg_mat, mean = z.mean, sd = z.sd)
    ## Second difference matrix
    D2 <- Matrix::bandSparse(n = K - 2, m = K, k = 0:2,
                             diagonals = list(rep(1, K),
                                              rep(-2, K),
                                              rep(1, K))
                             ) / (delta^2)

    list(x = x, f_Y = f_Y, C = C_mat, D2 = D2, delta = delta)
}

d <- make_data(Y = Y1, N = n, K = K, z.mean = 0, z.sd = sqrt(3.2))

f_X <- Variable(K)
lambda <- 0.011
objective <- Minimize(sum_squares(f_X - d$C %*% f_X) + lambda * p_norm(d$D2 %*% f_X))
constraints <- list(f_X >= 0,
                    sum(d$delta * f_X) == 1)
problem <- Problem(objective, constraints)
result <- solve(problem)





