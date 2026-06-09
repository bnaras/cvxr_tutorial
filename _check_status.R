check_solver_status <- function(object) {
    s <- status(object)
    if (!s %in% c("optimal", "optimal_inaccurate")) {
        cli::cli_warn(paste("Solver status:", s))
    }
    invisible(s)
}
