reset_phylosql_state <- function() {
  env <- phylosql:::.phylosql_state
  rm(list = ls(env), envir = env)
}
