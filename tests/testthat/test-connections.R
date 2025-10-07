test_that("credentials must exist and are cached", {
  reset_phylosql_state()

  creds <- data.frame(
    host = "localhost",
    dbname = "microbiome",
    port = 3306,
    user = "tester",
    stringsAsFactors = FALSE
  )

  cred_path <- tempfile(fileext = ".csv")
  on.exit(unlink(cred_path), add = TRUE)
  utils::write.csv(creds, cred_path, row.names = FALSE)

  fake_connection <- structure(list(id = 1L), class = "FakeConnection")

  with_mocked_bindings({
    con1 <- get_mtgn_connection(path = cred_path, key = "secret")
    expect_identical(con1, fake_connection)
    expect_identical(phylosql:::.phylosql_state$key, "secret")

    con2 <- get_mtgn_connection()
    expect_identical(con2, fake_connection)
  },
  DBI::dbConnect = function(...) fake_connection,
  DBI::dbDisconnect = function(...) NULL,
  DBI::dbIsValid = function(con) TRUE)
})


test_that("refresh forces reconnection", {
  reset_phylosql_state()

  creds <- data.frame(
    host = "localhost",
    dbname = "microbiome",
    port = 3306,
    user = "tester",
    stringsAsFactors = FALSE
  )

  cred_path <- tempfile(fileext = ".csv")
  on.exit(unlink(cred_path), add = TRUE)
  utils::write.csv(creds, cred_path, row.names = FALSE)

  calls <- 0
  fake_connection <- structure(list(id = 1L), class = "FakeConnection")

  with_mocked_bindings({
    get_mtgn_connection(path = cred_path, key = "secret")
    get_mtgn_connection(refresh = TRUE)
    expect_equal(calls, 2)
  },
  DBI::dbConnect = function(...) {
    calls <<- calls + 1
    fake_connection
  },
  DBI::dbDisconnect = function(...) NULL,
  DBI::dbIsValid = function(con) calls > 0)
})
