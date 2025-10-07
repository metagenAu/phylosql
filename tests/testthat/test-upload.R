test_that("upload functions validate required columns", {
  dummy_con <- structure(list(), class = "Dummy")

  with_mocked_bindings({
    expect_error(upload_lab_data(NULL, con = dummy_con), "must not be NULL")
    expect_error(upload_lab_data(data.frame(foo = 1), con = dummy_con), "Missing required column")
  },
  resolve_connection = function(con = NULL) dummy_con)
})


test_that("upload_lab_data skips existing rows", {
  dummy_con <- structure(list(), class = "Dummy")
  data <- data.frame(
    MetagenNumber = 1,
    variable = "depth",
    value = 42
  )

  with_mocked_bindings({
    expect_message(upload_lab_data(data, con = dummy_con), "No new records to upload.")
  },
  resolve_connection = function(con = NULL) dummy_con,
  collect_existing_keys = function(...) data[c("MetagenNumber", "variable")],
  compute_new_indices = function(existing, candidate, key_cols) integer())
})


test_that("uploadData executes bulk loader", {
  dummy_con <- structure(list(), class = "Dummy")
  executed <- FALSE
  temp_path <- file.path(tempdir(), "phylosql-test-upload.csv")
  if (file.exists(temp_path)) {
    file.remove(temp_path)
  }

  with_mocked_bindings({
    uploadData(data.frame(x = 1), "schema.table", con = dummy_con)
    expect_true(executed)
    expect_false(file.exists(temp_path))
  },
  resolve_connection = function(con = NULL) dummy_con,
  tempfile = function(...) temp_path,
  DBI::dbQuoteString = function(con, x) DBI::SQL(sprintf("'%s'", x)),
  DBI::dbQuoteIdentifier = function(con, x) DBI::SQL(x),
  DBI::dbExecute = function(con, query) {
    executed <<- TRUE
    invisible(0)
  },
  DBI::dbWithTransaction = function(con, code) code)
})
