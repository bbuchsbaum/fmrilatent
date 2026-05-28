# Tests for the bounded handle registry cache: growth bounds, LRU eviction,
# and strict fingerprint reuse (bd-01KSQQNHNXNBHBNGVQGA1RS9B2).

reg_register <- function(...) fmrilatent:::.latent_register_matrix(...)
reg_get      <- function(...) fmrilatent:::.latent_get_matrix(...)

test_that("registry growth is bounded by max_entries with LRU eviction", {
  fmrilatent_registry_clear()
  old <- getOption("fmrilatent.registry.max_entries")
  on.exit(options(fmrilatent.registry.max_entries = old), add = TRUE)
  options(fmrilatent.registry.max_entries = 5L)

  for (i in seq_len(20L)) {
    reg_register(paste0("g", i), matrix(i, 2, 2), type = "basis")
  }

  stats <- fmrilatent_registry_stats("basis")
  expect_lte(stats$basis$count, 5L)

  # The 5 most-recently-registered survive; older ones were evicted.
  ids <- fmrilatent_registry_list("basis")
  expect_setequal(unname(ids), paste0("g", 16:20))
  expect_null(reg_get("g1", "basis"))
  expect_false(is.null(reg_get("g20", "basis")))

  fmrilatent_registry_clear()
})

test_that("retrieval refreshes LRU recency so touched entries survive", {
  fmrilatent_registry_clear()
  old <- getOption("fmrilatent.registry.max_entries")
  on.exit(options(fmrilatent.registry.max_entries = old), add = TRUE)
  options(fmrilatent.registry.max_entries = 3L)

  reg_register("a", matrix(1, 2, 2), type = "basis")
  reg_register("b", matrix(2, 2, 2), type = "basis")
  reg_register("c", matrix(3, 2, 2), type = "basis")

  # Touch "a" so it is now most-recently-used; "b" becomes the eviction target.
  expect_false(is.null(reg_get("a", "basis")))
  reg_register("d", matrix(4, 2, 2), type = "basis")

  expect_setequal(unname(fmrilatent_registry_list("basis")), c("a", "c", "d"))
  expect_null(reg_get("b", "basis"))

  fmrilatent_registry_clear()
})

test_that("max_entries = Inf restores an unbounded cache", {
  fmrilatent_registry_clear()
  old <- getOption("fmrilatent.registry.max_entries")
  on.exit(options(fmrilatent.registry.max_entries = old), add = TRUE)
  options(fmrilatent.registry.max_entries = Inf)

  for (i in seq_len(50L)) {
    reg_register(paste0("u", i), matrix(i, 2, 2), type = "basis")
  }
  expect_equal(fmrilatent_registry_stats("basis")$basis$count, 50L)

  fmrilatent_registry_clear()
})

test_that("basis and loadings registries are capped independently", {
  fmrilatent_registry_clear()
  old <- getOption("fmrilatent.registry.max_entries")
  on.exit(options(fmrilatent.registry.max_entries = old), add = TRUE)
  options(fmrilatent.registry.max_entries = 2L)

  reg_register("b1", matrix(1, 2, 2), type = "basis")
  reg_register("b2", matrix(1, 2, 2), type = "basis")
  reg_register("b3", matrix(1, 2, 2), type = "basis")
  reg_register("l1", matrix(1, 3, 3), type = "loadings")
  reg_register("l2", matrix(1, 3, 3), type = "loadings")
  reg_register("l3", matrix(1, 3, 3), type = "loadings")

  expect_equal(fmrilatent_registry_stats("basis")$basis$count, 2L)
  expect_equal(fmrilatent_registry_stats("loadings")$loadings$count, 2L)

  fmrilatent_registry_clear()
})

test_that("clear resets LRU bookkeeping so a refilled cache stays bounded", {
  fmrilatent_registry_clear()
  old <- getOption("fmrilatent.registry.max_entries")
  on.exit(options(fmrilatent.registry.max_entries = old), add = TRUE)
  options(fmrilatent.registry.max_entries = 4L)

  for (i in seq_len(10L)) reg_register(paste0("p", i), matrix(i, 2, 2), type = "basis")
  expect_lte(fmrilatent_registry_stats("basis")$basis$count, 4L)

  fmrilatent_registry_clear()
  expect_equal(fmrilatent_registry_stats("basis")$basis$count, 0L)

  for (i in seq_len(10L)) reg_register(paste0("q", i), matrix(i, 2, 2), type = "basis")
  expect_lte(fmrilatent_registry_stats("basis")$basis$count, 4L)
  expect_setequal(unname(fmrilatent_registry_list("basis")), paste0("q", 7:10))

  fmrilatent_registry_clear()
})

test_that("fingerprinted reuse rejects an entry with a different fingerprint", {
  fmrilatent_registry_clear()

  reg_register("fp", matrix(1, 2, 2), type = "basis", fingerprint = "AAAA")
  expect_error(
    reg_register("fp", matrix(2, 2, 2), type = "basis", fingerprint = "BBBB"),
    "different fingerprint"
  )
  # Original entry is untouched (value and fingerprint preserved).
  cached <- reg_get("fp", "basis")
  expect_equal(as.vector(cached), rep(1, 4))
  expect_identical(attr(cached, "fmrilatent.handle_fingerprint"), "AAAA")

  fmrilatent_registry_clear()
})

test_that("fingerprinted reuse rejects an entry that lacks a fingerprint", {
  fmrilatent_registry_clear()

  # Legacy registration without a fingerprint.
  reg_register("legacy", matrix(1, 2, 2), type = "basis")
  expect_error(
    reg_register("legacy", matrix(1, 2, 2), type = "basis", fingerprint = "CCCC"),
    "without a fingerprint"
  )

  fmrilatent_registry_clear()
})

test_that("matching fingerprint reuse is a silent no-op (no warning, no error)", {
  fmrilatent_registry_clear()

  reg_register("same", matrix(1, 2, 2), type = "basis", fingerprint = "DDDD")
  expect_silent(
    result <- reg_register("same", matrix(1, 2, 2), type = "basis", fingerprint = "DDDD")
  )
  expect_false(result)

  fmrilatent_registry_clear()
})

test_that("non-fingerprinted collision still warns rather than errors", {
  fmrilatent_registry_clear()

  reg_register("legacy2", matrix(1, 2, 2), type = "basis")
  expect_warning(
    reg_register("legacy2", matrix(2, 2, 2), type = "basis"),
    "already registered"
  )

  fmrilatent_registry_clear()
})
