# Training dataset, cached after the first load. An environment is used
# because, unlike a plain variable, it can be modified from inside a function
# after the package namespace has been locked.
.brgdgt_cache <- new.env(parent = emptyenv())

load_reference_data <- function() {
  # reading the rds takes a moment, so only do it once per session
  if (is.null(.brgdgt_cache$reference)) {
    .brgdgt_cache$reference <- readRDS(extdata_path("reference_data.rds"))
  }
  .brgdgt_cache$reference
}

# Path to a file shipped in inst/extdata. system.file() returns "" rather than
# an error when the file is missing, so turn that into an informative error.
extdata_path <- function(file) {
  path <- system.file("extdata", file, package = "envirolipids")
  if (!nzchar(path)) {
    stop("File '", file, "' not found in the envirolipids installation",
         call. = FALSE)
  }
  path
}
