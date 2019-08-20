build_vignettes_to_inst <- function(rebuild = FALSE) {
  if (rebuild) {
    # Build vignettes to 'doc' and 'Meta'. Update '.gitignore'.
    devtools::build_vignettes() 
  }
  # Recreate empty directories
  unlink(c("inst/doc", "inst/Meta"), recursive = TRUE)
  dir.create("inst/doc")
  dir.create("inst/Meta") 
  # Copy files to 'inst' subfolders
  has_worked <- c(
    file.copy(list.files("doc", full.names = TRUE), to = "inst/doc"), 
    file.copy(list.files("Meta", full.names = TRUE), to = "inst/Meta")
  )
  # Return TRUE if everything worked OK
  return(all(has_worked)) 
}
