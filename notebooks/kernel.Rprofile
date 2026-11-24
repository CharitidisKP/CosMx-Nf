# The Positron kernel's R profile, which positron_kernel.sh sets as R_PROFILE_USER
# Giotto's Python for loadGiotto()
if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) options(giotto.py_path = Sys.getenv("RETICULATE_PYTHON"))
# The profile R reads when R_PROFILE_USER is unset
if (file.exists(".Rprofile")) {
  source(".Rprofile")
} else if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}
