# config.R — resolve dataset / output paths for the "complete rerun" (Option 2).
#
# Option 2 reads the locations of the raw CRISPR datasets and the output
# directory from a `~/.research_config` file that you create yourself (see the
# "Configure paths" step in the README). This helper looks up a single key from
# that file. It is sourced automatically by the analysis scripts, so the only
# thing you need to create is `~/.research_config`.
#
# Option 1 (figures from precomputed results) does not use this file at all.

.get_config_path <- function(dir_name) {
  cmd <- paste0("source ~/.research_config; echo $", dir_name)
  system(command = cmd, intern = TRUE)
}
