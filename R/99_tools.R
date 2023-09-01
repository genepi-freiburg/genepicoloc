split_params_df <- function(params_df) {
  split(params_df, seq(nrow(params_df)))
}
