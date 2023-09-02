params_df_to_list <- function(params_df) {
  split(params_df, seq(nrow(params_df)))
}
params_df_to_chunks <- function(params_df,
                                chunk_size = 1e4) {
  nrow_var <- nrow(params_df)
  N_chunks <- ceiling(nrow_var / chunk_size)
  params_chunks <- lapply(1:N_chunks, function(i) {
    i_end <- i*chunk_size
    if (i == N_chunks) { i_end <- nrow_var }
    params_df[((i-1)*chunk_size+1):i_end,]
  })
  return(params_chunks)
}
