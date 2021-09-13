.show_message_if_batch_size_1 <- function(batch_size) {
  if(batch_size == 1) {
    message("Single-threaded mode enabled. In order to speed up computations, increase defined batch_size or use a default value")
  }
}