#' @title Find Coalescence
#'


#
#
# # trace coalescent times back from the ancestry matrix
# ind1 <- 1
# haplo1 <- 1
# ind2 <- 1
# haplo2 <- 2
#
# row1 <- sum(coi[1:ind1]) - coi[ind1] + haplo1
# row2 <- sum(coi[1:ind2]) - coi[ind2] + haplo2
#
# pointer1 <- rep(row1, L)
# pointer2 <- rep(row2, L)
#
# has_coalesced <- rep(NA, L)
# for (g in length(anc):1) {
#   pointer1 <- anc[[g]][cbind(pointer1, 1:L)]
#   pointer2 <- anc[[g]][cbind(pointer2, 1:L)]
#
#   for (l in 1:L) {
#     if (is.na(has_coalesced[l])) {
#       if (pointer1[l] == pointer2[l]) {
#         has_coalesced[l] <- length(anc)-g + 1 # 1-based
#       }
#     }
#   }
#
# }
#
# plot(pos, has_coalesced, type = 's', ylim = c(1, (max(has_coalesced)+5)))
# has_coalesced
