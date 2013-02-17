bb.post4 <- bern.beta(c(20, 40), learn.sub$Score)
plot(bb.post4)
abline(v = mu.ci, col = "red")

mu2 <- (sum(learn.sub$Score) + 40)/120
n2 <- 120
mu2.ci <- c(mu2 - 1.96 * std.err.perc(mu2, n2), mu2 + 1.96 * std.err.perc(mu2, n2))
plot(bb.post3)
abline(v = mu2.ci, col = "red")


learn.sub.half1 <- learn.sub[1:(nrow(learn.sub)/2), ]
learn.sub.half2 <- learn.sub[(nrow(learn.sub)/2 + 1):nrow(learn.sub), ]

bb.h1.post1 <- bern.beta(c(1, 1), learn.sub.half1$Score)
plot(bb.h1.post1)
bb.h2.post1 <- bern.beta(unlist(bb.h1.post1$post.shape), learn.sub.half2$Score)
plot(bb.h2.post1)

bb.h2.post2 <- bern.beta(c(1, 1), learn.sub.half2$Score)
plot(bb.h2.post2)
bb.h1.post2 <- bern.beta(unlist(bb.h2.post2$post.shape), learn.sub.half1$Score)
plot(bb.h1.post2)

learn.sub.misses <- which(learn.sub$Score == 0)
learn.sub.hits <- which(learn.sub$Score == 1)
extra.hits <- length(learn.sub.hits) - nrow(learn.sub)/2
learn.sub.part1 <- learn.sub[c(learn.sub.misses, learn.sub.hits[1:extra.hits]), ]
learn.sub.part2 <- learn.sub[learn.sub.hits[(extra.hits+1):length(learn.sub.hits)], ]

bb.p1.post <- bern.beta(c(1, 1), learn.sub.part1$Score)
plot(bb.p1.post)
bb.p2.post <- bern.beta(unlist(bb.p1.post$post.shape), learn.sub.part2$Score)
plot(bb.p2.post)
