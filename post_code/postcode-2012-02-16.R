# Code from blog post here: 
source("../code/Chapter5/bernbetaShPsy.R")
source("../kruschke_code/HDIofICDF.R")
std.err.perc <- function(p, n) { sqrt(p*(1-p)/(n)) }
learn <- read.delim("../data/LangLearnPilot.txt")
learn.sub <- droplevels(learn[learn$Apt == "match.CaseB" &
                              learn$Test == 3 &
                              learn$Condition == "CaseB", ])

mu <- mean(learn.sub$Score)
n <- length(learn.sub$Score)
mu.ci <- c(mu - 1.96 * std.err.perc(mu, n), mu + 1.96 * std.err.perc(mu, n))

bb.post1 <- bern.beta(c(1, 1), learn.sub$Score)
plot(bb.post1)
abline(v = mu.ci, col = "red")

bb.post2 <- bern.beta(c(30, 30), learn.sub$Score)
plot(bb.post2)
abline(v = mu.ci, col = "red")

bb.post3 <- bern.beta(c(40, 20), learn.sub$Score)
plot(bb.post3)
abline(v = mu.ci, col = "red")
