library(data.table)
library(EdgeCount)

data("sample_ects")
t <- table(lengths(sample_ects@ecprob@adj[sample_ects@terms]))
u <- cumsum(rev(t))
print(u)
v <- round(seq(200, 3000, length.out = 10))
print(v)
