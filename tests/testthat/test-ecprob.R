library(testthat)
library(EdgeCount)

test_that("ECProb constructor", {
  
  toggle <- TRUE
  if (toggle){
    ecgraph <- ECGraph(test_path("network/toy.txt"))
    ecprob <- ECProb(ecgraph)
    # adj list
    expect_setequal(ecprob@adj$A, c("E"))
    expect_setequal(ecprob@adj$B, c("E","C"))
    expect_setequal(ecprob@adj$C, c("B","E","F"))
    expect_setequal(ecprob@adj$D, c("G"))
    expect_setequal(ecprob@adj$E, c("A", "B", "C", "F"))
    expect_setequal(ecprob@adj$F, c("E", "C"))
    expect_setequal(ecprob@adj$G, c("D", "H"))
    expect_setequal(ecprob@adj$H, c("G"))
    # degrees
    expect_equal(ecprob@degrees[["A"]], 1)
    expect_equal(ecprob@degrees[["B"]], 2)
    expect_equal(ecprob@degrees[["C"]], 3)
    expect_equal(ecprob@degrees[["D"]], 1)
    expect_equal(ecprob@degrees[["E"]], 4)
    expect_equal(ecprob@degrees[["F"]], 2)
    expect_equal(ecprob@degrees[["G"]], 2)
    expect_equal(ecprob@degrees[["H"]], 1)
    # names
    expect_setequal(ecprob@names, c("A", "B", "C", "D", "E", "F", "G", "H"))
    # graph size
    expect_equal(ecprob@graph_size, 8)
  }
})

test_that("ECProb get_lambda_between", {
  
  simulation <- TRUE
  if (simulation){
    ecgraph <- ECGraph(test_path("network/final_network.txt"))
    ecprob <- ECProb(ecgraph)
    nsim <- 1000
      sz <- 2^(c(1:9))
      sz1 <- NULL
      sz2 <- NULL
      sz1d <- NULL
      sz2d <- NULL
      lamb_n <- NULL
      lamb_o <- NULL
      time_n <- NULL
      time_o <- NULL
      outlier <- 0
      for (i in 1:length(sz)){
        for (j in 1:length(sz)){
          print(paste(i, j))
          for (k in 1:nsim){
            set1 <- as.character(sample(ecprob@names,sz[i]))
            t <- setdiff(ecprob@names,set1)
            set2 <- as.character(sample(t,sz[j]))
            set1d <- setdiff(set1,set2)
            set2d <- setdiff(set2,set1)
            time_start <- Sys.time()
            lambda_naive <- calculate_lambda_between_naive(ecprob, set1, set2)
            time_naive <- Sys.time() - time_start
            time_start <- Sys.time()
            lambda_optimal <- calculate_lambda_between(ecprob, set1, set2)
            time_optimal <- Sys.time() - time_start
            sz1 <- c(sz1, sz[i])
            sz2 <- c(sz2, sz[j])
            sz1d <- c(sz1d,length(set1d))
            sz2d <- c(sz2d,length(set2d))
            lamb_n <- c(lamb_n, lambda_naive)
            lamb_o <- c(lamb_o, lambda_optimal)
            time_n <- c(time_n, as.numeric(time_naive, units = "secs")*1000)
            time_o <- c(time_o, as.numeric(time_optimal, units = "secs")*1000)
          }
        }
      }
    lst <- list(sz1 = sz1, sz2 = sz2, sz1d = sz1d, sz2d = sz2d, lamb_n = lamb_n, lamb_o = lamb_o, time_n = time_n, time_o = time_o)
    df <- data.frame(lst)
    write.table(df,test_path("res/TestECPBetweenBig.txt"),sep="\t",row.names = F,quote = F)
  }

  df <- read.table(test_path("res/TestECPBetweenBig.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$sz1d_sz2d <- df$sz1d * df$sz2d
  df$sz1d_logsz1d_sz2d_logsz2d <- df$sz1d * log(pmax(1,df$sz1d)) + df$sz2d * log(pmax(1,df$sz2d))
  naive_summary <- aggregate(
    cbind(lamb_n, time_n, time_o, sz1d_logsz1d_sz2d_logsz2d) ~ sz1d_sz2d,
    data = df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  write.table(naive_summary,test_path("res/ECPBetweenNaive.txt"),sep="\t",row.names = F,quote = F)
  
  df <- read.table(test_path("res/TestECPBetweenBig.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$sz1d_sz2d <- df$sz1d * df$sz2d
  df$sz1d_logsz1d_sz2d_logsz2d <- df$sz1d * log(pmax(1,df$sz1d)) + df$sz2d * log(pmax(1,df$sz2d))
  optimized_summary <- aggregate(
    cbind(lamb_o, time_o, time_n, sz1d_sz2d) ~ sz1d_logsz1d_sz2d_logsz2d,
    data = df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  write.table(optimized_summary,test_path("res/ECPBetweenOptimized.txt"),sep="\t",row.names = F,quote = F)
  
  df <- read.table(test_path("res/TestECPBetweenBig.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  t <- df$lamb_o-df$lamb_n
  t <- t/((df$lamb_o+df$lamb_n)/2)
  
})

test_that("ECProb get_lambda_in", {
  
  simulation <- TRUE
  if (simulation){
    ecgraph <- ECGraph(test_path("network/final_network.txt"))
    ecprob <- ECProb(ecgraph)
    nsim <- 1000
    sz <- ceiling(1.37^(c(1:20))+0.5)
    lamb_n <- NULL
    lamb_o <- NULL
    lamb_f <- NULL
    time_n <- NULL
    time_o <- NULL
    time_f <- NULL
    size <- NULL
    for (i in 1:length(sz)){
      print(sz[i])
      for (sim in 1:nsim){
        set <- as.character(sample(ecprob@names,sz[i]))
        time_start <- Sys.time()
        lambda_naive <- calculate_lambda_in_naive(ecprob, set)
        time_naive <- Sys.time() - time_start
        time_start <- Sys.time()
        lambda_optimal <- calculate_lambda_in(ecprob, set)
        time_optimal <- Sys.time() - time_start
        time_start <- Sys.time()
        lambda_fast <- calculate_lambda_in_fast(ecprob, set)
        time_fast <- Sys.time() - time_start
        size <- c(size, sz[i])
        lamb_n <- c(lamb_n, lambda_naive)
        lamb_o <- c(lamb_o, lambda_optimal)
        lamb_f <- c(lamb_f,lambda_fast)
        time_n <- c(time_n, as.numeric(time_naive, units = "secs")*1000)
        time_o <- c(time_o, as.numeric(time_optimal, units = "secs")*1000)
        time_f <- c(time_f, as.numeric(time_fast, units = "secs")*1000)
      }
    }
    lst <- list(size = size, lamb_n = lamb_n, lamb_o = lamb_o, lamb_f = lamb_f, time_n = time_n, time_o = time_o, time_f = time_f)
    df <- data.frame(lst)
    write.table(df,test_path("res/TestECPInBig.txt"),sep="\t",row.names = F,quote = F)
  }
  
  df <- read.table(test_path("res/TestECPINBig.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$size_size <- df$size * df$size
  df$size_logsize <- df$size * log(df$size)
  naive_summary <- aggregate(
    cbind(lamb_n, time_n, time_o, time_f, size_logsize) ~ size_size,
    data = df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  write.table(naive_summary,test_path("res/ECPInNaive.txt"),sep="\t",row.names = F,quote = F)
  
  df <- read.table(test_path("res/TestECPInBig.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$size_size <- df$size * df$size
  df$size_logsize <- df$size * log(df$size)
  optimized_summary <- aggregate(
    cbind(lamb_o, time_o, time_n, time_f, size_size) ~ size_logsize,
    data = df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  write.table(optimized_summary,test_path("res/ECPInOptimized.txt"),sep="\t",row.names = F,quote = F)
  
  df <- read.table(test_path("res/TestECPInBig.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$size_size <- df$size * df$size
  lambda_summary <- aggregate(
    cbind(lamb_o, lamb_n, lamb_f, size_size) ~ size_size,
    data = df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  write.table(lambda_summary,test_path("res/ECPInLambda.txt"),sep="\t",row.names = F,quote = F)
})



test_that("ECProb p-value", {
  
  get_alpha <- function(lambda, m){
    
    t <- 1
    alpha <- 1
    for (i in 1:m){
      t <- t*(lambda/i)
      if (t < alpha / 1e12){
        break
      }
      alpha <- alpha + t
    }
    alpha <- alpha*exp(-lambda)
    return(as.numeric(alpha))
  }
  
  get_upper <- function(lambda, m, z){
    
    t <- lambda^z/factorial(z)
    p <- t
    for (i in (z+1):m){
      t <- t*(lambda/i)
      p <- p + t
    }
    p <- p*exp(-lambda)
    return(as.numeric(p))
  }
  
  get_PV <- function(lambda, m, z){
  
    p <- get_upper(lambda, m, z)
    p <- p/get_alpha(lambda, m)
    return(as.numeric(p))
  }
  
  sim <- FALSE
  if (sim){
    ecgraph <- ECGraph(test_path("network/final_network.txt"))
    ecprob <- ECProb(ecgraph)
    szs <- c(2:10)
    pn <- NULL
    p <- NULL
    sz <- NULL
    m <- NULL
    z <- NULL
    lambda <- NULL
    for (i in 1:length(szs)){
      set <- as.character(sample(ecprob@names,szs[i]))
      M <- length(set)*(length(set)-1)/2
      lamb <- calculate_lambda_in(ecprob, set)
      for (Z in 0:M){
        p1 <- get_PV(lamb, M, Z)
        lst <- edge_count_statistics(ecprob, Z, M, lamb)
        p2 <- lst$p_value
        pn <- c(pn, p1)
        p <- c(p, p2)
        sz <- c(sz, szs[i])
        m <- c(m, M)
        z <- c(z, Z)
      }
    }
    adiff = abs(pn-p)
    df <- data.frame(sz = sz, m = m, z = z, pn = pn, p = p, adiff = adiff)
    write.table(df, test_path("res/TestPValue.txt"),quote = F, sep = "\t", row.names = F)
    
    szs <- c(2:100)
    sz <- NULL
    m <- NULL
    z <- NULL
    tm <- NULL
    for (i in 1:length(szs)){
      print(szs[i])
      set <- as.character(sample(ecprob@names,szs[i]))
      M <- length(set)*(length(set)-1)/2
      lamb <- calculate_lambda_in(ecprob, set)
      for (Z in 0:M){
        sz <- c(sz, szs[i])
        m <- c(m, M)
        z <- c(z, Z)
        time_start <- Sys.time()
        lst <- edge_count_statistics(ecprob, Z, M, lamb)
        p <- lst$p_value
        time_p <- Sys.time() - time_start
        tm <- c(tm, as.numeric(time_p, units = "secs")*1000)
      }
    }
    df <- data.frame(sz = sz, m = m, z = z, tm = tm)
    write.table(df, test_path("res/TestPValueTime.txt"),quote = F, sep = "\t", row.names = F)
    sz_summary <- aggregate(
      tm ~ sz,
      data = df,
      FUN = function(x) c(mean = mean(x), sd = sd(x))
    )
    write.table(sz_summary, test_path("res/TestPValueTimeSz.txt"),quote = F, sep = "\t", row.names = F)
    m_summary <- aggregate(
      tm ~ m,
      data = df,
      FUN = function(x) c(mean = mean(x), sd = sd(x))
    )
    write.table(m_summary, test_path("res/TestPValueTimeM.txt"),quote = F, sep = "\t", row.names = F)
    z_summary <- aggregate(
      tm ~ z,
      data = df,
      FUN = function(x) c(mean = mean(x), sd = sd(x))
    )
    write.table(z_summary, test_path("res/TestPValueTimeZ.txt"),quote = F, sep = "\t", row.names = F)
  }
})

test_that("ECProb vectorize in", {
  
  find_kt_i <- function(i, k, m, M) {
    for (j in (i + 1):m) {
      if (k[i] * k[j] <= 2 * M) {
        return(j)
      }
    }
    return(m)
  }
  vfind_kt_i <- Vectorize(find_kt_i, vectorize.args = "i")
    
  ecgraph <- ECGraph(test_path("network/final_network.txt"))
  ecprob <- ECProb(ecgraph)
  M <- ecprob@graph_size
  n <- 1000
  set <- as.character(sample(ecprob@names,n))
  k <- as.numeric(unlist(ecprob@degrees[set]))
  kc <- numeric(n+1)
  kc[n+1] <- 0
  for (i in length(k):1){
    kc[i] <- kc[i+1]+k[i]
  }
  kc2 <- c(rev(cumsum(rev(k))), 0)
  expect_equal(kc, kc2)
  
  kt <- rep(n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (k[i]*k[j] <= 2*M) {
        kt[i] <- j
        break
      }
    }
  }
  kt2 <- vfind_kt_i(c(1:(n-1)), k, n, M)
  kt2 <- c(kt2, n)
  expect_equal(kt, kt2)
})

test_that("ECProb vectorize between", {
  
  k <- c(1:30)
  nsim <- 10
  flag <- 0
  for (sim in 1:nsim){
    k1 <- c(1:20)
    k2 <- c(5:30)
    M <- 200
    
    k1 <- sort(k1, decreasing = TRUE)
    k2 <- sort(k2, decreasing = TRUE)
    lambda <- 0
    
    lambda_n <- 0
    for (i in 1:length(k1)){
      for (j in 1:length(k2)){
        lambda_n <- lambda_n + min(1, k1[i]*k2[j]/(2*M))
      }
    }
    
    k1t <- rep(length(k2),length(k1))
    for (i in 1:length(k1)){
      for (j in 1:length(k2)){
        if (k1[i] * k2[j] <= 2*M){
          k1t[i] <- j
          break
        }
      }
    }
    find_k1t <- function(i, k1, k2, M){
      for (j in 1:length(k2)) {
        if (k1[i] * k2[j] <= 2 * M) {
          return(j)
        }
      }
      return(length(k1))
    }
    vfind_k1t <- Vectorize(find_k1t, vectorize.args = "i")
    k1tb <- vfind_k1t(c(1:length(k1)), k1, k2, M)
    if (!identical(k1t, k1tb)){
      flag <- 1
    }

    k2c <- numeric(length(k2)+1)
    k2c[length(k2)+1] <- 0
    for (i in (length(k2):1)){
      k2c[i] <- k2c[i+1]+k2[i]
    }
    k2cb <- c(rev(cumsum(rev(k2))), 0)
    if (!identical(k2c, k2cb)){
      flag <- 1
    }
    
    
    for (i in 1:length(k1)){
      lambda <- lambda + (k1t[i]-1) + k1[i]*k2c[k1t[i]]/(2*M)
    }
    
    lambda <- sum((k1t - 1) + (k1 * k2c[k1t]) / (2 * M))
    if (abs(lambda-lambda_n)/lambda_n > 1e-12){
     flag <- 1 
    }
  }
  expect_equal(flag, 0)
})

test_that("ECProb fast suitability", {
  
  ecgraph <- ECGraph(test_path("network/final_network.txt"))
  ecprob <- ECProb(ecgraph)
  res <- summarize_suitability_fast(ecprob)
  # print(res)
})

test_that("ECProb lambda_expected", {
  
  ecgraph <- ECGraph(test_path("network/final_network.txt"))
  ecprob <- ECProb(ecgraph)
  M <- ecprob@graph_size
  N <- ecprob@graph_order
  m <- c(2, 3, 4, 5, 10, 20, 50, 100)
  sim <- FALSE
  if (sim){
    nsim <- 100000
    m_lambda <- NULL
    sd_lambda <- NULL
    e_lambda <- NULL
    for (i in 1:length(m)){
      lambda <- numeric(nsim)
      for (j in 1:nsim){
        vset <- sample(ecprob@names, m[i])
        lambda[j] <- calculate_lambda_in_fast(ecprob, vset)
      }
      m_lambda <- c(m_lambda, mean(lambda))
      sd_lambda <- c(sd_lambda, sd(lambda))
      e_lambda <- c(e_lambda, m[i]*(m[i]-1)*M/(N*N))
    }
    df <- data.frame(m = m, m_lambda = m_lambda, e_lambda = e_lambda, sd_lambda = sd_lambda)
    write.table(df, test_path("res/TestExpectedLambdaIn.txt"), quote = F, sep = "\t", row.names = F)
  }
  # print(df)
})
