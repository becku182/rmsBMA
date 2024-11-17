test_that("ModelSpace builds correct ols_results table (1)", {
  x1<-rnorm(20, mean = 0, sd = 1)
  x2<-rnorm(20, mean = 0, sd = 2)
  x3<-rnorm(20, mean = 0, sd = 3)
  x4<-rnorm(20, mean = 0, sd = 1)
  x5<-rnorm(20, mean = 0, sd = 2)
  x6<-rnorm(20, mean = 0, sd = 4)
  e<-rnorm(20, mean = 0, sd = 0.5)
  y<-2+x1+2*x2+e
  data<-cbind(y,x1,x2,x3,x4,x5,x6)
  M<-6
  mSpace<-modelSpace(data,M)
  expect_equal(as.numeric(mSpace[3]),2^M)
})

test_that("ModelSpace builds correct ols_results table (2)", {
  x1<-rnorm(20, mean = 0, sd = 4)
  x2<-rnorm(20, mean = 0, sd = 2)
  x3<-rnorm(20, mean = 0, sd = 6)
  x4<-rnorm(20, mean = 0, sd = 1)
  x5<-rnorm(20, mean = 0, sd = 5)
  x6<-rnorm(20, mean = 0, sd = 4)
  e<-rnorm(20, mean = 0, sd = 2)
  y<-2+x1+2*x2+e
  data<-cbind(y,x1,x2,x3,x4,x5,x6)
  M<-6
  mSpace<-modelSpace(data,M)
  expect_equal(as.numeric(mSpace[3]),(2^M))
})

test_that("ModelSpace builds correct ols_results table (3)", {
  x1<-rnorm(20, mean = 0, sd = 1)
  x2<-rnorm(20, mean = 0, sd = 2)
  x3<-rnorm(20, mean = 0, sd = 3)
  x4<-rnorm(20, mean = 0, sd = 1)
  x5<-rnorm(20, mean = 0, sd = 2)
  x6<-rnorm(20, mean = 0, sd = 4)
  e<-rnorm(20, mean = 0, sd = 0.5)
  y<-2+x1+2*x2+e
  data<-cbind(y,x1,x2,x3,x4,x5,x6)
  M<-4
  mSpace<-modelSpace(data,M)
  expect_equal(as.numeric(mSpace[3]),1+choose(6,1)+choose(6,2)+choose(6,3)+choose(6,4))
})

test_that("ModelSpace builds correct ols_results table (4)", {
  x1<-rnorm(20, mean = 0, sd = 1)
  x2<-rnorm(20, mean = 0, sd = 2)
  x3<-rnorm(20, mean = 0, sd = 3)
  x4<-rnorm(20, mean = 0, sd = 1)
  x5<-rnorm(20, mean = 0, sd = 2)
  x6<-rnorm(20, mean = 0, sd = 4)
  e<-rnorm(20, mean = 0, sd = 0.5)
  y<-2+x1+2*x2+e
  data<-cbind(y,x1,x2,x3,x4,x5,x6)
  M<-3
  mSpace<-modelSpace(data,M)
  expect_equal(as.numeric(mSpace[3]),1+choose(6,1)+choose(6,2)+choose(6,3))
})


test_that("ModelSpace builds correct ols_results table (5)", {
  x1<-rnorm(20, mean = 0, sd = 1)
  x2<-rnorm(20, mean = 0, sd = 2)
  x3<-rnorm(20, mean = 0, sd = 3)
  x4<-rnorm(20, mean = 0, sd = 1)
  x5<-rnorm(20, mean = 0, sd = 2)
  x6<-rnorm(20, mean = 0, sd = 4)
  e<-rnorm(20, mean = 0, sd = 0.5)
  y<-2+x1+2*x2+e
  data<-cbind(y,x1,x2,x3,x4,x5,x6)
  mSpace<-modelSpace(data)
  expect_equal(as.numeric(mSpace[3]),1+choose(6,1)+choose(6,2)+choose(6,3)+choose(6,4)+choose(6,5)+choose(6,6))
})
