test_that("introduction of fixed effects leads to columns adding up to zero", {
  y<-matrix(1:20,nrow=20,ncol=1)
  x1<-matrix(21:40,nrow=20,ncol=1)
  x2<-matrix(41:60,nrow=20,ncol=1)
  data<-cbind(y,x1,x2)
  data
  new_data<-data_prep(data,FE=1,Time=5,Section=4,Time_FE=1,Section_FE=1,STD=0)
  new_data
  expect_equal(matrix(round(apply(new_data,2,mean),11),nrow=3,ncol=1),matrix(0,nrow=3,ncol=1))
})

test_that("introduction of fixed effects and standardization leads to columns adding up to zero and their standard deviations being equal to 1", {
  y<-rnorm(20, mean = 0, sd = 1)
  x1<-rnorm(20, mean = 0, sd = 1)
  x2<-rnorm(20, mean = 0, sd = 1)
  data<-cbind(y,x1,x2)
  data
  new_data<-data_prep(data,FE=1,Time=5,Section=4,Time_FE=1,Section_FE=1,STD=1)
  new_data
  expect_equal(matrix(round(apply(new_data,2,mean),11),nrow=3,ncol=1),matrix(0,nrow=3,ncol=1))
  expect_equal(matrix(round(apply(new_data,2,stats::sd),11),nrow=3,ncol=1),matrix(1,nrow=3,ncol=1))
})

test_that("standardization leads to columns adding up to zero and their standard deviations being equal to 1", {
  y<-rnorm(100, mean = 0, sd = 1)
  x1<-rnorm(100, mean = 0, sd = 1)
  x2<-rnorm(100, mean = 0, sd = 1)
  data<-cbind(y,x1,x2)
  data
  new_data<-data_prep(data,FE=1,Time=10,Section=10,Time_FE=1,Section_FE=1,STD=1)
  new_data
  expect_equal(matrix(round(apply(new_data,2,mean),11),nrow=3,ncol=1),matrix(0,nrow=3,ncol=1))
  expect_equal(matrix(round(apply(new_data,2,stats::sd),11),nrow=3,ncol=1),matrix(1,nrow=3,ncol=1))
})
