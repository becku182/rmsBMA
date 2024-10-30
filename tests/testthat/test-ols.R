test_that("ols gets the rights ols osbjects (1)", {
  x1<-rnorm(10, mean = 0, sd = 1)
  x2<-rnorm(10, mean = 0, sd = 2)
  x<-cbind(x1,x2)
  y<-2+x1+2*x2
  const<-1
  model<-ols(y,x,const)
  test<-lm(y~x)
  # Coefficients
  test1<-test[[1]][1:3]
  test1<-as.matrix(test1,nrow=3,ncol=1)
  test11<-matrix(0,nrow=3,ncol=1)
  for (i in 1:3){
    test11[i,1]<-test1[i,1]
  }
  expect_equal(round(as.matrix(model[[1]]),11),round(test11,11))
  # Standard errors
  test2<-summary(test)[[4]]
  test2<-as.matrix(test2,nrow=3,ncol=1)
  test2
  test21<-matrix(0,nrow=3,ncol=1)
  for (i in 1:3){
    test21[i,1]<-test2[i,2]
  }
  expect_equal(round(as.matrix(model[[2]]),11),round(test21,11))
  # R^2
  expect_equal(round(as.numeric(model[[4]],11)),round(as.numeric(summary(test)[[8]],11)))
  # Degrees of freedom
  expect_equal(as.numeric(model[5]),as.numeric(test[8]))
  # Dillution
  expect_equal(as.numeric(model[6]),det(stats::cor(x)))
})

test_that("ols gets the rights ols osbjects (2)", {
  x1<-rnorm(10, mean = 0, sd = 1)
  x2<-rnorm(10, mean = 0, sd = 2)
  e<-rnorm(10, mean = 0, sd =1)
  x<-cbind(x1,x2)
  y<-2+x1+2*x2+e
  const<-1
  model<-ols(y,x,const)
  test<-lm(y~x)
  # Coefficients
  test1<-test[[1]][1:3]
  test1<-as.matrix(test1,nrow=3,ncol=1)
  test11<-matrix(0,nrow=3,ncol=1)
  for (i in 1:3){
    test11[i,1]<-test1[i,1]
  }
  expect_equal(round(as.matrix(model[[1]]),11),round(test11,11))
  round(as.matrix(model[[1]]),11)==round(test11,11)
  # Standard errors
  test2<-summary(test)[[4]]
  test2<-as.matrix(test2,nrow=3,ncol=1)
  test2
  test21<-matrix(0,nrow=3,ncol=1)
  for (i in 1:3){
    test21[i,1]<-test2[i,2]
  }
  expect_equal(round(as.matrix(model[[2]]),11),round(test21,11))
  # R^2
  expect_equal(round(as.numeric(model[[4]],11)),round(as.numeric(summary(test)[[8]],11)))
  # Degrees of freedom
  expect_equal(as.numeric(model[5]),as.numeric(test[8]))
  # Dillution
  expect_equal(as.numeric(model[6]),det(stats::cor(x)))
})
