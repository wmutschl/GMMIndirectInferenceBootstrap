# Angrist and Krueger (1991)

library(foreign)
library(AER)
setwd("c:/temp")
x <- read.dta(file.choose())
dob <- x$yob+(x$qob-1)*0.25

# Replicate Figure I
Z <- matrix(NA,40,2)
Z[,1] <- seq(1930,1939.75,by=0.25)
for(i in 1:dim(Z)[1]) {
  Z[i,2] <- mean(x$educ[dob==Z[i,1]])
}
plot(Z,t="o",main="Figure I",xlab="Year of Birth",ylab="Years of Completed Education",ylim=c(12.2,13.2))

# Replicate Figure II
Z <- matrix(NA,40,2)
Z[,1] <- seq(1940,1949.75,by=0.25)
for(i in 1:dim(Z)[1]) {
  Z[i,2] <- mean(x$educ[dob==Z[i,1]])
}
plot(Z,t="o",main="Figure II",xlab="Year of Birth",ylab="Years of Completed Education",ylim=c(13,13.9))

# Replicate Figure V
Z <- matrix(NA,80,2)
Z[,1] <- seq(1930,1949.75,by=0.25)
for(i in 1:dim(Z)[1]) {
  Z[i,2] <- mean(x$lwklywge[dob==Z[i,1]])
}
plot(Z,t="o",main="Figure V",xlab="Year of Birth",ylab="Log Weekly Earnings")



# Replicate part of Table IV
# Drop all persons born after 1929Q4
x <- x[x$yob<1930,]
dob <- x$yob+(x$qob-1)*0.25

# Create yob dummies
Dyear <- as.factor(x$yob)

# Column (1)
regr <- lm(x$lwklywge~x$educ+Dyear)
summary(regr)

# Column (3)
age <- 1970-dob
regr <- lm(x$lwklywge~x$educ+Dyear+age+I(age^2))
summary(regr)

# Column (2)
Dq <- dob
Dq[Dq-floor(Dq)==0.75] <- 0
Dq <- factor(Dq)
regr <- ivreg(x$lwklywge~x$educ+Dyear|Dq+Dyear)
summary(regr)

# Column (4)
regr <- ivreg(x$lwklywge~x$educ+age+I(age^2)+Dyear|Dq+Dyear)
summary(regr)
