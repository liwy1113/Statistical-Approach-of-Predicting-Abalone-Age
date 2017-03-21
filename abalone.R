setwd("/Users/liwenyu/Documents/206/abalone/")
data=read.table('abalone.txt',sep=',')
names(data)=c("sex","length","diameter","height","whole_weight","shucked_weight","viscera_weight","shell_weight","rings")
head(data)
###missing values###
which(is.na(ab))

# check distributions of y, log transform y
par(mfrow=c(1,1))
hist(data$rings, main = "Distribution of rings", xlab="rings")
par(mfrow=c(2,2))
par(mar=rep(2,4))
boxcox(rings~.,data=data)
hist(1/(data$rings), main = "Distribution of 1/rings", xlab="rings")
hist(sqrt(data$rings), main = "Distribution of sqrt(rings)", xlab="rings")
hist(log(data$rings), main = "Distribution of log(rings)", xlab="rings")
data$rings = log(data$rings)

# check distribution of x, standardization
par(mfrow=c(3,3))
sapply(3:9, function(i) hist(data[,i], main=paste0("Distribution of ", names(data)[i])))
par(mfrow=c(1,1))
pie(table(data$sex), main="Distribution of sex")
for(i in c(1, 3:9)){
  data[,i] = scale(data[,i])*(1/sqrt(nrow(data)-1))
}
rm(i)

# check relationship between x and y
par(mfrow=c(1,1))
plot(data[,-2])
par(mar=rep(3,4))
boxplot(rings~sex, data=data, main="rings against sex")

# VIF
rxx=cor(data[,c(1,3:9)])[2:8,2:8]
rxy=cor(data[,c(1,3:9)])[2:8,1]
rxy=as.matrix(rxy)
VIF=diag(solve(rxx))
R=1-(1/VIF)
VIF; R
rm(VIF, R, rxx, rxy)

# split into train and test
set.seed(1)
index = sample(x = 1:nrow(data), size = round(nrow(data)/2))
train = data[index, ]
test = data[-index, ]
rm(index)

# first order model on train
fit.simple = lm(rings~., data=train)
fit.simple.mse = anova(fit.simple)['Residuals', 3]
library(leaps)
sub_set=regsubsets(rings~.,data=train,nbest=1,nvmax=10,method="exhaustive")
sum_sub=summary(sub_set)
n=nrow(train)
## number of coefficients in each model: p
p.m=as.integer(as.numeric(rownames(sum_sub$which))+1)
sse=sum_sub$rss
aic=n*log(sse/n)+2*p.m
bic=n*log(sse/n)+log(n)*p.m
res_sub=cbind(sum_sub$which,sse,sum_sub$rsq,sum_sub$adjr2,sum_sub$cp,
              aic, bic)
fit0=lm(rings~1,data=train) ##fit the model with only intercept
sse1=sum(fit0$residuals^2)
p=1
c1=sse1/fit.simple.mse-(n-2*p)
aic1=n*log(sse1/n)+2*p
bic1=n*log(sse1/n)+log(n)*p
none=c(1,rep(0,9),sse1,0,0,c1,bic1,aic1)
res_sub=rbind(none,res_sub) ##combine the results with other models
colnames(res_sub)=c(colnames(sum_sub$which),"sse", "R^2", "R^2_a", "Cp","aic", "bic")
res_sub
fit.final1 = lm(rings~., data=train)
summary(fit.final1)


# first order + interaction
fit.interaction = lm(rings~.^2, data=train)
fit0=lm(rings~1,data=train)
fit.final2=stepAIC(fit0,scope=list(upper=fit.interaction, lower=~1), direction="both", k=2)
summary(fit.final2)

# first order + interaction + second order
fit.full = lm(rings~.^2+I(length^2)+I(diameter^2)+I(height^2)+I(wholeWeight^2)+I(shuckedWeight^2)+I(visceraWeight^2)+I(shellWeight),data=train)
fit0=lm(rings~1,data=train)
fit.final3=stepAIC(fit0,scope=list(upper=fit.full, lower=~1), direction="both", k=2)
summary(fit.final3)

# internal validation
mse3= anova(fit.full)["Residuals",3] 

sse.fs1=anova(fit.final1)["Residuals",2] 
sse.fs2=anova(fit.final2)["Residuals",2]
sse.fs3=anova(fit.final3)["Residuals",2] 
sse.fs1; sse.fs2; sse.fs3
mse.fs1=anova(fit.final1)["Residuals",3] 
mse.fs2=anova(fit.final2)["Residuals",3] 
mse.fs3=anova(fit.final3)["Residuals",3] 
mse.fs1; mse.fs2; mse.fs3

p.fs1=length(fit.final1$coefficients) 
p.fs2=length(fit.final2$coefficients) 
p.fs3=length(fit.final3$coefficients) 
p.fs1; p.fs2; p.fs3


cp.fs1=sse.fs1/mse3-(n-2*p.fs1)
cp.fs2=sse.fs2/mse3-(n-2*p.fs2)
cp.fs3=sse.fs3/mse3-(n-2*p.fs3)
cp.fs1; cp.fs2; cp.fs3

press.fs1=sum(fit.final1$residuals^2/(1-influence(fit.final1)$hat)^2)
press.fs2=sum(fit.final2$residuals^2/(1-influence(fit.final2)$hat)^2)
press.fs3=sum(fit.final3$residuals^2/(1-influence(fit.final3)$hat)^2)
press.fs1; press.fs2; press.fs3

fit.final1; fit.final2; fit.final3

# external validate
pred.test1=predict.lm(fit.final1, test[,-1])
mspe.test1=mean((pred.test1-test[,1])^2)
pred.test2=predict.lm(fit.final2, test[,-1])
mspe.test2=mean((pred.test2-test[,1])^2)
pred.test3=predict.lm(fit.final3, test[,-1])
mspe.test3=mean((pred.test3-test[,1])^2)
mspe.test1; mspe.test2; mspe.test3

# use the third model as final model
fit.final=lm(fit.final3, data=data)
summary(fit.final)
anova(fit.final)

# model diagnostic
par(mfrow=c(1,2))
plot(fit.final, which=1)
plot(fit.final, which=2)

## check outliers in Y
res=residuals(fit.final)# residuals of the final model
n = nrow(data)
p = length(fit.final$coefficients)
h1 = influence(fit.final)$hat
d.res.std=studres(fit.final) #studentized deleted residuals
max(abs(d.res.std))
sort(abs(d.res.std),decreasing=T)
qt(1-0.1/(2*n),n-p-1) # bonferronis thresh hold
idx.Y = as.vector(which(abs(d.res.std)>=qt(1-0.1/(2*n),n-p-1)))
idx.Y ## outliers

## check outliers in X
idx.X = as.vector(which(h1>(2*p/n)))
idx.X ## two outliers
plot(h1,res,xlab="leverage",ylab="residuals")
idx.X ## outliers
length(idx.X)
## who is influencial?
plot(fit.final, which=4)
res = fit.final$residuals
mse = anova(fit.final)["Residuals", 3]
cook.d = res^2*h1/(p*mse*(1-h1)^2)
cook.d[which(h1>(2*p/n))]
head(round(sort(pf(cook.d[which(h1>(2*p/n))], p, n-p), decreasing = TRUE),3))
# only 1 influencial, 2052
fit.final2=lm(fit.final, data=data[-2052,])
f1=fitted(fit.final)
f2=fitted(fit.final2)
f1=f1[-2052]
f=f1-f2
summary(f/f1*100)
