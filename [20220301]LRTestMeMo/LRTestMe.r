# LRTestMe: Test the Mediation Effect By Likelihood Ratio Test                   #
#                                                                                #
# Author: Kazuha_Yae                                                             #
#                                                                                #
# This is Supplementary materials                                                #
#                                                                                #
# R Version: ≥ 3.4                                                               #
# Depends: lavaan                                                                #
#================================================================================#
# 0. import packages
require(mvtnorm)
require(boot)
require(lavaan)
require(sn)
require(moments)

# Section I. Latent variable model
# 1. Parameter Setting
# 参数设定部分参考了方杰和张敏强《参数和非参数方法的简单中介效应分析比较》一文，有改动（改动了样本量）
# 即样本量n设定为50, 100, 200, 500, 1000，6组参数设定为a = b = 0; a = sqrt(0.15), b = 0; a = 0, b = sqrt(0.35);
#     a = b = sqrt(0.02); a = b = sqrt(0.15); a = b = sqrt(0.35)
# 其余无关变量设定：intercept = 0; c' = 0; loadings = 0.7; indicator = 3 for each latent variable
# 1.1 设定样本量参数
N = c(50, 100, 200, 500, 1000)

# 1.2 根据参数生成协方差矩阵
# 本部分的基础功能：根据回归系数逆推协方差矩阵的部分
covMatGenerator = function(a, b, c){
	d   = c + a*b
	e   = a*c + b
	MAT = matrix(c(1, a, d, a, 1, e, d, e, 1), 3, 3)
	rownames(MAT) = colnames(MAT) = c('X', 'M', 'Y')
	return(MAT)
}
# 根据loadings和indicator指定载荷矩阵
loadingMatrix = matrix(0, nrow = 9, ncol = 3)
for(i in 1:9){loadingMatrix[i, floor((i-1)/3)+1] = 0.7}
# 根据不同组合的参数生成协方差矩阵
simCovMatrix = list(
a1 = loadingMatrix %*% covMatGenerator(a = 0         , b = 0         , c = 0) %*% t(loadingMatrix), 
a2 = loadingMatrix %*% covMatGenerator(a = sqrt(0.15), b = 0         , c = 0) %*% t(loadingMatrix), 
a3 = loadingMatrix %*% covMatGenerator(a = 0         , b = sqrt(0.35), c = 0) %*% t(loadingMatrix), 
b1 = loadingMatrix %*% covMatGenerator(a = sqrt(0.02), b = sqrt(0.02), c = 0) %*% t(loadingMatrix), 
b2 = loadingMatrix %*% covMatGenerator(a = sqrt(0.15), b = sqrt(0.15), c = 0) %*% t(loadingMatrix), 
b3 = loadingMatrix %*% covMatGenerator(a = sqrt(0.35), b = sqrt(0.35), c = 0) %*% t(loadingMatrix))
# 填入变量名，规范化对角线
for(i in 1:6){colnames(simCovMatrix[[i]]) = paste0(rep(c('X', 'M', 'Y'), each = 3), rep(1:3, 3)); diag(simCovMatrix[[i]]) = rep(1, 9)}

# 1.3 指定分析用模型
anaMod  = '
X =~ X1 + X2 + X3
M =~ M1 + M2 + M3
Y =~ Y1 + Y2 + Y3
Y ~  X + b * M
M ~  a * X
ab:= a*b'
# 特别地，指定LR检验用模型
anaModA = '
X =~ X1 + X2 + X3
M =~ M1 + M2 + M3
Y =~ Y1 + Y2 + Y3
Y ~  X + M'
anaModB = '
X =~ X1 + X2 + X3
M =~ M1 + M2 + M3
Y =~ Y1 + Y2 + Y3
Y ~  X
M ~  X'
anaModL = '
X =~ X1 + X2 + X3
M =~ M1 + M2 + M3
Y =~ Y1 + Y2 + Y3
Y ~  X + b * M
M ~  a * X
ab== 0'

# 2. 编写反馈结果的函数
# 特别需要注意的是，由于后续需要分析参数覆盖率，故这里直接给出了95% CI
# 2.1 Wald Method
# Wald检验（即Delta法），lavaan自带的默认检验方式
waldMethod = function(data, alpha = 0.95){
	fit  = sem(anaMod, data)
	res1 = parameterEstimates(fit, level = 0.95 )		# 第一个结果用来计算参数覆盖率
	res2 = parameterEstimates(fit, level = alpha)		# 第二个结果用来判定是否显著
	a    = res1[which(res1$label ==  'a'), c('ci.lower', 'ci.upper')]
	b    = res1[which(res1$label ==  'b'), c('ci.lower', 'ci.upper')]
	ab   = res1[which(res1$label == 'ab'), c('ci.lower', 'ci.upper')]
	test = (prod(res2[which(res2$label == 'ab'), c('ci.lower', 'ci.upper')])) > 0
	return(list(a, b, ab, test))
}
# 2.2 MC Method
# MC method取自Preacher和Selig（2012）的研究，该研究的思路是利用一阶泰勒展开近似的标准误（即Wald检验的标准误）不能很好反应乘积项的实际标准误/置信区间
#     故采用模拟的乘积分布替代近似得到的乘积分布
# 为了节约测试时间，这里仅抽样了1000个，可能对结果精度产生一定影响，但通常情况下这种影响不会很大
mcMethod = function(data, alpha = 0.95, seed = NA, REP = 1000L){
	if(!is.na(seed)) set.seed(seed)						# 设定随机数种子
	fit   = sem(anaMod, data)
	Sigma = vcov(fit)[c('a', 'b'), c('a', 'b')]
	Mu    = coef(fit)[c('a', 'b')]
	rnab  = rmvnorm(REP, mean = Mu, sigma = Sigma)		# 随机模拟	
	a     = quantile(rnab[, 1], c(0.025, 0.975))		# 参数覆盖
	b     = quantile(rnab[, 1], c(0.025, 0.975))
	rab   = rnab[, 1] * rnab[, 2]
	ab    = quantile(rab, c(0.025, 0.975))
	test  = (prod(quantile(rab, c(0.5-alpha/2, 0.5+alpha/2)))) > 0
	return(list(a, b, ab, test))
}
# 2.3 Bootstrap Method
# 没啥好说的，这里直接调用了lavaan集成的部分，设定重抽样次数为1000
bootMethod = function(data, alpha = 0.95, seed = NA, REP = 1000L, ci = 'perc'){
	if(!is.na(seed)) set.seed(seed)						# 设定随机数种子
	fit  = sem(anaMod, data, se = 'bootstrap', bootstrap = REP, parallel = 'snow', ncpus = 2)
	res1 = parameterEstimates(fit, level =  0.95, boot.ci.type = ci)
	res2 = parameterEstimates(fit, level = alpha, boot.ci.type = ci)
	a    = res1[which(res1$label ==  'a'), c('ci.lower', 'ci.upper')]
	b    = res1[which(res1$label ==  'b'), c('ci.lower', 'ci.upper')]
	ab   = res1[which(res1$label == 'ab'), c('ci.lower', 'ci.upper')]
	test = (prod(res2[which(res2$label == 'ab'), c('ci.lower', 'ci.upper')])) > 0	
	return(list(a, b, ab, test))
}
# 2.4 Likelihood Ratio
lrMethod = function(data, alpha = 0.95){
	fit  = sem(anaMod,  data)
	fita = sem(anaModA, data)
	fitb = sem(anaModB, data)
	res1 = parameterEstimates(fit, level =  0.95, boot.ci.type = ci)
	chia = fitmeasures(fita, 'chisq')
	chib = fitmeasures(fitb, 'chisq')
	a    = res1[which(res1$label ==  'a'), c('ci.lower', 'ci.upper')]
	b    = res1[which(res1$label ==  'b'), c('ci.lower', 'ci.upper')]
	ab   = res1[which(res1$label == 'ab'), c('ci.lower', 'ci.upper')]
	test = min(chia, chib) > qchisq(alpha, 1)
	return(list(a, b, ab, test))
}
# 3. 生成样本
simudataGenerator = function(n, COV, REP = 1000L, seed = NA){
	if(!is.na(seed)) set.seed(seed)
	dataSet = list()
	for(i in 1:REP){
		dataSet[[i]] = rmvnorm(n, mean = rep(0, 9), sigma = COV)
		colnames(dataSet[[i]]) = colnames(COV)
	}
	return(dataSet)
}
# 计算样本相关矩阵的参数估计值和参数覆盖率
parCorCov = function(dataSET, COV){
	R = length(dataSET)
	corEst = matrix(0, nrow = R, ncol = 36)
	corCov = matrix(0, nrow = R, ncol = 36)
	for(i in 1:R){
		temp = cor(dataSET[[i]])
		corEst[i, ] = temp[lower.tri(temp)]
		l = 1
		for(j in 1:8){
			for(k in (j+1):9){
				temp = as.numeric(cor.test(dataSET[[i]][, j], dataSET[[i]][, k])$conf)
				if(COV[k, j] < temp[1]){
					corEst[i, l] = -1
				}else if(COV[k, j] > temp[2]){
					corEst[i, l] = 1
				}
				l = l + 1
			}
		}
	} 
	corEst1 = apply(corEst, 2, mean)
	corEst2 = apply(corEst, 2, sd)
	corCov1 = matrix(0, nrow = 3, ncol = 36)
	for(i in 1:36){
		corCov1[1, i] = sum(corCov[, i] == -1)
		corCov1[2, i] = sum(corCov[, i] ==  0)
		corCov1[3, i] = sum(corCov[, i] ==  1)
	}
	return(rbind(corEst1, corEst2, corCov1))
}
R = length(dataSET)



# Test Section
dataSET = simudataGenerator(N[4], simCovMatrix$b1, seed = 12345678)
# 
parCorCov(dataSET, simCovMatrix$b1)


waldMethod(dataSET[[1000]])
mcMethod(dataSET[[1000]], seed = 87654321)
bootMethod(dataSET[[1000]], seed = 87654321)
lrMethod(dataSET[[1000]])














# Section II. Only explicit variables





# model <- ' 
#   # latent variable definitions
#      ind60 =~ x1 + x2 + x3
#      dem60 =~ y1 + a*y2 + b*y3 + c*y4
#      dem65 =~ y5 + a*y6 + b*y7 + c*y8
# 
#   # regressions
#     dem60 ~ ind60
#     dem65 ~ ind60 + dem60
# 
#   # residual correlations
#     y1 ~~ y5
#     y2 ~~ y4 + y6
#     y3 ~~ y7
#     y4 ~~ y8
#     y6 ~~ y8
# 	
# 	abc := a*b*c
# '
# 
# fit <- sem(model, data = PoliticalDemocracy)
# summary(fit, fit.measures = TRUE)
# 
# fit1 <- sem(model, data = PoliticalDemocracy, se = 'bootstrap')
# 
# fit2 <- sem(model, data = PoliticalDemocracy, se = 'bootstrap', test = 'boot')
# 
# parameterEstimates(fit2, level = 0.95, boot.ci.type = 'perc')
# parameterEstimates(fit2, level = 0.95, boot.ci.type = 'bca.simple')



#Import packages
if(!require(DoE.base)){install.packages('DoE.base')}
if(!require(boot)    ){install.packages('boot'    )}
if(!require(lavaan)  ){install.packages('lavaan'  )}
if(!require(mvtnorm) ){install.packages('mvtnorm' )}
#Section 0 Public Var
semmod = 'Y ~ b * M + c * X \n M ~ a * X \n ab := a * b'
semmo1 = 'Y ~ b * M + c * X \n M ~ a * X \n a == 0'
semmo2 = 'Y ~ b * M + c * X \n M ~ a * X \n b == 0'
lrmod1 = 'Y ~ M + X'
lrmod2 = 'M ~ X'
A = c(0, 0.1, -0.1, 0.3, -0.3, 0.5, -0.5)
B = c(0, 0.1, -0.1, 0.3, -0.3, 0.5, -0.5)
C = c(0, 0.1, -0.1, 0.3, -0.3, 0.5, -0.5)
N = c( 50, 100, 150, 250, 500, 750, 1000)
OA= oa.design(L49.7.8, randomize = FALSE)
seed = 72143767
#Section I Basic Method
#1.1 Joint Significance Method
JS_method = function(data, p = 0.05){
	fit1 = lm(lrmod1, data)
	fit2 = lm(lrmod2, data)
	pa   = 2 * pt(-abs(coef(fit1)[2] / sqrt(vcov(fit1)[2, 2])), df.residual(fit1))
	pb   = 2 * pt(-abs(coef(fit2)[2] / sqrt(vcov(fit2)[2, 2])), df.residual(fit2))
	return(max(pa, pb))
}
#1.2 Sobel Method (Wen et al., 2014)
Sobel_method = function(data, p = 0.05){
	fit1 = lm(lrmod1, data)
	fit2 = lm(lrmod2, data)
	a    = coef(fit1)[2]
	b    = coef(fit2)[2]
	va   = vcov(fit1)[2, 2]
	vb   = vcov(fit2)[2, 2]
	s    = sqrt(a*a*vb + b*b*va)
	pab  = pnorm(-abs(a*b/s))
	return(pab)
}
#1.3 Delta(Wald) Method
Wald_method = function(data, p = 0.05){
	fit  = sem(semmod, data)
	ab   = coef(fit)['a'] * coef(fit)['b']
	sab  = sqrt(sum(vcov(fit)[c('a', 'b'), c('a', 'b')]))
	pab  = pnorm(-abs(ab/sab))
	return(pab)
}
#1.4 Likelihood Ratio Method
LR_method = function(data, p = 0.05){
	fit  = sem(semmod, data)
	fit1 = sem(semmo1, data)
	fit2 = sem(semmo2, data)
	chi1 = fitmeasures(fit1, 'chisq')
	chi2 = fitmeasures(fit2, 'chisq')
	pab  = 1 - pchisq(min(chi1, chi2), df = 1)
	return(pab)
}
#1.5 Bootstrap #1 nonparametric / LReg
BootNparLReg_generate = function(data, p = 0.05){
	foo  = function(data, indices){
		D    = data[indices, ]
		fit1 = lm(lrmod1, D)
		fit2 = lm(lrmod2, D)
		return(coef(fit1)[2] * coef(fit2)[2])
	}
	set.seed(seed)
	res  = boot(data, statistic = foo, R = 5000L)
	return(res)
}
#1.6 Bootstrap #2 parametric / LReg
BootParLReg_perc = function(data, p = 0.05){
	fit1 = lm(lrmod1, data)
	fit2 = lm(lrmod2, data)
	a    = coef(fit1)[2]
	b    = coef(fit2)[2]
	sa   = sqrt(vcov(fit1)[2, 2])
	sb   = sqrt(vcov(fit2)[2, 2])
	set.seed(seed)
	rna  = rnorm(5000L, a, sa)
	rnb  = rnorm(5000L, b, sb)
	rab  = rna * rnb
	RAB  = ecdf(rab)
	return(1 - 2 * abs(RAB(0) - 0.5))
}
#1.7 Bootstrap #3 nonparametric / CB-SEM
BootNparSEM_generate = function(data, p = 0.05){
	foo  = function(data, indices){
		D    = data[indices, ]
		fit  = sem(semmod, data)
		ab   = coef(fit)['a'] * coef(fit)['b']
		return(ab)
	}
	set.seed(seed)
	res  = boot(data, statistic = foo, R = 5000L)
	return(res)
}
#1.8 Bootstrap #4 parametric / CB-SEM
BootParSEM_perc = function(data, p = 0.05){
	fit  = sem(semmod, data)
	Sigma= vcov(fit)[c('a', 'b'), c('a', 'b')]
	Mu   = coef(fit)[c('a', 'b')]
	set.seed(seed)
	rnab = rmvnorm(5000L, mean = Mu, sigma = Sigma)
	rab  = rnab[, 1] * rnab[, 2]
	RAB  = ecdf(rab)
	return(1 - 2 * abs(RAB(0) - 0.5))
}
#Section II Generate Data
#2.1 Generate Matrix
matrixGenerate = function(a, b, c){
	d   = c + a*b
	e   = a*c + b
	MAT = matrix(c(1, a, d, a, 1, e, d, e, 1), 3, 3)
	rownames(MAT) = colnames(MAT) = c('X', 'M', 'Y')
	return(MAT)
}
#2.2 Generate Data
dataGenerator = function(i){
	ind = as.numeric(OA[i])
	sig = matrixGenerate(A[ind[1]], B[ind[2]], C[ind[3]])
	dat = rmvnorm(N[ind[4]], sigma = sig)
	for(i in 1:3) dat[, i] = (dat[, i] - mean(dat[, i])) / sd(dat[, i])
	colnames(dat) = c('X', 'M', 'Y')
	return(data.frame(dat))
}
dataGenerator2= function(i, j, k, N){
	sig = matrixGenerate(A[i], B[j], C[k])
	dat = rmvnorm(N, sigma = sig)
	for(i in 1:3) dat[, i] = (dat[, i] - mean(dat[, i])) / sd(dat[, i])
	colnames(dat) = c('X', 'M', 'Y')
	return(data.frame(dat))
}
#Section III Simulate
#N   = 500
#Rep = 10
#p   = 0.05
#iL  = 4:4
#jL  = 4:4
#kL  = 4:4
Simufoo = function(N, Rep, p, iL, jL, kL = 1:7){
	RESULT = matrix(0, nrow = Rep, ncol = 6)
	colnames(RESULT) = c('JS', 'Sobel', 'Wald', 'LR', 'BPL', 'BPS')
	OUTCAL = numeric(15)
	if(!dir.exists('D:\\Simulate')){dir.create('D:\\Simulate')}
	if(!dir.exists('D:\\Simulate\\Logs')){dir.create('D:\\Simulate\\Logs')}
	WORKDIR= paste0('D:\\Simulate\\Logs\\N_', formatC(N, width = 4, flag = '0'), '_P_', formatC(p, digits = 3, flag = '0', format = 'f'))
	if(!dir.exists(WORKDIR)){dir.create(WORKDIR)}
	for(i in iL){
		for(j in jL){
			for(k in kL){
				#Generate the Testing Dataset
				DATASET = list()
				set.seed(seed)
				for(turn in 1:Rep){DATASET[[turn]] = dataGenerator2(i, j, k, N)}
				OUTTEMP = numeric(12)
				OUTTEMP[7]  = system.time({for(turn in 1:Rep){RESULT[turn, 1] <-        JS_method(DATASET[[turn]])}})[3]
				OUTTEMP[8]  = system.time({for(turn in 1:Rep){RESULT[turn, 2] <-     Sobel_method(DATASET[[turn]])}})[3]
				OUTTEMP[9]  = system.time({for(turn in 1:Rep){RESULT[turn, 3] <-      Wald_method(DATASET[[turn]])}})[3]
				OUTTEMP[10] = system.time({for(turn in 1:Rep){RESULT[turn, 4] <-        LR_method(DATASET[[turn]])}})[3]
				OUTTEMP[11] = system.time({for(turn in 1:Rep){RESULT[turn, 5] <- BootParLReg_perc(DATASET[[turn]])}})[3]
				OUTTEMP[12] = system.time({for(turn in 1:Rep){RESULT[turn, 6] <-  BootParSEM_perc(DATASET[[turn]])}})[3]
				OUTTEMP[1:6]= colSums(RESULT < p)
				write.csv(cbind(RESULT), paste0(WORKDIR, '\\', i, j, k, '.csv'))
				OUTCAL      = rbind(OUTCAL, c(A[i], B[j], C[k], OUTTEMP))
			}
		}
	}
	write.csv(OUTCAL, paste0('D:\\Simulate\\N_', formatC(N, width = 4, flag = '0'), '_P_', formatC(p, digits = 3, flag = '0', format = 'f'), '.csv'))
	OUTCAL
}
RR     = 1000L
i      = 1
ind    = as.numeric(OA[i])
n      = N[ind[4]]
RESOUT = c(n, Simufoo(n, Rep = RR, p = 0.05, iL = ind[1], jL = ind[2], kL = ind[3])[2, ])
for(i in 2:49){
	ind    = as.numeric(OA[i])
	n      = N[ind[4]]
	RESOUT = rbind(RESOUT, c(n, Simufoo(n, Rep = RR, p = 0.05, iL = ind[1], jL = ind[2], kL = ind[3])[2, ]))
}

#Simufoo(N = 50,   Rep = 1000, p = 0.05, iL = 1:7, jL = 1:7)
#Simufoo(N = 100,  Rep = 1000, p = 0.05, iL = 1:7, jL = 1:7)
#Simufoo(N = 400,  Rep = 1000, p = 0.05, iL = 1:7, jL = 1:7)
#Simufoo(N = 500,  Rep = 1000, p = 0.05, iL = 1:7, jL = 1:7)
#Simufoo(N = 1000, Rep = 1000, p = 0.05, iL = 1:7, jL = 1:7)

