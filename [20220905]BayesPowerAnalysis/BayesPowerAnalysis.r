# BayesPowerAnalysis: Power Analysis in Bayesian Data Analysis                   #
#                                                                                #
# Author: Kazuha_Yae                                                             #
# Ref: Kruschke, J. K. (2015)Doing Bayesian Data Analysis A Tutorial with R, JAGS#
#          , and Stan(2nd Ed.). Elsevier Inc.                                    #
#                                                                                #
# Program Version: Alpha 0.0.1                                                   #
# R Version: ≥ 4.1                                                               #
# Depends: rstan, hdrcde                                                         #
# License: GPL 3.0                                                               #
# Bug Report: sulfonamides@163.com                                               #
#================================================================================#
#-1. check packages
pkgs = c('rstan', 'hdrcde', 'bayestestR', 'BayesFactor')
for(i in 1:length(pkgs)){if(!(pkgs[i] %in% installed.packages()[, 'Package'])){install.packages(pkgs[i])}}

# 0. import packages
for(i in 1:length(pkgs)) library(pkgs[i])

# 1. base class MC
#--------------------------------------------------------------------------------
# 这部分给出整个模拟模块中的基础类MC
# 1.1 MC属性设定：
# 传入属性
#   generate.mod: 用以生成模拟数据的模型（以lavaan模型指定）
#   generate.cov: 用以生成模拟数据的协方差矩阵（直接以矩阵形式指定）
#   std.reg: 当以模型生成数据时，是否将输入模型的回归系数部分进行标准化处理
#   analysis.mod: 用以进行分析的模型
#   simulate.n: 样本容量
#   repeated: 重复次数
#   seed: 起始种子
# 模型内建属性
#   origin.mod: 原始模型
#   simulate.seed: 生成模拟数据用的种子
#   result.fitting: 模拟得到的模型拟合
#   result.paraest: 模拟得到的参数估计值
#   result.simcov: 模拟数据的协方差矩阵


foo = function(f, N){
	function(...){
		f(..., n = N)
	}
}
funif = foo(runif, N = 50)
set.seed(1234)
funif(max = 2, min = -1) -> test