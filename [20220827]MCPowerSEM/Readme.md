# 【MCPowerSEM: Structural Equation Modeling Power Analysis use Monte-Carlo Method】

基于Monte-Carlo法的SEM统计检验力分析

* 相关信息：

> Version: Alpha 0.0.2
> 
> R Version: ≥ 3.4
> 
> Depends: lavaan
> 
> Author: Kazuha_Yae
> 
> License: GPL ≥ 2

* 目的：基于R语言开发基于Monte-Carlo方法的SEM统计检验力分析

* 函数接口：

主函数：MCPowerSEM，该函数需要传入的参数包括genMod, genCov, anaMod, REP, N, std, seed，各参数含义如下

| 传入参数   | 说明                    |
|:------:|:---------------------:|
| genMod | 用以生成虚拟数据的理论模型         |
| genCov | 用以生成虚拟数据的协方差矩阵        |
| anaMod | 实际进行分析时所执行的模型         |
| REP    | 模拟的重复次数               |
| N      | 模拟样本量                 |
| std    | 是否需要根据“智能”还原标准化的协方差矩阵 |
| seed   | 随机数种子                 |

该函数构建一个MCPowerSEM类的实例，该实例中的各个属性记录了模型设置与模拟抽样的分析结果，该类的属性包括

| 属性             | 说明                        |
|:--------------:|:-------------------------:|
| generate.mod   | 用以生成模拟数据的模型（以lavaan模型指定）  |
| generate.cov   | 用以生成模拟数据的协方差矩阵（直接以矩阵形式指定） |
| std.reg        | 是否需要根据“智能”还原标准化的协方差矩阵     |
| analysis.mod   | 用以进行分析的模型                 |
| simulate.n     | 样本容量                      |
| repeated       | 重复次数                      |
| seed           | 起始种子                      |
| origin.mod     | 原始模型                      |
| simulate.seed  | 生成模拟数据用的种子                |
| result.fitting | 模拟得到的模型拟合                 |
| result.paraest | 模拟得到的参数估计值                |
| result.simcov  | 模拟数据的协方差矩阵                |

该类存在下述方法

| 方法       | 说明                                         |
|:--------:|:------------------------------------------:|
| parCover | 考察实际生成样本的协方差矩阵与生成数据的协方差矩阵的参数覆盖情况           |
| estCover | 考察实际生成样本的估计值与生成数据的模型之间的参数覆盖情况              |
| fitCover | 考察实际生成样本的拟合情况与生成数据的模型之间的参数覆盖情况             |
| estTest  | 给出指定决策标准下（默认是α为0.05和0.01）模型参数的双侧检验的统计检验力分析 |
| fitTest  | 给出指定决策标准下，模型拟合指数的单侧统计检验力分析                 |
| summary  | 打印综合结果                                     |

* 使用范例：

```
# 示例程序
# 模型来源于Sarah提出的问题
genMod = '
Y  ~ 0.1 * C + - 0.2 * X + 0.2 * M1 + 0.2 * M2
M1 ~ 0.1 * C + - 0.2 * X 
M2 ~ 0.1 * C + - 0.2 * X 
X  ~ 0.1 * C
'

anaMod = '
Y  ~ C + X + b1 * M1 + b2 * M2
M1 ~ C + a1 * X 
M2 ~ C + a2 * X 
X  ~ C
ab1 := a1 * b1
ab2 := a2 * b2
'

ana.test = MCPowerSEM(genMod = genMod, genCov = NULL, anaMod = anaMod, REP = 1000L, N = 401, std = TRUE, seed = 12345678)

summary(ana.test)

# 使用蒙特卡洛法对SEM模型的统计检验力进行分析
# 模型参数：
# 用以生成虚拟数据的协方差矩阵：
#             Y       M1       M2       X     C
# [1,]  1.00000  0.26496  0.26496 -0.2688 0.112
# [2,]  0.26496  1.00000  0.04640 -0.1920 0.080
# [3,]  0.26496  0.04640  1.00000 -0.1920 0.080
# [4,] -0.26880 -0.19200 -0.19200  1.0000 0.100
# [5,]  0.11200  0.08000  0.08000  0.1000 1.000
# 用以进行分析的模型：
# [1] "\nY  ~ C + X + b1 * M1 + b2 * M2\nM1 ~ C + a1 * X \nM2 ~ C + a2 * X \nX  ~ C\nab1 := a1 * b1\nab2 := a2 * b2\n"
# 模拟样本量： 401 
# 重复次数： 1000 
# 初始种子： 12345678 
# 
# 模型参数检验的统计检验力：
#                   est      se  z.value p.value power.at.0.05 power.at.0.01
# Y ~ C         0.09834 0.04611  2.13258 0.03296         0.577         0.338
# Y ~ X        -0.19826 0.04753 -4.17097 0.00003         0.989         0.948
# Y ~ M1        0.20932 0.04651  4.50044 0.00001         0.995         0.978
# Y ~ M2        0.20932 0.04651  4.50044 0.00001         0.990         0.961
# M1 ~ C        0.10020 0.04906  2.04236 0.04112         0.564         0.319
# M1 ~ X       -0.20202 0.04906 -4.11767 0.00004         0.987         0.936
# M2 ~ C        0.10020 0.04906  2.04236 0.04112         0.550         0.306
# M2 ~ X       -0.20202 0.04906 -4.11767 0.00004         0.977         0.927
# X ~ C         0.10000 0.04975  2.01008 0.04442         0.530         0.285
# Y ~~ Y        0.82477 0.05832 14.14214 0.00000         1.000         1.000
# M1 ~~ M1      0.95320 0.06740 14.14214 0.00000         1.000         1.000
# M2 ~~ M2      0.95320 0.06740 14.14214 0.00000         1.000         1.000
# X ~~ X        0.99000 0.07000 14.14214 0.00000         1.000         1.000
# C ~~ C        1.00000 0.00000      Inf 0.00000         1.000         1.000
# ab1 := a1*b1 -0.04229 0.01392 -3.03796 0.00238         0.957         0.744
# ab2 := a2*b2 -0.04229 0.01392 -3.03796 0.00238         0.942         0.726
# 
# 模型拟合检验的统计检验力：
#                crit power.Left power.Right
# rmsea          0.08      0.938       0.062
# rmsea          0.10      0.977       0.023
# rmsea          0.05      0.838       0.162
# rmsea.ci.upper 0.08      0.176       0.824
# cfi            0.90      0.000       1.000
# mfi            0.85      0.000       1.000
```

* 一些特别的备注：
  
  * **本项目为COSN Hackathon 2022的部分成果，COSN Hackathon的主项目地址见[这里](https://github.com/OpenSci-CN/SummerHackathon2022)**
  
  * **特别注意**：测试并不完全，实际使用过程中可能出现各类Bug，Bug反馈请联系sulfonamides@163.com
  
  * 目前只支持单层单组模型相关的统计检验力分析，分层/分组结构方程模型相关的功能尚需进一步开发
  
  * 理论上支持涉及调节项的模型，但不支持潜调节结构方程模型（lavaan包原生不支持）
    
    * 特别说明：使用涉及调节项的模型时，请将乘积项视为与其他变量地位对等的变量，指定蕴含所有变量的协方差矩阵或生成模型
      
      ```
      # 举一个简单的栗子，对于模型 Y ~ X + W + XW 来说，应当传入的是X, W, XW三个变量的协方差矩阵，而不仅仅是是X与W的协方差矩阵
      # 例如：
      # 以pwr2ppl包中的modmed14为例：
      # modmed14(rxw=.2, rxm=.25, rxww=.2, rxwy=-.2, rxxw=.35, rxy=.3, rwm=.4, rwy=.35, rmy=.3, n=200, rep=1000, alpha=.05)等价于下述代码：
      rxw=.2; rxm=.25; rxww=.2; rxwy=-.2; rxxw=.35; rxy=.3; rwm=.4; rwy=.35; rmy=.3; rxwm=0
      genCov = matrix(c(1.0, rxw, rxm, rxxw, rxy,
                      rxw, 1.0, rwm, rxww, rwy,
                      rxm, rwm, 1.0, rxwm, rmy,
                      rxxw, rxww, rxwm, 1.0, rxwy,
                      rxy, rwy, rmy, rxwy, 1.0), 5, 5)
      colnames(genCov) = c('x', 'w', 'm', 'xw', 'y')
      anaMod = '
      m ~ a * x
      y ~ x + m + w + b * xw
      ab:= a*b
      '
      ana.test = MCPowerSEM(genMod = NULL, genCov = genCov, anaMod = anaMod, REP = 1000L, N = 201, std = FALSE, seed = 12345678)
      
      summary(ana.test)
      
      #使用蒙特卡洛法对SEM模型的统计检验力进行分析
      #模型参数：
      #用以生成虚拟数据的协方差矩阵：
      #             x         w          m          xw          y
      #[1,] 1.0000000 0.1971191 0.25000000  0.33177706  0.2870921
      #[2,] 0.1971191 0.9985180 0.39392569  0.20735231  0.3170311
      #[3,] 0.2500000 0.3939257 1.00000000  0.01214862  0.3030372
      #[4,] 0.3317771 0.2073523 0.01214862  1.00000000 -0.1757028
      #[5,] 0.2870921 0.3170311 0.30303716 -0.17570275  1.0000000
      #用以进行分析的模型：
      #[1] "\n  m ~ a * x\n  y ~ x + m + w + b * xw\n  ab:= a*b\n  "
      #模拟样本量： 201 
      #重复次数： 1000 
      #初始种子： 12345678 
      #
      #模型参数检验的统计检验力：
      #               est      se  z.value p.value power.at.0.05 power.at.0.01
      #m ~ x      0.25000 0.06847  3.65148 0.00026         0.947         0.856
      #y ~ x      0.31513 0.06632  4.75154 0.00000         0.996         0.981
      #y ~ m      0.11844 0.06220  1.90413 0.05689         0.464         0.268
      #y ~ w      0.27908 0.06222  4.48534 0.00001         0.989         0.962
      #y ~ xw    -0.33956 0.06461 -5.25541 0.00000         0.999         0.990
      #m ~~ m     0.93750 0.09375 10.00000 0.00000         1.000         1.000
      #y ~~ y     0.72550 0.07255 10.00000 0.00000         1.000         1.000
      #x ~~ x     1.00000 0.00000      Inf 0.00000         1.000         1.000
      #x ~~ w     0.19712 0.00000      Inf 0.00000         1.000         1.000
      #x ~~ xw    0.33178 0.00000      Inf 0.00000         1.000         1.000
      #w ~~ w     0.99852 0.00000      Inf 0.00000         1.000         1.000
      #w ~~ xw    0.20735 0.00000      Inf 0.00000         1.000         1.000
      #xw ~~ xw   1.00000 0.00000      Inf 0.00000         1.000         1.000
      #ab := a*b -0.08489 0.02831 -2.99871 0.00271         0.919         0.707
      #
      #模型拟合检验的统计检验力：
      #               crit power.Left power.Right
      #rmsea          0.08      0.000       1.000
      #rmsea          0.10      0.000       1.000
      #rmsea          0.05      0.000       1.000
      #rmsea.ci.upper 0.08      0.000       1.000
      #cfi            0.90      0.996       0.004
      #mfi            0.85      0.008       0.992
      
      # 在这里，我们得到的结果是0.919的检验力，需要提醒注意的是，modmed采用了联合显著法，而我们默认使用了Wald检验，如果希望在这里换用联合显著法的话，则需要进行如下操作
      
      result.m = matrix(0, nrow = 1000, ncol = 2)
      result.sd= matrix(0, nrow = 1000, ncol = 2)
      for(i in 1:1000){
        result.m[i, ] = ana.test$result.paraest[[i]]$est[c(1, 5)]
        result.sd[i,] = ana.test$result.paraest[[i]]$se[c(1, 5)]
      }
      result.z = result.m / result.sd
      result   = apply(result.z, c(1, 2), function(x) ifelse(abs(x) > qnorm(0.975), 1, 0))
      (power   = mean(result[, 1] * result[, 2]))
      
      #前述操作实际上就是从ana.test的属性中提取出需要的数据进一步分析，该结果为0.946，与modmed14得到的结果（0.944）相似
      #另外在翻结果的时候里面好好像有一个说明错误，xw疑似应为mw
      ```
  
  * 请**严格注意**生成模型中的参数指定！该程序没有内置生成模型的检查功能，错误地指定生成模型（通常是遗漏关键参数或指定了冲突的参数）可能导致模型自动生成的理论协方差矩阵与您希望得到的协方差矩阵并不相同
  
  * ~~若使用该方法分析模型拟合指数（如RMSEA、CFI等）时，由于抽样误差的缘故，在使用根据前人研究或Meta分析生成的协方差矩阵或根据数据进行事后敏感性分析时，该方法会***严重***错估统计检验力，解决该问题的方法正在集成中，如果研究者希望*马上*使用该功能，可以参考Cheng和Wu（2017）的研究~~ 该方法已内嵌
  
  * 目前尚不支持指定生成模型的情况下对具有**任意**RMSEA特征的总体协方差矩阵进行模拟的功能，但相关功能已在开发中，如果研究者有紧迫需求，可以参考Cudeck与Browne（1992）的研究
  
  * 目前在对结构方程模型参数或参数函数的检验中仅支持Wald检验（Delta法），将来可能支持Lagrange乘子检验和似然比检验【但二者的主执行程序不同，可能需要重构类，但作者过懒\_(:3」∠)_，不好说啥时候能出来】，并且应该不会支持自举法【耗时过长，怎么想都不是一个划算的结果，并且作者不会写多线程（后者是重点）】；但一个好消息是，通常情况下Wald检验具有所有该类检验中最低的检验力，所以更换方法只会提升检验的有效性
  
  * 暂时不支持对估计时所选参数的更改（如拟合方法、似然函数类型等），原因与为何不支持Lagrange乘子检验和似然比检验一样【说白了还是懒】
  
  * 结果整理时默认使用了对称区间，没使用最高密度区间，该功能待完善
  
  * 在涉及潜变量模型时，事实上不是所有的生成数据都能良好拟合的，总会有可能出现一些“病态”的拟合（如迭代不收敛、残差方差为负值等），这些拟合理应从结果中排除，但这类结果的比例通常不高，对结果的影响通常不大，目前这部分内容也正在更新中（实际上使用try catch在一定程度上就可以解决此问题）
  
  * 只完成了简单的关于模型参数检验的部分，并且暂时只支持α = 0.05和0.01的双侧检验（power95和power99），但结果都在属性里了，研究者有需要可以自己定制功能
    
    * 例如：等效性检验。对于评价单个模型拟合的等效性检验，需要研究者参考王阳等人（2020）的研究，首先计算出等效性检验的chisq临界值，再传入计算（不过此类指标的统计检验力分析通常有更容易的参数方法，不一定非得使用模拟的方式）
      
      ```
      #例：以等效RMSEA_0 = 0.08的水平进行等效性检验，模型的设定如下
      genMod = '
      POS ~~ 0.23 * DES
      POS ~~ 0.40 * ANG
      DES ~~ 0.78 * ANG
      POS =~ 0.87 * X1 + 0.90 * X2 + 0.75 * X3 + 0.55 * X4
      DES =~ 0.59 * X5 + 0.72 * X6 + 0.58 * X7 + 0.71 * X8 + 0.59 * X9
      ANG =~ 0.73 * X10 + 0.78 * X11 + 0.57 * X12
      X1 ~~ 0.24 * X1
      X2 ~~ 0.19 * X2
      X3 ~~ 0.44 * X3
      X4 ~~ 0.70 * X4
      X5 ~~ 0.65 * X5
      X6 ~~ 0.47 * X6
      X7 ~~ 0.54 * X7
      X8 ~~ 0.49 * X8
      X9 ~~ 0.52 * X9
      X10 ~~ 0.45 * X10
      X11 ~~ 0.39 * X11
      X12 ~~ 0.55 * X12
      '
      
      anaMod = '
      POS =~  X1 + X2 + X3 + X4
      DES =~  X5 + X6 + X7 + X8
      ANG =~  X9 + X10 + X11 + X12
      '
      
      #假定样本为101人
      ana.test = MCPowerSEM(genMod = genMod, genCov = NULL, anaMod = anaMod, REP = 1000L, N = 101, std = FALSE, seed = 12345678)
      
      #确定临界值
      df = fitmeasures(ana.test$origin.mod, 'df')
      chisq.crit = qchisq(0.01, df, ncp = 0.08 * 0.08 * df * (101-1))
      
      #进行检验力分析
      fitTest(ana.test, 'chisq', chisq.crit)
      ```
  
  * 暂时没有更多了，想起来再补\_(:3」∠)_   