# Project Web Access Panel Sampling
# Author: Kazuha Yae
# Date: 2022-10-27
# Environment:
#   R version  4.1.3
#   sampling   2.9
# 只有一个地方需要用到sampling这个包, 而且还有替代方式, 所以装不装都OK啦【虽然确实建议装一下
#===========================================================================================
# 模拟部分1: 模拟的是先生成有限总体, 然后从有限总体里抽出网络访问固定样本和调查样本, 最后再从网络
#            访问固定样本和调查样本里抽样, 好处是似真, 坏处是不等概的入样概率太难调了, 费劲调半天
#            结果最高入样概率才是最低入样概率的两倍

# 指定参数
Parameter = list(
# 控制变量的随机数发生器
Controls = list(C1 = function(n){runif(n, 0, 1)},            # C1 ~ U(0, 1)                                                          
                C2 = function(n){
                    out = rpois(n, 1)
                    for(i in 1:n) if(out[i] > 6){out[i] = 6}
                    out},                                    # C2 ~ Poisson(1); C2 <= 6
                C3 = function(n){exp(rnorm(n, 1.5, 1))}),    # C3 ~ lognorm(1.5, 1^2)
# 样本量, totalN, netN, controlN, expN分别代表有限总体数量, 网络访问固定样本数量, 
#     随机样本数量和网络抽样研究样本数量
# 样本生成过程：有限总体中的一部分成为了网络访问固定样本，这部分是不等概的，现实
#     有两个研究，一个是严格从总体中按概率抽取的样本（随机样本），一个是我们常用
#     的从网络样本中抽样的样本（网络抽样研究样本）
N = list(totalN = 1e6+1e4, netN = 1e5, controlN = 1e4, expN = 1000),
# 参数, 包括用来生成X2, X1, Y以及不等概部分抽样权重的参数
para = list(beta =  c(1.3, 0.5, 0.3, 2, 0.35, 0.03), # 生成Y的参数
            alpha = c(1.1, 1.5, 0.2, 0.025),         # 生成X1的参数
            gamma = c(1.7, 1, 0.3, 0.05),            # 生成X2的参数
            theta = c(0.1, 5, 0.7, 1.5))             # 抽样权重参数
)
#===========================================================================================
# 求加权协方差矩阵
WeightedCOV = function(Data, weight){
    N       = nrow(Data)                                # 求行数
    Data    = as.matrix(Data)                           # 转化数据类型
    weight  = weight / sum(weight)                      # 归一化权向量
    Var.Bar = apply(Data * weight, 2, sum)              # 求均值
    Var.Div = Data - matrix(rep(Var.Bar, each = N), N)  # 求中心矩
    W.COV   = t(Var.Div * weight) %*% Var.Div * N/(N-1) # 求协方差矩阵
    rownames(W.COV) = colnames(W.COV) = colnames(Data)  # 贴上列名
    return(W.COV)
}
#===========================================================================================
# 正式模拟
simulation = function(seed = 1234, R = 1000L, Parameter){
    # 传入三个参数: seed是用来初始化随机数发生器的种子; R是后面模拟抽样的次数; Parameter是控
    #   制参数（关于parameter的规则参看前面给出的样例）

    # 0. 准备环节
    # Ps. 没有检查输入是否合法, 不过反正是自己人用无所谓啦
    # 又Ps. 其实写成OOP的形式可能更清楚一点, 不过R自带的OOP实际上跟这玩意儿也大差不差
    # 0.1 从Control中读取所需参数
    # 0.1.1 读取随机数发生器
    C1 = Parameter$Controls$C1
    C2 = Parameter$Controls$C2
    C3 = Parameter$Controls$C3
    # 0.1.2 读取样本量
    Nt = Parameter$N$totalN
    Nn = Parameter$N$netN
    Nc = Parameter$N$controlN
    Ne = Parameter$N$expN
    # 0.1.3 读取样本生成参数
    B  = Parameter$para$beta
    A  = Parameter$para$alpha
    G  = Parameter$para$gamma
    T  = Parameter$para$theta
    # 0.1.4 初始化随机数发生器
    set.seed(seed)

    #---------------------------------------------------------------------------------------
    # 1. 生成（有限）总体
    # 1.1 首先生成控制变量
    C1.sample = C1(Nt)
    C2.sample = C2(Nt)
    C3.sample = C3(Nt)
    # 1.2 然后根据需求生成三种不同情况下的模拟数据
    # 实际上生成的过程可以借助矩阵乘法写的更优雅一点, 但我懒.jpg
    # 1.2.1 独立的情况
    X1.indep  = rnorm(Nt)
    X2.indep  = rnorm(Nt)
    # 1.2.2 X2不独立的情况
    X2.depnd  = G[1] + G[2] * C1.sample + G[3] * C2.sample + G[4] * C3.sample + rnorm(Nt)    
    X2.mean   = G[1] + G[2] * 1/2 + G[3] + G[4] * exp(2)                                     # 理论均值
    X2.sd     = sqrt(G[2] * G[2] * 1/12 + G[3] * G[3] + G[4] * G[4] * (exp(5) - exp(4)) + 1) # 理论标准差
    X2.depnd  = (X2.depnd - X2.mean) / X2.sd                                                 # 标准化
    # 1.2.3 X1不独立的情况
    X1.depnd  = A[1] + A[2] * C1.sample + A[3] * C2.sample + A[4] * C3.sample + rnorm(Nt)
    X1.mean   = A[1] + A[2] * 1/2 + A[3] + A[4] * exp(2)
    X1.sd     = sqrt(A[2] * A[2] * 1/12 + A[3] * A[3] + A[4] * A[4] * (exp(5) - exp(4)) + 1)
    X1.depnd  = (X1.depnd - X1.mean) / X1.sd
    # 1.2.4 生成三种情况下的Y
    # X1/X2都独立于控制变量的情况
    Y1 = B[1] + B[2] * X1.indep + B[3] * X2.indep + B[4] * C1.sample + B[5] * C2.sample + B[6] * C3.sample + rnorm(Nt)
    # X2不独立于控制变量
    Y2 = B[1] + B[2] * X1.indep + B[3] * X2.depnd + B[4] * C1.sample + B[5] * C2.sample + B[6] * C3.sample + rnorm(Nt)
    # X1/X2都不独立于控制变量
    Y3 = B[1] + B[2] * X1.depnd + B[3] * X2.depnd + B[4] * C1.sample + B[5] * C2.sample + B[6] * C3.sample + rnorm(Nt)
    # 1.2.5 组装数据框, 这里为了省事儿就直接装一起了
    Population = data.frame(
    C1 = C1.sample, C2 = C2.sample, C3 = C3.sample,
    X1.ind = X1.indep, X2.ind = X2.indep, X1.dep = X1.depnd, X2.dep = X2.depnd,
    Y1 = Y1, Y2 = Y2, Y3 = Y3)
    
    #---------------------------------------------------------------------------------------
    # 2. 析出简单随机样本和不等概样本（简单随机样本模拟的是大规模社会调查, 不等概样本模拟的是网络访问固定样本）
    # 2.1 从有限总体（Nt）中析出简单随机样本, 样本量Nc
    inds = sample.int(Nt, size = Nc, replace = FALSE)
    Sample.control = Population[inds, ]
    Sample.comple  = Population[-inds,]
    # 2.2 从有限总体中析出不放回不等概随机样本, 模拟网络访问固定样本集, 样本量Nn
    # 生成权重
    W    = T[1] + T[2] * Sample.comple$C1 + T[3] * Sample.comple$C2 + T[4] * Sample.comple$C3
    P    = 1 / (1 + exp(-W))
    # 计算一阶入样概率, 在入样概率相差不是很悬殊的情况下, 这实际上等价于 pik = Nn * P / sum(P)
    # pik  = sampling::inclusionprobabilities(P, Nn)
    pik = Nn * P / sum(P)
    # 不等概不放回抽样（Brewer法, 虽然精度够但抽样效率太低了所以注释掉了）
    # inds = sampling::UPbrewer(pik)
    # 泊松抽样, 不保证样本数量但保证一阶入样概率
    inds = 1 * (runif(nrow(Sample.comple)) < pik)
    Sample.net = Sample.comple[which(inds == 1), ]
    # 因为使用的是泊松抽样, 所以实际的Nn发生了变化（相差不大就是啦）, 但在这里还需要更新一下Nn
    Nn = nrow(Sample.net)

    #---------------------------------------------------------------------------------------    
    # 3. 开始循环主体
    # 3.0 准备记录循环结果
    # 需要记录的东西太多了, 弄成一个向量太麻烦, 干脆直接整个List, 后面再整理
    # 记录list的格式是: list套list, 外层的list长度为Ne + 2, 包括Ne个模拟结果（序号从1一直到Ne）
    #   之后的Ne+1是对Population的计算结果, Ne+2是对Sample.net的计算结果; 外层list的每一个元素
    #   依然是list, 除最后两个外, 每个list的长度为4, 包含四种计算方法下的各种估计值, 方法一共四
    #   种: UW, 未加权; IW, 倾向值逆加权; WC, 组加权; PO, 事后分层调整加权。每个方法list中包含8
    #   个元素, 依次为: mean, 均值; cov, 协方差; ind, 假定自变量与协变量独立; hlf, 假定部分自变
    #   量与协变量独立; dep, 假定全部自变量与协变量相关; 含C表示这个回归方程加入了协变量
    out = list()
    # 3.1 进入循环
    for(i in 1:R){
        # 3.1.1 模拟抽样
        # 抽样规则: 假定研究者目前手上有两个数据集, 其一是按照人口学严格抽样的大型社会调查数据集
        #   Sample.control; 其二则是研究者在网络上发放问卷回收的结果。值得注意的是, 倾向值分析部
        #   分实际上是在做logistic回归（训练logistic分类器）, 对于该模型而言, 如果因变量的0-1配比
        #   严重失衡, 可能导致分类器过分学习0组或1组的特征导致对实际倾向值的误估, 所以我们这里从
        #   Sample.control中二次抽取了一个样本量与模拟调查样本量相当的样本
        # (1). 模拟在网络访问固定样本中抽样
        ind.exp = sample.int(Nn, size = Ne, replace = FALSE)
        Sample.exp = Sample.net[ind.exp, ]
        # (2). 在Sample.control中二次抽取样本
        ind.ctl = sample.int(Nc, size = Ne, replace = FALSE)
        Sample.ctl = Sample.control[ind.ctl, ]
        # (3). 捏合用于计算倾向值（Propensity Score）的混合样本, 该样本前Ne行为网络样本, 后Ne行为简单随机样本
        PSCalc = cbind(rbind(Sample.exp, Sample.ctl), Y = rep(c(1, 0), c(Ne, Ne)))
        #...................................................................................
        # 3.1.2 数据处理
        # (1). 计算倾向值得分
        PSfit  = glm(Y ~ C1 + C2 + C3, family = binomial(link = 'logit'), data = PSCalc) # 拟合logistic回归
        PS     = as.vector(1 / (1 + exp(-predict(PSfit))))                               # 计算倾向值得分
        # (2). 计算权重
        # 权重的计算方法来自刘展, 金勇进. (2016). 基于倾向得分匹配与加权调整的非概率抽样统计推断方法研究
        # Sample.exp等价于文中的匹配样本（记作M）, Sample.ctl等价于文中的目标样本（记作O）, 假定两样本的基础权数均相等
        # a. 逆加权调整
        IW = 1 / PS                      # 求倒数
        iw = IW[1:Ne] / sum(IW[1:Ne])    # 逆加权调整权重
        # b. 加权组调整
        G  = Ne/10                                                  # 设置每个组内的样本量为20/需要保证Ne是10的倍数
        WC = PS[order(PS)]                                          # 将倾向值按照大小排序
        WC = rep(tapply(WC, rep(1:G, each = 20), mean), each = 20)  # 对倾向值进行组内均值代换
        WC = 1 / WC[order(order(PS))]                               # 重排序并求倒数
        wc = WC[1:Ne] / sum(WC[1:Ne])                               # 加权组调整权重
        # c. 事后分层调整（直接用了b的分组）
        PO = PSCalc$Y                                               # 取出示性变量
        PO = PO[order(PS)]                                          # 对示性变量按照PS的顺序排序
        PO = rep(tapply(PO, rep(1:G, each = 20), mean), each = 20)  # 求出每个层中的个体来自匹配样本的概率
        PO = PO / (1 - PO)                                          # 求优势比（相当于换了个方式估计倾向值）
        PO = 1 / PO[order(order(PS))]                               # 重排序并求倒数
        po = PO[1:Ne] / sum(PO[1:Ne])                               # 事后分层调整权重
        #...................................................................................       
        # 3.1.3 各种计算
        # 未加权情况（最简单的情形）
        UW.mean = apply(Sample.exp, 2, mean)
        UW.cov  = cov(Sample.exp)
        UW.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp)
        UW.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp)
        UW.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp)
        UW.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp)
        UW.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp)
        UW.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp)
        # 逆加权调整
        IW.mean = apply(Sample.exp * iw, 2, sum)
        IW.cov  = WeightedCOV(Sample.exp, iw)
        IW.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp, weights = iw)
        IW.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp, weights = iw)
        IW.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp, weights = iw)
        IW.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp, weights = iw)
        IW.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp, weights = iw)
        IW.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp, weights = iw)        
        # 加权组调整
        WC.mean = apply(Sample.exp * wc, 2, sum)
        WC.cov  = WeightedCOV(Sample.exp, wc)
        WC.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp, weights = wc)
        WC.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp, weights = wc)
        WC.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp, weights = wc)
        WC.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp, weights = wc)
        WC.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp, weights = wc)
        WC.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp, weights = wc)         
        # 事后分层调整
        PO.mean = apply(Sample.exp * po, 2, sum)
        PO.cov  = WeightedCOV(Sample.exp, po)
        PO.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp, weights = po)
        PO.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp, weights = po)
        PO.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp, weights = po)
        PO.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp, weights = po)
        PO.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp, weights = po)
        PO.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp, weights = po)        
        #...................................................................................
        # 3.1.4 组装结果
        out[[i]] = list(UW = list(
        UW.mean = UW.mean, 
        UW.cov  = UW.cov , 
        UW.ind  = coef(UW.ind ),    # UW.ind  = UW.ind 
        UW.indC = coef(UW.indC),    # UW.indC = UW.indC
        UW.hlf  = coef(UW.hlf ),    # UW.hlf  = UW.hlf 
        UW.hlfC = coef(UW.hlfC),    # UW.hlfC = UW.hlfC
        UW.dep  = coef(UW.dep ),    # UW.dep  = UW.dep 
        UW.depC = coef(UW.depC)),   # UW.depC = UW.depC
        IW = list(
        IW.mean = IW.mean,
        IW.cov  = IW.cov ,
        IW.ind  = coef(IW.ind ),    # IW.ind  = IW.ind 
        IW.indC = coef(IW.indC),    # IW.indC = IW.indC
        IW.hlf  = coef(IW.hlf ),    # IW.hlf  = IW.hlf 
        IW.hlfC = coef(IW.hlfC),    # IW.hlfC = IW.hlfC
        IW.dep  = coef(IW.dep ),    # IW.dep  = IW.dep 
        IW.depC = coef(IW.depC)),   # IW.depC = IW.depC
        WC = list(
        WC.mean = WC.mean,
        WC.cov  = WC.cov ,
        WC.ind  = coef(WC.ind ),    # WC.ind  = WC.ind 
        WC.indC = coef(WC.indC),    # WC.indC = WC.indC
        WC.hlf  = coef(WC.hlf ),    # WC.hlf  = WC.hlf 
        WC.hlfC = coef(WC.hlfC),    # WC.hlfC = WC.hlfC
        WC.dep  = coef(WC.dep ),    # WC.dep  = WC.dep 
        WC.depC = coef(WC.depC)),   # WC.depC = WC.depC
        PO = list(
        PO.mean = PO.mean,
        PO.cov  = PO.cov ,
        PO.ind  = coef(PO.ind ),    # PO.ind  = PO.ind 
        PO.indC = coef(PO.indC),    # PO.indC = PO.indC
        PO.hlf  = coef(PO.hlf ),    # PO.hlf  = PO.hlf 
        PO.hlfC = coef(PO.hlfC),    # PO.hlfC = PO.hlfC
        PO.dep  = coef(PO.dep ),    # PO.dep  = PO.dep 
        PO.depC = coef(PO.depC)))   # PO.depC = PO.depC
    }
    # 3.2 贴上各种总体数据和基于总体数据的结果
    # 总体数据
    Tot.mean = apply(Sample.comple, 2, mean)
    Tot.cov  = cov(Sample.comple)
    Tot.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.comple)
    Tot.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.comple)
    Tot.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.comple)
    Tot.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.comple)
    Tot.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.comple)
    Tot.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.comple)    
    out[[R+1]] = list(
    # Tot.data = Sample.comple,
    Tot.mean = Tot.mean,
    Tot.cov  = Tot.cov ,
    Tot.ind  = coef(Tot.ind ),      # Tot.ind  = Tot.ind 
    Tot.indC = coef(Tot.indC),      # Tot.indC = Tot.indC
    Tot.hlf  = coef(Tot.hlf ),      # Tot.hlf  = Tot.hlf 
    Tot.hlfC = coef(Tot.hlfC),      # Tot.hlfC = Tot.hlfC
    Tot.dep  = coef(Tot.dep ),      # Tot.dep  = Tot.dep 
    Tot.depC = coef(Tot.depC))      # Tot.depC = Tot.depC
    # 网络总体数据
    Net.mean = apply(Sample.net, 2, mean)
    Net.cov  = cov(Sample.net)
    Net.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.net)
    Net.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.net)
    Net.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.net)
    Net.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.net)
    Net.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.net)
    Net.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.net)   
    out[[R+2]] = list(
    # Net.data = Sample.net,
    Net.mean = Net.mean,
    Net.cov  = Net.cov ,
    Net.ind  = coef(Net.ind ),      # Net.ind  = Net.ind 
    Net.indC = coef(Net.indC),      # Net.indC = Net.indC
    Net.hlf  = coef(Net.hlf ),      # Net.hlf  = Net.hlf 
    Net.hlfC = coef(Net.hlfC),      # Net.hlfC = Net.hlfC
    Net.dep  = coef(Net.dep ),      # Net.dep  = Net.dep 
    Net.depC = coef(Net.depC))      # Net.depC = Net.depC    
    # 返回结果
    return(out)
}

# 测试了一下, 差不多10s对标100个Trail, 就意外还挺快; 运存方面的话差不多1e6(百万)这个量级的话在10个G, 
#   当然实际运存用不了这么多, 但几乎会顶满; 所以对前面的部分改了一下, 删掉了没用的数据存储, 虽然增加
#   了后面整理的难度, 但这是值得的.jpg; 目前的情况是运行过程中大致需要700MB的运存(当然这里面包括了R自
#   己要占用的), 结束后的结果存储大致需要30MB。另, 务必记得回收内存, 模拟的部分写的过于随意, 虽然R也会
#   自动回收, 但不是即时回收的, 强迫症的话建议自己释放一下

# st = proc.time()
# result = simulation(seed = 1234, R = 1000L, Parameter)
# en = proc.time()
# gc()                      
# object.size(result)

#===========================================================================================
# 模拟部分2: 直接模拟从两个给定总体中抽样的过程, 相对而言没那么“拟真”, 但更容易测试一些更有压力的情况

# 指定参数
Parameter = list(
# 控制变量的随机数发生器（这种情况下需要传6个随机数发生器, 3个用来指定总体, 3个用来指定网络访问固定样本）
# 随机数拟似规则:
# 以均匀分布拟似年龄(目标是成年人样本, 不考虑金字塔结构, 实际上整体偏“老龄”)
# 以截断泊松分布拟似受教育水平, 0-6依次对标3+3n的受教育年限(主要是为了让这个变量不用编码成哑变量, 要不然太麻烦了)
# 以对数正态分布拟似收入水平, 设定位置参数为1.5, 尺度参数为1^2
# 以均匀分布的平方拟似网络群体年龄(低龄, 特别是18-30这个年龄段占比调的很大)
# 以二项分布拟似网络群体受教育水平(初中到研究生对应的概率质量依次是0.0256, 0.1536, 0.3456, 0.3456, 0.1296)
# 以对数正态分布拟似收入水平, 设定位置参数为2, 尺度参数为
Controls = list(C1.P = function(n){runif(n, 0, 1)},             # C1 ~ U(0, 1); E[C1] = 1/2, D[C1] = 1/12                                                    
                C2.P = function(n){
                       out = rpois(n, 1)
                       for(i in 1:n) if(out[i] > 6){out[i] = 6}
                       out},                                    # C2 ~ Poisson(1) & C2 <= 6; E[C2] = 0.9999, D[C2] = 0.9989 
                C3.P = function(n){exp(rnorm(n, 1.5, 1))},      # C3 ~ lognorm(1.5, 1^2); E[C3] = exp(2), D[C3] = exp(5) - exp(4)
                C1.N = function(n){runif(n, 0, 1) ** 2},        # sqrt(C1) ~ U(0, 1); E[C1] = 1/3, D[C1] = 4/45
                C2.N = function(n){rbinom(n, 4, 0.6) + 2},      # C2-2 ~ Binom(4, 0.6); E[C2] = 4.4, D[C2] = 0.96
                C3.N = function(n){exp(rnorm(n, 2.25, 0.5))}),  # C3 ~ lognorm(2.25, 0.5^2)
# 样本量, controlN, expN分别代表随机样本数量和网络抽样研究样本数量
# 样本生成过程：从总体中抽出一个较大的样本(样本量controlN)用来模拟以标准抽样框进行的大型社会调查数据集,
#     然后从网络访问固定样本中抽出一个样本量为expN的样本, 然后在大型社会调查数据集中抽出等大的子样本作为
#     配对样本用来计算倾向值加权的权数 
N = list(controlN = 5e4, expN = 1000),
# 用来生成X2, X1, Y的参数, 分布和参数都是给定的, 所以不再计算总体值, 只给出理论值
para = list(beta =  c(1.3, 0.5, 0.3, 2, 0.35, 0.03), # 生成Y的参数
            alpha = c(1.1, 1.0, 0.1, 0.015),         # 生成X1的参数
            gamma = c(1.7, 0.7, 0.2, 0.020))         # 生成X2的参数
)
#===========================================================================================
# 模拟部分: 除了抽样的部分和前面几乎一毛一样
simulation = function(seed = 1234, R = 1000L, Parameter){

    # 0. 准备环节
    # 0.1 从Control中读取所需参数
    # 0.1.1 读取总体的随机数发生器
    C1.P = Parameter$Controls$C1.P
    C2.P = Parameter$Controls$C2.P
    C3.P = Parameter$Controls$C3.P
    # 0.1.2 读取网络访问固定样本的随机数发生器
    C1   = Parameter$Controls$C1.N
    C2   = Parameter$Controls$C2.N
    C3   = Parameter$Controls$C3.N    
    # 0.1.3 读取样本量
    N = Parameter$N$controlN
    n = Parameter$N$expN
    # 0.1.4 读取样本生成参数
    B = Parameter$para$beta
    A = Parameter$para$alpha
    G = Parameter$para$gamma
    # 0.1.5 初始化随机数发生器
    set.seed(seed)

    #---------------------------------------------------------------------------------------
    # 1. 生成样本
    # 1.1 生成协变量
    # 1.1.1 生成总体中的简单随机样本, 模拟大规模社会调查
    C1.SRS = C1.P(N)
    C2.SRS = C2.P(N)
    C3.SRS = C3.P(N)
    # 1.1.2 生成网络访问固定样本的抽样池（其实不生成也行, 主要是懒得大改代码了, 这里设计池容量为1e7/千万级）
    Nt = 1e7
    C1.sample = C1(Nt)
    C2.sample = C2(Nt)
    C3.sample = C3(Nt)
    # 1.2 然后根据需求生成三种不同情况下的模拟数据
    # 实际上生成的过程可以借助矩阵乘法写的更优雅一点, 但我懒.jpg
    # 1.2.1 独立的情况
    X1.indep  = rnorm(Nt)
    X2.indep  = rnorm(Nt)
    # 1.2.2 X2不独立的情况
    X2.depnd  = G[1] + G[2] * C1.sample + G[3] * C2.sample + G[4] * C3.sample + rnorm(Nt)    
    # 给出C1-C3的均值和协差
    E = c(1/2,  0.9999053, exp(2))
    D = diag(c(1/12, 0.9989321, exp(5) - exp(4)))
    X2.mean   = as.numeric(G[1] + t(G[2:4]) %*% E)
    X2.sd     = as.numeric(sqrt(t(G[2:4]) %*% D %*% G[2:4] + 1))
    X2.depnd  = (X2.depnd - X2.mean) / X2.sd
    # 1.2.3 X1不独立的情况
    X1.depnd  = A[1] + A[2] * C1.sample + A[3] * C2.sample + A[4] * C3.sample + rnorm(Nt)
    X1.mean   = as.numeric(A[1] + t(A[2:4]) %*% E)
    X1.sd     = as.numeric(sqrt(t(A[2:4]) %*% D %*% A[2:4] + 1))
    X1.depnd  = (X1.depnd - X1.mean) / X1.sd
    # 1.2.4 生成三种情况下的Y
    # X1/X2都独立于控制变量的情况
    Y1 = B[1] + B[2] * X1.indep + B[3] * X2.indep + B[4] * C1.sample + B[5] * C2.sample + B[6] * C3.sample + rnorm(Nt)
    # X2不独立于控制变量
    Y2 = B[1] + B[2] * X1.indep + B[3] * X2.depnd + B[4] * C1.sample + B[5] * C2.sample + B[6] * C3.sample + rnorm(Nt)
    # X1/X2都不独立于控制变量
    Y3 = B[1] + B[2] * X1.depnd + B[3] * X2.depnd + B[4] * C1.sample + B[5] * C2.sample + B[6] * C3.sample + rnorm(Nt)
    # 1.2.5 组装数据框, 这里为了省事儿就直接装一起了
    # 先装网络样本池
    Population = data.frame(
    C1 = C1.sample, C2 = C2.sample, C3 = C3.sample,
    X1.ind = X1.indep, X2.ind = X2.indep, X1.dep = X1.depnd, X2.dep = X2.depnd,
    Y1 = Y1, Y2 = Y2, Y3 = Y3)
    # 再装模拟的SRS
    SRS = data.frame(C1 = C1.SRS, C2 = C2.SRS, C3 = C3.SRS)
    
    #---------------------------------------------------------------------------------------
    # 2. 没了
    # 惊不惊喜意不意外.jpg
    
    #---------------------------------------------------------------------------------------    
    # 3. 开始循环主体
    # 3.0 准备记录循环结果
    # 需要记录的东西太多了, 弄成一个向量太麻烦, 干脆直接整个List, 后面再整理
    # 记录list的格式是: list套list, 外层的list长度为Ne, 全是模拟结果（序号从1一直到Ne）外层list
    #   的每一个元素依然是list, 除最后两个外, 每个list的长度为4, 包含四种计算方法下的各种估计值,
    #   方法一共四种: UW, 未加权; IW, 倾向值逆加权; WC, 组加权; PO, 事后分层调整加权。每个方法
    #   list中包含8个元素, 依次为: mean, 均值; cov, 协方差; ind, 假定自变量与协变量独立; hlf, 假
    #   定部分自变量与协变量独立; dep, 假定全部自变量与协变量相关; 含C表示这个回归方程加入了协变量
    out = list()
    # 3.1 进入循环
    for(i in 1:R){
        # 3.1.1 模拟抽样
        # 抽样规则: 从网络抽样池里抽个样本量为n的样本, 然后从SRS里抽一个一样大小的样本做匹配样本
        # (1). 模拟在网络访问固定样本中抽样
        ind.exp = sample.int(Nt, size = n, replace = FALSE)
        Sample.exp = Population[ind.exp, ]
        # (2). 在SRS中二次抽取样本
        ind.ctl = sample.int(N, size = n, replace = FALSE)
        Sample.ctl = SRS[ind.ctl, ]
        # (3). 捏合用于计算倾向值（Propensity Score）的混合样本, 该样本前Ne行为网络样本, 后Ne行为简单随机样本
        PSCalc = cbind(rbind(Sample.exp[, c('C1', 'C2', 'C3')], Sample.ctl), Y = rep(c(1, 0), c(n, n)))
        #...................................................................................
        # 3.1.2 数据处理
        # (1). 计算倾向值得分
        PSfit  = glm(Y ~ C1 + C2 + C3, family = binomial(link = 'logit'), data = PSCalc) # 拟合logistic回归
        PS     = as.vector(1 / (1 + exp(-predict(PSfit))))                               # 计算倾向值得分
        # (2). 计算权重
        # 权重的计算方法来自刘展, 金勇进. (2016). 基于倾向得分匹配与加权调整的非概率抽样统计推断方法研究
        # Sample.exp等价于文中的匹配样本（记作M）, Sample.ctl等价于文中的目标样本（记作O）, 假定两样本的基础权数均相等
        # a. 逆加权调整
        IW = 1 / PS                      # 求倒数
        iw = IW[1:n] / sum(IW[1:n])      # 逆加权调整权重
        # b. 加权组调整
        G  = n/10                                                   # 设置每个组内的样本量为20/需要保证Ne是10的倍数
        WC = PS[order(PS)]                                          # 将倾向值按照大小排序
        WC = rep(tapply(WC, rep(1:G, each = 20), mean), each = 20)  # 对倾向值进行组内均值代换
        WC = 1 / WC[order(order(PS))]                               # 重排序并求倒数
        wc = WC[1:n] / sum(WC[1:n])                                 # 加权组调整权重
        # c. 事后分层调整（直接用了b的分组）
        PO = PSCalc$Y                                               # 取出示性变量
        PO = PO[order(PS)]                                          # 对示性变量按照PS的顺序排序
        PO = rep(tapply(PO, rep(1:G, each = 20), mean), each = 20)  # 求出每个层中的个体来自匹配样本的概率
        PO = PO / (1 - PO)                                          # 求优势比（相当于换了个方式估计倾向值）
        PO = 1 / PO[order(order(PS))]                               # 重排序并求倒数
        po = PO[1:n] / sum(PO[1:n])                                 # 事后分层调整权重
        #...................................................................................       
        # 3.1.3 各种计算
        # 未加权情况（最简单的情形）
        UW.mean = apply(Sample.exp, 2, mean)
        UW.cov  = cov(Sample.exp)
        UW.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp)
        UW.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp)
        UW.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp)
        UW.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp)
        UW.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp)
        UW.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp)
        # 逆加权调整
        IW.mean = apply(Sample.exp * iw, 2, sum)
        IW.cov  = WeightedCOV(Sample.exp, iw)
        IW.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp, weights = iw)
        IW.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp, weights = iw)
        IW.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp, weights = iw)
        IW.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp, weights = iw)
        IW.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp, weights = iw)
        IW.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp, weights = iw)        
        # 加权组调整
        WC.mean = apply(Sample.exp * wc, 2, sum)
        WC.cov  = WeightedCOV(Sample.exp, wc)
        WC.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp, weights = wc)
        WC.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp, weights = wc)
        WC.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp, weights = wc)
        WC.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp, weights = wc)
        WC.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp, weights = wc)
        WC.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp, weights = wc)         
        # 事后分层调整
        PO.mean = apply(Sample.exp * po, 2, sum)
        PO.cov  = WeightedCOV(Sample.exp, po)
        PO.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp, weights = po)
        PO.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp, weights = po)
        PO.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp, weights = po)
        PO.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp, weights = po)
        PO.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp, weights = po)
        PO.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp, weights = po)        
        #...................................................................................
        # 3.1.4 组装结果
        out[[i]] = list(UW = list(
        UW.mean = UW.mean, 
        UW.cov  = UW.cov , 
        UW.ind  = coef(UW.ind ),    # UW.ind  = UW.ind 
        UW.indC = coef(UW.indC),    # UW.indC = UW.indC
        UW.hlf  = coef(UW.hlf ),    # UW.hlf  = UW.hlf 
        UW.hlfC = coef(UW.hlfC),    # UW.hlfC = UW.hlfC
        UW.dep  = coef(UW.dep ),    # UW.dep  = UW.dep 
        UW.depC = coef(UW.depC)),   # UW.depC = UW.depC
        IW = list(
        IW.mean = IW.mean,
        IW.cov  = IW.cov ,
        IW.ind  = coef(IW.ind ),    # IW.ind  = IW.ind 
        IW.indC = coef(IW.indC),    # IW.indC = IW.indC
        IW.hlf  = coef(IW.hlf ),    # IW.hlf  = IW.hlf 
        IW.hlfC = coef(IW.hlfC),    # IW.hlfC = IW.hlfC
        IW.dep  = coef(IW.dep ),    # IW.dep  = IW.dep 
        IW.depC = coef(IW.depC)),   # IW.depC = IW.depC
        WC = list(
        WC.mean = WC.mean,
        WC.cov  = WC.cov ,
        WC.ind  = coef(WC.ind ),    # WC.ind  = WC.ind 
        WC.indC = coef(WC.indC),    # WC.indC = WC.indC
        WC.hlf  = coef(WC.hlf ),    # WC.hlf  = WC.hlf 
        WC.hlfC = coef(WC.hlfC),    # WC.hlfC = WC.hlfC
        WC.dep  = coef(WC.dep ),    # WC.dep  = WC.dep 
        WC.depC = coef(WC.depC)),   # WC.depC = WC.depC
        PO = list(
        PO.mean = PO.mean,
        PO.cov  = PO.cov ,
        PO.ind  = coef(PO.ind ),    # PO.ind  = PO.ind 
        PO.indC = coef(PO.indC),    # PO.indC = PO.indC
        PO.hlf  = coef(PO.hlf ),    # PO.hlf  = PO.hlf 
        PO.hlfC = coef(PO.hlfC),    # PO.hlfC = PO.hlfC
        PO.dep  = coef(PO.dep ),    # PO.dep  = PO.dep 
        PO.depC = coef(PO.depC)))   # PO.depC = PO.depC
    }
    ## 3.2 贴上各种总体数据和基于总体数据的结果
    ## 贴上一个网络总体数据的均值和协方差
    #Net.mean = apply(Population, 2, mean)
    #Net.cov  = cov(Population)  
    #out[[R+1]] = list(
    ## Net.data = Sample.net,
    #Net.mean = Net.mean,
    #Net.cov  = Net.cov)   
    ## 返回结果
    return(out)
}
#===========================================================================================



# st = proc.time()
# result = simulation(seed = 1234, R = 1000L, Parameter)
# en = proc.time()
# gc()                      
# object.size(result)
# 
# RESULT = matrix(0, 1000, 84)
# 
# for(i in 1:1000){
#     TMP = numeric(0)
#     for(j in 1:4){
#         for(k in 3:8){
#             TMP = c(TMP, result[[i]][[j]][[k]][-1])
#         }
#     }
#     RESULT[i, ] = TMP
# }
# 
# THETA = rep(c(IND, INDC, HLF, HLFC, DEP, DEPC), 4)






#===========================================================================================
# 最麻烦的部分要来了, 根据既有的参数还原总体参数、网络样本对应的总体参数, 以及网络样本对应的回归标准误
# 这部分是用来计算bias, MSE以及参数覆盖率的, 此外, 为了省事儿, 没有考虑截距
# 总体数据:
#   C1 ~ U(0, 1);              E[C1] = 1/2,        D[C1] = 1/12 
#   C2 ~ Poisson(1) & C2 <= 6; E[C2] = 0.9999053,  D[C2] = 0.9989321 
#   C3 ~ lognorm(1.5, 1^2);    E[C3] = exp(2),     D[C3] = exp(5) - exp(4)
# 模拟网络调查样本:
#   sqrt(C1) ~ U(0, 1);        E[C1] = 1/3,        D[C1] = 4/45
#   C2-2 ~ Binom(4, 0.6);      E[C2] = 4.4,        D[C2] = 0.96
#   C3 ~ lognorm(2.25, 0.5^2); E[C3] = exp(2.375), D[C3] = exp(5) - exp(4.75)
# 指定参数:
B = c(1.3, 0.5, 0.3, 2, 0.35, 0.03)
A = c(1.1, 1.0, 0.1, 0.015)
G = c(1.7, 0.7, 0.2, 0.020)
C = c(1/12, 0.9989321, exp(5) - exp(4))
D = diag(C)
k2= as.numeric(sqrt(t(G[2:4]) %*% D %*% G[2:4] + 1))
k1= as.numeric(sqrt(t(A[2:4]) %*% D %*% A[2:4] + 1))
# 装配回归矩阵
Beta = matrix(0, 10, 10)
rownames(Beta) = colnames(Beta) = c('C1', 'C2', 'C3', 'X1i', 'X2i', 'X1d', 'X2d', 'Y1', 'Y2', 'Y3')
Beta['X1d', c('C1', 'C2', 'C3')] = A[2:4] / k1
Beta['X2d', c('C1', 'C2', 'C3')] = G[2:4] / k2
Beta['Y1' , c('X1i', 'X2i', 'C1', 'C2', 'C3')] = B[2:6]
Beta['Y2' , c('X1i', 'X2d', 'C1', 'C2', 'C3')] = B[2:6]
Beta['Y3' , c('X1d', 'X2d', 'C1', 'C2', 'C3')] = B[2:6]
# 装配总体的残差矩阵
Psi.P = diag(c(C, 1, 1, 1/k1/k1, 1/k2/k2, 1, 1, 1))
# 计算理论总体协方差
Cov.P = solve(diag(10) - Beta) %*% Psi.P %*% t(solve(diag(10) - Beta))
# 装配网络样本的残差矩阵
Psi.N = diag(c(4/45, 0.96, exp(5) - exp(4.75), 1, 1, 1/k1/k1, 1/k2/k2, 1, 1, 1))
# 计算网络样本理论协方差
Cov.N = solve(diag(10) - Beta) %*% Psi.N %*% t(solve(diag(10) - Beta))
# 计算总体参数
# 计算六个回归的总体参数
IND  = as.vector(solve(Cov.P[c('X1i', 'X2i')                  , c('X1i', 'X2i')])                   %*% Cov.P[c('X1i', 'X2i')                  , 'Y1'])
INDC = as.vector(solve(Cov.P[c('X1i', 'X2i', 'C1', 'C2', 'C3'), c('X1i', 'X2i', 'C1', 'C2', 'C3')]) %*% Cov.P[c('X1i', 'X2i', 'C1', 'C2', 'C3'), 'Y1'])
HLF  = as.vector(solve(Cov.P[c('X1i', 'X2d')                  , c('X1i', 'X2d')])                   %*% Cov.P[c('X1i', 'X2d')                  , 'Y1'])
HLFC = as.vector(solve(Cov.P[c('X1i', 'X2d', 'C1', 'C2', 'C3'), c('X1i', 'X2d', 'C1', 'C2', 'C3')]) %*% Cov.P[c('X1i', 'X2d', 'C1', 'C2', 'C3'), 'Y1'])
DEP  = as.vector(solve(Cov.P[c('X1d', 'X2d')                  , c('X1d', 'X2d')])                   %*% Cov.P[c('X1d', 'X2d')                  , 'Y1'])
DEPC = as.vector(solve(Cov.P[c('X1d', 'X2d', 'C1', 'C2', 'C3'), c('X1d', 'X2d', 'C1', 'C2', 'C3')]) %*% Cov.P[c('X1d', 'X2d', 'C1', 'C2', 'C3'), 'Y1'])
# 计算六个回归的网络样本对应参数
NIND  = as.vector(solve(Cov.N[c('X1i', 'X2i')                  , c('X1i', 'X2i')])                   %*% Cov.N[c('X1i', 'X2i')                  , 'Y1'])
NINDC = as.vector(solve(Cov.N[c('X1i', 'X2i', 'C1', 'C2', 'C3'), c('X1i', 'X2i', 'C1', 'C2', 'C3')]) %*% Cov.N[c('X1i', 'X2i', 'C1', 'C2', 'C3'), 'Y1'])
NHLF  = as.vector(solve(Cov.N[c('X1i', 'X2d')                  , c('X1i', 'X2d')])                   %*% Cov.N[c('X1i', 'X2d')                  , 'Y1'])
NHLFC = as.vector(solve(Cov.N[c('X1i', 'X2d', 'C1', 'C2', 'C3'), c('X1i', 'X2d', 'C1', 'C2', 'C3')]) %*% Cov.N[c('X1i', 'X2d', 'C1', 'C2', 'C3'), 'Y1'])
NDEP  = as.vector(solve(Cov.N[c('X1d', 'X2d')                  , c('X1d', 'X2d')])                   %*% Cov.N[c('X1d', 'X2d')                  , 'Y1'])
NDEPC = as.vector(solve(Cov.N[c('X1d', 'X2d', 'C1', 'C2', 'C3'), c('X1d', 'X2d', 'C1', 'C2', 'C3')]) %*% Cov.N[c('X1d', 'X2d', 'C1', 'C2', 'C3'), 'Y1'])


# UW.ind  = lm(Y1 ~ X1.ind + X2.ind               , Sample.exp)
# UW.indC = lm(Y1 ~ X1.ind + X2.ind + C1 + C2 + C3, Sample.exp)
# UW.hlf  = lm(Y1 ~ X1.ind + X2.dep               , Sample.exp)
# UW.hlfC = lm(Y1 ~ X1.ind + X2.dep + C1 + C2 + C3, Sample.exp)
# UW.dep  = lm(Y1 ~ X1.dep + X2.dep               , Sample.exp)
# UW.depC = lm(Y1 ~ X1.dep + X2.dep + C1 + C2 + C3, Sample.exp)










