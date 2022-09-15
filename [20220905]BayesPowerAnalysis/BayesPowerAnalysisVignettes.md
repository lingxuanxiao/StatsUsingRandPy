# 贝叶斯估计与假设检验中的样本量规划问题

Kazuha_Yae

2022/09/15

主要参考文献: Kruschke, J. K. (2015) Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan (2nd Ed.). Elsevier Inc.

> 相关信息：
> 
> R版本：≥ 4.1
> 
> 依赖包：rstan，bayestestR，BayesFactor
> 
> ___特别说明：这部分主要是关于基础理论、所用包与示例的介绍，方便的检验工具尚在开发中___
> 
> 联系作者：[sulfonamides@163.com](mailto:sulfonamides@163.com)

## 一、贝叶斯估计与贝叶斯假设检验

#### &emsp;&ensp;1. 贝叶斯定理、贝叶斯估计

&emsp;&emsp;首先，我们默认读者已经具备了基础的概统知识。根据条件概率与联合概率之间的关系，容易有：

$$
P(A|B)*P(B)=P(AB)=P(B|A)*P(A)
$$

据此，容易有贝叶斯公式：

$$
P(A|B)=\frac{P(AB)}{P(B)}=P(B|A)\ast\frac{P(A)}{P(B)}
$$

现假定样本 $\mathbf{x}=(x_1,x_2,...,x_n)$ 来自于参数为 $\theta$ 的某一总体，参数 $\theta$ 属于参数空间 $\Theta$。根据贝叶斯公式，当我们已知样本 $\mathbf{x}$ 时，参数 $\theta$ 的分布为 $f_{\theta|X}(\theta|\mathbf{x})=f_{X|\theta}(\mathbf{x}|\theta)\ast \frac{f_\theta(\theta)}{f_X(\mathbf{x})}$ 。又注意到，在给定样本的情况下，$f_X(\mathbf{x})=\int\limits_{t\in \Theta}^{} f_{X|\theta}(\mathbf{x}|t)f_\theta(t)\mathrm{d}t$ 是一个定值，故可将其视为对分子乘积项 $f_{X|\theta}(\mathbf{x}|\theta)\ast f_\theta(\theta)$ 的归一化常数，实际决定后验分布 $f_{\theta|X}(\theta|\mathbf{x})$ 的只有 似然 $f_{X|\theta}(\mathbf{x}|\theta)$ 和先验分布 $f_\theta(\theta)$ ，即 $posterior \propto likelihood \ast prior$。由于后验分布的估计依赖于先验分布，故在贝叶斯估计中，根据某些准则选择合适的先验分布是必要的步骤。~~【诚然，根据“力大砖飞”的自然规律，只要样本量足够大，什么样的先验分布（当然，不能是那些病态的先验分布）都能被样本信息修正成应该有的样子，但考虑到实验经费和审稿人，最好还是不要选择一些奇奇怪怪的分布】~~ 通常情况下，可以考虑下述先验分布：

- 扁平先验分布/扩散先验分布/非主观先验分布
  
  此类先验分布可能是最接近对“无知”的模拟，这类分布的特点是概率质量/概率密度为一常数，值得注意的是，当选择此类分布作为先验分布时，分布可以是“不恰当”（即积分为 $\infty$）的，只需保证后验分布是恰当的即可。
- Jeffreys先验/无信息先验分布/不变先验分布
  
  此类先验分布的特点是在单调变换中具有不变性，对此类先验分布而言，不管我们选择了哪种参数化的方法，只要参数之间的变换具有“一对一”的性质，那么这些参数的先验分布总是相通的。换言之，Jeffreys先验保障了我们选择不同参数模型分析相同数据时所得结果的稳定性。
  
  例如：对于二项分布 $B(n,p)$ 来说，其参数 $p$ （事件发生几率）的Jeffreys先验是 $Beta(0.5, 0.5)$，而当我们以对其做logit变换后 $\mu=\ln{\frac{p}{1-p}}$ 的量作为参数重构该分布时（这个过程通常也叫重参数化/re-parameterization），新参数 $\mu$ 的Jeffreys先验与原参数 $p$ 的Jeffreys先验具有下述关系：$f_\mu(\mu)=f_{Beta}(\frac{e^\mu}{1+e^\mu},0.5,0.5)\ast\frac{e^\mu}{(1+e^\mu)^2} $
  即新参数的Jeffreys先验可由原参数的Jeffreys先验进行密度函数变换得来。
- 模糊先验分布，弱信息先验分布
  
  此类分布与其说是一类具体的分布，不如说是指定尺度参数时的策略。通常情况下，一个方差很大的分布总能更广泛地涵盖各种参数可能，并给予“远处”的值更可观的概率质量/概率密度。在实践中，此类分布多以方差巨大的正态分布给出，例如 $N(0,10^2)$，$N(0,10000^2)$ 等。
- 一般弱信息先验分布，具体信息先验分布
  
  当有足够的理论、先例或学术习惯支持时，我们也可以选择一些信息更多的分布。例如在心理学中，对t检验标准化效应量（即两组样本的均值差除以联合标准差，其实就是Cohen's d啦）的先验分布往往选择位置参数为0，尺度参数为 $\frac{\sqrt[]{2}}{2}$ 的柯西分布。
- 共轭先验分布族
  
  共轭先验分布族可能不是所有先验分布中最好理解的，但一定是所有先验分布中最好计算的，共轭先验分布来自于这样一类有趣的性质：
  
  > 对某一分布 $f(x|\theta)$ ，当其参数的先验分布为 $p(\theta)$ 时，后验分布 $p(\theta|x)$ 与先验分布从属于相同的分布族（直白一点，相当于后验分布和先验分布属于同一类分布，只不过参数不同）。
  
  最有名的一对共轭分布莫过于Beta-二项共轭，即当某二项分布随机变量的参数的先验分布为Beta分布时，则该参数的后验分布也是Beta分布。共轭分布理论大大降低了贝叶斯估计中的计算量，将复杂的积分计算变为小学加减法。不过，不是所有分布的共轭先验分布都能显式求得，但对于我们最常见的一类分布族——指数分布族（日常应用中常见的分布，如正态分布、指数分布、Beta分布、Gamma分布、二项分布、几何分布等都属于指数分布族）而言，所有指数分布族的成员都有共轭先验分布。关于共轭先验分布的具体信息（共轭形式、后验超参数的估计等）可以参看[维基百科：共轭先验分布](https://en.wikipedia.org/wiki/Conjugate_prior)。

&emsp;&emsp;当我们已知先验分布与似然时，我们就可以通过某些方式求得后验分布，这些方式包括~~口算法【Bushi~~借助共轭先验分布的性质计算后验分布、马尔科夫链蒙特卡罗（Markov Chain Monte Carlo, MCMC）法与变分贝叶斯近似/变分贝叶斯推断法。特别地，如果我们只希望得到后验分布最大密度处的参数取值，在后验分布的概率密度容易求导的情况下，可以使用共轭梯度法或牛顿法（不仅需要一阶导数，还需要二阶导数）等求取后验概率密度的极大值点。如果连求导都成为了一种奢侈的时候，我们还有期望最大化方法可以用，[这里](https://www.bilibili.com/read/cv18588785)对期望最大化方法进行了一个简单的介绍。求出后验分布后，我们可以根据后验分布对参数 $\theta$ 进行估计。常用的估计方式包括：

| 类属   | 方法名称     | 估计方式                                                |
|:----:|:--------:|:---------------------------------------------------:|
| 点估计  | 后验均值估计   | 后验分布的期望                                             |
| 点估计  | 后验中位数估计  | 后验分布的中位数                                            |
| 点估计  | 最大后验估计   | 后验分布的众数/后验密度最大时对应的参数                                |
| 区间估计 | 最大后验密度区间 | 满足置信度 $1-\alpha$ 的最小（短）区间                           |
| 区间估计 | 等尾区间     | 后验分布的 $\frac{\alpha}{2}$ 与 $1-\frac{\alpha}{2}$ 分位点 |

#### &emsp;&ensp;2. 贝叶斯假设检验

&emsp;&emsp;考虑来自总体密度函数为 $f(x;\theta)$ 的样本 $\mathbf{x}=(x_1,x_2,...,x_n)$，$\theta \in \Theta$，$\Theta_0$ 和 $\Theta_1$ 是对参数空间 $\Theta$ 的一个划分，待检验假设为：

$$
H_0: \theta \in \Theta_0 \qquad vs \qquad H_1: \theta \in \Theta_1 
$$

显然，在我们知晓后验分布 $p(\theta|\mathbf{x})$ 的情况下，我们可以分别计算原假设和备择假设成立的概率：

$$
P(H_0|\mathbf{x})=\int\limits_{\theta \in \Theta_0}^{} p(\theta|\mathbf{x}) \mathrm{d}\theta \qquad P(H_1|\mathbf{x})=\int\limits_{\theta \in \Theta_1}^{} p(\theta|\mathbf{x}) \mathrm{d}\theta
$$

考虑二者的比值，容易有：

$$
\frac{P(H_0|\mathbf{x})}{P(H_1|\mathbf{x})} = \frac{\int\limits_{\theta \in \Theta_0}^{} p(\mathbf{x}|\theta) g(\theta) \mathrm{d}\theta}{\int\limits_{\theta \in \Theta_1}^{} p(\mathbf{x}|\theta) g(\theta) \mathrm{d}\theta}
$$

式中 $g(\theta)$ 表示参数的先验分布。此时，仿照频率学派中的“似然比”概念，定义：

$$
BF_{01}=\frac{P(\mathbf{x}|H_0)}{P(\mathbf{x}|H_1)}=\frac{P(H_0|\mathbf{x})}{P(H_1|\mathbf{x})}\frac{P(H_1)}{P(H_0)} = \frac{\int\limits_{\theta \in \Theta_0}^{} p(\mathbf{x}|\theta) g(\theta) \mathrm{d}\theta}{\int\limits_{\theta \in \Theta_1}^{} p(\mathbf{x}|\theta) g(\theta) \mathrm{d}\theta} \frac{\int\limits_{\theta \in \Theta_1}^{} g(\theta) \mathrm{d}\theta}{\int\limits_{\theta \in \Theta_0}^{} g(\theta) \mathrm{d}\theta}
$$

特别地，若原假设为简单假设，备择假设为复合假设时，如$H_0:\theta=0\quad vs \quad H_1:\theta \ne 0$。此时 $P(H_0|\mathbf{x})$ 与 $P(H_0)$ 均为0，$P(H_1|\mathbf{x})$ 与 $P(H_1)$ 均为1，有研究者建议当遇到这种情形时，以Savage-Dickey密度比（即 $\frac{p(H_0|\mathbf{x})}{p(H_0)}$）替代原式。$BF_{01}$ 意味着“相较于备择假设，实际数据在多大程度上更加支持原假设”。虽然 $BF_{01}$ 不像似然比那样具有明确的可讨论的分布，但其本身就蕴含有一定的决策性质。有学者给出了如下的决策标准：

| 取值范围           | 如果 $BF_{01}$ 落在取值范围内 | 如果 $BF_{01}$ 的倒数落在取值范围内 |
|:--------------:|:--------------------:|:-----------------------:|
| $(100,\infty)$ | 极强证据支持原假设            | 极强证据支持备择假设              |
| $(30,100]$     | 非常强证据支持原假设           | 非常强证据支持备择假设             |
| $(10,30]$      | 强证据支持原假设             | 强证据支持备择假设               |
| $(3,10]$       | 中等强度证据支持原假设          | 中等强度证据支持备择假设            |
| $(1,3]$        | 弱证据支持原假设             | 弱证据支持备择假设               |
| $1$            | 没有证据支持任何假设           | 没有证据支持任何假设              |

&emsp;&emsp;当然，除采用 $BF_{01}$ 作为假设检验的决策标准之外，与频率学派类似，在贝叶斯学派中，研究者也可以使用后验区间估计作为假设检验的决策标准。例如：后验区间的95%最大后验密度区间不包含0；后验区间的95%等尾区间不包含0等。比较特殊的是，在贝叶斯学派中，还可以通过实际对等域（Region of Practical Equivalence，ROPE）来进行决策，实际对等域

SSSS

#### &emsp;&ensp;3. 对bayestestR 与 BayesFactor 的简单介绍

#### &emsp;&ensp;4. rstan与概率编程

## 二、贝叶斯估计与假设检验中的样本量规划问题

#### &emsp;&ensp;1. 统计检验力的扩展

https://easystats.github.io/bayestestR/

https://richarddmorey.github.io/BayesFactor/

10.3758/PBR.16.2.225