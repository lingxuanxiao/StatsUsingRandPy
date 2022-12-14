# 调查结果是否会误导我们的认识？——非概率样本中的统计推断问题

Kazuha_Yae

2022/10/30

主要参考文献: 

[1] 刘展,金勇进.网络访问固定样本调查的统计推断研究[J].统计与信息论坛,2017,32(02):3-10.

[2] 刘展,金勇进.基于倾向得分匹配与加权调整的非概率抽样统计推断方法研究[J].统计与决策,2016(21):4-8.DOI:10.13546/j.cnki.tjyjc.2016.21.001.

[3] 刘展,潘莹丽,涂朝凤,张梦.基于倾向得分的伪权数构造与混合样本推断[J].统计与决策,2021,37(02):20-24.DOI:10.13546/j.cnki.tjyjc.2021.02.004.

> 相关信息：
> 
> R版本：≥ 4.1
> 
> 依赖包：sampling
> 
> ___特别说明：这部分主要是大概的摘要思路和简单的实验报告___
> 
> 联系作者：[sulfonamides@163.com](mailto:sulfonamides@163.com)

## 一、调查结果是否误导了我们的认识？

+ 关于抽样偏差带来的统计推断问题，不同的学科似乎都形成了相似的共识：非概率样本似乎总有着这样那样的问题，以至基于非概率样本的统计推断总是“不可靠”的

+ 但随着技术的发展，网络时代产生的各种数据集中很少有实质上的严格概率样本（事实上哪怕是非网络时代，概率样本也是很难的，毕竟谁能保证上门调查的拒绝率是0%呢）

+ 但在回归系数为常数的假定下，当我们观测到了应该观测的控制变量时，我们得到的对偏回归系数的估计总是无偏的

+ 但话又说回来了，如果否认回归系数为常数的这个假定，那为什么还要用OLS回归呢。

+ 如果我们没有采集所有的控制变量进行控制，而这些控制变量恰恰既影响了自变量又影响了因变量，这种情况下如果目标是“知晓X与Y之间的关系”，那么无论如何我们都不能得到“干净”的结果，在这种情况下，甚至样本越有偏，反而控制的效果越好。特别地，在所有偏回归系数同向的情况下，样本方差越小（越不随机），反而得到的偏回归系数越“干净”。

+ 一个可能的攻讦是偏回归系数受到了未观测到的变量的影响（调节效应模型），但话又又说回来了，这种情况下抽样有没有偏得到的估计结果总归是有偏的，这又有什么意义呢？而反过来说，如果没有任何的有偏样本——特殊样本的加入，我们甚至永远无法知道是否有这种潜在的调节变量存在。

+ 有偏样本真正影响的是它的外推，这也正是为什么很多目标是均值的调查不能用非随机样本，你的目的都是外推了，怎么能不控制抽样；但对析因设计而言有偏无偏影响也就那回事儿

+ 至少在析因设计里，抽样框的设计远没有研究过程中的严格执行重要，特别是实验中的随机化处理以及问卷中的测量方法偏差控制（对风笑天的那个回应就是要区分出调查中的测量问题 vs 调查中的抽样问题）

+ 另外，说道方法控制，不是说“我给被试重赏”就一定能收到高质量数据，而且现代统计学的发展本身就是效率和准确不断相互妥协的结果，统计学是高度模型相依的，随着时代的变化模型的发展千变万化，80年代的抽样思想怎么能指导当今的研究，呵呵.jpg

+ 有偏样本可能导致对误差的过学习，确实会影响外推，但可能以某些方法纠正

+ 事实上，一个良好的态度是对研究的开放态度，即严格执行 + 标注全部信息 + 开放数据共享，Meta分析在很大程度上能够汇总前人研究，让我们向真理更加迈进

+ 总之，对于配额抽样而言，在直接能把配额的变量采为控制变量的情况下单纯从析因的角度完全没必要