# 使用方法

​			进入src文件夹, 在确保安装了gfortran的前提下, 输入:

```shell
gfortran main.f90 -o target
```

​			那么就会在src目录下面生成一个名为target的可执行文件, 然后, 我们需要修改在同一个目录下的名为parameter.nml的文件, 里面的具体内容如下:

```
&myList
length = 500.d0		! 多晶长度
gNumber = 5			! 第一层晶粒个数
del = 1				! 每一层的晶粒变化个数
layer = 20			! 一个多少个晶粒层
isAtom = .false.	! 是否写出ATOMSK格式的文件
Period = .false.	! 是否让晶粒周期性分布
cFactor = 0.9d0		! 扰动因子的大小, 越小偏离平衡位置越小
/

```

​			梯度生成方法和思想在补充材料中可以看到, 请同学们自行观看学习.