program main
  use iso_fortran_env, only: stdout => output_unit, wp => real64
  implicit none
  real(wp), allocatable :: dK(:)             !< 某一层晶粒的尺寸
  real(wp), allocatable :: coordinate(:, :)  !< 含有所有控制点坐标的数组
  real(wp) :: L        !< 二维模型的宽度
  real(wp) :: xlim, ylim   !< 两个方向上的的坐标极限
  real(wp) :: chaosFactor  !< 一个在(0, 1)之间的实数, 决定了控制点偏离平衡位置的程度
  integer :: index = 1 !< 指向某一个晶粒的坐标
  integer :: K        !< 梯度晶粒的层数
  integer :: first    !< 第一层晶粒数
  integer :: delta    !< 每一层减少的晶粒数
  integer :: gNum     !< 每一层的晶粒数
  integer :: layers   !< 一共多少个晶粒
  integer :: i, j     !< 循环变量
  logical :: isAtomsk !< 判断是否写出atomsk的文件
  logical :: isPeriod !< 是否考虑让晶粒对称周期性分布

  call ReadParameter('parameter.nml', L, first, delta, K, isAtomsk, isPeriod, chaosFactor)
  layers = K * first + (K * (K - 1)) / 2 * delta !! 利用等差数列计算需要分配多少内存

  !! 分配内存
  allocate(dk(K))                                
  allocate(coordinate(layers, 3))   !! 有3列数据, 分别是控制节点的xy坐标, 和属于第几层的晶粒

  !! 初始化所有晶粒尺寸为0
  dk = 0.d0   
  call random_seed  !! 初始化随机数种子
  
  do i = 1, K
    gNum = first + (i - 1 ) * delta         !! 第 i 层的晶粒个数
    do j = 1, gNum      
      dK(i) = L / (gNum)                    !! 计算第 i 层晶粒的尺寸
      if ( i == 1) then
        coordinate(index, 1) = (0.5 * dK(i)) + (j - 1) * dK(i)  !! 计算晶粒节点的x坐标
        coordinate(index, 2) = dK(i) * 0.5d0                    !! 计算晶粒节点的y坐标
        coordinate(index, 3) = real(i, kind=wp)                 !! 声明是第几层的晶粒, 方便后续调用
        index = index + 1
      else
        coordinate(index, 1) = (0.5 * dK(i)) + (j - 1) * dK(i) 
        coordinate(index, 2) = (sum(dK(1 : i - 1)) + dk(i) * 0.5d0)  
        coordinate(index, 3) = real(i, kind=wp)
        index = index + 1
      end if
    end do
  end do
  !! 循环结束, 此时的index - 1 也正好等于晶粒个数
  xlim = L
  if (isPeriod) then  !! 计算x和y区域的上下限
    ylim = 2.d0 * sum(dK) - dk(1)
  else
    ylim = sum(dK)      
  end if
  call WriteData(coordinate, 'seed.txt')
  !! 生成完之后将相关信息输出到屏幕上面
  if (isPeriod) then
    write(stdout, '(A,I0,A)') '控制点生成结束, 一共生成了', 2*(index - 1)-first,'个控制节点.'
  else
    write(stdout, '(A,I0,A)') '控制点生成结束, 一共生成了', index - 1,'个控制节点.'
  end if

  write(stdout, '(A,I0,A)') '其中第一层的控制节点个数为: ', first, '个.'
  write(stdout, '(A,I0,A)') '每层变化的控制节点个数为: ', delta, '个.'

  if (isPeriod) then
    write(stdout, '(A,I0,A)') '总计生成', 2*K-1, '层梯度层.'
  else
    write(stdout, '(A,I0,A)') '总计生成', K, '层梯度层.'
  end if
  write(stdout, '(A, G0, A, G0)') 'x方向上长度: ', xlim, ', y方向上长度: ', ylim
  if (isAtomsk) call WriteAtomsk(coordinate, 'gradient.txt', xlim, ylim)    
contains 
  subroutine ReadParameter(nmlName, length, gNumber, del, layer, isAtom, period, cFactor)
    character(len=*), intent(in) :: nmlName    !< 含有nml格式的文件的名称
    real(wp), intent(out) :: length            !< 多晶模型中的长度
    real(wp), intent(out) :: cFactor           !< 一个在(0, 1)之间的实数, 决定了控制点偏离平衡位置的程度
    integer, intent(out) :: gNumber            !< 第一层晶粒的数量
    integer, intent(out) :: del                !< 每一层减少的晶粒数量
    integer, intent(out) :: layer              !< 一共有多少层
    logical, intent(out) :: isAtom             !< 是否写出atomsk格式的控制节点文件
    logical, intent(out) :: period             !< 是否考虑让晶粒对称周期性分布
    integer :: nmlID                           !< 文件ID通道
    integer :: ioflag                          !< 看看文件读取是否成功
    namelist/myList/ length, gNumber, del, layer, isAtom, period, cFactor

    open(newunit=nmlID, file=nmlName, action='read')    !! 打开文件
    read(unit=nmlID, nml=myList, iostat=ioflag)         !! 按namelist方式读取
    close(nmlID)                                        !! 关闭文件
    !! 读取之后应该对参数进行一些基本的判断
    if (del == 0 ) then
      write(stdout , '(A)') 'You must make variable del not equal to zero.'
      stop
    else if (del < 0 ) then
      if ((gNumber + (layer - 1) * del) < 0) then
        write(stdout , '(A)') 'Negative number of grain, please check your parameter.'
        stop
      end if
    end if
    if (gNumber < 0.d0) stop 'Please give parameter gNumer a positive number.'
    if (length < 0.0d0) stop 'Please give parameter length a positive number.'
    if (cFactor < 0.d0 .or. cFactor > 1.d0) stop 'chaosFactor must between 0 and 1, please check your input.'

  end subroutine ReadParameter

  subroutine WriteData(array, dataname)
    !! 将计算出来的结果写入文件中方便MATLAB用来绘图
    real(wp), intent(in) :: array(:, :)      !< 包含所有控制节点坐标的数组
    character(len=*), intent(in) :: dataname !< 写出去的数据文件的名称
    real(wp) :: chaosx   !< 晶粒的扰动, 在x方向上
    real(wp) :: chaosy   !< 晶粒的扰动, 在y方向上
    integer :: dataID !< data文件的通道ID         
    integer :: x, y   !< 循环变量
    integer :: ndim   !< 晶粒的个数
    integer :: grainNum   !< 每一层的晶粒个数
    integer :: floor  !< 第几层

    ndim = size(array, dim=1)   !! 获取第一个维度的尺寸
    open(newunit=dataID, file=dataname, action='write')
    if (isPeriod) then  !! 如果要求周期性的晶粒分布
      do x = 1, ndim
        floor = nint(coordinate(x, 3))        !! 查询第x个晶粒是第几层
        call random_number(chaosx)            !! 生成一个(0, 1)之间的随机数
        call random_number(chaosy)
        chaosx = (chaosx - 0.5d0) * dK(floor) * chaosFactor     !! 保证扰动在正负的晶粒尺寸一半以内
        chaosy = (chaosy - 0.5d0) * dK(floor) * chaosFactor     !! 0.95是扰动因子, 越小说明晶粒扰动越小
        write(dataID, *) array(x, 1) + chaosx, array(x, 2) + chaosy + sum(dK(2:K))
      end do
      do x = first + 1, index - 1      !! 将相关的晶粒输出一遍
        floor = nint(coordinate(x, 3))
        call random_number(chaosx)            !! 生成一个(0, 1)之间的随机数
        call random_number(chaosy)
        chaosx = (chaosx - 0.5d0) * dK(floor) * chaosFactor     !! 保证扰动在正负的晶粒尺寸一半以内
        chaosy = (chaosy - 0.5d0) * dK(floor) * chaosFactor     !! chaosFactor是扰动因子, 越小说明晶粒扰动越小
        write(dataID, *) coordinate(x, 1) + chaosx, 2.d0 * coordinate(1, 2) - coordinate(x, 2) + chaosy + sum(dK(2:K))
      end do
    else
      do x = 1, ndim
        floor = nint(coordinate(x, 3))        !! 查询第x个晶粒是第几层
        call random_number(chaosx)            !! 生成一个(0, 1)之间的随机数
        call random_number(chaosy)
        chaosx = (chaosx - 0.5d0) * dK(floor) * chaosFactor     
        chaosy = (chaosy - 0.5d0) * dK(floor) * chaosFactor     
        write(dataID, *) array(x, 1) + chaosx, array(x, 2) + chaosy
      end do
    end if
    close(dataID)

  end subroutine WriteData

  subroutine WriteAtomsk(array, modelname, xdim, ydim)
    !! 以atomsk可以读取的格式写入文本文件之中
    real(wp), intent(in) :: array(:, :)       !< 包含所有控制节点坐标的数组
    real(wp), intent(in) :: xdim, ydim        !< 仿真盒子在xy方向上的长度
    character(len=*), intent(in) :: modelname !< 写出去的数据文件的名称
    real(wp) :: chaosx   !< 晶粒的扰动, 在x方向上
    real(wp) :: chaosy   !< 晶粒的扰动, 在y方向上
    integer :: modelID                     
    integer :: x
    integer :: ndim 
    integer :: floor  !< 第几层

    ndim = size(array, dim=1)   !! 获取第一个维度的尺寸
    open(newunit=modelID, file=modelname, action='write')
    write(modelID, *) 'box', xdim, ydim, 0.d0
    if (isPeriod) then  !! 如果要求周期性的晶粒分布
      do x = 1, ndim
        floor = nint(coordinate(x, 3))        !! 查询第x个晶粒是第几层
        call random_number(chaosx)            !! 生成一个(0, 1)之间的随机数
        call random_number(chaosy)
        chaosx = (chaosx - 0.5d0) * dK(floor) * chaosFactor     !! 保证扰动在正负的晶粒尺寸一半以内
        chaosy = (chaosy - 0.5d0) * dK(floor) * chaosFactor     !! 0.95是扰动因子, 越小说明晶粒扰动越小
        write(modelID, *) 'node', array(x, 1) + chaosx, array(x, 2) + chaosy + sum(dK(2:K)), 0.d0, 'random'
      end do
      do x = first + 1, index - 1      !! 将相关的晶粒输出一遍
        floor = nint(coordinate(x, 3))        !! 查询第x个晶粒是第几层
        call random_number(chaosx)            !! 生成一个(0, 1)之间的随机数
        call random_number(chaosy)
        chaosx = (chaosx - 0.5d0) * dK(floor) * chaosFactor     !! 保证扰动在正负的晶粒尺寸一半以内
        chaosy = (chaosy - 0.5d0) * dK(floor) * chaosFactor     !! 0.95是扰动因子, 越小说明晶粒扰动越小
        write(modelID, *) 'node', coordinate(x, 1) + chaosx, 2.d0 * coordinate(1, 2) - coordinate(x, 2) &
         + chaosy + sum(dK(2:K)), 0.d0, 'random'
      end do
    else 
      do x = 1, ndim
        floor = nint(coordinate(x, 3))        !! 查询第x个晶粒是第几层
        call random_number(chaosx)            !! 生成一个(0, 1)之间的随机数
        call random_number(chaosy)
        chaosx = (chaosx - 0.5d0) * dK(floor) * chaosFactor     !! 保证扰动在正负的晶粒尺寸一半以内
        chaosy = (chaosy - 0.5d0) * dK(floor) * chaosFactor     !! 0.95是扰动因子, 越小说明晶粒扰动越小
        write(modelID, *) 'node', array(x, 1) + chaosx, array(x, 2) + chaosy, 0.d0, 'random'
      end do
    end if
    close(modelID)
  end subroutine WriteAtomsk

end program main