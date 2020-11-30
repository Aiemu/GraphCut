# 基于 GraphCut 的纹理生成
2017011438 曾正

### 一、开发环境及编译运行方式
#### 1. 开发环境
- `Clion 2020.2`
- `C++ 14`
- `cmake version 3.18.2`

#### 2. 编译及运行方式
在 `src` 目录下执行：
``` shell
cmake .
make
```
即可生成可执行文件 `GraphCut`，接着执行：
```
./GraphCut <model> <iteration_time> <patch_path> 
           <output_path> <output_width> <output_height>
```
参数说明如下：  
- `model`：块偏移生成算法，0: 随机生成，1: 匹配整个块的最佳位置，2: 匹配子块的最佳位置
- `iteration_time`：迭代次数，在一定范围内，次数越大，效果越好。建议设置为 `70`
- `patch_path`：纹理块输入路径
- `output_path`：合成结果输出路径
- `output_width`、`output_height`：合成结果高度和宽度

如：
```
./GraphCut 2 70 img/green.gif green.png 100 100
```

执行后，将在 `src/output` 中生成结果图片，包含原图和标出缝的 seams 图。

</br>
</br>

### 二、实现的功能
#### 1. 基础功能
实现了Project1要求的基本算法，包括三种块偏移生成算法
1. 随机生成 *Random Replacement*
2. 匹配整个块的最佳位置 *Entire patch matching*
3. 匹配子块的最佳位置 *Sub-patch matching*

#### 2. 附加功能
- `old cuts`，标记出了 `Seam nodes`（输出结果 `seams-*.png`）
- 能量函数引入梯度
- `FFT` 加速：在 *Entire patch matching* 算法中调用了 FFT 算法，具体实现在 `GraphCutTexture.fftAcc()` 中

</br>
</br>

### 三、代码实现思路
#### 1. 算法基本思路
1. 在开始阶段，将原始纹理图片放入画布（大小即目标图像的大小）上某一随机位置，作为起始图像。
2. 根据三种块偏移算法计算出下一纹理图放入的位置 
3. 对于新放入图和原图的重叠部分的像素点，靠近原图的像素点与 source 相连，靠近新放入图的点与 sink 相连
4. 计算出图的最小割（使用最大流、最小割算法），分别从新放入图和原图中取出合适的像素点作为重叠部分
5. 重复上述步骤，直到画布填满
6. 对于 Old cuts，只需要保存每个像素的来源（原图还是新放入图），若一对像素点有不同的来源，则它们是一对 Old cuts。然后将这一对点与中间的 seam 点连接即可

#### 2. 块偏移生成算法
1. 随机生成 *Random Replacement*  
   对随机生成的 offset，计算重叠部分像素的比例是否超过阈值，若是则作为下一张图像放入的位置

2. 匹配整个块的最佳位置 *Entire patch matching*  
    对随机生成的 offset，计算重叠部分像素的比例是否超过阈值，若是则生成一个 `[0, 1]` 的浮点数与全体 cost 比较，若 cost 较大则将 offset 作为下一张图像放入的位置

3. 匹配子块的最佳位置 *Sub-patch matching*  
   对随机生成的 offset，计算重叠部分像素的比例是否超过阈值，若是则生成一个 `[0, 1]` 的浮点数与所有可能的子块的 cost 比较，若 cost 较大则将 offset 作为下一张图像放入的位置

#### 3. 最小割算法
使用基于 BFS 和 DFS 的 Dinitz 算法，具体实现在 `MacFlow.dinitz()` 中，由于不是此次 Proj 的重点，此处不进行详细描述。

#### 4. 最大流算法
具体实现在 `MacFlow` 类中，由于不是此次 Proj 的重点，此处不进行详细描述。

</br>
</br>

### 四、实验结果及分析
实验结果如下：

<div style="display: flex;">
    <div style="width: 100%;">
        <img src="1.png">
        <div style="text-align: center;">图4.1.1</div>
    </div>
    <div style="width: 100%; margin-left:5%;">
        <img src="o1.png">
        <div style="text-align: center;">图4.1.2</div>
    </div>
</div> 

</br>
</br>
<div style="display: flex;">
    <div style="width: 100%;">
        <img src="2.png">
        <div style="text-align: center;">图4.2.1</div>
    </div>
    <div style="width: 100%; margin-left:5%;">
        <img src="o2.png">
        <div style="text-align: center;">图4.2.2</div>
    </div>
</div> 
</br>
</br>

<div style="display: flex;">
    <div style="width: 100%;">
        <img src="3.png">
        <div style="text-align: center;">图4.3.1</div>
    </div>
    <div style="width: 100%; margin-left:5%;">
        <img src="o3.png">
        <div style="text-align: center;">图4.3.2 len = 8</div>
    </div>
</div> 

可以看到其中草莓的效果较好，这也得益于本身纹理的不规则性。