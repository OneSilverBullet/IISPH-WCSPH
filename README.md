# IISPH-WCSPH


### 0. 项目初衷

该项目是实验室SPH方向的基础代码框架。SPH领域的开源代码实在太少了，因此我开发了这样的SPH流体仿真代码框架，供SPH领域的相关同学可以参考借鉴。

本项目有两个目的：

（1）应实验室的要求，构建一个完善的SPH流体仿真平台。

（2）开源给SPH流体仿真领域的相关同学，希望对相关领域研究的同学有所帮助。


### 1. 项目基本信息
* 编程语言：C++
* 基于框架：Windows
* 其他引用库：QT
* 实验环境：Visual Studio 2017
* 说明：这个项目框架使用的是QT，并且引用QT当中的openGL API进行绘制。为了更加专注于SPH的实现，本项目90%都是纯C++开发。涉及QT及渲染部分的代码非常少。这意味着，只需要你具有C++知识就可以进行阅读学习。


### 2. 项目实现的功能
* 基础的数学库：包含MATRIX,VECTOR等。
* 空间分割框架的实现：这是SPH粒子流实现的基础。
* WCSPH流体仿真：SPH的经典方法之一。
* IISPH流体仿真：SPH的经典方法之一。
* 基于距离场的边界判断：这是实现固液耦合的基础。
* Marching Cube算法：用于将粒子流体转换为有表面的流体。
* 完善的QT可视化界面：便于更改方法，可视化流体仿真效果的QT界面。

### 3. 目前实现的效果

![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/1.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/2.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/3.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/4.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/5.png)
![Image](https://github.com/OneSilverBullet/IISPH-WCSPH/blob/master/DEMO/6.png)
