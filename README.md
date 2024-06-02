# Numerical Linear Algebra

C++实现数值线性代数经典算法。

参考教材：*《数值线性代数（第2版）》——徐树方、高立、张平文*

### **0.【代码结构】**

程序包含两个主要文件 `Funcion.h` 和 `Ecercise.h` 。 

`Funcion.h` 中实现矩阵类（支持各种基本运算、矩阵转置、LU 分解、 Cholesky 分解、QR分解、上Hessenberg化、双重步位移QR迭代、二对角化），基本方程组求解方法（上三角、下三角、Guass、全主元Guass、列主元Guass、Cholesky、Cholesky改进），范数计算方法（1范数、无穷范数），方程组古典迭代解法（Jacobi、G-S、JOR），实用共轭梯度法，幂法求模最大根，隐式QR算法，过关Jacobi法，二分法求第K大特征值，反幂法，SVD迭代。

`Ecercise.h` 中构建矩阵并求解。

### **1.【线性方程组直接解法】**

不选主元、全主元、列主元三种Guass消去法，Cholesky分解及其改进版。

[【report】](./code/homework1/repoort.md)

### **2.【方程组解误差分析】**

矩阵范数计算、方程求解误差分析。

[【report】](./code/homework2/repoort.md)

### **3.【最小二乘】**

QR分解算法求解线性方程组、最小二乘问题。

[【report】](./code/homework3/repoort.md)

### **4.【线性方程组古典迭代解法】**

Jacobi迭代法、G-S迭代法、SOR迭代法求解方程组。

[【report】](./code/homework4/repoort.md)

### **5.【共轭梯度法】**

实用共轭梯度法。

[【report】](./code/homework5/repoort.md)

### **6.【非对称特征值】**

幂法求模特征根、QR方法（上Hessenberg分解、双重步位移QR迭代、隐式QR法）

[【report】](./code/homework6/repoort.md)

### **7.【对称特征值】**

过关Jacobi法、二分法、反幂法。

[【report】](./code/homework7/repoort.md)

### **8.【对称特征值】**

矩阵二对角化、SVD迭代。

[【report】](./code/homework8/repoort.md)