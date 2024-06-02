

## 【Report】Homework 1**

### **一.【问题描述】**

实现不选主元、全主元、列主元三种Guass消去法，Cholesky分解及其改进版。并分别用其求解方程组。

## **二.【程序介绍】**

程序包含两个主要文件 `Funcion.h` 和 `Ecercise.h` 。

`Funcion.h` 中实现矩阵类（支持各种基本运算、矩阵转置、LU 分解和 Cholesky 分解）和基本方程组求解方法。

```cpp
class Matrix {
public:
    /*分解成功返回1，否则返回0*/
    bool LU();//LU分解
	bool FPLU(vector<int>&P,vector<int>&Q);//全主元LU分解
	bool CPLU(vector<int>&P);//列主元LU分解
	bool LL();//LL分解
	bool LDL();//LDL分解
private:
};

/*求解成功返回1，否则返回0*/
bool Sovle_Lxb(Matrix& L,Matrix& b,bool flag=0);//前代法(flag=1时表示单位下三角)
bool Sovle_Uxb(Matrix& U,Matrix& b);//回代法
/*求解成功返回答案矩阵，否则返回空矩阵*/
Matrix Solve_LU(Matrix A, Matrix b);//Guass
Matrix Solve_FPLU(Matrix A, Matrix b);//全主元Guass
Matrix Solve_CPLU(Matrix A, Matrix b);//列主元Guass
Matrix Solve_LL(Matrix A,Matrix b);//Cholesky
Matrix Solve_LDL(Matrix A,Matrix b);//Cholesky+
```

`Ecercise.h` 中分别构造矩阵 `A` 和 `b` 并进行求解：

```cpp
void LU(Matrix A,Matrix b){//三种LU分解法
	...
    Solve_LU(A,b);
    ...
    Solve_FPLU(A,b);
    ...
    Solve_CPLU(A,b);
	...
}
void Cholesky(Matrix A,Matrix b){//两种Cholesky解法
	...
    Solve_LL(A,b);
    ...
    Solve_LDL(A,b);
    ...
}
```

<div STYLE="page-break-after: always;"></div>
## **三.【实验结果】**

### **1.三种LU分解法比较**

对于同一方程组的求解，LU分解法运行速度最快，但精度较低。全主元LU分解法和列主元LU分解法精度相差无几，但列主元LU分解法运行速度相对更快（与LU分解法相差不大）。

可见列主元LU分解法为综合较优的选择。若精度要求不高、时限较严，可采用LU分解法。

![](./_1.png)

<div STYLE="page-break-after: always;"></div>
### **2.两种Cholesky分解比较**

改进版精度略优（`sqrt` 函数精度差），但会消耗更多的时间。

对于部分方程组，未改进版无法求解，可能是精度问题所致。

综合下来选择改进版更优。

![](./_2.png)

![](./_4.png)

<div STYLE="page-break-after: always;"></div>
### **3.两类分解法比较**

改进版 Cholesky 分解法在时间消耗和精度上都略逊于列主元 LU 分解法。

综合考量下来，列主元 LU 分解法为最佳选择。

![](./_3.png)

![](./_5.png)