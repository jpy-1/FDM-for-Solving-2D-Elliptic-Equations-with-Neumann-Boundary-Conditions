# FDM-for-Solving-2D-Elliptic-Equations-with-Neumann-Boundary-Conditions
使用有限差分法中的Crank-Nicolson格式，数值求解以下二维抛物型方程：

## 求解问题
### 方程与求解区域
求解二维区域 $(x, y) \in[0,1] \times[0,1]$ 上的抛物型方程：
$$u_{t}-\nabla(p(\vec{x}) \nabla u)+q(\vec{x}) u=f$$

### 系数与初边值条件
- 系数函数：
  $$p(\vec{x})=2 x$$，
  $$q(\vec{x})=3 y^{2}$$
- 初始条件：
  $u(x, y, 0)=0$
- 边界条件：
  1. Dirichlet边界：
     $u(0, y, t)=u(x, 0, t)=0$
  2. Neumann边界：
     $u_{x}(1, y, t)=u_{y}(x, 1, t)=0$
- 解析解取为：
$u(x,y,t) = t^2 * sin(πx/2) * sin(πy/2)$

右端源项f由解析解代入原方程推导得到。

---

## 文件说明


| 文件名 | 核心功能 |
|--------|----------|
| fdm_class.m | 定义求解问题基本参数 |
| fun_fdm.m | 有限差分核心求解函数 |
| error_compute.m | 误差计算函数 |
| converge_order.m | 收敛阶计算程序 |
| fdm_plot.m | 绘图程序|


---

## 实现原理与结论
见有限差分程序实验报告.pdf




