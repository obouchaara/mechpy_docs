---
title: "Example - Material"
author: omar
date: 2023-12-08 10:00:00 +0000
categories: example
pin: false
---

```python
import sys
import numpy as np
import sympy as sp

sys.path.append("../src")

from core.numeric.material import IsotropicMaterial, TransverseIsotropicMaterial
from core.symbolic.material import SymbolicIsotropicMaterial, SymbolicTransverseIsotropicMaterial,SymbolicOrthotropicMaterial

np.set_printoptions(formatter={"float": "{:0.2e}".format})
```


```python
isotropic_material = IsotropicMaterial(210e9, 0.3)
stiffness_tensor = isotropic_material.stiffness_tensor()
compliance_tensor = isotropic_material.compliance_tensor()
display(stiffness_tensor)
display(compliance_tensor)
```


    StiffnessTensor(
    [[2.83e+11 1.21e+11 1.21e+11 0.00e+00 0.00e+00 0.00e+00]
     [1.21e+11 2.83e+11 1.21e+11 0.00e+00 0.00e+00 0.00e+00]
     [1.21e+11 1.21e+11 2.83e+11 0.00e+00 0.00e+00 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 8.08e+10 0.00e+00 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 0.00e+00 8.08e+10 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 8.08e+10]]
    )



    ComplianceTensor(
    [[4.76e-12 -1.43e-12 -1.43e-12 0.00e+00 0.00e+00 0.00e+00]
     [-1.43e-12 4.76e-12 -1.43e-12 0.00e+00 0.00e+00 0.00e+00]
     [-1.43e-12 -1.43e-12 4.76e-12 0.00e+00 0.00e+00 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 1.24e-11 0.00e+00 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.24e-11 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.24e-11]]
    )



```python
transverse_isotropic_material = TransverseIsotropicMaterial(
    150e12, 10e12, 0.3, 70e12, 5e12
)
stiffness_tensor = transverse_isotropic_material.stiffness_tensor()
compliance_tensor = transverse_isotropic_material.compliance_tensor()
display(stiffness_tensor)
display(compliance_tensor)
```


    StiffnessTensor(
    [[1.65e+14 6.43e+13 6.43e+13 0.00e+00 0.00e+00 0.00e+00]
     [6.43e+13 1.65e+14 6.43e+13 0.00e+00 0.00e+00 0.00e+00]
     [6.43e+13 6.43e+13 1.00e+13 0.00e+00 0.00e+00 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 7.00e+13 0.00e+00 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 0.00e+00 7.00e+13 0.00e+00]
     [0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 5.00e+12]]
    )



    None



```python
E, nu = sp.symbols("E nu")
lamda, mu = sp.symbols("lamda mu", cls=sp.Function)
lamda = (E * nu) / ((1 + nu) * (1 - 2 * nu))
mu = E / (2 * (1 + nu))
display(lamda, mu)
E, nu = sp.symbols("E nu", cls=sp.Function)
E = mu * (3 * lamda + 2 * mu) / (lamda + mu)
nu = lamda / (2 * (lamda + mu))
display(E, nu)
display(sp.simplify((E)), sp.simplify(nu))
```


$\displaystyle \frac{E \nu}{\left(1 - 2 \nu\right) \left(\nu + 1\right)}$



$\displaystyle \frac{E}{2 \nu + 2}$



$\displaystyle \frac{E \left(\frac{3 E \nu}{\left(1 - 2 \nu\right) \left(\nu + 1\right)} + \frac{2 E}{2 \nu + 2}\right)}{\left(2 \nu + 2\right) \left(\frac{E \nu}{\left(1 - 2 \nu\right) \left(\nu + 1\right)} + \frac{E}{2 \nu + 2}\right)}$



$\displaystyle \frac{E \nu}{\left(1 - 2 \nu\right) \left(\nu + 1\right) \left(\frac{2 E \nu}{\left(1 - 2 \nu\right) \left(\nu + 1\right)} + \frac{2 E}{2 \nu + 2}\right)}$



$\displaystyle E$



$\displaystyle \nu$



```python
symbolic_isotropic_material = SymbolicIsotropicMaterial()
display(symbolic_isotropic_material)
display(symbolic_isotropic_material.stiffness_tensor().data)
display(symbolic_isotropic_material.compliance_tensor().data)
```


    SymbolicIsotropicMaterial(E, nu)



$\displaystyle \left[\begin{matrix}\lambda + 2 \mu & \lambda & \lambda & 0 & 0 & 0\\\lambda & \lambda + 2 \mu & \lambda & 0 & 0 & 0\\\lambda & \lambda & \lambda + 2 \mu & 0 & 0 & 0\\0 & 0 & 0 & \mu & 0 & 0\\0 & 0 & 0 & 0 & \mu & 0\\0 & 0 & 0 & 0 & 0 & \mu\end{matrix}\right]$



$\displaystyle \left[\begin{matrix}\frac{1}{E} & - \frac{\nu}{E} & - \frac{\nu}{E} & 0 & 0 & 0\\- \frac{\nu}{E} & \frac{1}{E} & - \frac{\nu}{E} & 0 & 0 & 0\\- \frac{\nu}{E} & - \frac{\nu}{E} & \frac{1}{E} & 0 & 0 & 0\\0 & 0 & 0 & \frac{2 \left(\nu + 1\right)}{E} & 0 & 0\\0 & 0 & 0 & 0 & \frac{2 \left(\nu + 1\right)}{E} & 0\\0 & 0 & 0 & 0 & 0 & \frac{2 \left(\nu + 1\right)}{E}\end{matrix}\right]$



```python
E_a, nu_b = sp.symbols("E_a nu_b")
symbolic_isotropic_material = SymbolicIsotropicMaterial(E_a, nu_b)
display(symbolic_isotropic_material)
display(symbolic_isotropic_material.stiffness_tensor().data)
display(symbolic_isotropic_material.stiffness_tensor(lames_param=False).data)
display(symbolic_isotropic_material.compliance_tensor().data)
```


    SymbolicIsotropicMaterial(E_a, nu_b)



$\displaystyle \left[\begin{matrix}\lambda + 2 \mu & \lambda & \lambda & 0 & 0 & 0\\\lambda & \lambda + 2 \mu & \lambda & 0 & 0 & 0\\\lambda & \lambda & \lambda + 2 \mu & 0 & 0 & 0\\0 & 0 & 0 & \mu & 0 & 0\\0 & 0 & 0 & 0 & \mu & 0\\0 & 0 & 0 & 0 & 0 & \mu\end{matrix}\right]$



$\displaystyle \left[\begin{matrix}\frac{E_{a} \left(\nu_{b} - 1\right)}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & - \frac{E_{a} \nu_{b}}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & - \frac{E_{a} \nu_{b}}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & 0 & 0 & 0\\- \frac{E_{a} \nu_{b}}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & \frac{E_{a} \left(\nu_{b} - 1\right)}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & - \frac{E_{a} \nu_{b}}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & 0 & 0 & 0\\- \frac{E_{a} \nu_{b}}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & - \frac{E_{a} \nu_{b}}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & \frac{E_{a} \left(\nu_{b} - 1\right)}{\left(\nu_{b} + 1\right) \left(2 \nu_{b} - 1\right)} & 0 & 0 & 0\\0 & 0 & 0 & \frac{E_{a}}{2 \left(\nu_{b} + 1\right)} & 0 & 0\\0 & 0 & 0 & 0 & \frac{E_{a}}{2 \left(\nu_{b} + 1\right)} & 0\\0 & 0 & 0 & 0 & 0 & \frac{E_{a}}{2 \left(\nu_{b} + 1\right)}\end{matrix}\right]$



$\displaystyle \left[\begin{matrix}\frac{1}{E_{a}} & - \frac{\nu_{b}}{E_{a}} & - \frac{\nu_{b}}{E_{a}} & 0 & 0 & 0\\- \frac{\nu_{b}}{E_{a}} & \frac{1}{E_{a}} & - \frac{\nu_{b}}{E_{a}} & 0 & 0 & 0\\- \frac{\nu_{b}}{E_{a}} & - \frac{\nu_{b}}{E_{a}} & \frac{1}{E_{a}} & 0 & 0 & 0\\0 & 0 & 0 & \frac{2 \left(\nu_{b} + 1\right)}{E_{a}} & 0 & 0\\0 & 0 & 0 & 0 & \frac{2 \left(\nu_{b} + 1\right)}{E_{a}} & 0\\0 & 0 & 0 & 0 & 0 & \frac{2 \left(\nu_{b} + 1\right)}{E_{a}}\end{matrix}\right]$



```python
lamda, mu = sp.symbols("lamda mu")
symbolic_isotropic_material = SymbolicIsotropicMaterial(lamda=lamda, mu=mu)
display(symbolic_isotropic_material)
display(symbolic_isotropic_material.stiffness_tensor().data)
display(symbolic_isotropic_material.compliance_tensor().data)

```


    SymbolicIsotropicMaterial(mu*(3*lamda + 2*mu)/(lamda + mu), lamda/(2*lamda + 2*mu))



$\displaystyle \left[\begin{matrix}\lambda + 2 \mu & \lambda & \lambda & 0 & 0 & 0\\\lambda & \lambda + 2 \mu & \lambda & 0 & 0 & 0\\\lambda & \lambda & \lambda + 2 \mu & 0 & 0 & 0\\0 & 0 & 0 & \mu & 0 & 0\\0 & 0 & 0 & 0 & \mu & 0\\0 & 0 & 0 & 0 & 0 & \mu\end{matrix}\right]$



$\displaystyle \left[\begin{matrix}\frac{\lambda + \mu}{\mu \left(3 \lambda + 2 \mu\right)} & - \frac{\lambda}{2 \mu \left(3 \lambda + 2 \mu\right)} & - \frac{\lambda}{2 \mu \left(3 \lambda + 2 \mu\right)} & 0 & 0 & 0\\- \frac{\lambda}{2 \mu \left(3 \lambda + 2 \mu\right)} & \frac{\lambda + \mu}{\mu \left(3 \lambda + 2 \mu\right)} & - \frac{\lambda}{2 \mu \left(3 \lambda + 2 \mu\right)} & 0 & 0 & 0\\- \frac{\lambda}{2 \mu \left(3 \lambda + 2 \mu\right)} & - \frac{\lambda}{2 \mu \left(3 \lambda + 2 \mu\right)} & \frac{\lambda + \mu}{\mu \left(3 \lambda + 2 \mu\right)} & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{\mu} & 0 & 0\\0 & 0 & 0 & 0 & \frac{1}{\mu} & 0\\0 & 0 & 0 & 0 & 0 & \frac{1}{\mu}\end{matrix}\right]$



```python
lamda_a, mu_b = sp.symbols("lamda_a mu_b")
symbolic_isotropic_material = SymbolicIsotropicMaterial(lamda=lamda_a, mu=mu_b)
display(symbolic_isotropic_material)
display(symbolic_isotropic_material.stiffness_tensor().data)
display(symbolic_isotropic_material.compliance_tensor().data)
```


    SymbolicIsotropicMaterial(mu_b*(3*lamda_a + 2*mu_b)/(lamda_a + mu_b), lamda_a/(2*lamda_a + 2*mu_b))



$\displaystyle \left[\begin{matrix}\lambda_{a} + 2 \mu_{b} & \lambda_{a} & \lambda_{a} & 0 & 0 & 0\\\lambda_{a} & \lambda_{a} + 2 \mu_{b} & \lambda_{a} & 0 & 0 & 0\\\lambda_{a} & \lambda_{a} & \lambda_{a} + 2 \mu_{b} & 0 & 0 & 0\\0 & 0 & 0 & \mu_{b} & 0 & 0\\0 & 0 & 0 & 0 & \mu_{b} & 0\\0 & 0 & 0 & 0 & 0 & \mu_{b}\end{matrix}\right]$



$\displaystyle \left[\begin{matrix}\frac{\lambda_{a} + \mu_{b}}{\mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & - \frac{\lambda_{a}}{2 \mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & - \frac{\lambda_{a}}{2 \mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & 0 & 0 & 0\\- \frac{\lambda_{a}}{2 \mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & \frac{\lambda_{a} + \mu_{b}}{\mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & - \frac{\lambda_{a}}{2 \mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & 0 & 0 & 0\\- \frac{\lambda_{a}}{2 \mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & - \frac{\lambda_{a}}{2 \mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & \frac{\lambda_{a} + \mu_{b}}{\mu_{b} \left(3 \lambda_{a} + 2 \mu_{b}\right)} & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{\mu_{b}} & 0 & 0\\0 & 0 & 0 & 0 & \frac{1}{\mu_{b}} & 0\\0 & 0 & 0 & 0 & 0 & \frac{1}{\mu_{b}}\end{matrix}\right]$



```python
transverse_isotropic_material = SymbolicTransverseIsotropicMaterial()
display(transverse_isotropic_material.stiffness_tensor().data)
# display(transverse_isotropic_material.compliance_tensor().data)
```


$\displaystyle \left[\begin{matrix}- \frac{E_{L}}{\nu^{2} - 1} & - \frac{E_{L} \nu}{\nu - 1} & - \frac{E_{L} \nu}{\nu - 1} & 0 & 0 & 0\\- \frac{E_{L} \nu}{\nu - 1} & - \frac{E_{L}}{\nu^{2} - 1} & - \frac{E_{L} \nu}{\nu - 1} & 0 & 0 & 0\\- \frac{E_{L} \nu}{\nu - 1} & - \frac{E_{L} \nu}{\nu - 1} & E_{T} & 0 & 0 & 0\\0 & 0 & 0 & G_{L} & 0 & 0\\0 & 0 & 0 & 0 & G_{L} & 0\\0 & 0 & 0 & 0 & 0 & G_{T}\end{matrix}\right]$



```python
E1 = sp.symbols("E1")
E2 = E1 * sp.symbols("n")
E3 = E2
G12 = sp.symbols("G12")
G23 = G12 * sp.symbols("m")
G31 = G12 * sp.symbols("p")
orthotropic_material = SymbolicOrthotropicMaterial(E1=E1, E2=E2, E3=E3, G12=G12, G23=G23, G31=G31, nu12=0.3, nu23=0.3, nu31=0.3)
display(orthotropic_material)
display(orthotropic_material.mechanical_props)
display(orthotropic_material.stiffness_tensor().data)
display(orthotropic_material.compliance_tensor().data)
```


    SymbolicOrthotropicMaterial({'E1': E1, 'E2': E1*n, 'E3': E1*n, 'G12': G12, 'G23': G12*m, 'G31': G12*p, 'nu12': 0.3, 'nu23': 0.3, 'nu31': 0.3})



    {'E1': E1,
     'E2': E1*n,
     'E3': E1*n,
     'G12': G12,
     'G23': G12*m,
     'G31': G12*p,
     'nu12': 0.3,
     'nu23': 0.3,
     'nu31': 0.3}



$\displaystyle \left[\begin{matrix}E_{1} & 0.3 E_{1} n & 0.3 E_{1} & 0 & 0 & 0\\0.3 E_{1} n & E_{1} n & 0.3 E_{1} n & 0 & 0 & 0\\0.3 E_{1} & 0.3 E_{1} n & E_{1} n & 0 & 0 & 0\\0 & 0 & 0 & G_{12} m & 0 & 0\\0 & 0 & 0 & 0 & G_{12} p & 0\\0 & 0 & 0 & 0 & 0 & G_{12}\end{matrix}\right]$



$\displaystyle \left[\begin{matrix}\frac{n \left(0.91 - 0.0819 n\right)}{E_{1} \cdot \left(0.0081 n^{3} - 0.17676 n^{2} + 0.9721 n - 0.09\right)} & \frac{0.027 n^{2} - 0.3081 n + 0.09}{E_{1} \cdot \left(0.0081 n^{3} - 0.17676 n^{2} + 0.9721 n - 0.09\right)} & \frac{- 0.0081 n^{2} + 0.117 n - 0.3}{E_{1} \cdot \left(0.0081 n^{3} - 0.17676 n^{2} + 0.9721 n - 0.09\right)} & 0 & 0 & 0\\\frac{0.027 n^{2} - 0.3081 n + 0.09}{E_{1} \cdot \left(0.0081 n^{3} - 0.17676 n^{2} + 0.9721 n - 0.09\right)} & \frac{- 0.09 n^{2} + 1.0081 n - 0.09}{E_{1} n \left(0.0081 n^{3} - 0.17676 n^{2} + 0.9721 n - 0.09\right)} & \frac{0.0189 n - 0.21}{E_{1} \cdot \left(0.0081 n^{3} - 0.17676 n^{2} + 0.9721 n - 0.09\right)} & 0 & 0 & 0\\\frac{0.3 - 0.09 n}{E_{1} \cdot \left(0.09 n^{2} - 0.964 n + 0.09\right)} & \frac{0.21}{E_{1} \cdot \left(0.09 n^{2} - 0.964 n + 0.09\right)} & \frac{0.09 n - 1.0}{E_{1} \cdot \left(0.09 n^{2} - 0.964 n + 0.09\right)} & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{G_{12} m} & 0 & 0\\0 & 0 & 0 & 0 & \frac{1}{G_{12} p} & 0\\0 & 0 & 0 & 0 & 0 & \frac{1}{G_{12}}\end{matrix}\right]$



```python

```
