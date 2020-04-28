from sympy import *

init_printing(use_unicode=True)

a,b,q,s,v,x,y,z,K,K0 = symbols('a b q s v x y z K K0', real=True)
T_0xy = exp(-v*x/(2*a))/(2*pi*K)*K0*(1/(2*a)*sqrt((x**2+y**2)*(4*a*b+v**2)))*q

dT_0xy = simplify(diff(T_0xy,v))
latex(dT_0xy)

###################### 1d
T_0x = a*exp(-(v*x + sqrt(x**2*(4*a*b + v**2))/(2*a)))/(K*sqrt(4*a*b+v**2))
dT_0x = simplify(diff(T_0x,v))

ftv = q*exp(v*x/(2*a))/(2*s*K)*(exp(-sqrt(x**2*(4*a*b+v**2))/(2*a))-(v*exp(-sqrt(x**2*(4*a*b+v**2))/(2*a)))/(sqrt(4*a*b + v**2))\
                                -exp(-sqrt(x**2*(4*a*s+4*a*b+v**2))/(2*a))+(v*exp(-sqrt(x**2*(4*a*s+4*a*b+v**2))/(2*a)))/(sqrt(4*a*s+4*a*b + v**2)))
ftv0 = simplify(limit(ftv, s, 0))
ftv1 = q*exp(v*x/(2*a))/(2*K)*(exp(-sqrt(x**2*(4*a*b+v**2))/(2*a))-(v*exp(-sqrt(x**2*(4*a*b+v**2))/(2*a)))/(sqrt(4*a*b + v**2))\
                                -exp(-sqrt(x**2*(4*a*s+4*a*b+v**2))/(2*a))+(v*exp(-sqrt(x**2*(4*a*s+4*a*b+v**2))/(2*a)))/(sqrt(4*a*s+4*a*b + v**2)))
dftv1 = diff(ftv1, s)
ftv0 = dftv1.subs(s, 0)
ftv0 = simplify(ftv0)

simplify(dT_0x - ftv0)
dT_0x.equals(ftv0)

###################### 3d
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
fdtTQ = exp(v*x/(2*a))/(2*pi*K*R)*exp(-(D*R)/(2*a))

T_0xyz = q*exp(v*(x-sqrt(x**2+y**2+z**2))/(2*a))/(2*pi*a*sqrt(x**2+y**2+z**2))
dT_0xyz = simplify(diff(T_0xyz,v))
latex(dT_0xyz)

a,q,v,x,y,z,s = symbols('a q v x y z s', real=True)
T_0xyzs = q*exp(-(v*x+sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*pi*a*sqrt(x**2+y**2+z**2))
dT_0xyzs = simplify(diff(T_0xyzs,v))
latex(dT_0xyzs)

Tv0xyz = limit(dT_0xyzs, s, 0)

R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
XX = (v*R - x*D)*exp(-R*D/(2*a))/D
XX0 = XX.subs(s,0)
ft_vel2 = -q*exp(v*x/(2*a))/(4*s*pi*a**2*R)*(XX-XX0)
latex(ft_vel2)
Tv0xyz2 = simplify(limit(ft_vel2, s, 0))

# tf => -q*(v*sqrt(x**2 + y**2 + z**2) + x*sqrt(4*a*s + v**2))*exp(-(v*x + sqrt(4*a*s
# + v**2)*sqrt(x**2 + y**2 + z**2))/(2*a))/(4*pi*a**2*sqrt(4*a*s + v**2)*sqrt(x*
# *2 + y**2 + z**2))

R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
XX = exp(-R*D/(2*a))
XX0 = XX.subs(s,0)
ft_vel3 = exp(v*x/(2*a))/(2*s*pi*a*K*R)*(XX-XX0)
latex(ft_vel3)
Tv0xyz2 = simplify(limit(ft_vel3, s, 0))

#### 2d with K0 as function
a,b,q,v,x,y,K = symbols('a b q v x y K', real=True)
K0=Function('K0')
T_0xy = exp(-v*x/(2*a))/(2*pi*K)*K0((1/(2*a)*sqrt((x**2+y**2)*(4*a*b+v**2))))*q

dT_0xy = simplify(diff(T_0xy,v))
latex(dT_0xy)

#### 2d with K0 as besselk
a,b,q,v,x,y,K = symbols('a b q v x y K', real=True)
T_0xy = exp(-v*x/(2*a))/(2*pi*K)*besselk(0,(1/(2*a)*sqrt((x**2+y**2)*(4*a*b+v**2))))*q

dT_0xy = simplify(diff(T_0xy,v))
latex(dT_0xy)