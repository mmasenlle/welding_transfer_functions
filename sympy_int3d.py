from sympy import *
# x,y,z,x1,y1,z1,b,c=symbols('x y z x1 y1 z1 b c',real=True)
# dT=(c*x**2-x+c*y**2+c*z**2-c*x*sqrt(x**2+y**2+z**2))/(sqrt(x**2+y**2+z**2)**3)
# G=exp(c*(x1-x)-((x1-x)**2+(y1-y)**2)/b)*(exp(-(z1+z)**2/b)+exp(-(z1-z)**2/b))
# f=G*dT
# #init_printing(use_unicode=True)
# print('f=',f)
# F=integrate(f,(x,-oo,oo),(y,-oo,oo),(z,0,oo))
# print('F=',F)
#
# # Fx = ((e^(-(z_1+z)^2/b)+e^(-(z_1-z)^2/b))*(-c*x*sqrt(x^2+z^2+y^2)+c*x^2-x+c*z^2+c*y^2)*e^(c*(x_1-x)-((x_1-x)^2+(y_1-y)^2)/b))/(x^2+z^2+y^2)^(3/2)