from sympy import *

init_printing(use_unicode=True)

s,x,y,z = symbols('s x y z', real=True)
a,q,v,K = symbols('a q v K', real=True, positive=True)

T_0xyz = q*exp(v*(x-sqrt(x**2+y**2+z**2))/(2*a))/(2*pi*a*sqrt(x**2+y**2+z**2))
dT_0xyz = simplify(diff(T_0xyz,v))
fgood = lambdify((a,q,v,x,y,z), dT_0xyz)
# ftest(6.6e-06,1000,.002,.3,.2,.1)

values = ((6.6e-06,1000,.002,.3,.2,.1),(6.6e-06,3000,.002,.3,.2,.1),(6.6e-06,1000,.003,.3,.2,.1),(6.6e-06,1000,.002,.3,.0,.1),(6.6e-06,1000,.003,.2,.0,.0))
def check(ftest):
    res = 0
    for v in values:
        res += (ftest(*v)-fgood(*v))**2
    return res


class BIdx:
    def __init__(self):
        self.i = -1
    def b(self,j):
        self.i += 1
        return not (j & (1<<self.i))


R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)

def create(idx):
    b = BIdx()
    func = 1

    for e in (1/q,4,pi,a,a,s,K,R):
        if (b.b(idx)):
            func = func / e

    sums = 0
    sum1 = 1
    (v * R - x * D) * exp(-R * D / (2 * a)) / D
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum1 = sum1 * e
    sums = sums + sum1
    sum2 = 1
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum2 = sum2 * e
    sums = sums + sum2
    sum3 = 1
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum3 = sum3 * e
    sums = sums + sum3
    sum4 = 1
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum4 = sum4 * e
    sums = sums + sum4
    func = func * sums

    exp1 = 1
    for e in (v, x, 1/(a*2)):
        if (b.b(idx)):
            exp1 = exp1 * e
    func = func * exp(exp1)

    exp2 = 1
    for e in (v, -R, 1/(a*2)):
        if (b.b(idx)):
            exp2 = exp2 * e
    func = func * exp(exp2)

    # print('func:',func)
    return func

check(fgood)

goodness = []
mingn = 1e20
mingni = -1
for i in range(2):
    ffs = create(i)
    ff = simplify(limit(ffs, s, 0))
    ff_ = simplify(ff - dT_0xyz)
    print('************',i,'*************')
    print('ffs:', ffs)
    print('ff:', ff)
    print('ff_:',ff_)
    if (ff.equals(dT_0xyz)):
        print('equals !!!!!!!!')
        # break
    ffl = lambdify((a, q, v, x, y, z), ff)
    gn = check(ffl)
    print('gn:', gn)
    goodness.append(gn)
    if (gn < mingn):
        mingn = gn
        mingni = i

print('migni:',mingni,mingn)



R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)

def create2(idx):
    b = BIdx()
    func = 1

    for e in (1/q,4,pi,a,a,s,K,R):
        if (b.b(idx)):
            func = func / e

    sums = 0
    sum1 = 1
    (v * R - x * D) * exp(-R * D / (2 * a)) / D
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum1 = sum1 * e
    sums = sums + sum1
    sum2 = 1
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum2 = sum2 * e
    sums = sums + sum2
    sum3 = 1
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum3 = sum3 * e
    sums = sums + sum3
    sum4 = 1
    for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
        if (b.b(idx)):
            sum4 = sum4 * e
    sums = sums + sum4
    func = func * sums

    exp1 = 1
    for e in (v, x, 1/(a*2)):
        if (b.b(idx)):
            exp1 = exp1 * e
    func = func * exp(exp1)

    print('func:',func)
    return func


def create3(idx):
    b = BIdx()
    func = q*exp(v*x/(2*a))/(4*pi*s*a*a*R)

    sums = 0
    for i in range(4):
        sum1 = 1
        for e in (-1,x,v,R,D,exp(-R * D / (2 * a)), exp(-R * v / (2 * a)), 1/D, 1/v):
            if (b.b(idx)):
                sum1 = sum1 * e
        sums = sums + sum1
    func = func * sums

    return func

ffs = create3(513)
ffs = create3(1557)
ffs = create3(1049093)
ff = simplify(limit(ffs, s, 0))
ff_ = simplify(ff - dT_0xyz)
ff.equals(dT_0xyz)

# test 1
T_0xyzs = q*exp(-(v*x+sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*pi*a*sqrt(x**2+y**2+z**2))
dT_0xyzs = simplify(diff(T_0xyzs,v))
Tv0xyz1 = simplify(limit(dT_0xyzs, s, 0))
f1 = lambdify((a,q,v,x,y,z), Tv0xyz1)
check(f1)

# test 2
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
XX = (v*R - x*D)*exp(-R*D/(2*a))/D
XX0 = XX.subs(s,0)
ft_vel2 = -q*exp(v*x/(2*a))/(4*s*pi*a**2*R)*(XX-XX0)
Tv0xyz2 = simplify(limit(ft_vel2, s, 0))
ff_ = simplify(Tv0xyz2 - dT_0xyz)
Tv0xyz2.equals(dT_0xyz)
f2 = lambdify((a,q,v,x,y,z,K), Tv0xyz2)
check(f2)

# test 3
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
XX = exp(-R*D/(2*a))
XX0 = XX.subs(s,0)
ft_vel3 = exp(v*x/(2*a))/(2*s*pi*a*24*R)*(XX-XX0)
# latex(ft_vel3)
Tv0xyz3 = simplify(limit(ft_vel3, s, 0))
ff_ = simplify(Tv0xyz3 - dT_0xyz)
Tv0xyz3.equals(dT_0xyz)
f3 = lambdify((a,q,v,x,y,z,K), Tv0xyz3)
check(f3)

# test 4
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
ft_vel4 = exp(v*x/(2*a))/(2*s*pi*a*K*R)*(exp(-R*v/(2*a))-v*exp(-R*v/(2*a))/v-exp(-R*D/(2*a))+v*exp(-R*D/(2*a))/D)
# latex(ft_vel3)
Tv0xyz4 = simplify(limit(ft_vel4, s, 0))
ff_ = simplify(Tv0xyz4 - dT_0xyz)
Tv0xyz4.equals(dT_0xyz)
f4 = lambdify((a,q,v,x,y,z,K), Tv0xyz4)
check(f4)


# test 5
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
ft_vel5 = -q*(x-R)*v*v*exp(v*x/(2*a))/(8*s*pi*a*a*a*R)*(exp(-R*v/(2*a))-v*exp(-R*v/(2*a))/v-exp(-R*D/(2*a))+v*exp(-R*D/(2*a))/D)
# latex(ft_vel3)
Tv0xyz5 = simplify(limit(ft_vel5, s, 0))
ff_ = simplify(Tv0xyz5 - dT_0xyz)
Tv0xyz5.equals(dT_0xyz)
f5 = lambdify((a,q,v,x,y,z,K), Tv0xyz5)
check(f5)

# test 6
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
ft_vel6 = -q*(x-R)*v*v*exp(v*x/(2*a))/(8*s*pi*a*a*a*R)*(-exp(-R*D/(2*a))+v*exp(-R*D/(2*a))/D)
nft_vel6 = -exp(-R*D/(2*a))+v*exp(-R*D/(2*a))/D
nft_vel6.subs(s,0)
dnft_vel6 = simplify(diff(nft_vel6,s))
dnft_vel6.subs(s,0)

Tv0xyz6 = simplify(limit(ft_vel6, s, 0))
ff_ = simplify(Tv0xyz6 - dT_0xyz)
Tv0xyz6.equals(dT_0xyz)
f6 = lambdify((a,q,v,x,y,z,K), Tv0xyz6)
check(f6)

# test 7
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
ft_vel7 = -q*(x-R)*v*exp(v*x/(2*a))/(4*s*pi*a*a*R*R)*(exp(-R*D/(2*a))-exp(-R*v/(2*a)))
nft_vel7 = exp(-R*D/(2*a))-exp(-R*v/(2*a))
nft_vel7.subs(s,0)
dnft_vel7 = simplify(diff(nft_vel7,s))
dnft_vel7_2 = simplify(diff(exp(-R*D/(2*a)),s))
dnft_vel7.subs(s,0)

Tv0xyz7 = simplify(limit(ft_vel7, s, 0))
ff_ = simplify(Tv0xyz7 - dT_0xyz)
Tv0xyz7.equals(dT_0xyz)
f7 = lambdify((a,q,v,x,y,z,K), Tv0xyz7)
check(f7)


# test 8
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
ft_vel8 = q*v*v*exp(v*x/(2*a))/(4*s*pi*a*a*R)*(-exp(-R*D/(2*a))/sqrt(-2*s*x*v + v**2)+exp(-R*v/(2*a))/v)
nft_vel8 = exp(R*D/(2*a))/sqrt(2*s*x*v + v**2)-exp(R*v/(2*a))/v
simplify(nft_vel8.subs(s,0))
dnft_vel8 = simplify(diff(nft_vel8,s))
dnft_vel8.subs(s,0)

df=diff(v*v*exp((-R)*D/(2*a))/sqrt(-4*s*x*x*v + v**2)+v*exp((-R)*D/(2*a)),s)
df.subs(s,0)
simplify(diff(v*exp((-R)*D/(2*a))/D,s))
simplify(diff(-R*exp(-R*D/(2*a)),s))

Tv0xyz8 = simplify(limit(ft_vel8, s, 0))
ff_ = simplify(Tv0xyz8 - dT_0xyz)
Tv0xyz8.equals(dT_0xyz)


# test 9
R = sqrt(x**2+y**2+z**2)
D = sqrt(4*a*s + v**2)
ft_vel9 = q*v*exp(v*x/(2*a))/(8*s*pi*a*a*R)*(v*exp(-R*D/(2*a))/sqrt(4*s*x*v + v**2)-exp(-R*v/(2*a))+exp(-R*D/(2*a))-exp(-R*v/(2*a)))
nft_vel9 = v*exp(-R*D/(2*a))/sqrt(-2*s*x*v + v**2)-exp(-R*v/(2*a)+exp((-R)*D/(2*a))-exp((-R)*v/(2*a)))
simplify(nft_vel9.subs(s,0))
dnft_vel9 = simplify(diff(nft_vel9,s))
dnft_vel9.subs(s,0)

Tv0xyz9 = simplify(limit(ft_vel9, s, 0))
ff_ = simplify(Tv0xyz9 - dT_0xyz)
Tv0xyz9.equals(dT_0xyz)