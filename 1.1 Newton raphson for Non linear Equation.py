import numpy as np 
import sympy as smp 
from sympy import diff,solve,symbols

# Newton Non- linear
x=symbols('x')
y=symbols('y')
z = symbols('z')

u= 3*x-smp.cos(y*z)-(1/2) 
v= x**2-81*(y+0.1)**2+smp.sin(z)+1.06
w = np.e**(-x*y)+20*z+((10*np.pi-3)/3)

# x = 2
# y = 1
# n = 5
# x=int(input('enter the initial point x'))
# y=int(input('enter the initial point y'))
# n=int(input('enter number of iterations'))

dudx= smp.lambdify([x,y,z],diff(u,x))
dudy=smp.lambdify([x,y,z],diff(u,y))
dudz =smp.lambdify([x,y,z],diff(u,z)) 
dvdx=smp.lambdify([x,y,z],diff(v,x))
dvdy=smp.lambdify([x,y,z],diff(v,y))
dvdz =smp.lambdify([x,y,z],diff(v,z)) 
dwdx=smp.lambdify([x,y,z],diff(w,x))
dwdy=smp.lambdify([x,y,z],diff(w,y))
dwdz =smp.lambdify([x,y,z],diff(w,z))

p = 0.1
q = 0.1
r = -0.1
n = 9
# a = np.array([p,q,r])
# print(a)
u_val=smp.lambdify([x,y,z],u)
v_val=smp.lambdify([x,y,z],v)
w_val = smp.lambdify([x,y,z],w)
i=0
while i<n:
    J= np.array([[dudx(p,q,r),dudy(p,q,r),dudz(p,q,r)],
                  [dvdx(p,q,r),dvdy(p,q,r),dvdz(p,q,r)],
                  [dwdx(p,q,r),dwdy(p,q,r),dwdz(p,q,r)]])
    # print(J)
    P=np.array([u_val(p,q,r),v_val(p,q,r),w_val(p,q,r)])
    inv_J=np.linalg.inv(J)
    ans=np.dot(inv_J,-P)
    # ans = np.array(ans)
    # print(ans)
    p+=ans[0]
    q+=ans[1]
    r+=ans[2]
    a = np.array([p,q,r])
    i+=1
    m = np.round(a+ans,4)
    print(m)
    