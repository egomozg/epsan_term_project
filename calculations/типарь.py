import numpy as np
from sympy import *
import matplotlib.pyplot as plt
init_printing()

k0u, k1u = symbols('k0u k1u')
Te = 1
Td0 = 14.4
Td = 0.6
xd = 2.5
Tu1 = 0.9
p, w = symbols('p w')

W_sw = 1 / (Te * p + 1)
W_g = 1 / (Td0 * p + 1)
W_s = xd / (Td * p + 1)
W_gvv = W_s / (W_sw * W_g)


pprint(expand(factor(W_gvv)))

W_raz = k0u * W_sw * W_g
W = W_raz / (W_raz + 1)

pprint(W_raz)
pprint(simplify(W).expand(denom=True))

k0u = 60
k1u = 20

W_raz = (k0u * (Tu1*p+1) + k1u*p)/((Tu1*p+1) * (Te*p+1) * (Td0*p+1))
w_1 = solveset(Eq(Abs(W_raz.subs(p, 1j*w)), 1), w, domain=S.Reals)
print(max(w_1))
print('угол при r=1 = ', arg(W_raz.subs(p, 1j*max(w_1))).evalf(4)* 180/ pi.evalf(4))
print('угол запаса = ', arg(W_raz.subs(p, 1j*max(w_1))).evalf(3) * 180 / pi.evalf(3) + 180)
W_raz2 = lambdify(p, W_raz, "numpy")
w = np.linspace(0,1000,1000000)
W_raz_num_real = np.zeros((5, np.size(w)))
W_raz_num_imag = np.zeros((5, np.size(w)))
i = 0
for k1u in range(0, 25, 5):
    W_raz = (k0u * (Tu1*p+1) + k1u*p)/((Tu1*p+1) * (Te*p+1) * (Td0*p+1))
    W_raz2 = lambdify(p, W_raz, "numpy")
    W_raz_num_real[i, :] = np.real(W_raz2(1j * w))
    W_raz_num_imag[i, :] = np.imag(W_raz2(1j * w))
    i += 1

#print(W_raz_num_real[1, :])

#print(solveset(Eq(Abs(W_raz.subs(p, 1j*w)), 1), w, domain=S.Reals))
#print(W_raz.evalf(subs={p: 1j*1e30}))

phi = np.linspace(0, 2*np.pi, 200)
r = 1
x = r*np.cos(phi)
y = r*np.sin(phi)

#for i in range(0, 5):
#    plt.plot(W_raz_num_real[i, :], W_raz_num_imag[i, :])

plt.plot(W_raz_num_real[0, :], W_raz_num_imag[0, :], label='k1u=0')
plt.plot(W_raz_num_real[1, :], W_raz_num_imag[1, :], label='k1u=5')
plt.plot(W_raz_num_real[2, :], W_raz_num_imag[2, :], label='k1u=10')
#plt.plot(W_raz_num_real[3, :], W_raz_num_imag[3, :], label='k1u=15')
plt.plot(W_raz_num_real[4, :], W_raz_num_imag[4, :], '--', label='k1u=20')
plt.plot(x,y, '-.')
plt.plot(-1, 0, 'ro')
plt.annotate('(-1; 0j)', xy=(-1, 0), xytext=(-2, 0.25),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )

plt.annotate(r'$\omega = \infty$', xy=(0, 0), xytext=(-0.3, 0.1)
             #arrowprops=dict(facecolor='black', shrink=0.05),
             )

plt.grid()
plt.title(r'Амплитудно-фазовая характеристика при $k_{U0} = %i$' % k0u)
plt.xlim(-3, 2)
plt.ylim(-2, 1.5)
plt.xlabel(r'$Re(W_{раз})$')
plt.ylabel(r'$Im(W_{раз})$')
plt.legend(loc='upper right')
#plt.axis("equal")
plt.show()
