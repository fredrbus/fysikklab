# TFY41xx Fysikk vaaren 2021.
#
# Programmet tar utgangspunkt i hoeyden til de 8 festepunktene.
# Deretter beregnes baneformen y(x) ved hjelp av 7 tredjegradspolynomer,
# et for hvert intervall mellom to festepunkter, slik at baade banen y,
# dens stigningstall y' = dy/dx og dens andrederiverte
# y'' = d2y/dx2 er kontinuerlige i de 6 indre festepunktene.
# I tillegg velges null krumning (andrederivert)
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# De ulike banene er satt opp med tanke paa at kula skal
# (1) fullfoere hele banen selv om den taper noe mekanisk energi underveis;
# (2) rulle rent, uten aa gli ("slure").

# Importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

# Horisontal avstand mellom festepunktene er 0.200 m
#h = 0.200
# xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
# xfast = np.asarray([-2.373E-2, 0.175, 0.373, 0.569,
#                    0.768, 0.969, 1.171, 1.370])

f_av_x = open('tracked_data/y_fast.txt')
x_arr = []
y_arr = []
for line in f_av_x.readlines()[2:]:
    x, y = line[:-1].split('\t')
    x_arr.append(float(x))
    y_arr.append(float(y))
f_av_x.close()
x_0 = x_arr[0]
x_arr_adjusted = []
for x in x_arr:
    x_arr_adjusted.append(x - x_0)
xfast = np.asarray(x_arr_adjusted)
yfast = np.asarray(y_arr)


# Skriv inn y-verdiene til banens 8 festepunkter i tabellen yfast.
# Et vilkaarlig eksempel:
#yfast = np.asarray([0.391, 0.295, 0.272, 0.233, 0.130, 0.103, 0.182, 0.204])
# Erstatt med egne tallverdier avlest i tracker.
# Programmet beregner de 7 tredjegradspolynomene, et
# for hvert intervall mellom to festepunkter,
# med funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')
# Funksjonen cs kan naa brukes til aa regne ut y(x), y'(x) og y''(x)
# for en vilkaarlig horisontal posisjon x, eventuelt for mange horisontale
# posisjoner lagret i en tabell:
# cs(x)   tilsvarer y(x)
# cs(x,1) tilsvarer y'(x)
# cs(x,2) tilsvarer y''(x)
# Her lager vi en tabell med x-verdier mellom 0 og 1.4 m
xmin = 0.000
xmax = max(xfast)
dx = 0.001
x = np.arange(xmin, xmax, dx)
# Funksjonen arange returnerer verdier paa det "halvaapne" intervallet
# [xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
# dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400,
# dvs x blir en tabell med 1401 elementer
Nx = len(x)
y = cs(x)  # y=tabell med 1401 verdier for y(x)
dy = cs(x, 1)  # dy=tabell med 1401 verdier for y'(x)
d2y = cs(x, 2)  # d2y=tabell med 1401 verdier for y''(x)

g = 9.81
c = 2/5
M = 0.031


def v_av_x(x):
    return np.sqrt(((2*g)*(cs(0)-cs(x)))/(1+c))


def k_av_x(x):
    return cs(x, 2)/((1+(cs(x, 1)**2))**(3/2))


def n_av_x(x):
    return M*(g*np.cos(beta_av_x(x))+a_sentripetal_av_x(x))


def a_sentripetal_av_x(x):
    return (v_av_x(x)**2)*k_av_x(x)


def beta_av_x(x):
    return np.arctan(cs(x, 1))


def a_av_x(x):
    return -((5*g*np.sin(beta_av_x(x)))/7)


def f_av_x(x):
    return (2*M*g*np.sin(beta_av_x(x)))/7


def forward_euler_v(t):
    if t == 0:
        return 0
    return forward_euler_v


def forward_euler_x(t):
    if t == 0:
        return 0


# Tidsutvikling
v_x = list(map(v_av_x, x))
t_n = [0]
t = 0
for n in range(1, len(x)):
    dx = x[n] - x[n-1]
    v_x1_x2 = v_x[n-1] + v_x[n]
    dt = 2*dx / v_x1_x2
    t_n.append((dt + t))
    t += dt


cs_t = CubicSpline(t_n, x, bc_type='natural')


def x_av_t(t):
    return cs_t(t)

# # Plotteeksempel: Banen y(x)
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, y, xfast, yfast, '*')
# plt.title('Banen')
# plt.xlabel('x (m)', fontsize=20)
# plt.ylabel('y(x) (m)', fontsize=20)
# plt.grid()
# plt.show()

# # Plotteeksempel: Kulens hastighet v(x)
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, v_av_x(x))
# plt.title('Kulens hastighet')
# plt.xlabel('$x$ (m)', fontsize=20)
# plt.ylabel('$v$ (m/s)', fontsize=20)
# plt.grid()
# plt.show()


# # Plotteeksempel: Kulens sentripetalakselerasjonen a(x)
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, a_sentripetal_av_x(x))
# plt.title('Kulens sentripetalakselerasjonen')
# plt.xlabel('$x$ (m)', fontsize=20)
# plt.ylabel('$a(x)$ (m/s^2)', fontsize=20)
# plt.grid()
# plt.show()


# # Plot Beta
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, beta_av_x(x))
# plt.title('Beta')
# plt.xlabel('x(m)', fontsize=20)
# plt.ylabel('β (radianer)', fontsize=20)
# plt.grid()
# plt.show()

# # Plot krumning
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, k_av_x(x))
# plt.title('Banens krumning')
# plt.xlabel('x (m)', fontsize=20)
# plt.ylabel('k(x) (1/m)', fontsize=20)
# plt.grid()
# plt.show()


# # Plot normalkraft
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, n_av_x(x))
# plt.title('Normalkraft')
# plt.xlabel('x (m)', fontsize=20)
# plt.ylabel('N/Mg', fontsize=20)
# plt.grid()
# plt.show()

# # Plot Forholdet mellom f og N
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, abs(f_av_x(x)/n_av_x(x)))
# plt.title('Forholdet mellom f og n')
# plt.xlabel('x (m)', fontsize=20)
# plt.ylabel('|f/N|', fontsize=20)
# plt.grid()
# plt.show()


# # Plot a
# baneform = plt.figure('y(x)', figsize=(12, 6))
# plt.plot(x, a_av_x(x))
# plt.title('Kulens akselerasjon')
# plt.xlabel('x (m)', fontsize=20)
# plt.ylabel('$a(x)$ (m/s^2)', fontsize=20)
# plt.grid()
# plt.show()
