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
from scipy.optimize import curve_fit


# Horisontal avstand mellom festepunktene er 0.200 m
# h = 0.200
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

for n in range(0, len(xfast)):
    print(str(np.round((xfast[n]), 3)) + " & " + str(np.round((yfast[n]), 3)))


# def realverdier():
#     t_vals = []
#     x_vals = [0 for x in range(0, 37)]
#     for n in range(1, 11):
#         path = "./tracked_data/txy_" + str(n) + ".txt"
#         f = open(path)
#         i = 0
#         for ln in f.readlines()[2:]:
#             if n < 2:
#                 t_vals.append(np.round(float(ln.split('\t')[0]), 3))
#             print(ln.split('\t'))
#             x_vals[i] += np.round(float(ln.split('\t')[1]), 3)
#             i += 1
#         f.close()
#     x_vals_mean = list(map(lambda x: x/10, x_vals))
#     x_vals_mean_adjusted = list(
#         map(lambda x: np.round(x - x_vals_mean[0], 3), x_vals_mean))
#     return x_vals_mean_adjusted, t_vals


# x_real, t_real = realverdier()

# print(x_real)
# print(t_real)


# Skriv inn y-verdiene til banens 8 festepunkter i tabellen yfast.
# Et vilkaarlig eksempel:
# yfast = np.asarray([0.391, 0.295, 0.272, 0.233, 0.130, 0.103, 0.182, 0.204])
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


def v_x_av_x(x):
    return v_av_x(x)*np.cos(np.cos(cs()))


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
    v_x1_x2 = v_x[n-1]*np.cos(beta_av_x(x[n-1])) + \
        v_x[n]*np.cos(beta_av_x(x[n]))
    dt = 2*dx / v_x1_x2
    t_n.append((dt + t))
    t += dt


cs_t = CubicSpline(t_n, x, bc_type='natural')


def x_av_t(t):
    return cs_t(t)


print(t_n[-1])

sluttfarter = []

t_real = [0 for x in range(0, 36)]
x_real = [0 for x in range(0, 36)]
y_real = [0 for x in range(0, 36)]
for n in range(1, 11):
    path = "./tracked_data/txy_" + str(n) + ".txt"
    f = open(path)
    lines = f.readlines()
    t0, x0, y0 = [float(p) for p in lines[2].split('\t')]
    t_real_l = []
    x_real_l = []
    y_real_l = []
    for i in range(2, len(lines)):
        t_real[i-2] += float(lines[i].split('\t')[0]) - t0
        x_real[i-2] += float(lines[i].split('\t')[1]) - x0
        y_real[i-2] += float(lines[i].split('\t')[2]) - y0
        t_real_l.append(float(lines[i].split('\t')[0]) - t0)
        x_real_l.append(float(lines[i].split('\t')[1]) - x0)
        y_real_l.append(float(lines[i].split('\t')[2]) - y0)
    print(t_real_l)
    print(t_real_l)
    csr_l_x = CubicSpline(t_real_l, x_real_l)
    csr_l_y = CubicSpline(t_real_l, y_real_l)
    der_x_l = csr_l_x.derivative()
    der_y_l = csr_l_y.derivative()
    sluttfarter.append(np.sqrt(csr_l_x(t_n[-1])**2+csr_l_y(t_n[-1])**2))
    f.close()

t_real = [t/10 for t in t_real]
x_real = [x/10 for x in x_real]
y_real = [y/10 for y in y_real]

print('')
print('')
print(sluttfarter)
print('')
print('')
csrx = CubicSpline(t_real, x_real)
csry = CubicSpline(t_real, y_real)


def v_real(t):
    der_x = csrx.derivative()
    der_y = csry.derivative()
    return np.sqrt(der_x(t)**2+der_y(t)**2)


def a_real(t):
    der_x = csrx.derivative(2)  # csrx.derivative()
    der_y = csry.derivative(2)
    return np.sqrt(der_x(t)**2+der_y(t)**2)

    # # Plotteeksempel: Banen y(x)
    # baneform = plt.figure('y(x)', figsize=(12, 6))
    # plt.plot(x, y, xfast, yfast, '*')
    # plt.title('Banen')
    # plt.xlabel('x (m)', fontsize=20)
    # plt.ylabel('y(x) (m)', fontsize=20)
    # plt.grid()
    # plt.show()
    # Plotteeksempel: Kulens hastighet v(x)
baneform = plt.figure('y(x)', figsize=(12, 6))
plt.plot(x_av_t(t_n), v_av_x(x_av_t(t_n)))
plt.plot(csrx(t_n), v_real(t_n))
plt.title('Kulens hastighet')
plt.xlabel('$x$ (m)', fontsize=20)
plt.ylabel('$v$ (m/s)', fontsize=20)
plt.grid()
plt.show()


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

# Plot x av t
baneform = plt.figure('y(x)', figsize=(12, 6))
plt.plot(t_n, x_av_t(t_n))
plt.plot(t_n, csrx(t_n))
plt.title('Kulens x-posisjon')
plt.xlabel('t (s)', fontsize=20)
plt.ylabel('$x$ (m)', fontsize=20)
plt.grid()
plt.show()

# print("Beregnet: ")
# print(v_av_x(x[-1]))
# print("Målt: ")
# print(np.mean(sluttfarter))
# print("Diff: ")
# print(np.mean(sluttfarter)/v_av_x(x[-1]))

sluttfarter_snitt = np.mean(sluttfarter)

sluttfarter_sd = np.std(sluttfarter)
sluttfarter_sd = np.std(sluttfarter)

sluttfarter_varians = np.var(sluttfarter)

print("\nForsøk \t\t Sluttfart")
print("-"*27)
for n in range(1, 11):
    print(str(n) + " \t\t " + str(np.round(sluttfarter[n-1], 3)))


print("")
print("Snittfart:" + "\t" + str(sluttfarter_snitt))
print("Standardavvik:" + "\t" + str(sluttfarter_sd))
print("Varians:" + "\t" + str(sluttfarter_varians))
print("Standardfeil:" + "\t" + str(np.round(np.sqrt(sluttfarter_varians /
      len(sluttfarter))*100 / sluttfarter_snitt, 3)) + "%")
print("Beregnet fart:" + "\t" + str(np.round(v_av_x(x[-1]), 3)))
print("Avvik:" + "\t\t" +
      str(np.round((np.mean(sluttfarter)/v_av_x(x[-1])-1)*100, 3)) + "%")
