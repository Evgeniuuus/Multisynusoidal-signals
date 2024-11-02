import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from math import sin, sqrt, atan, cos
from sympy.functions.elementary.complexes import arg

mu = 1  # Задаем все параметры
w = 1
phi = 0
delta_t = 0.25
n = 100
x = [0]

print("----------------------------Промежуточные действия----------------------------")
func = sp.sympify("1+sin(x)+sin(2*x)+sin(3*x)")  # Задаем и проверяем функции с одной и тремя гармониками
# 1+sin(1*x+1)      и     1+sin(x)+sin(2*x)+sin(3*x)

for i in range(n):
    x.append(x[i] + delta_t)

y = [func.subs(sp.Symbol('x'), x[i]) for i in range(len(x))]

plt.plot(x, y, label="Функция")
plt.legend()
plt.grid()
plt.show()

matrixa = np.zeros((n + 1, n + 1))  # Создаем нулевую матрицу

matrixa[0][0] = 1  # Заполняем по правилу
for i in range(n):
    matrixa[1][i] = y[i]
eps = 10 ** -3
m = 0
for j in range(2, n):
    for k in range(0, n + 1 - j):
        matrixa[j][k] = (matrixa[j - 2][k + 1] / matrixa[j - 2][0] - matrixa[j - 1][k + 1] / matrixa[j - 1][0])
    print(sum(abs(matrixa[j])))
    if sum(abs(matrixa[j])) < eps:
        print("l=", j)
        m = j + 1
        break

if m == 0:
    print("Гармоник нет.")
    exit()

print("----------------------------Матрица----------------------------")
print("Номер нулевой строки: ", m)
n = int((m - 4) / 4)
print("Количество гармоник: ", n)

fraction = str(y[0]) + "/"  # Делаем дробь
for i in range(2, m - 1):
    fraction += "(1+("
    fraction += str(matrixa[i][0]) + "*z**-1)"
    fraction += "/"
fraction = fraction[:-1]

for i in range(2, m - 1):
    fraction += ")"
c = sp.sympify(fraction)
a = sp.simplify(fraction)  # Упрощаем дробь

b = str(a)  # Вытаскиваем знаменатель этой дроби
index = b.find('/') + 1
denominator = b[index:]
print("----------------------------Работа с корнями----------------------------")
print("Полученный знаменатель: ", denominator)

roots = sp.solve(denominator, sp.Symbol('z'))  # Находим корни нашего знаменателя
print("\nВсе корни: ", roots)

complex_nuber_type = "5 + 2*I"
complex_nuber_type = sp.sympify(complex_nuber_type)
args = []
for i in range(len(roots)):
    if type(roots[i]) is type(complex_nuber_type):  # Если корень комплексное число добавляем в список
        args.append(arg(roots[i]))
print("\nАргументы комплексных: ", args)

wi = []
args.sort()
args_New = []

for i in range(len(args)):  # Только положительные корни
    if args[i] > 0:
        args_New.append(args[i])

for i in range(n):
    wi.append(args_New[i] / delta_t)

print("\nПревращенные в омегу", wi, "\n")

print("По формулам: ")

element1 = "mu+"
for i in range(n):
    element1 += "b" + str(i) + "+"
element1 = element1[:-1]
element1 += "==" + str(y[0])

equations = []  # Дальше все по формулам
for i in range(2 * n):
    y_in_delta = func.subs(sp.Symbol('x'), (i + 1) * delta_t)
    our_sum = "mu+"
    for j in range(n):
        our_sum += "a" + str(j) + "*sin(" + str(i + 1) + "*" + str(wi[j]) + "*" + str(delta_t) + ")+b" + str(
            j) + "*cos(" + str(i + 1) + "*" + str(wi[j]) + "*" + str(delta_t) + ")+"
    our_sum = our_sum[:-1]
    our_sum += "==" + str(y_in_delta)
    equations.append(our_sum)

for i in range(2 * n):
    print(equations[i])

element1 = sp.N(sp.sympify(element1, evaluate=False))
equations = [sp.N(sp.sympify(equations[i], evaluate=False)) for i in range(len(equations))]
equations.append(element1)

print("\n")
for i in range(len(equations)):
    print(equations[i])
unknowns = []

for i in range(n):
    unknowns.append(sp.Symbol('a' + str(i)))
for i in range(n):
    unknowns.append(sp.Symbol('b' + str(i)))
unknowns.append(sp.Symbol('mu'))
res = sp.linsolve(equations, unknowns)
print("\nРезультат:", res)

c = [sqrt(pow(res.args[0][i], 2) + pow(res.args[0][i + n], 2)) for i in range(n)]
phi = [atan(res.args[0][i + n] / res.args[0][i]) for i in range(n)]
a = [c[i] * cos(phi[i]) for i in range(n)]
b = [c[i] * sin(phi[i]) for i in range(n)]
mu = res.args[0][2 * n]

harmonics = []
harmonics_y = []
for i in range(n):  # Для n штук
    harmonics.append(sp.sympify(
        str(mu) + "+" + str(a[i]) + "*sin(" + str(wi[i]) + "*x)+" + str(b[i]) + "*cos(" + str(wi[i]) + "*" + "x)"))
    harmonics_y.append([])
    harmonics_y[i] = [harmonics[i].subs(sp.Symbol('x'), x[j]) for j in range(len(x))]  # Подставляем
print("\nГармоники:")

for i in range(n):
    print("Гармоника " + str(i) + ":", end="")
    print(harmonics[i])

for i in range(n):
    plt.plot(x, harmonics_y[i], '--', label=str("Гармоника " + str(i + 1)))

plt.plot(x, y, label="Наша функция")
plt.legend()
plt.grid()
plt.show()

y_polyharmony = []  # Складываем гармоники
for i in range(len(x)):
    y_polyharmony.append(0)
for i in range(len(x)):
    for j in range(n):
        y_polyharmony[i] += harmonics_y[j][i]
    y_polyharmony[i] -= mu * (n - 1)  # Сделаем через костыль * (n - 1)
plt.plot(x, y_polyharmony, label='Полигармоника')
plt.plot(x, y, '--', label='Наша функция', color="green", )
plt.grid()
plt.legend()
plt.show()
