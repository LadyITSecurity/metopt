import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def func(t: float):
    return np.log(t ** 2 - 4 * t + 5)


def method_dichotomy(left, right, epsilon, delta):
    count_iter = 0
    count_calc = 0
    y1 = 0
    y2 = 0
    while abs(right - left) > 2*delta:
        x1 = (right + left - epsilon)/2
        x2 = (right + left + epsilon)/2
        y1 = func(x1)
        y2 = func(x2)
        count_calc += 2
        if y1 > y2:
            left = x1
        else:
            right = x2
        count_iter += 1
    return count_calc, count_iter, y1, x1


def method_golden_ratio(left, right, delta):
    last_side, count_calc, count_iter = 0, 0, 0
    const = 0.618
    y0 = 0
    while abs(right - left) > 2 * delta:
        if count_iter == 0:
            x2 = (left + const * (right - left))
            x1 = (left + (1 - const) * (right - left))
            y1, y2 = func(x1), func(x2)
            count_calc += 1
        elif last_side == 0:
            x1, y1 = x2, y2
            x2 = (left + const * (right - left))
            y2 = func(x2)
        else:
            x2, y2 = x1, y1
            x1 = left + (1 - const) * (right - left)
            y1 = func(x1)
        count_calc += 1
        if y1 > y2:
            left = x1
            y0 = y2
            last_side = 0
        else:
            right = x2
            y0 = y1
            last_side = 1
        count_iter += 1
        #print(left, right)
    result = (left + right)/2
    return count_calc, count_iter, y0, result


def method_fibonacci(left, right, count_iter):
    fib = []
    p, c = 0, 1
    for i in range(count_iter + 1):
        fib += [c]
        c, p = c+p, c
    #print(fib)
    count_calc = 0
    l = (right - left) / fib[count_iter]

    x1 = left + l * fib[count_iter - 2]
    x2 = right - l * fib[count_iter - 2]
    y1 = func(x1)
    y2 = func(x2)
    count_calc = 2
    k = 1
    while k < count_iter - 1:
        if y1 < y2:
            right = x2
            k += 1
            x2, y2 = x1, y1
            x1 = left + l * fib[count_iter - 1 - k]
            y1 = func(x1)
        else:
            left = x1
            k += 1
            x1, y1 = x2, y2
            x2 = right - l * fib[count_iter - 1 - k]
            y2 = func(x2)
        count_calc += 1
        #print(right - left)
        if y1 > y2:
            left = x1
        else:
            right = x2
        x0 = (left + right) / 2
        y0 = y1
        #print(right - left)
    result = (right - left) / 2
    return count_calc, count_iter, y0, x0, result


from sympy import symbols
x = symbols('x')
# f = sp.sympify(log(x ** 2 - 4 * x + 5))
a = lambda x: sp.ln(x ** 2 - 4 * x + 5)
# f = sp.ln(x ** 2 - 4 * x + 5)

count_iter, count_calc = 0, 0
x0, y0 = 0, 0
left = 0.99
right = 4.0
delta = 0.001
epsilon = 0.0001
t = np.linspace(1, 4, 100)
f = np.log(t ** 2 - 4 * t + 5)



'''
data = []
col_names = ['Название метода', ' итераций ', 'вычислений', '    x0    ', 'y0']
count_calc, count_iter, y0, x0 = method_dichotomy(left, right, epsilon, delta)
data.append(['дихотомия', count_iter, count_calc, x0, y0])
x0, y0 = 0.0, 0.0
count_calc, count_iter, y0, x0 = method_golden_ratio(left, right, delta)
data.append(['золотое сечение', str(count_iter), str(count_calc), str(x0), str(y0)])
x0, y0 = 0.0, 0.0
count_iter = 16
count_calc, count_iter, y0, x0, result = method_fibonacci(left, right, count_iter)
data.append(['Фибоначчи', str(count_iter), str(count_calc), str(x0), str(y0)])

from tabulate import tabulate
print(tabulate(data, headers=col_names, tablefmt="fancy_grid", numalign='center'))


print(result)
plt.plot(t, f)
plt.show()

'''
'''
count_calc, count_iter, y0, x0 = method_dichotomy(left, right, epsilon, delta)
print('result of method of dichotomy')
print('count of calculation ', count_calc)
print('count of iteration ', count_iter)
print('x0 =  ', x0)
print('y0 =  ', y0)
print()
x0, y0 = 0, 0

count_calc, count_iter, y0, x0 = method_golden_ratio(left, right, delta)
print('result of method of golden ratio')
print('count of calculation ', count_calc)
print('count of iteration ', count_iter)
print('x0 =  ', x0)
print('y0 =  ', y0)
print()
x0, y0 = 0, 0

count_iter = 50
count_calc, count_iter, y0, x0 = method_fibonacci(left, right, count_iter)
print('result of method of Fibonacci')
print('count of calculation ', count_calc)
print('count of iteration ', count_iter)
print('x0 =  ', x0)
print('y0 =  ', y0)
print()
x0, y0 = 0, 0
'''