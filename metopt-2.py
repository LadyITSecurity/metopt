import numpy as np
from sympy import symbols


def func(x, y, z, alpha, x0, z0):
    fx = (x - x0) * np.cos(alpha) + (z - z0) * np.sin(alpha)  # x`
    fz = (z - z0) * np.cos(alpha) - (x - x0) * np.sin(alpha)  # z`
    return fx ** 2 + np.cosh(y) + fz ** 2  # f(...)


def dfdx(x, y, z, alpha, x0, z0):
    fx = (x - x0) * np.cos(alpha) + (z - z0) * np.sin(alpha)
    fz = (z - z0) * np.cos(alpha) - (x - x0) * np.sin(alpha)
    return 2 * fx * np.cos(alpha) - 2 * fz * np.sin(alpha)


def dfdy(x, y, z, alpha, x0, z0):
    return np.sinh(y)


def dfdz(x, y, z, alpha, x0, z0):
    fx = (x - x0) * np.cos(alpha) + (z - z0) * np.sin(alpha)
    fz = (z - z0) * np.cos(alpha) - (x - x0) * np.sin(alpha)
    return 2 * fx * np.sin(alpha) + 2 * fz * np.cos(alpha)


def dfdydy(yc):
    return np.cosh(yc)


def get_min_coordinate(x0, z0, xc, yc, zc, alpha, epsilon, delta, count):
    print('Координатный метод: \n')
    this, prev = 0, 0
    orts = np.array([xc, yc, zc])
    step = np.array([[epsilon, 0, 0], [0, epsilon, 0], [0, 0, epsilon]])
    distance = 2 * delta

    while distance ** 0.5 > delta:
        first = np.array([orts[0], orts[1], orts[2]])
        for i in step:
            this = func(orts[0], orts[1], orts[2], alpha, x0, z0)
            orts += i
            next = func(orts[0], orts[1], orts[2], alpha, x0, z0)
            if next > this:
                i *= -1
                orts += i
            prev, this = this, func(orts[0], orts[1], orts[2], alpha, x0, z0)
            while this <= prev:
                orts += i
                prev, this = this, func(orts[0], orts[1], orts[2], alpha, x0, z0)
                print("\t", orts[0].round(3), orts[1].round(3), orts[2].round(3), "\t", prev)
            orts -= i
            print(orts[0].round(3), orts[1].round(3), orts[2].round(3), prev)
        second = np.array([orts[0], orts[1], orts[2]])
        distance = 0
        for i in range(count):
            distance += (second[i] - first[i]) ** 2
        print()
    return orts, prev


def get_min_gradient(x0, z0, xc, yc, zc, alpha, epsilon, delta, count):
    print('Метод градиента: \n')
    this, prev = 0, 0
    orts = np.array([xc, yc, zc])
    distance = 2 * delta
    iter = 1
    while distance ** 0.5 > delta:
        first = np.array([orts[0], orts[1], orts[2]])
        step = np.array([epsilon * dfdx(orts[0], orts[1], orts[2], alpha, x0, z0),
                         epsilon * dfdy(orts[0], orts[1], orts[2], alpha, x0, z0),
                         epsilon * dfdz(orts[0], orts[1], orts[2], alpha, x0, z0)])
        orts -= step
        prev, this = this, func(orts[0], orts[1], orts[2], alpha, x0, z0)
        second = np.array([orts[0], orts[1], orts[2]])
        distance = 0
        for i in range(count):
            distance += (second[i] - first[i]) ** 2
        print(iter, '-ая итерация: \t', end='')
        iter += 1
        print(orts[0].round(3), orts[1].round(3), orts[2].round(3), '\tf(...) = ', prev)
    return orts, prev


def get_grad(x0, z0, xc, yc, zc, alpha):
    result = [[dfdx(xc, yc, zc, alpha, x0, z0)],
              [dfdy(xc, yc, zc, alpha, x0, z0)],
              [dfdz(xc, yc, zc, alpha, x0, z0)]]
    return result


def method_Newton(x0, z0, xc, yc, zc, alpha, delta, count):
    print('Метод Ньютона: \n')
    x_next = [[xc],
              [yc],
              [zc]]
    distance = 2 * delta
    print('начальные координаты')
    print(x_next, '\n')
    iter = 1
    while distance ** 0.5 > delta:
        xk = x_next
        H = [[2.0, 0.0, 0.0],
             [0.0, dfdydy(xk[1][0]), 0.0],
             [0.0, 0.0, 2.0]]
        grad = get_grad(x0, z0, xk[0][0], xk[1][0], xk[2][0], alpha)
        inv = np.linalg.inv(H)  # обратная матрица
        x_next = xk - np.dot(inv, grad)
        print('---------', iter, '-ая итерация -----------')
        iter += 1
        print('матрица Н')
        print(H[0])
        print(H[1])
        print(H[2], '\n')
        print('матрица градиента')
        print(grad, '\n')
        print('результат')
        print(x_next, '\n')
        distance = 0
        for i in range(count):
            distance += (x_next[i][0] - xk[i][0]) ** 2
    return xk


def main():
    x0, z0 = 1, 3
    xc, yc, zc = 3.0, 3.0, 3.0
    alpha = np.pi/6
    epsilon = 0.5
    delta = 0.01

    orts, f = get_min_gradient(x0, z0, xc, yc, zc, alpha, epsilon, delta, 3)    
    print("Минимум функции: ")
    print(orts[0].round(3), orts[1].round(3), orts[2].round(3))
    print(f, '\n')

    result = method_Newton(x0, z0, xc, yc, zc, alpha, delta, 3)
    print('Минимум функции (альфа = ', alpha, ' ):')
    print(result)
    print(func(result[0][0], result[1][0], result[2][0], alpha, x0, z0))


if __name__ == "__main__":
    main()
