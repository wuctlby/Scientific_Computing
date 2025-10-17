import numpy as np
from scipy.integrate import quad

expected_result = (np.exp(np.pi/2) - 1) / 2

def f(x):
    return np.exp(x) * np.cos(x)

def relative_error(approx, exact):
    return abs(approx - exact) / abs(exact)

if __name__ == "__main__":
    result, error = quad(f, 0, np.pi/2)
    print(f"Integral result: {result:.15f}, Estimated error: {error}")
    print(f"Expected result: {expected_result:.15f}")
    print(f"Relative error: {relative_error(result, expected_result)}")

    # integrate from .dat file
    with open("data_points.dat", "r") as file:
        data = np.loadtxt(file, skiprows=1)
    x_data, y_data = data[:, 0], data[:, 1]

    # Perform numerical integration using the trapezoidal rule
    result = np.trapezoid(y_data, x_data)
    print(f"Integral result from data: {result:.15f}")
    print(f"Relative error from data: {relative_error(result, expected_result)}")