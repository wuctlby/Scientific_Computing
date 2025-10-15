import time
import numpy as np

# scalar
a = None

# vector
x = []
y = []

#results
result_d = []

# dimension
N = [10, 10**6, 10**8]


def compute_hello_world(a_value, x_value, y_value, dimension, equation, result_d_value=None):
    global a, x, y, result_d
    a = a_value
    x = [x_value] * dimension
    y = [y_value] * dimension

    for i in range(dimension):
        result_d.append(eval(equation, {'a': a, 'x': x[i], 'y': y[i]}))
        if result_d_value is not None and result_d[i] != result_d_value:
            raise ValueError(f"Computation error at index {i}: expected {result_d_value}, got {result_d[i]}")
    print(f"Computation completed for dimension {dimension}. All results match expected value: {result_d_value}")
    return result_d, a, x, y

if __name__ == "__main__":
    a_value = 3
    x_value = 0.1
    y_value = 7.1

    start_time = time.time()
    print("Starting computations with normal Python lists...")
    equation = 'a * x + y'
    result_d_value = eval(equation, {'a': a_value, 'x': x_value, 'y': y_value})
    with open("results.txt", "w") as f:
        f.write(f"Equation: {equation}\n")
        f.write(f"Expected result for all elements: {result_d_value}\n")
    print(f"Expected result for equation '{equation}' for all elements: {result_d_value}")

    for dimension in N:
        results, a, x, y = compute_hello_world(a_value, x_value, y_value, dimension, equation, result_d_value)
        
        middle_time = time.time()
        print(f"Time taken for dimension {dimension}: {middle_time - start_time:.2f} seconds")
        with open("results.txt", "a") as f:
            f.write(f"Dimension: {dimension}, All results match expected value: {result_d_value}\n")
            f.write(f"Time taken: {middle_time - start_time:.2f} seconds for dimension {dimension}\n")
        result_d = []

    end_time = time.time()
    print(f"All computations completed in {end_time - start_time:.2f} seconds \n")

    start_time = time.time()
    print("Starting computations with NumPy arrays...")
    a_value = 3
    x_value = 0.1
    y_value = 7.1
    equation = 'a * x + y'
    result_d_value = eval(equation, {'a': a_value, 'x': x_value, 'y': y_value})
    with open("results.txt", "a") as f:
        f.write(f"\nUsing NumPy arrays:\n")
        f.write(f"Equation: {equation}\n")
        f.write(f"Expected result for all elements: {result_d_value}\n")
    print(f"Expected result for equation '{equation}' for all elements: {result_d_value}")
    for dimension in N:
        a = a_value
        x = np.full(dimension, x_value)
        y = np.full(dimension, y_value)
        result_d = eval(equation, {'a': a, 'x': x, 'y': y})
        if not np.all(result_d == result_d_value):
            raise ValueError(f"Computation error: not all results match expected value {result_d_value}")
        print(f"Computation completed for dimension {dimension}. All results match expected value: {result_d_value}")
        middle_time = time.time()
        print(f"Time taken for dimension {dimension}: {middle_time - start_time:.2f} seconds")
        with open("results.txt", "a") as f:
            f.write(f"Dimension: {dimension}, All results match expected value: {result_d_value}\n")
            f.write(f"Time taken: {middle_time - start_time:.2f} seconds for dimension {dimension}\n")
        result_d = []
    end_time = time.time()
    print(f"All computations completed in {end_time - start_time:.2f} seconds")