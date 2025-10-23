# daxpy.py
import numpy as np
import pytest

def test_daxpy():
    # case 1
    a = 2.0
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    expected = np.array([6.0, 9.0, 12.0])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()

    # case 2
    a = 3.0
    x = np.array([0.0, -1.0, 1.0])
    y = np.array([1.0, 1.0, 1.0])
    expected = np.array([1.0, -2.0, 4.0])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()

    # case 3
    a = 0.0
    x = np.array([1.5, 2.5])
    y = np.array([3.5, 4.5])
    expected = np.array([3.5, 4.5])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()

    # case 4
    a = 5.0
    x = np.array([])
    y = np.array([])
    expected = np.array([])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()
    print("All test cases passed!")

    # case 5 wrong onput
    a = -1.0
    x = np.array([1.0])
    y = np.array([2.0])
    expected = np.array([2.0])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()


def daxpy(a, x, y):
    if x.shape != y.shape:
        raise ValueError("Input vectors must have the same shape.")
    return a * x + y

if __name__ == "__main__":
    N = 1000
    # seed generator from device
    rng = np.random.default_rng()
    # generate random numbers from standard gaussian distribution
    x = rng.standard_normal((0, 1, N))
    y = rng.standard_normal((0, 1, N))
    a = rng.uniform(-10, 10)
    result = daxpy(a, x, y)
    print("DAXPY computation completed.")
    # Run tests
    test_daxpy()