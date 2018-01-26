import src._hilbert as hilbert

bits = 5
ndim = 3

def test_integer_to_transpose():
    """Assert that a 15 bit Hilbert integer is correctly transposed into a 3-d vector
                  ABCDEFGHIJKLMNO
         10590 (0b010100101011110)
                  ADGJM
         X[0] = 0b01101 = 13
                  BEHKN
         X[1] = 0b10011 = 19
                  CFILO
         X[2] = 0b00110 = 6
    """
    assert hilbert._hilbert_integer_to_transpose(10590, bits, ndim) == [13, 19, 6]

def test_transpose_to_integer():
    """Assert that a 15 bit Hilbert integer is correctly recovered from its transposed 3-d vector
                  ABCDEFGHIJKLMNO
         10590 (0b010100101011110)
                  ADGJM
         X[0] = 0b01101 = 13
                  BEHKN
         X[1] = 0b10011 = 19
                  CFILO
         X[2] = 0b00110 = 6
    """ 
    assert hilbert._transpose_to_hilbert_integer([13, 19, 6], bits, ndim) == 10590

def test_reversibility():
    """Assert coordinates_from_distance and distance_from_coordinates are inverse operations."""
    n_h = 2**(ndim * bits)
    for h in range(n_h):
        x = hilbert.coordinates_from_distance(h, bits, ndim)
        h_test = hilbert.distance_from_coordinates(x, bits, ndim)
        assert h == h_test