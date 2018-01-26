def _binary_repr(num, width):
    """Return a binary string representation of `num` with zero padded to `width` bits."""
    return format(num, 'b').zfill(width)

def _hilbert_integer_to_transpose(l, bits, ndim):
    """Store a Hilbert integer (`l`) as its transpose (`x`).
    :param l: integer distance along Hilbert curve
    :type l: ``int``
    :param bits: number of iterations in Hilbert curve
    :type bits: ``int``
    :param ndim: number of dimensions
    :type ndim: ``int``
    """
    l_bit_str = _binary_repr(l, bits * ndim)
    x = [int(l_bit_str[i::ndim], 2) for i in range(ndim)]
    return x

def _transpose_to_hilbert_integer(x, bits, ndim):
    """Restore a Hilbert integer (`l`) from its transpose (`x`).
    :param x: the transpose of a Hilbert integer (ndim components of length bits)
    :type x: ``list`` of ``int``
    :param bits: number of iterations in Hilbert curve
    :type bits: ``int``
    :param ndim: number of dimensions
    :type ndim: ``int``
    """
    x_bit_str = [_binary_repr(x[i], bits) for i in range(ndim)]
    l = int(''.join([y[i] for i in range(bits) for y in x_bit_str]), 2)
    return l

def coordinates_from_distance(l, bits, ndim):
    """Return the coordinates for a given Hilbert distance.
    :param l: integer distance along the curve
    :type l: ``int``
    :param bits: side length of hyper-cube is 2^bits
    :type bits: ``int``
    :param ndim: number of dimensions
    :type ndim: ``int``
    """
    x = _hilbert_integer_to_transpose(l, bits, ndim)
    Z = 2 << (bits-1)

    # Gray decode by H ^ (H/2)
    t = x[ndim-1] >> 1
    for i in range(ndim-1, 0, -1):
        x[i] ^= x[i-1]
    x[0] ^= t

    # Undo excess work
    Q = 2
    while Q != Z:
        P = Q - 1
        for i in range(ndim-1, -1, -1):
            if x[i] & Q:
                # invert
                x[0] ^= P
            else:
                # exchange
                t = (x[0] ^ x[i]) & P
                x[0] ^= t
                x[i] ^= t
        Q <<= 1

    # done
    return x

def distance_from_coordinates(x, bits, ndim):
    """Return the Hilbert distance for a given set of coordinates.
    :param x: coordinates len(x) = ndim
    :type x: ``list`` of ``int``
    :param bits: side length of hyper-cube is 2^bits
    :type bits: ``int``
    :param ndim: number of dimensions
    :type ndim: ``int``
    """
    M = 1 << (bits - 1)

    # Inverse undo excess work
    Q = M
    while Q > 1:
        P = Q - 1
        for i in range(ndim):
            if x[i] & Q:
                x[0] ^= P
            else:
                t = (x[0] ^ x[i]) & P
                x[0] ^= t
                x[i] ^= t
        Q >>= 1

    # Gray encode
    for i in range(1, ndim):
        x[i] ^= x[i-1]
    t = 0
    Q = M
    while Q > 1:
        if x[ndim-1] & Q:
            t ^= Q - 1
        Q >>= 1
    for i in range(ndim):
        x[i] ^= t

    l = _transpose_to_hilbert_integer(x, bits, ndim)
    return l