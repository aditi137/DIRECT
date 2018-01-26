import Hilbert

def test_Hilbert_import():
    assert dir(Hilbert) == ['__doc__', '__file__', '__loader__', '__name__', '__package__', '__spec__', 'h2r', 'r2h']

def test_Hilbert_curve():
    X = Hilbert.r2h(5, 3, [5, 10, 20])
    assert X == [7, 21, 25]
    assert Hilbert.h2r(5, 3, X) == [5, 10, 20]