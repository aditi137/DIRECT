//TODO: Unhandled exception, most probably a memory leak somewhere
#include "hilbert.h"

/*
void LinetoAxes (
coord_t* Axes,	// multidimensional geometrical axes
coord_t* Line,	// linear serial #
int b,			// # bits
int n)			// dimension
{
if (n <= 1)		// trivial case
*Axes = *Line;
else {
LinetoTranspose(Axes, Line, b, n);
TransposetoAxes(Axes, b, n);
}
}


void AxestoLine (
coord_t* Line,	// linear serial #
coord_t* Axes,	// multidimensional geometrical axes
int b,			// # bits
int n)			// dimension
{
coord_t store[1024];	// avoid overwriting Axes
int i;

// trivial case
if (n <= 1)
*Line = *Axes;
// surely the usual case
else if (n <= 1024) {
for (i = 0; i < n; i++)
store[i] = Axes[i];
AxestoTranspose(store, b, n);
TransposetoLine(Line, store, b, n);
}
// must do in place at greater cost
else {
AxestoTranspose(Axes, b, n);
TransposetoLine(Line, Axes, b, n);
TransposetoAxes(Axes, b, n);
}
}


static void LinetoTranspose (
coord_t* X,		// transpose
coord_t* Line,	// Hilbert integer
int b,			// # bits
int n)			// dimension
{
coord_t M = 1 << (b - 1), p, q;
int i, j = 0;
p = M;
for (i = 0; i < n; i++)
X[i] = 0;
for (i = 0; i < n; i++)
for (q = M; q > 0 ; q >>= 1) {
if (Line[i] & q)
X[j] |= p;
if (++j == n) {
j = 0;
p >>= 1;
}
}
}


static void TransposetoLine(
coord_t* Line,	// Hilbert integer
coord_t* X,		// transpose
int b,			// # bits
int n)			// dimension
{
coord_t M = 1 << (b - 1), p, q;
int i, j = 0;
p = M;

for (i = 0; i < n; i++) {
Line[i] = 0;
for (q = M; q > 0; q >>= 1) {
if (X[j] & p)
Line[i] |= q;
if (++j == n) {
j = 0;
p >>= 1;
}
}
}
}
*/

static PyObject* TransposetoAxes(PyObject* self, PyObject* args)
{
	coord_t* X = new coord_t;	// position
	int b, n, i;				// # bits, dimension, counter
	PyObject *pX, *pItem;

	// assigning python input to c++ types
	if (!PyArg_ParseTuple(args, "iiO!", &b, &n, &PyList_Type, &pX)) {
		PyErr_SetString(PyExc_TypeError, "Parameters must be an int, an int, and a list.");
		return NULL;
	}

	for (i = 0; i < PyList_Size(pX); i++) {
		pItem = PyList_GetItem(pX, i);
		if (!PyLong_Check(pItem)) {
			PyErr_SetString(PyExc_TypeError, "List items must be integers.");
			return NULL;
		}
		X[i] = PyLong_AsUnsignedLong(pItem);
	}

	// TransposetoAxes ...
	coord_t M = 2 << (b - 1), p, q, t;
	// Gray decode by H ^ (H/2)
	t = X[n - 1] >> 1;
	for (i = n - 1; i > 0; i--)
		X[i] ^= X[i - 1];
	X[0] ^= t;

	// Undo excess work
	for (q = 2; q != M; q <<= 1) {
		p = q - 1;
		for (i = n - 1; i >= 0; i--) {
			// invert
			if (X[i] & q)
				X[0] ^= p;
			// exchange
			else {
				t = (X[0] ^ X[i]) & p;
				X[0] ^= t;
				X[i] ^= t;
			}
		}
	}

	// return X as python list
	for (i = 0; i < PyList_Size(pX); i++) {
		PyList_SetItem(pX, i, PyLong_FromLong(X[i]));
	}
	delete X;
	return pX;
}


static PyObject* AxestoTranspose(PyObject* self, PyObject* args)
{
	coord_t* X = new coord_t;	// position
	int b, n, i;				// # bits, dimension, counter
	PyObject *pX, *pItem;

	// assigning python input to c++ types
	if (!PyArg_ParseTuple(args, "iiO!", &b, &n, &PyList_Type, &pX)) {
		PyErr_SetString(PyExc_TypeError, "Parameters must be an int, an int, and a list.");
		return NULL;
	}

	for (i = 0; i < PyList_Size(pX); i++) {
		pItem = PyList_GetItem(pX, i);
		if (!PyLong_Check(pItem)) {
			PyErr_SetString(PyExc_TypeError, "List items must be integers.");
			return NULL;
		}
		X[i] = PyLong_AsUnsignedLong(pItem);
	}

	// AxestoTranspose ...
	coord_t M = 1 << (b - 1), p, q, t;
	// Inverse undo
	for (q = M; q > 1; q >>= 1) {
		p = q - 1;
		for (i = 0; i < n; i++) {
			// invert
			if (X[i] & q)
				X[0] ^= p;
			// exchange
			else {
				t = (X[0] ^ X[i]) & p;
				X[0] ^= t;
				X[i] ^= t;
			}
		}
	}
	// Gray encode (inverse of decode)
	for (i = 1; i < n; i++)
		X[i] ^= X[i - 1];
	t = 0;
	for (q = M; q > 1; q >>= 1)
	if (X[n - 1] & q)
		t ^= q - 1;
	for (i = 0; i < n; i++)
		X[i] ^= t;

	// return X as python list
	for (i = 0; i < PyList_Size(pX); i++) {
		PyList_SetItem(pX, i, PyLong_FromLong(X[i]));
	}
	delete X;
	return pX;
}