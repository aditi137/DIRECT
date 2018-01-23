#include "hilbert.h"
#include <iostream>
using namespace std;


static PyObject* LinetoAxes(PyObject* self, PyObject* args) {
	coord_t *Axes = new coord_t, *Line = new coord_t;	// multidimensional geometrical axes, linear serial # 
	int b, n;	// # bits, dimension
	PyObject *pX;

	// assigning python input to c++ types
	if (!PyArg_ParseTuple(args, "iiO!", &b, &n, &PyList_Type, &pX)) {
		PyErr_SetString(PyExc_TypeError, "Parameters must be an int, an int, and a list.");
		return NULL;
	}

	// copying PyList to Line array
	for (int i = 0; i < PyList_Size(pX); i++) {
		PyObject* pItem = PyList_GetItem(pX, i);
		if (!PyLong_Check(pItem)) {
			PyErr_SetString(PyExc_TypeError, "List items must be integers.");
		}
		Line[i] = PyLong_AsUnsignedLong(pItem);
	}

	// LinetoAxes ...
	if (n <= 1)		// trivial case
		*Axes = *Line;
	else {
		LinetoTranspose(Axes, Line, b, n);
		TransposetoAxes(Axes, b, n);
	}

	// return Axes as python list
	for (int i = 0; i < PyList_Size(pX); i++) {
		PyList_SetItem(pX, i, PyLong_FromLong(Axes[i]));
	}
	delete Axes;
	delete Line;
	return pX;
}


static PyObject* AxestoLine(PyObject* self, PyObject* args) {
	coord_t *Line = new coord_t, *Axes = new coord_t;	// linear serial #, multidimensional geometrical axes,
	coord_t store[1024];	// avoid overwriting Axes
	int b, n, i;			// # bits, dimension, counter
	PyObject *pX;

	// assigning python input to c++ types
	if (!PyArg_ParseTuple(args, "iiO!", &b, &n, &PyList_Type, &pX)) {
		PyErr_SetString(PyExc_TypeError, "Parameters must be an int, an int, and a list.");
		return NULL;
	}

	// copying PyList to Axes array
	for (int i = 0; i < PyList_Size(pX); i++) {
		PyObject* pItem = PyList_GetItem(pX, i);
		if (!PyLong_Check(pItem)) {
			PyErr_SetString(PyExc_TypeError, "List items must be integers.");
		}
		Axes[i] = PyLong_AsUnsignedLong(pItem);
	}

	// AxestoLine ...
	if (n <= 1)	// trivial case
		*Line = *Axes;
	else if (n <= 1024) {	// surely the usual case
		for (i = 0; i < n; i++)
			store[i] = Axes[i];
		AxestoTranspose(store, b, n);
		TransposetoLine(Line, store, b, n);
	}
	else {		// must do in place at greater cost
		AxestoTranspose(Axes, b, n);
		TransposetoLine(Line, Axes, b, n);
		TransposetoAxes(Axes, b, n);
	}

	// return Line as python list
	for (int i = 0; i < PyList_Size(pX); i++) {
		Py_INCREF(pX);
		PyList_SetItem(pX, i, PyLong_FromLong(Line[i]));
	}
	delete Line;
	delete Axes;
	return pX;
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


static void TransposetoAxes(
	coord_t* X,	// position: the Hilbert index stored in transposed form -> coordinate vector
	int b,		// # bits per coordinate
	int n		// dimension
)
/* Convert the Hilbert index into an N-dimensional point expressed as a vector of u-ints.*/
{
	coord_t M = 2 << (b - 1), p, q, t;
	int i;
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
}


static void AxestoTranspose(
	coord_t* X,	// position: point in n-space -> the Hilbert distance (or index) as a transposed Hilbert index
	int b,		// # bits: depth of the Hilbert curve
	int n		// dimension, 
)
/*	Given the axes (coordinates) of a point in n-dimensional space, find the distance to that point along the Hilbert curve.
	That distance will be transposed; broken into pieces and distributed into an array.
	The number of dimensions is the length of the Axes array.
*/
{
	coord_t M = 1 << (b - 1), p, q, t;
	int i;
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
}