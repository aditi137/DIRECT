//=============================================================================
//              Hilbert-curve (a space-filling Peano curve) library
//=============================================================================
#include <Python.h>


#ifndef HILBERT_H
#define HILBERT_H

typedef unsigned long coord_t;	// char, short, int for upto 8, 16, 32 bits per word

#ifdef __cplusplus
extern "C" {
#endif
	// Functions: TransposetoAxes
	//            AxestoTranspose
	// Purpose:   Transform in-place between Hilbert transpose and geometrical axes
	// Example:   b=5 bits for each of n=3 coordinates.
	//            15-bit Hilbert integer = A B C D E F G H I J K L M N O  is stored as its Transpose
	//                   X[0] = A D G J M                X[2]|
	//                   X[1] = B E H K N    <------->       | /X[1]
	//                   X[2] = C F I L O               axes |/
	//                          high  low                    0------ X[0]
	//            Axes are stored conventially as b-bit integers.
	//-----------------------------------------------------------------------------
	static PyObject* TransposetoAxes(PyObject*, PyObject*);	// # bits, dimension, position
	static PyObject* AxestoTranspose(PyObject*, PyObject*);	// # bits, dimension, position
#ifdef __cplusplus
};
#endif

#endif