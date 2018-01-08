#include "hilbert.cpp"
#include <iostream>

using namespace std;


static PyObject* some_func(PyObject* self, PyObject* args)
{
	__int64 input_value;
	if (!PyArg_ParseTuple(args, "L", &input_value))
		return 0;
	return PyLong_FromLongLong(input_value + 1);
}


PyMODINIT_FUNC PyInit_Hilbert(void)
{
	static PyMethodDef HilbertMethods[] = {
		{ "add_one", some_func, METH_VARARGS , NULL },
		{ "t2r", TransposetoAxes, METH_VARARGS, "Transpose to Real" },
		{ "r2t", AxestoTranspose, METH_VARARGS, "Real to Transpose" },
		{ NULL, NULL, 0/* flag telling the interpreter the calling convention to be used for the C function */, NULL }		/* Sentinel */
	};

	static struct PyModuleDef hilbertmodule = {
		PyModuleDef_HEAD_INIT,
		"Hilbert",	/* name of module */
		NULL,		/* module documentation, may be NULL */
		-1,			/* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
		HilbertMethods	/* array of exported functions */
	};

	PyObject *m = PyModule_Create(&hilbertmodule);
	if (m == NULL)
		return NULL;
	return m;
}