#include "hilbert.cpp"


PyMODINIT_FUNC PyInit_Hilbert(void)
{
	static PyMethodDef HilbertMethods[] = {
		{ "h2r", LinetoAxes, METH_VARARGS, "Hilbert line to real coordinates" },
		{ "r2h", AxestoLine, METH_VARARGS, "Real coordinates to Hilbert line" },
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