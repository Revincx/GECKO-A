/* Python bindings for GECKO-A C modules
 * This module provides Python access to keyparameter and keyflag modules
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "keyparameter.h"
#include "keyflag.h"
#include "minidict.h"
#include "toolbox.h"

/* Python wrapper for define_defaults */
static PyObject* py_define_defaults(PyObject* self, PyObject* args) {
    (void)self;  /* Unused */
    (void)args;  /* No arguments */
    
    define_defaults();
    Py_RETURN_NONE;
}

/* Python wrapper for getting configuration values */
static PyObject* py_get_config(PyObject* self, PyObject* args) {
    (void)self;  /* Unused */
    (void)args;  /* No arguments */
    
    PyObject* dict = PyDict_New();
    if (!dict) return NULL;
    
    PyDict_SetItemString(dict, "dirgecko", PyUnicode_FromString(dirgecko));
    PyDict_SetItemString(dict, "dirout", PyUnicode_FromString(dirout));
    PyDict_SetItemString(dict, "maxgen", PyLong_FromLong(maxgen));
    PyDict_SetItemString(dict, "TK", PyFloat_FromDouble(TK));
    PyDict_SetItemString(dict, "critvp", PyFloat_FromDouble(critvp));
    PyDict_SetItemString(dict, "brcut", PyFloat_FromDouble(brcut));
    PyDict_SetItemString(dict, "yldcut", PyFloat_FromDouble(yldcut));
    
    return dict;
}

/* Python wrapper for minidict functions */
static PyObject* py_clean_minid(PyObject* self, PyObject* args) {
    (void)self;
    (void)args;
    
    clean_minid();
    Py_RETURN_NONE;
}

static PyObject* py_add_fo(PyObject* self, PyObject* args) {
    const char *nam, *formula;
    
    if (!PyArg_ParseTuple(args, "ss", &nam, &formula)) {
        return NULL;
    }
    
    add_fo(nam, formula);
    Py_RETURN_NONE;
}

static PyObject* py_get_fo(PyObject* self, PyObject* args) {
    const char *nam;
    char formula[MXLFO+1];
    
    if (!PyArg_ParseTuple(args, "s", &nam)) {
        return NULL;
    }
    
    get_fo(nam, formula);
    return PyUnicode_FromString(formula);
}

/* Python wrapper for countstring */
static PyObject* py_countstring(PyObject* self, PyObject* args) {
    const char *line, *str;
    
    if (!PyArg_ParseTuple(args, "ss", &line, &str)) {
        return NULL;
    }
    
    int count = countstring(line, str);
    return PyLong_FromLong(count);
}

/* Python wrapper for kval (Arrhenius rate constant) */
static PyObject* py_kval(PyObject* self, PyObject* args) {
    PyObject *list_obj;
    double T;
    double arrh[3];
    
    if (!PyArg_ParseTuple(args, "Od", &list_obj, &T)) {
        return NULL;
    }
    
    if (!PyList_Check(list_obj) || PyList_Size(list_obj) != 3) {
        PyErr_SetString(PyExc_ValueError, "arrh must be a list of 3 floats");
        return NULL;
    }
    
    for (int i = 0; i < 3; i++) {
        PyObject *item = PyList_GetItem(list_obj, i);
        arrh[i] = PyFloat_AsDouble(item);
        if (PyErr_Occurred()) return NULL;
    }
    
    double k = kval(arrh, T);
    return PyFloat_FromDouble(k);
}

/* Method definitions */
static PyMethodDef GeckoMethods[] = {
    {"define_defaults", py_define_defaults, METH_NOARGS, 
     "Initialize default parameters"},
    {"get_config", py_get_config, METH_NOARGS,
     "Get configuration dictionary"},
    {"clean_minid", py_clean_minid, METH_NOARGS,
     "Clean the mini dictionary"},
    {"add_fo", py_add_fo, METH_VARARGS,
     "Add a formula to mini dictionary: add_fo(name, formula)"},
    {"get_fo", py_get_fo, METH_VARARGS,
     "Get formula from mini dictionary: get_fo(name)"},
    {"countstring", py_countstring, METH_VARARGS,
     "Count substring occurrences: countstring(line, substring)"},
    {"kval", py_kval, METH_VARARGS,
     "Calculate Arrhenius rate constant: kval([A, n, Ea_R], T)"},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef geckomodule = {
    PyModuleDef_HEAD_INIT,
    "geckoa",
    "GECKO-A C extension module",
    -1,
    GeckoMethods
};

/* Module initialization */
PyMODINIT_FUNC PyInit_geckoa(void) {
    PyObject *m = PyModule_Create(&geckomodule);
    if (m == NULL) return NULL;
    
    /* Add constants */
    PyModule_AddIntConstant(m, "MXLCO", MXLCO);
    PyModule_AddIntConstant(m, "MXLFO", MXLFO);
    PyModule_AddIntConstant(m, "MXNODE", MXNODE);
    PyModule_AddIntConstant(m, "MXPS", MXPS);
    PyModule_AddIntConstant(m, "MXPD", MXPD);
    
    return m;
}
