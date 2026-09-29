/* Test-only NumPy allocator: its returned pointers cannot be passed to free(). */
#define PY_SSIZE_T_CLEAN
#define NPY_TARGET_VERSION NPY_1_22_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static size_t outstanding = 0;
static size_t allocations = 0;

static void *shift_malloc(void *ctx, size_t size) {
    char *base;
    if (size > SIZE_MAX - 64) return NULL;
    base = (char *)malloc(size + 64);
    if (base == NULL) return NULL;
    ++outstanding;
    ++allocations;
    return base + 64;
}

static void *shift_calloc(void *ctx, size_t n, size_t size) {
    void *ptr;
    if (size != 0 && n > SIZE_MAX / size) return NULL;
    ptr = shift_malloc(ctx, n * size);
    if (ptr != NULL) memset(ptr, 0, n * size);
    return ptr;
}

static void shift_free(void *ctx, void *ptr, size_t size) {
    if (ptr != NULL) {
        --outstanding;
        free((char *)ptr - 64);
    }
}

static void *shift_realloc(void *ctx, void *ptr, size_t size) {
    char *base;
    if (ptr == NULL) return shift_malloc(ctx, size);
    if (size > SIZE_MAX - 64) return NULL;
    base = (char *)realloc((char *)ptr - 64, size + 64);
    return base == NULL ? NULL : base + 64;
}

static PyDataMem_Handler handler = {
    "simcoon_test_shifted", 1,
    {NULL, shift_malloc, shift_calloc, shift_realloc, shift_free}
};

static PyObject *install(PyObject *self, PyObject *unused) {
    PyObject *capsule = PyCapsule_New(&handler, "mem_handler", NULL);
    PyObject *previous;
    if (capsule == NULL) return NULL;
    previous = PyDataMem_SetHandler(capsule);
    Py_DECREF(capsule);
    return previous;
}

static PyObject *restore(PyObject *self, PyObject *previous) {
    return PyDataMem_SetHandler(previous);
}

static PyObject *count(PyObject *self, PyObject *unused) {
    return PyLong_FromSize_t(outstanding);
}

static PyObject *total(PyObject *self, PyObject *unused) {
    return PyLong_FromSize_t(allocations);
}

static PyMethodDef methods[] = {
    {"install", install, METH_NOARGS, NULL},
    {"restore", restore, METH_O, NULL},
    {"count", count, METH_NOARGS, NULL},
    {"total", total, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL}
};
static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_numpy_allocator_probe", NULL, -1, methods
};
PyMODINIT_FUNC PyInit__numpy_allocator_probe(void) {
    import_array();
    return PyModule_Create(&module);
}
