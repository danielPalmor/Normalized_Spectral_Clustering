#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdlib.h>
#include "spkmeans.h"
#include "matrixOperations.h"

#define EPSILON 0
#define MAX_ITER 300

/**
 * CAPI of Spectral Kmeans
 */
static PyObject *fit(PyObject *self, PyObject *args)
{
    double **xMatrix;
    double **result;

    int k, goal, numOfPoints, dim;
    int i, j, rowResult, colResult;

    PyObject * py_x_matrix;
    PyObject * py_result;


    if (!PyArg_ParseTuple(args, "Oiiii", &py_x_matrix, &goal, &numOfPoints,
                          &dim, &k))
        return NULL;

    xMatrix = matrixAllocation(numOfPoints, dim);

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim; j++)
        {
            xMatrix[i][j] = PyFloat_AsDouble(
                    PyList_GetItem(py_x_matrix, j + (dim * i)));
        }
    }

    if (goal == 1)/* spk */
    {
        result = spectralKmeans(xMatrix, numOfPoints, dim, &k);
    }
    else if (goal == 2)/* wam */
    {
        result = weightedAdjacencyMatrix(xMatrix, numOfPoints, dim);
    }
    else if (goal == 3)/* ddg */
    {
        result = diagonalDegreeMatrix(xMatrix, numOfPoints, dim);
    }
    else if (goal == 4)/* lnorm */
    {
        result = lNorm(xMatrix, numOfPoints, dim);
    }
    else if (goal == 5)/* jacobi */
    {
        result = jacobi(xMatrix, dim);
    }
    else
    {
        freeMatrix(xMatrix);
        printf("An Error Has Occurred");
        exit(1);
    }

    if (goal == 1)/* spk */
    {
        rowResult = numOfPoints;
        colResult = k;
        py_result = PyList_New((rowResult * colResult) + 1);
        PyList_SetItem(py_result, (rowResult * colResult),
                       PyFloat_FromDouble(k));
    }
    else if (goal == 5)/* jacobi */
    {
        rowResult = numOfPoints + 1;
        colResult = numOfPoints;
        py_result = PyList_New(rowResult * colResult);
    }
    else
    {
        rowResult = numOfPoints;
        colResult = numOfPoints;
        py_result = PyList_New(rowResult * colResult);
    }


    for (i = 0; i < rowResult; i++)
    {
        for (j = 0; j < colResult; j++)
        {
            PyList_SetItem(py_result, (i * colResult) + j,
                           PyFloat_FromDouble(result[i][j]));
        }
    }

    freeMatrix(result);
    freeMatrix(xMatrix);
    return Py_BuildValue("O", py_result);
}

/**
 * CAPI of kmeans
 */
static PyObject *kmeans_fit(PyObject *self, PyObject *args)
{

    PyObject * matrix_points;
    PyObject * initial_centroids;

    int i, j, k, dim, numOfPoints;
    double **cMatrix, **cMatrixInitialCentroids;

    double **result;
    PyObject * centroids;


    if (!PyArg_ParseTuple(args, "iOOii", &k, &matrix_points, &initial_centroids,
                          &dim, &numOfPoints))
        return NULL;

    cMatrix = matrixAllocation(numOfPoints, dim);
    cMatrixInitialCentroids = matrixAllocation(k, dim);

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim; j++)
        {
            cMatrix[i][j] = PyFloat_AsDouble(
                    PyList_GetItem(matrix_points, j + (dim * i)));
        }
    }

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dim; j++)
        {
            cMatrixInitialCentroids[i][j] = PyFloat_AsDouble(
                    PyList_GetItem(initial_centroids, j + (dim * i)));
        }
    }

    result = kmeans(k, MAX_ITER, EPSILON, cMatrix, cMatrixInitialCentroids,
                    numOfPoints, dim);

    centroids = PyList_New(k * dim);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dim; j++)
        {
            PyList_SetItem(centroids, (i * dim) + j,
                           PyFloat_FromDouble(result[i][j]));
        }
    }

    freeMatrix(cMatrixInitialCentroids);
    return Py_BuildValue("O", centroids);
}

static PyMethodDef
kmeansMethods[] = {{"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("")},
                   {"kmeans_fit", (PyCFunction) kmeans_fit,METH_VARARGS,
                    PyDoc_STR("")},
                   {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduleDef = {PyModuleDef_HEAD_INIT, "spectralkmeans",
                                       NULL, -1, kmeansMethods};

PyMODINIT_FUNC PyInit_spectralkmeans(void)
{
    PyObject * m;
    m = PyModule_Create(&moduleDef);
    if (!m)
        return NULL;
    return m;
}

