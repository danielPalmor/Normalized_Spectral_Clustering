#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdlib.h>
#include "spkmeans.h"

/**
  NEW CHANGES HERE
  IDEAS: include MatrixOperations to prevent mult code
*/


static PyObject* fit(PyObject *self, PyObject *args)
{
    double **xMatrix, *xArray;
    double **result;

    int k, goal, numOfPoints, dim;
    int i, j, rowResult, colResult;

    PyObject* py_x_matrix;
    PyObject* py_result;


    if (!PyArg_ParseTuple(args, "Oiiii", &py_x_matrix, &goal, &numOfPoints, &dim, &k))
        return NULL;

    xArray = calloc(numOfPoints * dim, sizeof(double));
    if(xArray == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    xMatrix = calloc(numOfPoints, sizeof(double*));
    if (xMatrix == NULL)
    {
        printf("An Error Has Occurred");
        free(xArray);
        exit(1);
    }

    for (i = 0; i < numOfPoints; i++)
    {
        xMatrix[i] = xArray + i * dim;
    }

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim; j++){
            xMatrix[i][j] = PyFloat_AsDouble(PyList_GetItem(py_x_matrix, j+(dim*i)));
        }
    }

    if (goal == 1)
    {
        result = spectralKmeans(xMatrix, numOfPoints, dim, &k);
    }
    else if (goal == 2)
    {
        result = weightedAdjacencyMatrix(xMatrix, numOfPoints, dim);
    }
    else if (goal == 3)
    {
        result = diagonalDegreeMatrix(xMatrix, numOfPoints, dim);
    }
    else if (goal == 4)
    {
        result = lNorm(xMatrix, numOfPoints, dim);
    }
    else if (goal == 5)
    {
        result = jacobi(xMatrix, dim);
    }


    if (goal == 1)
    {
        rowResult = numOfPoints;
        colResult = k;
        py_result = PyList_New((rowResult * colResult) + 1);
        PyList_SetItem(py_result, (rowResult * colResult),PyFloat_FromDouble(k));
    }
    else if (goal == 5)
    {
        rowResult = numOfPoints+1;
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
            PyList_SetItem(py_result, (i*colResult)+j,PyFloat_FromDouble(result[i][j]));
        }
    }

    free(result[0]);
    free(result);
    free(xArray);
    free(xMatrix);
    return Py_BuildValue("O",py_result);
}

static PyObject* kmeans_fit(PyObject *self, PyObject *args)
{
    const int MAX_ITER = 300;
    const double EPS = 0;

    PyObject* matrix_points;
    PyObject* initial_centroids;

    int i, k, dim, numOfPoints;
    double **cMatrix, **cMatrixInitialCentroids;
    double *cArray, *cArrayInitialCentroids;

    double **result;
    PyObject *centroids;


    if (!PyArg_ParseTuple(args, "iOOii", &k, &matrix_points, &initial_centroids, &dim, &numOfPoints))
        return NULL;

    cArray = calloc(numOfPoints * dim, sizeof(double));
    if(cArray == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    cMatrix = calloc(numOfPoints, sizeof(double*));
    if(cMatrix == NULL)
    {
        free(cArray);
        printf("An Error Has Occurred");
        exit(1);
    }

    for (i = 0; i < numOfPoints; i++) {
        cMatrix[i] = cArray + i * dim;
    }

    cArrayInitialCentroids = calloc(k * dim, sizeof(double));
    if(cArrayInitialCentroids == NULL)
    {
        free(cArray);
        free(cMatrix);
        printf("An Error Has Occurred");
        exit(1);
    }
    cMatrixInitialCentroids = calloc(k, sizeof(double*));
    if(cMatrixInitialCentroids == NULL)
    {
        free(cArrayInitialCentroids);
        free(cArray);
        free(cMatrix);
        printf("An Error Has Occurred");
        exit(1);
    }

    for (i = 0; i < k; i++) {
        cMatrixInitialCentroids[i] = cArrayInitialCentroids + i * dim;
    }

    for (i = 0; i < numOfPoints; i++)
    {
        int j;
        for (j = 0; j < dim; j++){
            cMatrix[i][j] = PyFloat_AsDouble(PyList_GetItem(matrix_points, j+(dim*i)));
        }
    }

    for (i = 0; i < k; i++)
    {
        int j;
        for (j = 0; j < dim; j++){
            cMatrixInitialCentroids[i][j] = PyFloat_AsDouble(PyList_GetItem(initial_centroids, j+(dim*i)));
        }
    }

    result = kmeans(k, MAX_ITER, EPS, cMatrix, cMatrixInitialCentroids,
                    dim, numOfPoints);

    centroids =  PyList_New(k*dim);
    for (i = 0; i < k; i++)
    {
        int j;
        for (j = 0; j < dim; j++)
        {
            PyList_SetItem(centroids, (i*dim)+j,PyFloat_FromDouble(result[i][j]));
        }
    }
    free(cMatrixInitialCentroids);
    free(cArrayInitialCentroids);
    return Py_BuildValue("O",centroids);
}

static PyMethodDef kmeansMethods[] = {{"fit", (PyCFunction) fit, METH_VARARGS,PyDoc_STR("")},
                                      {"kmeans_fit", (PyCFunction) kmeans_fit, METH_VARARGS,PyDoc_STR("")},
                                      {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduleDef = {
        PyModuleDef_HEAD_INIT,
        "spectralkmeans",
        NULL,
        -1,
        kmeansMethods
};

PyMODINIT_FUNC
PyInit_spectralkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduleDef);
    if (!m)
        return NULL;
    return m;
}

