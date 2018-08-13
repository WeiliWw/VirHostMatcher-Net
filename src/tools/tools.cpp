#include "Python.h"
#include <numpy/arrayobject.h>
#include "kmer_count_multithreads.h" 
#include "WIsH.h"
#include <atomic>
#include <string>

static PyObject *wish(PyObject *self, PyObject *args)
{
   char *genomeDir, *modelDir, *resultDir, *command;
   unsigned int thread;
   int res;
  if (!PyArg_ParseTuple(args, "ssssi", &genomeDir, &modelDir, &resultDir, &command, &thread))
     return NULL;
   res = wish(genomeDir,modelDir,resultDir,command,thread);
    return PyLong_FromLong(res);
}


static PyObject *kmer_count(PyObject *self, PyObject *args)
{
    char* filename;
    int K;
    int Num_Threads;
    bool Reverse;
    if (!PyArg_ParseTuple(args, "siii", &filename, &Num_Threads, &Reverse, &K))
        return NULL;
    //npy_intp SIZE = pow(4, K);
    //std::vector<int> count_res;
    std::vector<std::atomic<int>> count_array;
    //std::cout<<filename<<std::endl;
    //count_res = count_vector(filename, K, Reverse);
    count_array = count(filename, K, Reverse, Num_Threads);
    //std::cout << SIZE << count_array.size() << std::endl;
    /*
    for (int i=0;i<10;i++){
        std::cout<<count_array[i]<<std::endl;
    }
    */
    int SIZE = pow(4, K);
    PyObject *result = PyList_New(SIZE);
    for (int i=0;i<SIZE;i++){
	  PyList_SetItem(result, i, PyLong_FromLong(count_array[i]));
    }
    //result = PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_res.data()));
    return result;
}

static PyMethodDef module_methods[] = {
    {"kmer_count", kmer_count, METH_VARARGS, ""},
    {"wish", wish, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef toolsmodule = {
    PyModuleDef_HEAD_INIT,
    "tools",   /* name of module */
    "", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_tools(void)
{
    import_array();
    return PyModule_Create(&toolsmodule);
}

