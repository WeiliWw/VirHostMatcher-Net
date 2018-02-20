/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#ifndef _Utils_H
#define	_Utils_H

#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <iomanip> 
#include <time.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <tuple>
#include <map>
#include <unordered_map>


#define BASE 4
#define MAX_ORDER 10
#define LOG2 log(2)

// if the char char_arg_nt violates ACGTOnly condition, return -1; Otherwise, 
// map A/a, C/c, G/g, T/t to 0, 1, 2, 3 respectively.
int nt2int(char char_arg_nt);

// if the char char_arg_nt violates ACGTOnly condition, return N; Otherwise, 
// map A/a, C/c, G/g, T/t to T, G, C, A respectively.
char nt2ComplementNt(char char_arg_nt);

// if the seq violates ACGTOnly condition, return -1; Otherwise, map the seq
// into an integer ranging from 0 to 4^{length(seq)}-1, i.e. the corresponding
// index in the Markov transition matrix or k-mer vector. For example, 
// nt2index("AAAA")=0, nt2index("TAAA")=64, nt2index("TTTT")=255, etc.
// 
// str_arg_seq: the m-mer in m-order Markov model or k-mer seq
unsigned long long nt2index(const std::string & str_arg_seq);

unsigned long long index2revCompleIdx(unsigned long long i_arg_index, int i_arg_seqLength);

// the inverse operation of nt2index
// 
// i_arg_index: the index in the Markov transition matrix or k-mer vector
// i_arg_seqLength: the length of seq
std::string index2nt(unsigned long long i_arg_index, int i_arg_seqLength);

// reverse complement of the input str, expect to satisfy ACGTOnly condition
std::string revComplementStr(const std::string & currStr);

//complement of the input str, expect to satisfy ACGTOnly condition
std::string complementStr(const std::string & currStr);

// reverse of the input str
std::string revStr(const std::string & currStr);

std::string trim(std::string currStr);

void split(const std::string& currStr, std::string delim, std::vector<std::string > & ret);

std::string toLowerCase(std::string currStr);

bool endsWith(std::string const & currStr, std::string const & ending);

/*unsigned long long updateIdx(unsigned long long i_arg_currIdx, const std::string & str_arg_appendStr)
{
	unsigned long long i_tmp_newIdx = i_arg_currIdx;

	for (unsigned int i = 0; i < str_arg_appendStr.size(); ++i)
	{
		i_tmp_newIdx = (i_tmp_newIdx << 2) + nt2int(str_arg_appendStr.at(i));

		if (i_tmp_newIdx < i_arg_currIdx)
			throw std::runtime_error(" left shifting cause overflow! ");
	}
	return i_tmp_newIdx;
}*/

bool file_exists(const std::string & str_arg_filename);
bool dir_exists(const std::string & str_arg_directory);
std::string getFileName(const std::string & str_arg_path);

template <typename T> void free_vec_ptr(std::vector <T* > & v)
{
	int size = v.size();
	T* p = NULL;
	for (int i = 0; i < size; i++)
	{
		p = v[i];
		delete[] p;
	}
	v.clear();
}

// find the max and argmax in an array
template <typename T> T max(const T * x, int n, int* argmax)
{
	*argmax = 0;
	T max_val = x[0];
	for (int i = 1; i < n; i++)
		if (x[i] > max_val)
		{
		max_val = x[i];
		*argmax = i;
		}
	return max_val;
}

// find the max and argmax in an vector
template <typename T> T max_vec(std::vector<T> & v, int n, int* argmax)
{
	*argmax = 0;
	T max_val = v[0];
	for (int i = 1; i < n; i++)
		if (v[i] > max_val)
		{
		max_val = v[i];
		*argmax = i;
		}
	return max_val;
}


// find the min and argmin in an vector
template <typename T> T min_vec(std::vector<T> & v, int n, int* argmin)
{
	*argmin = 0;
	T min_val = v[0];
	for (int i = 1; i < n; i++)
		if (v[i] < min_val)
		{
		min_val = v[i];
		*argmin = i;
		}
	return min_val;
}

//given log(a) and log(b), return log(a + b)
double log_sum(double log_a, double log_b);

// give a_1, ..., a_n,
// return log(exp(a_1)+...+exp(a_n))
double log_normalize(double * array, int nlen);

// the vector version
double log_normalize(std::vector<double> & vec, int nlen);

//given log(a) and log(b), return log(a - b) a>b
double log_subtract(double log_a, double log_b);

bool almostEquals(double a, double b);

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

#endif