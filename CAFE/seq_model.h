/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#ifndef _SEQ_MODEL_H
#define	_SEQ_MODEL_H

#include "utils.h"

class MarkovModel
{
public:
	MarkovModel(int i_arg_order);
	~MarkovModel();
	void normalize();
	void print();
	//void constructMarg(unsigned long long* vec_arg_orderDim);
	void constructMarg(unsigned long long* vec_arg_orderkmerIdx, unsigned long long* vec_arg_orderkmerCnt, unsigned long long l_arg_dim);
	void addMargProb(unsigned long long l_arg_currKmerIdx, double d_value);
	double getMargProb(unsigned long long l_arg_currKmerIdx);
	void addTransProb(unsigned long long l_arg_currKmerIdx, unsigned long long i_arg_route, double d_value);
	double getTransProb(unsigned long long l_arg_currKmerIdx, unsigned long long i_arg_route);
	int getOrder() { return i_order; }

private:
	int i_order;
	unsigned long long i_rowDim;
	double** arr_obvTransProb; //non-empty rows by 4 matrix
	double* vec_obvMargProb; //non-empty rows by 1 matrix
	std::unordered_map<unsigned long long, unsigned long long>* kmerIdxRowIdxTable;
};


#endif