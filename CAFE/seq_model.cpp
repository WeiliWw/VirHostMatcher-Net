#include "seq_model.h"
#include <cmath>
#include <cstring>

MarkovModel::MarkovModel(int i_arg_order) { i_order = i_arg_order; i_rowDim = 0; }

MarkovModel::~MarkovModel()
{
	for (unsigned long long i = 0; i < i_rowDim; ++i) delete[] arr_obvTransProb[i];
	delete[] arr_obvTransProb; delete[] vec_obvMargProb; delete kmerIdxRowIdxTable;
}

void MarkovModel::constructMarg(unsigned long long* vec_arg_orderkmerIdx, unsigned long long* vec_arg_orderkmerCnt, unsigned long long l_arg_dim)
{
	kmerIdxRowIdxTable = new std::unordered_map<unsigned long long, unsigned long long>();

	i_rowDim = l_arg_dim;
	vec_obvMargProb = new double[i_rowDim]; memset(vec_obvMargProb, 0, sizeof(double) * i_rowDim);

	for (unsigned long long rowIdx = 0; rowIdx < l_arg_dim; ++rowIdx)
	{
		unsigned long long currKmerIdx = vec_arg_orderkmerIdx[rowIdx];
		(*kmerIdxRowIdxTable)[currKmerIdx] = rowIdx;
		vec_obvMargProb[rowIdx] = vec_arg_orderkmerCnt[rowIdx];
	}

	arr_obvTransProb = new double*[i_rowDim];
	for (unsigned long long i = 0; i < i_rowDim; ++i)
	{
		arr_obvTransProb[i] = new double[BASE];
		memset(arr_obvTransProb[i], 0, sizeof(double) * BASE);
	}
}

void MarkovModel::normalize()
{
	double d_tmp_totalCnt = 0;
	for (std::unordered_map<unsigned long long, unsigned long long>::iterator iter = kmerIdxRowIdxTable->begin(); iter != kmerIdxRowIdxTable->end(); iter++)
	{
		unsigned long long currKmerIdx = iter->first;
		unsigned long long newIdx = (*kmerIdxRowIdxTable)[currKmerIdx];

		d_tmp_totalCnt += vec_obvMargProb[newIdx];
		double d_tmp_localCnt = 0;
		for (unsigned long long j = 0; j < BASE; ++j) d_tmp_localCnt += arr_obvTransProb[newIdx][j];

		if (d_tmp_localCnt > 0)
		{
			//for (unsigned long long j = 0; j < BASE; ++j) arr_obvTransProb[newIdx][j] /= d_tmp_localCnt;
			for (unsigned long long j = 0; j < BASE; ++j)
			{
				if (arr_obvTransProb[newIdx][j] > 0) arr_obvTransProb[newIdx][j] = log(arr_obvTransProb[newIdx][j]) - log(d_tmp_localCnt);
			}
		}
	}

	if (d_tmp_totalCnt > 0)
	{
		double d_tmp_totalCnt_log = log(d_tmp_totalCnt);
		for (std::unordered_map<unsigned long long, unsigned long long>::iterator iter = kmerIdxRowIdxTable->begin(); iter != kmerIdxRowIdxTable->end(); iter++)
		{
			unsigned long long currKmerIdx = iter->first;
			unsigned long long newIdx = (*kmerIdxRowIdxTable)[currKmerIdx];
			//vec_obvMargProb[newIdx] /= d_tmp_totalCnt;
			if (vec_obvMargProb[newIdx] > 0) vec_obvMargProb[newIdx] = log(vec_obvMargProb[newIdx]) - d_tmp_totalCnt_log;
		}
	}

}

void MarkovModel::print()
{
	std::vector<unsigned long long> idxVec;
	for (std::unordered_map<unsigned long long, unsigned long long>::iterator iter = kmerIdxRowIdxTable->begin(); iter != kmerIdxRowIdxTable->end(); iter++) idxVec.push_back(iter->first);
	std::sort(idxVec.begin(), idxVec.end());

	for (std::vector<unsigned long long>::iterator iter = idxVec.begin(); iter != idxVec.end(); iter++)
	{
		unsigned long long currKmerIdx = *iter;
		unsigned long long newIdx = (*kmerIdxRowIdxTable)[currKmerIdx];
		std::cout << currKmerIdx << "\t" << vec_obvMargProb[newIdx] << "\t";
		std::cout << arr_obvTransProb[newIdx][0] << " " << arr_obvTransProb[newIdx][1] << " " << arr_obvTransProb[newIdx][2] << " " << arr_obvTransProb[newIdx][3] << std::endl;
	}

}

void MarkovModel::addMargProb(unsigned long long l_arg_currKmerIdx, double d_value)
{
	if (kmerIdxRowIdxTable->find(l_arg_currKmerIdx) == kmerIdxRowIdxTable->end()) return;
	vec_obvMargProb[(*kmerIdxRowIdxTable)[l_arg_currKmerIdx]] += d_value;
}

double MarkovModel::getMargProb(unsigned long long l_arg_currKmerIdx)
{
	if (kmerIdxRowIdxTable->find(l_arg_currKmerIdx) == kmerIdxRowIdxTable->end()) return 0;
	return vec_obvMargProb[(*kmerIdxRowIdxTable)[l_arg_currKmerIdx]];
}

void MarkovModel::addTransProb(unsigned long long l_arg_currKmerIdx, unsigned long long i_arg_route, double d_value)
{
	if (kmerIdxRowIdxTable->find(l_arg_currKmerIdx) == kmerIdxRowIdxTable->end()) return;
	arr_obvTransProb[(*kmerIdxRowIdxTable)[l_arg_currKmerIdx]][i_arg_route] += d_value;
}

double MarkovModel::getTransProb(unsigned long long l_arg_currKmerIdx, unsigned long long i_arg_route)
{
	if (kmerIdxRowIdxTable->find(l_arg_currKmerIdx) == kmerIdxRowIdxTable->end()) return 0;
	return arr_obvTransProb[(*kmerIdxRowIdxTable)[l_arg_currKmerIdx]][i_arg_route];
}
