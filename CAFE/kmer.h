/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#ifndef _KMER_H
#define	_KMER_H

#include "seq_model.h"

class AbsIter
{
public:
	AbsIter(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap) { i_k = i_arg_k; currKmer = 0; kmerCntUnorderMap = arg_kmerCntUnorderMap; };
	virtual ~AbsIter() {}

	virtual void operator++() = 0;
	virtual bool hasNext() = 0;
	virtual double operator*() = 0;

	void operator++(int) { ++(*this); };
	unsigned long long getCurrKmer() { return currKmer; }

public:
	int i_k;
	unsigned long long currKmer;
	std::unordered_map<unsigned long long, unsigned long>* kmerCntUnorderMap;
};

class KmerCntTraverseIter : public AbsIter
{
public:
	KmerCntTraverseIter(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap) : AbsIter(i_arg_k, arg_kmerCntUnorderMap) { maxAllowedIdx = (unsigned long long)pow(BASE, i_arg_k) - 1; }
	void operator++() { currKmer++; };
	bool hasNext() { return (currKmer <= maxAllowedIdx); };
	double operator*() { return (*kmerCntUnorderMap)[currKmer]; };

public:
	unsigned long long maxAllowedIdx;
};

class AbsHashIter : public AbsIter
{
public:
	AbsHashIter(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec) : AbsIter(i_arg_k, arg_kmerCntUnorderMap) { currKmerVecIdx = 0; kmerVec = arg_kmerVec; currKmer = kmerVec->at(currKmerVecIdx); }
	void operator++() { currKmerVecIdx++; if (currKmerVecIdx < kmerVec->size()) currKmer = kmerVec->at(currKmerVecIdx); };
	bool hasNext() { return (currKmerVecIdx < kmerVec->size()); };
	virtual double operator*() = 0;

public:
	unsigned long long currKmerVecIdx;
	std::vector<unsigned long long>* kmerVec;
};

class KmerCntHashIter : public AbsHashIter
{
public:
	KmerCntHashIter(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec) : AbsHashIter(i_arg_k, arg_kmerCntUnorderMap, arg_kmerVec) {}
	double operator*() { return (*kmerCntUnorderMap)[currKmer]; }
};

class KmerFreqHashIter : public AbsHashIter
{
public:
	KmerFreqHashIter(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec, unsigned long l_arg_totalKmer) : AbsHashIter(i_arg_k, arg_kmerCntUnorderMap, arg_kmerVec) { totalKmerInv = 1.0 / l_arg_totalKmer; }
	double operator*() { return totalKmerInv * (*kmerCntUnorderMap)[currKmer]; }

private:
	double totalKmerInv;
};

class KmerProbDelegate
{
public:
	KmerProbDelegate(int i_arg_k, MarkovModel* arg_mrkvModel, bool b_arg_isRevCompl = false);

	void init();
	double getKmerlogProb(unsigned long long queryNextKmerIdx);

private:
	bool push(unsigned long long idx);
	bool push(std::stack<unsigned long long> & indices);
	bool increment();
	unsigned long long pop();

private:
	MarkovModel* markovModel;
	bool isEnd, isRevCompl;
	int i_k;
	unsigned long long nextPosKmerIdx, maxAllowedIdx;

	std::stack<unsigned long long> kmer_traceStack;
	std::stack<unsigned long long> orderIdx_traceStack;
	std::stack<double> logProbProd_traceStack;
	std::stack<unsigned long long> lowerIdx_traceStack, upperIdx_traceStack;
};

class KmerProbEnsembDelegate
{
public:
	KmerProbEnsembDelegate(int i_arg_k, MarkovModel* arg_mrkvModel, bool b_arg_singleStrain);
	~KmerProbEnsembDelegate();

	void init();
	double getKmerlogProb(unsigned long long queryNextKmerIdx);

private:
	bool b_singleStrain;
	int i_k;
	MarkovModel* markovModel;
	KmerProbDelegate *kmerProbDelegate, *revComplKmerProbDelegate;
};

class AbsDistStrategy
{
public:
	AbsDistStrategy(int i_arg_k, bool b_arg_singleStrain) { i_k = i_arg_k; b_singleStrain = b_arg_singleStrain; }
	virtual double getDist() = 0;

public:
	int i_k;
	bool b_singleStrain;
};


class AbsTupleDistStrategy : public AbsDistStrategy
{
public:
	AbsTupleDistStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsDistStrategy(i_arg_k, b_arg_singleStrain){}
	virtual void dealWithTuple(double src_X_w, double trgt_X_w) = 0;
};


class AbsQuadStrategy : public AbsDistStrategy
{
public:
	AbsQuadStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsDistStrategy(i_arg_k, b_arg_singleStrain){}
	virtual void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w) = 0;
};

class AbsMrkvStrategy : public AbsDistStrategy
{
public:
	AbsMrkvStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsDistStrategy(i_arg_k, b_arg_singleStrain){}
	virtual void dealWithMrkv(MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel) = 0;
};


class IterFactory
{
public:
	static IterFactory *getInstance();

	double getFreqDist(AbsTupleDistStrategy* distStrategy, 
		std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec, unsigned long src_totalKmer,
		std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec, unsigned long trgt_totalKmer);

	double getCntDist(AbsTupleDistStrategy* distStrategy, int i_arg_lowerCnt,
		std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec,
		std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec);

	//specific to ChiSq Distance
	double getCntDist(AbsQuadStrategy* distStrategy,
		std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec,
		std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec);

	double getCntExpDist(AbsQuadStrategy* distStrategy, int i_arg_lowerCnt, 
		std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec, MarkovModel* src_mrkvModel, unsigned long src_totalKmer,
		std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec, MarkovModel* trgt_mrkvModel, unsigned long trgt_totalKmer);

	double getMrkvDist(AbsMrkvStrategy* distStrategy, MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel);

	double getCoPhylogDist(int i_arg_k, 
		std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec,
		std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec);

	AbsIter* getKmerCntIterator(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec, int i_arg_lowerCnt);
	AbsIter* getKmerFreqIterator(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec, unsigned long l_arg_totalKmer);
	KmerProbEnsembDelegate* getKmerProbDelegate(int i_arg_k, bool b_arg_singleStrain, MarkovModel* arg_mrkvModel);

private:
	IterFactory(){}
	static IterFactory* instance;
};


class KmerModel
{
public:
	KmerModel(int i_arg_k, bool b_arg_singleStrain) { i_k = i_arg_k; b_singleStrain = b_arg_singleStrain; kmerCntUnorderMap = new std::unordered_map<unsigned long long, unsigned long>(); kmerVec = new std::vector<unsigned long long>(); }
	~KmerModel() { delete kmerCntUnorderMap; delete kmerVec; }

	bool load(int i_arg_k, std::string str_arg_inputURL);
	bool saveFromLargerK(int i_arg_k, int i_arg_larger_k, std::string str_arg_inputURL, std::string str_arg_outputURL);
	bool saveFromJellyFish(std::string str_arg_jfTxtURL, std::string str_arg_outputURL);
	bool saveFromFasta(int i_arg_k, std::string str_arg_fastaFileURL, std::string str_arg_outputURL);

	unsigned long totalKmer();
	MarkovModel* getMarkovModel(int i_arg_order, std::string str_arg_saveURLPrefix);

private:
	bool save(std::map<unsigned long long, unsigned long> *kmerCntMap, std::string str_arg_outputURL);

public:
	int i_k;
	bool b_singleStrain;
	std::unordered_map<unsigned long long, unsigned long> *kmerCntUnorderMap;
	std::vector<unsigned long long>* kmerVec;
};


int getEstMarkovOrder(int i_arg_k, std::string str_arg_saveURLPrefix, std::string str_arg_seqName);

#endif