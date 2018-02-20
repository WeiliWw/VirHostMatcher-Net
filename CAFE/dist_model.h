/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#ifndef _DIST_MODEL_H
#define	_DIST_MODEL_H

#include "kmer.h"

enum dist{
	D2, D2STAR, D2SHEPP, CVtree, 
	Ch, Eu, Ma, FFP,
	CHISQ, JS, Co_Phylog,
	COSINE, PEARSON, CANBERRA, HAMMING,
	MATCHING, JACCARD, TANIMOTO, DICE, ANTIDICE, SNEATH, HAMMAN, PHI, ANDERBERG, GOWER, RUSSEL, YULE, OCHIAI, KULCZYNSKI
};

class L1FreqStrategy : public AbsTupleDistStrategy
{
public:
	L1FreqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff); d_result += d_tmp_diff; }
	double getDist(){ return d_result; }

private:
	double d_result;
};

class L2FreqStrategy : public AbsTupleDistStrategy
{
public:
	L2FreqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; d_tmp_diff = d_tmp_diff*d_tmp_diff; d_result += d_tmp_diff; }
	double getDist(){ return sqrt(d_result); }

private:
	double d_result;
};

class FFPStrategy : public AbsTupleDistStrategy
{
public:
	FFPStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) 
	{ 
		if(0 == src_X_w || 0 == trgt_X_w) return;
		d_result += src_X_w*(log(src_X_w)-log(trgt_X_w)) / LOG2; d_result += trgt_X_w*(log(trgt_X_w)-log(src_X_w)) / LOG2;
	}
	double getDist(){ return 0.5*d_result; }

private:
	double d_result;
};

class LInfFreqStrategy : public AbsTupleDistStrategy
{
public:
	LInfFreqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff); d_result = std::max(d_result, d_tmp_diff); }
	double getDist(){ return d_result; }

private:
	double d_result;
};

class PearsonStrategy : public AbsTupleDistStrategy
{
public:
	PearsonStrategy(int i_arg_k, bool b_arg_singleStrain) :AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ sum_src_X_w_trgt_X_w = 0; sum_src_X_w_sq = 0; sum_trgt_X_w_sq = 0; sum_src_X_w = 0; sum_trgt_X_w = 0; sum_count=0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { sum_src_X_w_trgt_X_w += (src_X_w * trgt_X_w); sum_src_X_w_sq += (src_X_w * src_X_w); sum_trgt_X_w_sq += (trgt_X_w * trgt_X_w); sum_src_X_w += src_X_w; sum_trgt_X_w += trgt_X_w; sum_count++; }
	double getDist(){ return 1.0-(sum_src_X_w_trgt_X_w-sum_src_X_w*sum_trgt_X_w/sum_count)/(sqrt(sum_src_X_w_sq-sum_src_X_w*sum_src_X_w/sum_count)*sqrt(sum_trgt_X_w_sq-sum_trgt_X_w*sum_trgt_X_w/sum_count)); }

private:
	double sum_src_X_w_trgt_X_w, sum_src_X_w_sq, sum_trgt_X_w_sq, sum_src_X_w, sum_trgt_X_w, sum_count;
};

class CanberraStrategy : public AbsTupleDistStrategy
{
public:
	CanberraStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff); d_result += (d_tmp_diff/(src_X_w + trgt_X_w)); }
	double getDist(){ return d_result; }

private:
	double d_result;
};

class HammingStrategy : public AbsTupleDistStrategy
{
public:
	HammingStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; sum = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { sum++; if((src_X_w>0 && 0==trgt_X_w) || (trgt_X_w>0 && 0==src_X_w)) d_result++;}
	double getDist(){ return d_result/sum; }

private:
	double d_result, sum;
};

class AbsBinaryTupleStrategy : public AbsTupleDistStrategy
{
public:
	AbsBinaryTupleStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ A = 0; B = 0; C = 0; D = 0; N = 0;}
	virtual void dealWithTuple(double src_X_w, double trgt_X_w) { N++; if(src_X_w>0 && trgt_X_w>0) A++; if(src_X_w>0 && 0==trgt_X_w) B++; if(0==src_X_w && trgt_X_w>0) C++; if(0==src_X_w && 0==trgt_X_w) D++;}

public:
	double A,B,C,D,N;
};

class MatchingStrategy : public AbsBinaryTupleStrategy
{
public:
	MatchingStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-(A+D)/N; }
};

class JaccardStrategy : public AbsBinaryTupleStrategy
{
public:
	JaccardStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-A/(N-D); }
};

class TanimotoStrategy : public AbsBinaryTupleStrategy
{
public:
	TanimotoStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-(A+D)/((A+D)+2*(B+C)); }
};

class DiceStrategy : public AbsBinaryTupleStrategy
{
public:
	DiceStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-2*A/(2*A+B+C); }
};

class AntidiceStrategy : public AbsBinaryTupleStrategy
{
public:
	AntidiceStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-A/(A+B+C); }
};

class SneathStrategy : public AbsBinaryTupleStrategy
{
public:
	SneathStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-2*(A+D)/(2*(A+D)+(B+C)); }
};

class HammanStrategy : public AbsBinaryTupleStrategy
{
public:
	HammanStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-(((A+D)-(B+C))/N)*(((A+D)-(B+C))/N); }
};

class PhiStrategy : public AbsBinaryTupleStrategy
{
public:
	PhiStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-1-((A*D-B*C)/sqrt((A+B)*(A+C)*(D+B)*(D+C)))*((A*D-B*C)/sqrt((A+B)*(A+C)*(D+B)*(D+C))); }
};

class AnderbergStrategy : public AbsBinaryTupleStrategy
{
public:
	AnderbergStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-(A/(A+B)+A/(A+C)+D/(C+D)+D/(B+D))/4; }
};

class GowerStrategy : public AbsBinaryTupleStrategy
{
public:
	GowerStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-A*D/sqrt((A+B)*(A+C)*(D+B*(D+C))); }
};

class RusselStrategy : public AbsBinaryTupleStrategy
{
public:
	RusselStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-A/N; }
};

class YuleStrategy : public AbsBinaryTupleStrategy
{
public:
	YuleStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-((A*D-B*C)/(A*D+B*C))*((A*D-B*C)/(A*D+B*C)); }
};

class OchiaiStrategy : public AbsBinaryTupleStrategy
{
public:
	OchiaiStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-A/sqrt((A+B)*(A+C)); }
};

class KulczynskiStrategy : public AbsBinaryTupleStrategy
{
public:
	KulczynskiStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsBinaryTupleStrategy(i_arg_k, b_arg_singleStrain){ }
	double getDist(){ return 1-(A/(A+B)+A/(A+C))/2; }
};

class D2Strategy : public AbsTupleDistStrategy
{
public:
	D2Strategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ sum_src_X_w_trgt_X_w = 0; sum_src_X_w_sq = 0; sum_trgt_X_w_sq = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { sum_src_X_w_trgt_X_w += (src_X_w * trgt_X_w); sum_src_X_w_sq += (src_X_w * src_X_w); sum_trgt_X_w_sq += (trgt_X_w * trgt_X_w); }
	double getDist(){ return 1.0 - sum_src_X_w_trgt_X_w / (sqrt(sum_src_X_w_sq)*sqrt(sum_trgt_X_w_sq)); }

private:
	double sum_src_X_w_trgt_X_w, sum_src_X_w_sq, sum_trgt_X_w_sq;
};

class D2starStrategy : public AbsQuadStrategy
{
public:
	D2starStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum_numerator = 0; sum_src_X_w_tilde_sq_div_EX_w = 0; sum_trgt_X_w_tilde_sq_div_EX_w = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return 0.5*(1.0 - sum_numerator / (sqrt(sum_src_X_w_tilde_sq_div_EX_w)*sqrt(sum_trgt_X_w_tilde_sq_div_EX_w))); }

private:
	double sum_numerator, sum_src_X_w_tilde_sq_div_EX_w, sum_trgt_X_w_tilde_sq_div_EX_w;
};

class D2sheppStrategy : public AbsQuadStrategy
{
public:
	D2sheppStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum_numerator = 0; sum_src_X_w_tilde_sq_div_sqr_sum = 0; sum_trgt_X_w_tilde_sq_div_sqr_sum = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return 0.5*(1.0 - sum_numerator / (sqrt(sum_src_X_w_tilde_sq_div_sqr_sum)*sqrt(sum_trgt_X_w_tilde_sq_div_sqr_sum))); }

private:
	double sum_numerator, sum_src_X_w_tilde_sq_div_sqr_sum, sum_trgt_X_w_tilde_sq_div_sqr_sum;
};

class HaoStrategy : public AbsQuadStrategy
{
public:
	HaoStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum_numerator = 0; sum_sq_src_X_w_tilde_div_EX_w = 0; sum_sq_trgt_X_w_tilde_div_EX_w = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return 0.5*(1.0 - sum_numerator / (sqrt(sum_sq_src_X_w_tilde_div_EX_w)*sqrt(sum_sq_trgt_X_w_tilde_div_EX_w))); }

private:
	double sum_numerator, sum_sq_src_X_w_tilde_div_EX_w , sum_sq_trgt_X_w_tilde_div_EX_w;
};

class ChiSqStrategy : public AbsQuadStrategy
{
public:
	ChiSqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return sum; }
private:
	double sum ;
};

class JensenShannonStrategy : public AbsMrkvStrategy
{
public:
	JensenShannonStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsMrkvStrategy(i_arg_k, b_arg_singleStrain){ sum_entropy = 0; src_entropy = 0; trgt_entropy = 0; }
	void dealWithMrkv(MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel);
	double getDist(){ return sqrt( - sum_entropy + (src_entropy + trgt_entropy) / 2); }

private:
	double sum_entropy, src_entropy, trgt_entropy;
};


class DistFactory
{
public:
	static DistFactory *getInstance();

	double getL1dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getL2dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getLInfdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getChiSqdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getFFPdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getCoPhylogdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getPearsondist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getCanberradist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getHammingdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getMatchingdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getJaccarddist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getTanimotodist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getDicedist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getAntidicedist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getSneathdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getHammandist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getPhidist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getAnderbergdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getGowerdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getRusseldist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getYuledist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getOchiaidist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getKulczynskidist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);


	double getD2dist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getD2stardist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, 
		int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getD2sheppdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
		int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getHaodist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
		std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getJensenShannondist(int i_arg_k, bool b_arg_singleStrain, 
		int i_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);


private:
	DistFactory(){}
	static DistFactory* instance;
};




#endif