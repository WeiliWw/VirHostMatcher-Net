#include "dist_model.h"

DistFactory* DistFactory::instance = 0;

DistFactory* DistFactory::getInstance()
{
	if (!instance) instance = new DistFactory();
	return instance;
}

void D2starStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	if (0 == src_EX_w || 0 == trgt_EX_w) return;

	double src_X_w_tilde = src_X_w - src_EX_w;
	double trgt_X_w_tilde = trgt_X_w - trgt_EX_w;
	double numerator = src_X_w_tilde * trgt_X_w_tilde / sqrt(src_EX_w*trgt_EX_w);
	double src_X_w_tilde_sq_div_EX_w = src_X_w_tilde*src_X_w_tilde / src_EX_w;
	double trgt_X_w_tilde_sq_div_EX_w = trgt_X_w_tilde*trgt_X_w_tilde / trgt_EX_w;

	sum_numerator += numerator;
	sum_src_X_w_tilde_sq_div_EX_w += src_X_w_tilde_sq_div_EX_w;
	sum_trgt_X_w_tilde_sq_div_EX_w += trgt_X_w_tilde_sq_div_EX_w;
}

void D2sheppStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	double src_X_w_tilde = src_X_w - src_EX_w;
	double trgt_X_w_tilde = trgt_X_w - trgt_EX_w;
	double src_X_w_tilde_sq = src_X_w_tilde*src_X_w_tilde;
	double trgt_X_w_tilde_sq = trgt_X_w_tilde*trgt_X_w_tilde;
	double denominator = sqrt(src_X_w_tilde_sq + trgt_X_w_tilde_sq);

	if (0 == denominator) return;

	sum_numerator += src_X_w_tilde * trgt_X_w_tilde / denominator;
	sum_src_X_w_tilde_sq_div_sqr_sum += src_X_w_tilde_sq / denominator;
	sum_trgt_X_w_tilde_sq_div_sqr_sum += trgt_X_w_tilde_sq / denominator;
}

void HaoStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	if (0 == src_EX_w || 0 == trgt_EX_w) return;

	double src_X_w_tilde_div_EX_w = ((double)src_X_w - src_EX_w) / src_EX_w;
	double trgt_X_w_tilde_div_EX_w = ((double)trgt_X_w - trgt_EX_w) / trgt_EX_w;

	sum_numerator += src_X_w_tilde_div_EX_w*trgt_X_w_tilde_div_EX_w;
	sum_sq_src_X_w_tilde_div_EX_w += src_X_w_tilde_div_EX_w*src_X_w_tilde_div_EX_w;
	sum_sq_trgt_X_w_tilde_div_EX_w += trgt_X_w_tilde_div_EX_w*trgt_X_w_tilde_div_EX_w;
}

void ChiSqStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	//double tmp1 = src_X_w_1 - src_X_w*all_X_w_1 / all_X_w; sum += (tmp1*tmp1*all_X_w / (src_X_w*all_X_w_1));
	double tmp1 = src_EX_w - src_X_w*trgt_EX_w / trgt_X_w; sum += (tmp1*tmp1*trgt_X_w / (src_X_w*trgt_EX_w));
}

void JensenShannonStrategy::dealWithMrkv(MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel)
{
	int i_src_order = src_mrkvModel->getOrder();

	for (unsigned long long i = 0; i < (unsigned long long)pow(BASE, i_src_order); ++i)
	{
		double sum_entropyOverCol = 0, src_entropyOverCol = 0, trgt_entropyOverCol = 0;

		for (unsigned int j = 0; j < BASE; ++j)
		{
			if (0 != src_mrkvModel->getTransProb(i, j)) src_entropyOverCol += exp(src_mrkvModel->getTransProb(i, j)) * src_mrkvModel->getTransProb(i, j) / LOG2;
			if (0 != trgt_mrkvModel->getTransProb(i, j)) trgt_entropyOverCol += exp(trgt_mrkvModel->getTransProb(i, j)) * trgt_mrkvModel->getTransProb(i, j) / LOG2;
			double d_tmp_sumTransProb = (exp(src_mrkvModel->getTransProb(i, j)) + exp(trgt_mrkvModel->getTransProb(i, j))) / 2;

			if (0 != d_tmp_sumTransProb) sum_entropyOverCol += d_tmp_sumTransProb * log(d_tmp_sumTransProb) / LOG2;
		}
		src_entropy += exp(src_mrkvModel->getMargProb(i)) * src_entropyOverCol;
		trgt_entropy += exp(trgt_mrkvModel->getMargProb(i)) * trgt_entropyOverCol;
		sum_entropy += (exp(src_mrkvModel->getMargProb(i)) + exp(trgt_mrkvModel->getMargProb(i))) / 2 * sum_entropyOverCol;
	}
}


double DistFactory::getL1dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	L1FreqStrategy* strategy = new L1FreqStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getFreqDist(strategy, 
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getL2dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	L2FreqStrategy* strategy = new L2FreqStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getLInfdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	LInfFreqStrategy* strategy = new LInfFreqStrategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getFFPdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	FFPStrategy* strategy = new FFPStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getFreqDist(strategy, 
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getCoPhylogdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	double dist = IterFactory::getInstance()->getCoPhylogDist(i_arg_k, 
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);
	return dist;
}

double DistFactory::getPearsondist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	PearsonStrategy* strategy = new PearsonStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getCanberradist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	CanberraStrategy* strategy = new CanberraStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getHammingdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	HammingStrategy* strategy = new HammingStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getMatchingdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	MatchingStrategy* strategy = new MatchingStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getJaccarddist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	JaccardStrategy* strategy = new JaccardStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getTanimotodist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	TanimotoStrategy* strategy = new TanimotoStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getDicedist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	DiceStrategy* strategy = new DiceStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getAntidicedist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	AntidiceStrategy* strategy = new AntidiceStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getSneathdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	SneathStrategy* strategy = new SneathStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getHammandist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	HammanStrategy* strategy = new HammanStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getPhidist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	PhiStrategy* strategy = new PhiStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getAnderbergdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	AnderbergStrategy* strategy = new AnderbergStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getGowerdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	GowerStrategy* strategy = new GowerStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getRusseldist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	RusselStrategy* strategy = new RusselStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getYuledist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	YuleStrategy* strategy = new YuleStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getOchiaidist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	OchiaiStrategy* strategy = new OchiaiStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getKulczynskidist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	KulczynskiStrategy* strategy = new KulczynskiStrategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getD2dist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	D2Strategy* strategy = new D2Strategy(i_arg_k, b_arg_singleStrain);
	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getChiSqdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	ChiSqStrategy* strategy = new ChiSqStrategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getCntDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec,
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getD2stardist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
	int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	D2starStrategy* strategy = new D2starStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_src_order, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_trgt_order, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getCntExpDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, src_markovModel, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, trgt_markovModel, arg_trgtKmerModel->totalKmer());

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

double DistFactory::getD2sheppdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
	int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	D2sheppStrategy* strategy = new D2sheppStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_src_order, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_trgt_order, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getCntExpDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, src_markovModel, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, trgt_markovModel, arg_trgtKmerModel->totalKmer());

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

double DistFactory::getHaodist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
	std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	HaoStrategy* strategy = new HaoStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_arg_k - 2, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_arg_k - 2, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getCntExpDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, src_markovModel, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, trgt_markovModel, arg_trgtKmerModel->totalKmer());

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

double DistFactory::getJensenShannondist(int i_arg_k, bool b_arg_singleStrain,
	int i_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	JensenShannonStrategy* strategy = new JensenShannonStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_order, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_order, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getMrkvDist(strategy, src_markovModel, trgt_markovModel);

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

