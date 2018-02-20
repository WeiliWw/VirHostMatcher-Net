#include "kmer.h"
#include <cstring>

IterFactory* IterFactory::instance = 0;

IterFactory* IterFactory::getInstance()
{
	if (!instance) instance = new IterFactory();
	return instance;
}

AbsIter* IterFactory::getKmerCntIterator(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec, int i_arg_lowerCnt)
{
	if (i_arg_lowerCnt <= 0) return new KmerCntTraverseIter(i_arg_k, arg_kmerCntUnorderMap);
	else return new KmerCntHashIter(i_arg_k, arg_kmerCntUnorderMap, arg_kmerVec);
}

AbsIter* IterFactory::getKmerFreqIterator(int i_arg_k, std::unordered_map<unsigned long long, unsigned long>* arg_kmerCntUnorderMap, std::vector<unsigned long long>* arg_kmerVec, unsigned long l_arg_totalKmer)
{
	return new KmerFreqHashIter(i_arg_k, arg_kmerCntUnorderMap, arg_kmerVec, l_arg_totalKmer);
}

KmerProbEnsembDelegate* IterFactory::getKmerProbDelegate(int i_arg_k, bool b_arg_singleStrain, MarkovModel* arg_mrkvModel)
{
	KmerProbEnsembDelegate* kmerProbDelegate = new KmerProbEnsembDelegate(i_arg_k, arg_mrkvModel, b_arg_singleStrain);
	kmerProbDelegate->init();
	return kmerProbDelegate;
}

double IterFactory::getFreqDist(AbsTupleDistStrategy* distStrategy,
	std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec, unsigned long src_totalKmer,
	std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec, unsigned long trgt_totalKmer)
{
	AbsIter* src_freqIter = getKmerFreqIterator(distStrategy->i_k, src_kmerCntUnorderMap, src_kmerVec, src_totalKmer);
	AbsIter* trgt_freqIter = getKmerFreqIterator(	distStrategy->i_k, trgt_kmerCntUnorderMap, trgt_kmerVec, trgt_totalKmer);

	std::queue<unsigned long long> src_kmer_queue, trgt_kmer_queue;
	std::queue<double> src_freq_queue, trgt_freq_queue;

	while (src_freqIter->hasNext() || trgt_freqIter->hasNext())
	{
		if (src_freqIter->hasNext())
		{
			double src_X_w = *(*src_freqIter); double src_kmer = src_freqIter->getCurrKmer(); (*src_freqIter)++;
			src_kmer_queue.push(src_kmer); src_freq_queue.push(src_X_w);
		}

		if (trgt_freqIter->hasNext())
		{
			double trgt_X_w = *(*trgt_freqIter); double trgt_kmer = trgt_freqIter->getCurrKmer(); (*trgt_freqIter)++;
			trgt_kmer_queue.push(trgt_kmer); trgt_freq_queue.push(trgt_X_w);
		}

		while (!src_kmer_queue.empty() && !trgt_kmer_queue.empty())
		{
			if (src_kmer_queue.front() < trgt_kmer_queue.front())
			{
				distStrategy->dealWithTuple(src_freq_queue.front(), 0);
				src_kmer_queue.pop(); src_freq_queue.pop();
			}
			else if (src_kmer_queue.front() > trgt_kmer_queue.front())
			{
				distStrategy->dealWithTuple(0, trgt_freq_queue.front());
				trgt_kmer_queue.pop(); trgt_freq_queue.pop();
			}
			else
			{
				distStrategy->dealWithTuple(src_freq_queue.front(), trgt_freq_queue.front());
				src_kmer_queue.pop(); src_freq_queue.pop(); trgt_kmer_queue.pop(); trgt_freq_queue.pop();
			}
		}
	}
	while (!src_kmer_queue.empty())
	{
		distStrategy->dealWithTuple(src_freq_queue.front(), 0);
		src_kmer_queue.pop(); src_freq_queue.pop();
	}
	while (!trgt_kmer_queue.empty())
	{
		distStrategy->dealWithTuple(0, trgt_freq_queue.front());
		trgt_kmer_queue.pop(); trgt_freq_queue.pop();
	}

	delete src_freqIter; delete trgt_freqIter;
	return distStrategy->getDist();
}

double IterFactory::getCoPhylogDist(int i_arg_k, 
		std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec,
		std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec)
{
	int maskLen = i_arg_k; maskLen = (maskLen >> 1); maskLen = (maskLen << 1);
	unsigned long long mask_low = ((unsigned long long)1 << maskLen) - 1;

	std::unordered_map<unsigned long long, int>* src_context_dup = new std::unordered_map<unsigned long long, int>();
	std::unordered_map<unsigned long long, int>* src_context_obj = new std::unordered_map<unsigned long long, int>();
	//std::unordered_map<unsigned long long, double>* src_context_cnt = new std::unordered_map<unsigned long long, double>();
	std::unordered_map<unsigned long long, int>* trgt_context_dup = new std::unordered_map<unsigned long long, int>();
	std::unordered_map<unsigned long long, int>* trgt_context_obj = new std::unordered_map<unsigned long long, int>();
	//std::unordered_map<unsigned long long, double>* trgt_context_cnt = new std::unordered_map<unsigned long long, double>();

	AbsIter* src_CntIter = getKmerCntIterator(i_arg_k, src_kmerCntUnorderMap, src_kmerVec, 0);
	AbsIter* trgt_CntIter = getKmerCntIterator(i_arg_k, trgt_kmerCntUnorderMap, trgt_kmerVec, 0);

	while (src_CntIter->hasNext())
	{
		double src_X_w = *(*src_CntIter); unsigned long long src_kmer = src_CntIter->getCurrKmer(); (*src_CntIter)++; if (src_X_w <= 0) continue;
		
		unsigned long long kmer_low  = (src_kmer & mask_low);
		unsigned long long kmer_high = (src_kmer >> maskLen);
		int objVal = (kmer_high % BASE); kmer_high = (kmer_high >> 2);
		unsigned long long kmer_new = ((kmer_high << maskLen) + kmer_low);

		(*src_context_dup)[kmer_new] ++;
		(*src_context_obj)[kmer_new] = objVal;
		//(*src_context_cnt)[kmer_new] = src_X_w;
	}

	while (trgt_CntIter->hasNext())
	{
		double trgt_X_w = *(*trgt_CntIter); unsigned long long trgt_kmer = trgt_CntIter->getCurrKmer(); (*trgt_CntIter)++; if (trgt_X_w <= 0) continue;
		
		unsigned long long kmer_low  = (trgt_kmer & mask_low);
		unsigned long long kmer_high = (trgt_kmer >> maskLen);
		int objVal = (kmer_high % BASE); kmer_high = (kmer_high >> 2);
		unsigned long long kmer_new = ((kmer_high << maskLen) + kmer_low);

		(*trgt_context_dup)[kmer_new] ++;
		(*trgt_context_obj)[kmer_new] = objVal;
		//(*trgt_context_cnt)[kmer_new] = trgt_X_w;
	}

	double sum=0, hit=0;
	for (std::unordered_map<unsigned long long, int>::iterator iter = src_context_dup->begin(); iter != src_context_dup->end(); iter++)
	{
		unsigned long long currKmerIdx = iter->first;
		int objDup_src = (*src_context_dup)[currKmerIdx]; int objDup_trgt = (*trgt_context_dup)[currKmerIdx];
		
		if(1 == objDup_src && 1 == objDup_trgt) 
		{
			sum ++;
			int obj_src = (*src_context_obj)[currKmerIdx]; int obj_trgt = (*trgt_context_obj)[currKmerIdx];
			
			if(obj_src != obj_trgt) hit ++;
		}
	}

	delete src_CntIter; delete trgt_CntIter;
	delete src_context_dup; delete src_context_obj; //delete src_context_cnt;
	delete trgt_context_dup; delete trgt_context_obj; //delete trgt_context_cnt;
	
	if(0 == sum) return 0;
	else return hit/sum;
}

double IterFactory::getCntDist(AbsTupleDistStrategy* distStrategy, int i_arg_lowerCnt,
	std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec,
	std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec)
{
	AbsIter* src_CntIter = getKmerCntIterator(distStrategy->i_k, src_kmerCntUnorderMap, src_kmerVec, i_arg_lowerCnt);
	AbsIter* trgt_CntIter = getKmerCntIterator(distStrategy->i_k, trgt_kmerCntUnorderMap, trgt_kmerVec, i_arg_lowerCnt);

	std::queue<unsigned long long> src_kmer_queue, trgt_kmer_queue;
	std::queue<double> src_cnt_queue, trgt_cnt_queue;

	while (src_CntIter->hasNext() || trgt_CntIter->hasNext())
	{
		if (src_CntIter->hasNext())
		{
			double src_X_w = *(*src_CntIter); double src_kmer = src_CntIter->getCurrKmer(); (*src_CntIter)++;
			if (src_X_w >= i_arg_lowerCnt) { src_kmer_queue.push(src_kmer); src_cnt_queue.push(src_X_w); }
		}

		if (trgt_CntIter->hasNext())
		{
			double trgt_X_w = *(*trgt_CntIter); double trgt_kmer = trgt_CntIter->getCurrKmer(); (*trgt_CntIter)++;
			if (trgt_X_w >= i_arg_lowerCnt) { trgt_kmer_queue.push(trgt_kmer); trgt_cnt_queue.push(trgt_X_w); }
		}

		while (!src_kmer_queue.empty() && !trgt_kmer_queue.empty())
		{
			if (src_kmer_queue.front() < trgt_kmer_queue.front())
			{
				src_kmer_queue.pop(); src_cnt_queue.pop();
			}
			else if (src_kmer_queue.front() > trgt_kmer_queue.front())
			{
				trgt_kmer_queue.pop(); trgt_cnt_queue.pop();
			}
			else
			{
				distStrategy->dealWithTuple(src_cnt_queue.front(), trgt_cnt_queue.front());
				src_kmer_queue.pop(); src_cnt_queue.pop(); trgt_kmer_queue.pop(); trgt_cnt_queue.pop();
			}
		}
	}
	delete src_CntIter; delete trgt_CntIter;
	return distStrategy->getDist();
}

double IterFactory::getCntDist(AbsQuadStrategy* distStrategy,
	std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec,
	std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec)
{
	AbsIter* src_kmerCntIter = getKmerCntIterator(distStrategy->i_k, src_kmerCntUnorderMap, src_kmerVec, 0);
	AbsIter* trgt_kmerCntIter = getKmerCntIterator(distStrategy->i_k, trgt_kmerCntUnorderMap, trgt_kmerVec, 0);

	double src_X_w_1, src_X_w_2, src_X_w_3, src_X_w_4, src_X_w;
	double trgt_X_w_1, trgt_X_w_2, trgt_X_w_3, trgt_X_w_4, trgt_X_w;
	double all_X_w_1, all_X_w_2, all_X_w_3, all_X_w_4, all_X_w;

	while (true)
	{
		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_1 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_1 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_2 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_2 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_3 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_3 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_4 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_4 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
		
		src_X_w = src_X_w_1 + src_X_w_2 + src_X_w_3 + src_X_w_4;
		trgt_X_w = trgt_X_w_1 + trgt_X_w_2 + trgt_X_w_3 + trgt_X_w_4;
		all_X_w_1 = src_X_w_1 + trgt_X_w_1;
		all_X_w_2 = src_X_w_2 + trgt_X_w_2;
		all_X_w_3 = src_X_w_3 + trgt_X_w_3;
		all_X_w_4 = src_X_w_4 + trgt_X_w_4;
		all_X_w = src_X_w + trgt_X_w;

		if (0 == all_X_w) continue;

		if (all_X_w_1 > 0)
		{
			if (src_X_w > 0){ distStrategy->dealWithQuad(src_X_w, src_X_w_1, all_X_w, all_X_w_1); }
			if (trgt_X_w > 0){ distStrategy->dealWithQuad(trgt_X_w, trgt_X_w_1, all_X_w, all_X_w_1); }
		}
		if (all_X_w_2 > 0)
		{
			if (src_X_w > 0){ distStrategy->dealWithQuad(src_X_w, src_X_w_2, all_X_w, all_X_w_2); }
			if (trgt_X_w > 0){ distStrategy->dealWithQuad(trgt_X_w, trgt_X_w_2, all_X_w, all_X_w_2); }
		}
		if (all_X_w_3 > 0)
		{
			if (src_X_w > 0){ distStrategy->dealWithQuad(src_X_w, src_X_w_3, all_X_w, all_X_w_3); }
			if (trgt_X_w > 0){ distStrategy->dealWithQuad(trgt_X_w, trgt_X_w_3, all_X_w, all_X_w_3); }
		}
		if (all_X_w_4 > 0)
		{
			if (src_X_w > 0){ distStrategy->dealWithQuad(src_X_w, src_X_w_4, all_X_w, all_X_w_4); }
			if (trgt_X_w > 0){ distStrategy->dealWithQuad(trgt_X_w, trgt_X_w_4, all_X_w, all_X_w_4); }
		}
	}

	delete src_kmerCntIter; delete trgt_kmerCntIter;
	return distStrategy->getDist();
}

double IterFactory::getCntExpDist(AbsQuadStrategy* distStrategy, int i_arg_lowerCnt,
	std::unordered_map<unsigned long long, unsigned long>* src_kmerCntUnorderMap, std::vector<unsigned long long>* src_kmerVec, MarkovModel* src_mrkvModel, unsigned long src_totalKmer,
	std::unordered_map<unsigned long long, unsigned long>* trgt_kmerCntUnorderMap, std::vector<unsigned long long>* trgt_kmerVec, MarkovModel* trgt_mrkvModel, unsigned long trgt_totalKmer)
{
	AbsIter* src_CntIter = getKmerCntIterator(distStrategy->i_k, src_kmerCntUnorderMap, src_kmerVec, i_arg_lowerCnt);
	AbsIter* trgt_CntIter = getKmerCntIterator(distStrategy->i_k, trgt_kmerCntUnorderMap, trgt_kmerVec, i_arg_lowerCnt);

	KmerProbEnsembDelegate* src_kmerProbDelegate = getKmerProbDelegate(distStrategy->i_k, distStrategy->b_singleStrain, src_mrkvModel);
	KmerProbEnsembDelegate* trgt_kmerProbDelegate = getKmerProbDelegate(distStrategy->i_k, distStrategy->b_singleStrain, trgt_mrkvModel);

	double log_src_totalKmerLen = log(src_totalKmer), log_trgt_totalKmerLen = log(trgt_totalKmer);

	std::queue<unsigned long long> src_kmer_queue, trgt_kmer_queue;
	std::queue<double> src_cnt_queue, trgt_cnt_queue;

	while (src_CntIter->hasNext() || trgt_CntIter->hasNext())
	{
		if (src_CntIter->hasNext())
		{
			double src_X_w = *(*src_CntIter); double src_kmer = src_CntIter->getCurrKmer(); (*src_CntIter)++;
			if (src_X_w >= i_arg_lowerCnt) { src_kmer_queue.push(src_kmer); src_cnt_queue.push(src_X_w); }
		}

		if (trgt_CntIter->hasNext())
		{
			double trgt_X_w = *(*trgt_CntIter); double trgt_kmer = trgt_CntIter->getCurrKmer(); (*trgt_CntIter)++;
			if (trgt_X_w >= i_arg_lowerCnt) { trgt_kmer_queue.push(trgt_kmer); trgt_cnt_queue.push(trgt_X_w); }
		}

		while (!src_kmer_queue.empty() && !trgt_kmer_queue.empty())
		{
			if (src_kmer_queue.front() < trgt_kmer_queue.front())
			{
				src_kmer_queue.pop(); src_cnt_queue.pop();
			}
			else if (src_kmer_queue.front() > trgt_kmer_queue.front())
			{
				trgt_kmer_queue.pop(); trgt_cnt_queue.pop();
			}
			else
			{
				double src_X_w = src_cnt_queue.front(); double src_prob_X_w = src_kmerProbDelegate->getKmerlogProb(src_kmer_queue.front());
				double trgt_X_w = trgt_cnt_queue.front(); double trgt_prob_X_w = trgt_kmerProbDelegate->getKmerlogProb(trgt_kmer_queue.front());

				src_kmer_queue.pop(); src_cnt_queue.pop(); trgt_kmer_queue.pop(); trgt_cnt_queue.pop();
				
				double src_EX_w = ((0 == src_prob_X_w) ? 0 : exp(log_src_totalKmerLen + src_prob_X_w));
				double trgt_EX_w = ((0 == trgt_prob_X_w) ? 0 : exp(log_trgt_totalKmerLen + trgt_prob_X_w));

				distStrategy->dealWithQuad(src_X_w, src_EX_w, trgt_X_w, trgt_EX_w);
			}
		}
	}

	delete src_CntIter; delete trgt_CntIter; delete src_kmerProbDelegate; delete trgt_kmerProbDelegate;
	return distStrategy->getDist();
}

double IterFactory::getMrkvDist(AbsMrkvStrategy* distStrategy, MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel)
{
	distStrategy->dealWithMrkv(src_mrkvModel, trgt_mrkvModel);
	return distStrategy->getDist();
}

KmerProbDelegate::KmerProbDelegate(int i_arg_k, MarkovModel* arg_mrkvModel, bool b_arg_isRevCompl)
{
	i_k = i_arg_k;
	maxAllowedIdx = (unsigned long long)pow(BASE, i_arg_k) - 1;
	markovModel = arg_mrkvModel;
	isRevCompl = b_arg_isRevCompl;

	//markovModel->print();
}

void KmerProbDelegate::init()
{
	nextPosKmerIdx = 0;
	isEnd = false;

	while (!kmer_traceStack.empty()) kmer_traceStack.pop();
	while (!orderIdx_traceStack.empty()) orderIdx_traceStack.pop();
	while (!logProbProd_traceStack.empty()) logProbProd_traceStack.pop();

	while (!lowerIdx_traceStack.empty()) lowerIdx_traceStack.pop();
	while (!upperIdx_traceStack.empty()) upperIdx_traceStack.pop();

	//while (kmer_traceStack.size() < i_k) push(0);
	std::stack<unsigned long long> route_traceStack;
	for (int i = 0; i < i_k; ++i) route_traceStack.push(0);
	push(route_traceStack);
}

bool KmerProbDelegate::push(unsigned long long idx)
{
	int i_order = markovModel->getOrder();

	unsigned long long l_tmp_currOrderIdx = 0;
	if (!orderIdx_traceStack.empty()) l_tmp_currOrderIdx = orderIdx_traceStack.top();

	if (!isRevCompl) //regular case
	{
		if (kmer_traceStack.size() == i_order)
		{
			if (0 == markovModel->getMargProb(l_tmp_currOrderIdx) && 0 != i_order) return false;
		}

		if (kmer_traceStack.size() >= i_order)
		{
			if (0 == markovModel->getTransProb(l_tmp_currOrderIdx, idx)) return false;
		}

		unsigned long long l_tmp_newOrderIdx = 0; // when i_order = 0
		if (1 == i_order)
			l_tmp_newOrderIdx = idx;
		else if (i_order > 1) l_tmp_newOrderIdx = ((~((unsigned long long)3 << ((i_order - 1) * 2)) & l_tmp_currOrderIdx) << 2) + idx;
		orderIdx_traceStack.push(l_tmp_newOrderIdx);

		double d_tmp_newLogProbProd = 0;
		if (!logProbProd_traceStack.empty()) d_tmp_newLogProbProd = logProbProd_traceStack.top();

		//if (kmer_traceStack.size() == i_order) d_tmp_newLogProbProd += log(markovModel->getMargProb(l_tmp_currOrderIdx));
		//if (kmer_traceStack.size() >= i_order) d_tmp_newLogProbProd += log(markovModel->getTransProb(l_tmp_currOrderIdx, idx));
		if (kmer_traceStack.size() == i_order) d_tmp_newLogProbProd += markovModel->getMargProb(l_tmp_currOrderIdx);
		if (kmer_traceStack.size() >= i_order) d_tmp_newLogProbProd += markovModel->getTransProb(l_tmp_currOrderIdx, idx);

		logProbProd_traceStack.push(d_tmp_newLogProbProd);
	}
	else //reverse complement case
	{
		unsigned long long compl_idx = BASE - 1 - idx;
		unsigned long long compl_idx_toRemove = l_tmp_currOrderIdx % BASE;

		unsigned long long l_tmp_newOrderIdx = 0; // when i_order = 0
		if (0 == i_order)
			compl_idx_toRemove = compl_idx;
		else if (1 == i_order)
			l_tmp_newOrderIdx = compl_idx;
		else if (i_order > 1)
		{
			if (kmer_traceStack.size() < i_order)
				l_tmp_newOrderIdx = (compl_idx << (kmer_traceStack.size() * 2)) + l_tmp_currOrderIdx;
			else
				l_tmp_newOrderIdx = (compl_idx << ((i_order - 1) * 2)) + (l_tmp_currOrderIdx >> 2);
		}

		if (kmer_traceStack.size() >= i_order)
		{
			if (0 == markovModel->getTransProb(l_tmp_newOrderIdx, compl_idx_toRemove)) return false;
		}

		if (kmer_traceStack.size() == i_k - 1)
		{
			if (0 == markovModel->getMargProb(l_tmp_newOrderIdx) && 0 != i_order) return false;
		}

		orderIdx_traceStack.push(l_tmp_newOrderIdx);
		double d_tmp_newLogProbProd = 0;
		if (!logProbProd_traceStack.empty()) d_tmp_newLogProbProd = logProbProd_traceStack.top();

		//if (kmer_traceStack.size() >= i_order) d_tmp_newLogProbProd += log(markovModel->getTransProb(l_tmp_newOrderIdx, compl_idx_toRemove));
		//if (kmer_traceStack.size() == i_k - 1) d_tmp_newLogProbProd += log(markovModel->getMargProb(l_tmp_newOrderIdx));
		if (kmer_traceStack.size() >= i_order) d_tmp_newLogProbProd += markovModel->getTransProb(l_tmp_newOrderIdx, compl_idx_toRemove);
		if (kmer_traceStack.size() == i_k - 1) d_tmp_newLogProbProd += markovModel->getMargProb(l_tmp_newOrderIdx);
		logProbProd_traceStack.push(d_tmp_newLogProbProd);
	}


	kmer_traceStack.push(idx);
	nextPosKmerIdx = (nextPosKmerIdx << 2) + idx;

	unsigned long long idxLowBound = 0;
	unsigned long long gap = (unsigned long long)pow(BASE, i_k - lowerIdx_traceStack.size() - 1);
	if (!lowerIdx_traceStack.empty()) idxLowBound = lowerIdx_traceStack.top();
	lowerIdx_traceStack.push(idxLowBound + idx*gap);
	upperIdx_traceStack.push(idxLowBound + (idx + 1)*gap - 1);

	return true;
}

unsigned long long KmerProbDelegate::pop()
{
	if (kmer_traceStack.empty()) throw std::runtime_error(" cannot pop from empty stack! ");

	unsigned long long i_tmp_curr = kmer_traceStack.top();
	nextPosKmerIdx = (nextPosKmerIdx >> 2);

	kmer_traceStack.pop();
	orderIdx_traceStack.pop();
	logProbProd_traceStack.pop();
	lowerIdx_traceStack.pop();
	upperIdx_traceStack.pop();

	return i_tmp_curr;
}

bool KmerProbDelegate::push(std::stack<unsigned long long> & indices)
{
	while (kmer_traceStack.size() < i_k)
	{
		unsigned long long i_tmp_nextIdx = 0;
		if (!indices.empty()) { i_tmp_nextIdx = indices.top(); indices.pop(); }

		bool b_tmp_succeed = push(i_tmp_nextIdx);
		while (!b_tmp_succeed)
		{
			while (!indices.empty()) indices.pop();
			if (i_tmp_nextIdx >= BASE - 1) break;
			i_tmp_nextIdx++;
			b_tmp_succeed = push(i_tmp_nextIdx);
		}
		if (!b_tmp_succeed)
		{
			bool b_tmp_succeed1 = increment();
			if (!b_tmp_succeed1) return false;
		}
	}
	return true;
}

bool KmerProbDelegate::increment()
{
	while (!kmer_traceStack.empty() && kmer_traceStack.top() >= BASE - 1) pop();
	if (kmer_traceStack.empty()) return false;

	unsigned long long i_tmp_currIdx = pop();
	unsigned long long i_tmp_nextIdx = i_tmp_currIdx + 1;

	bool b_tmp_succeed = push(i_tmp_nextIdx);
	while (!b_tmp_succeed)
	{
		if (i_tmp_nextIdx >= BASE - 1) break;

		i_tmp_nextIdx++;
		b_tmp_succeed = push(i_tmp_nextIdx);
	}

	if (!b_tmp_succeed) return increment();
	return true;
}

double KmerProbDelegate::getKmerlogProb(unsigned long long queryNextKmerIdx)
{
	if (isEnd) return 0;

	if (nextPosKmerIdx < queryNextKmerIdx)
	{
		unsigned long long currUpperIdx = upperIdx_traceStack.top(), currLowerIdx = lowerIdx_traceStack.top();
		while (currUpperIdx < queryNextKmerIdx)
		{
			pop();
			if (!upperIdx_traceStack.empty()) { currUpperIdx = upperIdx_traceStack.top(); currLowerIdx = lowerIdx_traceStack.top(); }
			else { currUpperIdx = maxAllowedIdx; currLowerIdx = 0; }
		}

		unsigned long long offset = queryNextKmerIdx - currLowerIdx;
		std::stack<unsigned long long> route_traceStack;
		for (int i = 0; i < i_k - lowerIdx_traceStack.size(); ++i)
		{
			unsigned long long i_tmp_route = offset % BASE;
			route_traceStack.push(i_tmp_route);
			offset = (offset >> 2);
		}

		bool b_tmp_succeed = push(route_traceStack);
		if (!b_tmp_succeed) { isEnd = true; return 0; }
	}
	//std::cout << "nextPosKmerIdx\t" << nextPosKmerIdx << std::endl;

	if (nextPosKmerIdx > queryNextKmerIdx) return 0;
	else if (queryNextKmerIdx == nextPosKmerIdx) { return logProbProd_traceStack.empty()?0:logProbProd_traceStack.top(); }
	else throw std::runtime_error("nextPosKmerIdx < queryNextKmerIdx cannot be true! ");
}

KmerProbEnsembDelegate::KmerProbEnsembDelegate(int i_arg_k, MarkovModel* arg_mrkvModel, bool b_arg_singleStrain)
{
	i_k = i_arg_k;
	markovModel = arg_mrkvModel;
	b_singleStrain = b_arg_singleStrain;

	kmerProbDelegate = new KmerProbDelegate(i_k, markovModel, false);
	if (!b_singleStrain) revComplKmerProbDelegate = new KmerProbDelegate(i_k, markovModel, true);
}

KmerProbEnsembDelegate::~KmerProbEnsembDelegate()
{
	delete kmerProbDelegate;
	if (!b_singleStrain) delete revComplKmerProbDelegate;
	//delete markovModel;
}

void KmerProbEnsembDelegate::init()
{
	kmerProbDelegate->init();
	if (!b_singleStrain) revComplKmerProbDelegate->init();
}

double KmerProbEnsembDelegate::getKmerlogProb(unsigned long long queryNextKmerIdx)
{
	double d_log_kmerProb = kmerProbDelegate->getKmerlogProb(queryNextKmerIdx);
	if (!b_singleStrain)
	{
		double d_log_revComplKmerProb = revComplKmerProbDelegate->getKmerlogProb(queryNextKmerIdx);

		if (0 == d_log_kmerProb && 0 == d_log_revComplKmerProb) return 0;
		else if (0 == d_log_kmerProb || 0 == d_log_revComplKmerProb) return d_log_revComplKmerProb + d_log_kmerProb - LOG2;
		else return log_sum(d_log_kmerProb, d_log_revComplKmerProb) - LOG2;
	}
	return d_log_kmerProb;
}

bool KmerModel::load(int i_arg_k, std::string str_arg_inputURL)
{
	kmerCntUnorderMap->clear(); kmerVec->clear();

	unsigned long long i_tmp_rowDim = 0;
	std::ifstream tmp_finPipe(str_arg_inputURL.c_str(), std::ios::in | std::ios::binary);
	tmp_finPipe.read((char*)&i_tmp_rowDim, sizeof(unsigned long long));

	unsigned long long* vec_kmerIdx = new unsigned long long[i_tmp_rowDim]; memset(vec_kmerIdx, 0, sizeof(unsigned long long) * i_tmp_rowDim);
	unsigned long* vec_kmerCnt = new unsigned long[i_tmp_rowDim]; memset(vec_kmerCnt, 0, sizeof(unsigned long) * i_tmp_rowDim);

	tmp_finPipe.read((char*)vec_kmerIdx, sizeof(unsigned long long) * i_tmp_rowDim);
	tmp_finPipe.read((char*)vec_kmerCnt, sizeof(unsigned long) * i_tmp_rowDim);
	tmp_finPipe.close();

	for (unsigned long long rowIdx = 0; rowIdx < i_tmp_rowDim; ++rowIdx)
	{
		unsigned long long currKmerIdx = vec_kmerIdx[rowIdx];
		unsigned long currKmerCnt = vec_kmerCnt[rowIdx];

		(*kmerCntUnorderMap)[currKmerIdx] += currKmerCnt;
		if (!b_singleStrain) (*kmerCntUnorderMap)[index2revCompleIdx(currKmerIdx, i_arg_k)] += currKmerCnt;

		kmerVec->push_back(currKmerIdx);

		//std::cout << rowIdx << " = " << currKmerIdx << " = " << currKmerCnt << std::endl;
	}

	//std::cout << i_tmp_rowDim << " ; "<<kmerCntUnorderMap->size() << " ; " << kmerVec->size() << std::endl;

	delete[] vec_kmerIdx; delete[] vec_kmerCnt;
	return true;
}

bool KmerModel::save(std::map<unsigned long long, unsigned long> *kmerCntMap, std::string str_arg_outputURL)
{
	unsigned long long i_tmp_rowDim = kmerCntMap->size();
	unsigned long long* vec_kmerIdx = new unsigned long long[i_tmp_rowDim]; memset(vec_kmerIdx, 0, sizeof(unsigned long long) * i_tmp_rowDim);
	unsigned long* vec_kmerCnt = new unsigned long[i_tmp_rowDim]; memset(vec_kmerCnt, 0, sizeof(unsigned long) * i_tmp_rowDim);

	unsigned long long currIdx = 0;
	for (std::map<unsigned long long, unsigned long>::iterator iter = kmerCntMap->begin(); iter != kmerCntMap->end(); iter++)
	{
		vec_kmerIdx[currIdx] = iter->first;
		vec_kmerCnt[currIdx] = iter->second;
		currIdx++;
	}

	std::ofstream tmp_ofsPipe(str_arg_outputURL.c_str(), std::ios::out | std::ios::binary);
	tmp_ofsPipe.write((char*)&i_tmp_rowDim, sizeof(unsigned long long));
	tmp_ofsPipe.write((char*)vec_kmerIdx, sizeof(unsigned long long)*i_tmp_rowDim);
	tmp_ofsPipe.write((char*)vec_kmerCnt, sizeof(unsigned long)*i_tmp_rowDim);
	tmp_ofsPipe.flush();
	tmp_ofsPipe.close();

	delete[] vec_kmerIdx; delete[] vec_kmerCnt;
	return true;
}

bool KmerModel::saveFromLargerK(int i_arg_k, int i_arg_larger_k, std::string str_arg_inputURL, std::string str_arg_outputURL)
{
	if (i_arg_larger_k < i_arg_k)
	{
		std::cout << "[error]: i_arg_larger_k < i_arg_k in KmerModel::loadFromLargerK. " << std::endl;
		return false;
	}
	if (i_arg_larger_k == i_arg_k) return load(i_arg_k, str_arg_inputURL);

	std::map<unsigned long long, unsigned long> *kmerCntMap = new std::map<unsigned long long, unsigned long>();

	unsigned long long i_tmp_rowDim = 0;
	std::ifstream tmp_finPipe(str_arg_inputURL.c_str(), std::ios::in | std::ios::binary);
	tmp_finPipe.read((char*)&i_tmp_rowDim, sizeof(unsigned long long));

	unsigned long long* vec_kmerIdx = new unsigned long long[i_tmp_rowDim]; memset(vec_kmerIdx, 0, sizeof(unsigned long long) * i_tmp_rowDim);
	unsigned long* vec_kmerCnt = new unsigned long[i_tmp_rowDim]; memset(vec_kmerCnt, 0, sizeof(unsigned long) * i_tmp_rowDim);

	tmp_finPipe.read((char*)vec_kmerIdx, sizeof(unsigned long long) * i_tmp_rowDim);
	tmp_finPipe.read((char*)vec_kmerCnt, sizeof(unsigned long) * i_tmp_rowDim);
	tmp_finPipe.close();

	for (unsigned long long rowIdx = 0; rowIdx < i_tmp_rowDim; ++rowIdx)
	{
		unsigned long long currKmerIdx = vec_kmerIdx[rowIdx];
		unsigned long currKmerCnt = vec_kmerCnt[rowIdx];

		unsigned long long newKmerIdx = (currKmerIdx >> (2 * (i_arg_larger_k - i_arg_k)));
		(*kmerCntMap)[newKmerIdx] += currKmerCnt;
		//if (!b_singleStrain) (*kmerCntMap)[index2revCompleIdx(newKmerIdx, i_arg_k)] += currKmerCnt;
	}

	delete[] vec_kmerIdx; delete[] vec_kmerCnt;
	save(kmerCntMap, str_arg_outputURL);
	delete kmerCntMap;

	return true;
}

bool KmerModel::saveFromJellyFish(std::string str_arg_jfTxtURL, std::string str_arg_outputURL)
{
	std::map<unsigned long long, unsigned long> *kmerCntMap = new std::map<unsigned long long, unsigned long>();

	std::ifstream fastaStream(str_arg_jfTxtURL.c_str());
	if (!fastaStream.is_open()) return false;

	std::string str_line1 = "", currKmer = "";
	while (getline(fastaStream, str_line1))
	{
		if (str_line1.empty() || "" == str_line1) continue;
		getline(fastaStream, currKmer);

		str_line1.erase(str_line1.begin());
		long currKmerCnt = atol(str_line1.c_str());

		unsigned long long currKmerIdx = nt2index(currKmer);
		(*kmerCntMap)[currKmerIdx] += currKmerCnt;
		//if (!b_singleStrain) (*kmerCntTable)[index2revCompleIdx(currKmerIdx, currKmer.length())] += currKmerCnt;
	}
	this->i_k = currKmer.length();
	fastaStream.close();

	save(kmerCntMap, str_arg_outputURL);
	delete kmerCntMap;
	return true;
}

bool KmerModel::saveFromFasta(int i_arg_k, std::string str_arg_fastaFileURL, std::string str_arg_outputURL)
{
	std::map<unsigned long long, unsigned long> *kmerCntMap = new std::map<unsigned long long, unsigned long>();

	std::ifstream fastaStream(str_arg_fastaFileURL.c_str());
	if (fastaStream.is_open())
	{
		std::string str_tmp_line = ""; 
		int totalCharCnt = 0;
		unsigned long long currKmerIdx = 0;
			
		while (getline(fastaStream, str_tmp_line))
		{
			if (str_tmp_line.empty() || "" == str_tmp_line) continue;
			if (str_tmp_line.substr(0, 1) == ">")
			{
				totalCharCnt = 0; currKmerIdx = 0;
				continue;
			}
	
			for (int i = 0; i < str_tmp_line.length(); ++i)
			{
				int idx = nt2int(str_tmp_line.at(i)); if (idx < 0) continue;
				if (totalCharCnt < i_arg_k) totalCharCnt++;
				currKmerIdx = ((~((unsigned long)3 << ((i_arg_k - 1) * 2)) & currKmerIdx) << 2) + idx;
					
				if (totalCharCnt < i_arg_k) continue;
				(*kmerCntMap)[currKmerIdx] ++;
			}
		}
		fastaStream.close();
	}
	this->i_k = i_arg_k;

	save(kmerCntMap, str_arg_outputURL);
	delete kmerCntMap;
	return true;
}

unsigned long KmerModel::totalKmer()
{
	unsigned long totalKmerCnt = 0;
	for (std::unordered_map<unsigned long long, unsigned long>::iterator iter = kmerCntUnorderMap->begin(); iter != kmerCntUnorderMap->end(); iter++) totalKmerCnt += (*kmerCntUnorderMap)[iter->first];
	return totalKmerCnt;
}

MarkovModel* KmerModel::getMarkovModel(int i_arg_order, std::string str_arg_saveURLPrefix)
{
	MarkovModel* tmp_mrkvModel = new MarkovModel(i_arg_order);

	std::string str_orderURL = str_arg_saveURLPrefix + patch::to_string(i_arg_order);
	std::string str_orderPlusOneURL = str_arg_saveURLPrefix + patch::to_string(i_arg_order + 1);

	if ((!file_exists(str_orderURL) && i_arg_order > 0) || !file_exists(str_orderPlusOneURL))
	{
		std::cout << "[error]: file not exists! in KmerModel::getMarkovModel !!!" << std::endl;
		return NULL;
	}

	unsigned long long i_orderDim = 0;
	unsigned long long *vec_orderkmerIdx, *vec_orderkmerCnt;
	if (i_arg_order <= 0)
	{
		vec_orderkmerIdx = new unsigned long long[1]; vec_orderkmerIdx[0] = 0;
		vec_orderkmerCnt = new unsigned long long[1]; vec_orderkmerCnt[0] = 1;
		i_orderDim = 1;
	}
	else
	{
		std::ifstream tmp_finPipe(str_orderURL.c_str(), std::ios::in | std::ios::binary);
		tmp_finPipe.read((char*)&i_orderDim, sizeof(unsigned long long));

		vec_orderkmerIdx = new unsigned long long[i_orderDim]; memset(vec_orderkmerIdx, 0, sizeof(unsigned long long) * i_orderDim);
		vec_orderkmerCnt = new unsigned long long[i_orderDim]; memset(vec_orderkmerCnt, 0, sizeof(unsigned long long) * i_orderDim);

		tmp_finPipe.read((char*)vec_orderkmerIdx, sizeof(unsigned long long) * i_orderDim);
		tmp_finPipe.read((char*)vec_orderkmerCnt, sizeof(unsigned long long) * i_orderDim);
		tmp_finPipe.close();
	}
	tmp_mrkvModel->constructMarg(vec_orderkmerIdx, vec_orderkmerCnt, i_orderDim);
	delete[] vec_orderkmerIdx; delete[] vec_orderkmerCnt;


	unsigned long long i_orderPlusOneDim = 0;
	std::ifstream tmp_finPipePlusOne(str_orderPlusOneURL.c_str(), std::ios::in | std::ios::binary);
	tmp_finPipePlusOne.read((char*)&i_orderPlusOneDim, sizeof(unsigned long long));

	unsigned long long *vec_orderPlusOnekmerIdx = new unsigned long long[i_orderPlusOneDim]; memset(vec_orderPlusOnekmerIdx, 0, sizeof(unsigned long long) * i_orderPlusOneDim);
	unsigned long long *vec_orderPlusOnekmerCnt = new unsigned long long[i_orderPlusOneDim]; memset(vec_orderPlusOnekmerCnt, 0, sizeof(unsigned long long) * i_orderPlusOneDim);

	tmp_finPipePlusOne.read((char*)vec_orderPlusOnekmerIdx, sizeof(unsigned long long) * i_orderPlusOneDim);
	tmp_finPipePlusOne.read((char*)vec_orderPlusOnekmerCnt, sizeof(unsigned long long) * i_orderPlusOneDim);
	tmp_finPipePlusOne.close();

	for (unsigned long long rowIdx = 0; rowIdx < i_orderPlusOneDim; ++rowIdx)
	{
		unsigned long long currKmerIdx = vec_orderPlusOnekmerIdx[rowIdx];
		unsigned long long i_tmp_route = currKmerIdx % BASE;
		unsigned long long newKmerIdx = (currKmerIdx >> 2);
		tmp_mrkvModel->addTransProb(newKmerIdx, i_tmp_route, vec_orderPlusOnekmerCnt[rowIdx]);
	}
	delete[] vec_orderPlusOnekmerIdx; delete[] vec_orderPlusOnekmerCnt;

	tmp_mrkvModel->normalize();
	return tmp_mrkvModel;
}


int getEstMarkovOrder(int i_arg_k, std::string str_arg_saveURLPrefix, std::string str_arg_seqName)
{
	std::cout << "Now estimating markov order for " << str_arg_seqName << " ..." << std::endl;

	KmerModel* kmerModel = new KmerModel(i_arg_k, true);
	kmerModel->load(i_arg_k, str_arg_saveURLPrefix + patch::to_string(i_arg_k));
	unsigned long l_totalKmer = kmerModel->totalKmer();
	delete kmerModel;
	kmerModel = new KmerModel(i_arg_k, true);

	int maxOrder = i_arg_k; if(maxOrder > MAX_ORDER) maxOrder = MAX_ORDER;
	std::vector<double> result_traceVec;
	for (int currOrder = 0; currOrder < maxOrder; ++currOrder)
	{
		double BIC = (BASE - 1) * pow(BASE, currOrder) * log(l_totalKmer + i_arg_k - currOrder);

		MarkovModel* markovModel = kmerModel->getMarkovModel(currOrder, str_arg_saveURLPrefix);

		KmerModel* tmpKmerModel = new KmerModel(currOrder + 1, true);
		tmpKmerModel->load(currOrder + 1, str_arg_saveURLPrefix + patch::to_string(currOrder + 1));

		AbsIter* tmpkmerCntIter = IterFactory::getInstance()->getKmerCntIterator(i_arg_k, tmpKmerModel->kmerCntUnorderMap, tmpKmerModel->kmerVec, 1);

		double logLH = 0;
		while (tmpkmerCntIter->hasNext())
		{
			double src_X_w = *(*tmpkmerCntIter); (*tmpkmerCntIter)++; if (0 == src_X_w) continue;
			unsigned long long kmerIdx = tmpkmerCntIter->getCurrKmer();
			
			int route = kmerIdx % BASE;
			unsigned long long idx_new = (kmerIdx >> 2);

			//if (0 != markovModel->getTransProb(idx_new, route)) logLH += src_X_w*log(markovModel->getTransProb(idx_new, route));
			if (0 != markovModel->getTransProb(idx_new, route)) logLH += src_X_w*markovModel->getTransProb(idx_new, route);
		}
		BIC = -2 * logLH + BIC;
		result_traceVec.push_back(BIC);

		delete tmpkmerCntIter; delete tmpKmerModel; delete markovModel;
	}
	delete kmerModel;


	if (result_traceVec.size() <= 1) return 0;

	int min_idx = 0;
	min_vec<double>(result_traceVec, result_traceVec.size(), &min_idx);
	for (int i = 0; i < result_traceVec.size(); ++i) std::cout << "order = " << i << "  BIC = " << result_traceVec.at(i) << std::endl;
	std::cout << "The selected order = " << min_idx << std::endl;
	return min_idx;
}