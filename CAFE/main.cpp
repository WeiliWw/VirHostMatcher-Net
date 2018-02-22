/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#include "utils.h"
#include "seq_model.h"
#include "kmer.h"
#include "dist_model.h"
#include "output.h"


#if defined(__unix__) || defined(__APPLE__) || defined(__MACH__)
#include <dirent.h>
int getdir(std::string dir, std::vector<std::string> &files)
{
	std::string str_currDir = dir;
	if (!str_currDir.empty() && !endsWith(str_currDir, "/")) str_currDir.append("/");

	DIR *dp;
	struct dirent *dirp;
	if ((dp = opendir(dir.c_str())) == NULL) {
		std::cout << "Error(" << errno << ") opening " << dir << std::endl;
		return errno;
	}

	while ((dirp = readdir(dp)) != NULL)
	{
		if (endsWith(dirp->d_name, ".fasta") || endsWith(dirp->d_name, ".fa") || endsWith(dirp->d_name, ".fna"))
		{
			std::string fastaURL = str_currDir;
			fastaURL.append(std::string(dirp->d_name));
			files.push_back(fastaURL);
		}
	}
	closedir(dp);
	return 0;
}
#elif defined(_WIN32) || defined(WIN32)
#include <windows.h>
std::string wchar_t2string(const wchar_t *wchar)
{
	std::string str = "";
	int index = 0;
	while (wchar[index] != 0)
	{
		str += (char)wchar[index];
		++index;
	}
	return str;
}

wchar_t *string2wchar_t(const std::string &str)
{
	wchar_t wchar[260];
	int index = 0;
	while (index < str.size())
	{
		wchar[index] = (wchar_t)str[index];
		++index;
	}
	wchar[index] = 0;
	return wchar;
}
int getdir(std::string dir, std::vector<std::string> &files)
{
	std::string str_currDir = dir;
	if (!str_currDir.empty() && !endsWith(str_currDir, "/")) str_currDir.append("/");

	std::string searchDir = str_currDir; searchDir.append("*.*");
	WIN32_FIND_DATA FindFileData;
	wchar_t * FileName = string2wchar_t(searchDir);
	HANDLE hFind = FindFirstFile(FileName, &FindFileData);

	while (FindNextFile(hFind, &FindFileData))
	{
		std::string file = wchar_t2string(FindFileData.cFileName);
		if (endsWith(file, ".fasta") || endsWith(file, ".fa") || endsWith(file, ".fna"))
		{
			std::string fastaURL = str_currDir;
			fastaURL.append(std::string(file));
			files.push_back(fastaURL);
		}
	}
	return 0;
}
#endif

void print_usage_and_exit()
{
	printf("CAFE:\t aCcelerated Alignment-FrEe sequence analysis\n");
	printf("Description:\t The program provides 29 alignment-free sequence distance measures.\n");
	printf("Authors:\t Yang Lu and Prof. Fengzhu Sun, Computational and Molecular Biology, University of Southern California.\n");
	printf("\nusage:\n");
	printf("./cafe [options]* -D <dist> -I <fa_files> -K <intK>\n");
	printf("\nMain arguments\n");
	printf("\t-D <dist>\tComma-separated list of distance measurements, E.g. -D D2star,Ma,CVtree. The options include: \n");
	printf("\t Conventional measures based on k-mer counts : \n");

	printf("\t\t Ch: Chebyshev distance \n");
	printf("\t\t Canberra: Canberra distance \n");
	printf("\t\t Chisq: Chi-Square distance \n");
	printf("\t\t Cosine: Cosine distance \n");
	printf("\t\t Co-phylog: Co-phylog distance with the seed C_{(k-1)/2,(k-1)/2}O_{1} when k is odd or C_{k/2-1,k/2}O_{1} when k is even \n");
	printf("\t\t D2: D2 distance \n");
	printf("\t\t Eu: Euclidean distance \n");
	printf("\t\t FFP: Feature frequency profiles (FFP) \n");
	printf("\t\t JS: Jensen-Shannon divergence \n");
	printf("\t\t Ma: Manhattan distance \n");
	printf("\t\t Pearson: Pearson distance \n");

	printf("\t Newly developed measures based on background adjusted k-mer counts: \n");
	printf("\t\t CVtree: CVtree distance \n");
	printf("\t\t D2shepp: D2shepp distance \n");
	printf("\t\t D2star: D2star distance \n");

	printf("\t Measures based on presence/absence of k-mers: \n");
	printf("\t\t Anderberg: Anderberg distance \n");
	printf("\t\t Antidice: anti-Dice distance \n");
	printf("\t\t Dice: Dice distance \n");
	printf("\t\t Gower: Gower distance \n");
	printf("\t\t Hamman: Hamman distance \n");
	printf("\t\t Hamming: Hamming distance \n");
	printf("\t\t Jaccard: Jaccard distance \n");
	printf("\t\t Kulczynski: Kulczynski distance \n");
	printf("\t\t Matching: Matching distance \n");
	printf("\t\t Ochiai: Ochiai distance \n");
	printf("\t\t Phi: Pearson Phi distance \n");
	printf("\t\t Russel: Russel-Rao distance \n");
	printf("\t\t Sneath: Sneath-Sokal distance \n");
	printf("\t\t Tanimoto: Rogers-Tanimoto distance \n");
	printf("\t\t Yule: Yule distance \n");
	
	printf("\t-F1 <fa_Dir>\tFolder for virus containing only fasta files with extension '.fasta', '.fa', and '.fna'. \n");
	printf("\t-F2 <fa_Dir>\tFolder for bacteria containing only fasta files with extension '.fasta', '.fa', and '.fna'. \n");
	//printf("\t-I <fa_files>\tComma-separated list of sequence fasta files, e.g. -I speciesA.fa,speciesB.fa,speciesC.fa. Pairwise similarity is calculated based upon the sequences specified with this option. \n");
	printf("\t-K <intK>\tKmer Length\n");
	printf("\nOptions\n");
	printf("\t-J <jfexe_path>\tUse jellyfish to accelerate kmer counting. <jfexe_path> denotes the file path of jellyfish executable file, e.g. jellyfish-2.2.4/bin/./jellyfish \n");
	printf("\t-L <lower>\tOnly consider k-mer with occurrence >= <lower>. The default value is 0. \n");
	printf("\t-M <order>\tMarkov Order involved in D2star and D2shepp. There are two possible options. The first option is one single value indicating that all the sequences use the same order. The second option is comma-separated list of orders. Notice that the length of the list should match the number of fasta files. The order value could be non-negative integer but less than Kmer length or \"-1\" with the special intention to automatically infer the suitable order (not suitable for JS). The default Markov Order is -1 (Automaticcaly determine by BIC).\n");
	printf("\t-R \t\tConsider Reverse Complement in kmer counting. \n");
	printf("\t-S <dir>\tSave/Load calculated k-mer count binary files to the folder <dir>. Each input fasta file corresponds to particular model. \n");
	printf("\t-O <path>\tOutput results to file at <path> \n");
	printf("\t-T <type>\tThe output type as the input to downstream analysis, including: plain, phylip (as hierarchical clustering), cytoscape (as network analysis) and mds (Multidimensional Scaling as 2D plotting). E.g. -T mds. The default type is plain. \n");
	//printf("\t-V <dir>\tSave visualization result to the folder <dir>. \n");
	printf("\nExamples:\n");
	printf("\t./cafe -M 0 -O v1_h1.txt -S hash/ -T plain -F1 virus1/ -F2 host1/ -K 10 -D D2star,Ma\n");
	//printf("\t./cafe -M 0 -S model_dir -I speciesA.fa,speciesB.fa -J jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star,Ma\n");
	//printf("\t./cafe -M 0 -L 2 -I speciesA.fa,speciesB.fa -J jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star,Ma -R\n");
	printf("\n");

	exit(0);
}




int cafe(int argc, char* argv[])
{
	clock_t startTime, endTime;
	startTime = clock();

	std::vector<std::string> vec_distStr;
	std::vector<std::string> vec_orderStr;
	std::vector<std::string> vec_fastaFiles, vec_fastaFiles1, vec_fastaFiles2;
	int fasta1Idx = 0, fasta2Idx = 0;

	std::vector<int> vec_order;
	std::vector<dist> vec_dist;
	std::vector<std::string> vec_namelist, vec_namelist1, vec_namelist2;
	std::vector<std::string> vec_saveURLlist;

	int i_k = 0, i_lowerCnt = 0;
	OUTPUT_TYPE outputType = PHYLIP;
    // set jellyfishValid to false
	bool singleStrain = true, containChiSq = false, containCvtree = false, jellyfishValid = false;
	std::string str_save_modelDir = "", str_save_vizDir = "", str_outputFileURL = "", str_jellyfishExeURL = "";
	
	printf("Start parsing the arguments... \n");
	if (argc < 2 || !strcmp(argv[1], "-help") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--usage")) print_usage_and_exit();
	for (int i = 1; i < argc; ++i)
	{
		//if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "-i")) split(std::string(argv[++i]), ",", vec_fastaFiles);
		if (!strcmp(argv[i], "-F1") || !strcmp(argv[i], "-f1"))
		{
			std::string fastaDir = std::string(argv[++i]);
			std::cout << fastaDir << std::endl;
			getdir(fastaDir, vec_fastaFiles1);

			fasta1Idx = vec_fastaFiles.size();
			for (int idx = 0; idx < vec_fastaFiles1.size(); ++idx)
			{
				//std::cout << "Append: " << vec_fastaFiles1[idx] << std::endl;
				vec_fastaFiles.push_back(vec_fastaFiles1[idx]);
			}
		}
		else if (!strcmp(argv[i], "-F2") || !strcmp(argv[i], "-f2"))
		{
			std::string fastaDir = std::string(argv[++i]);
			std::cout << fastaDir << std::endl;
			getdir(fastaDir, vec_fastaFiles2);

			fasta2Idx = vec_fastaFiles.size();
			for (int idx = 0; idx < vec_fastaFiles2.size(); ++idx)
			{
				//std::cout << "Append: " << vec_fastaFiles2[idx] << std::endl;
				vec_fastaFiles.push_back(vec_fastaFiles2[idx]);
			}
		}
		else if (!strcmp(argv[i], "-K") || !strcmp(argv[i], "-k")) i_k = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-M") || !strcmp(argv[i], "-m")) split(std::string(argv[++i]), ",", vec_orderStr); 
		else if (!strcmp(argv[i], "-L") || !strcmp(argv[i], "-l")) i_lowerCnt = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "-s")) str_save_modelDir = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-O") || !strcmp(argv[i], "-o")) str_outputFileURL = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-V") || !strcmp(argv[i], "-v")) str_save_vizDir = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-J") || !strcmp(argv[i], "-j"))
		{
			str_jellyfishExeURL = std::string(argv[++i]);
			if (!file_exists(str_jellyfishExeURL)) 
			{ 
				jellyfishValid = false; 
				std::cout << "[error]: Jellyfish executable file not exist at " << str_jellyfishExeURL << std::endl; 
				//print_usage_and_exit();
			}
		}
		else if (!strcmp(argv[i], "-T") || !strcmp(argv[i], "-t"))
		{
			++i;
			if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "plain")) outputType = PLAIN;
			else if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "phylip")) outputType = PHYLIP;
			else if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "cytoscape")) outputType = CYTOSCAPE;
			else if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "mds")) outputType = MDS;
			else printf("[warning]: The output type is unrecognized! \n");
		}
		else if (!strcmp(argv[i], "-D") || !strcmp(argv[i], "-d"))
		{
			split(std::string(argv[++i]), ",", vec_distStr);
		}
		else if (!strcmp(argv[i], "-R") || !strcmp(argv[i], "-r")) { singleStrain = false; }
	}

	for (int j = 0; j < vec_distStr.size(); ++j)
	{
		if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2")) vec_dist.push_back(D2);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2star")) vec_dist.push_back(D2STAR);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2shepp")) vec_dist.push_back(D2SHEPP);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "ma")) vec_dist.push_back(Ma);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "eu")) vec_dist.push_back(Eu);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "ch")) vec_dist.push_back(Ch);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "ffp")) vec_dist.push_back(FFP);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "co-phylog")) vec_dist.push_back(Co_Phylog);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "cvtree")) { containCvtree = true; vec_dist.push_back(CVtree); }
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "js")) vec_dist.push_back(JS);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "chisq")) { containChiSq = true; vec_dist.push_back(CHISQ); }
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "cosine")) vec_dist.push_back(COSINE);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "pearson")) vec_dist.push_back(PEARSON);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "canberra")) vec_dist.push_back(CANBERRA);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "hamming")) vec_dist.push_back(HAMMING);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "matching")) vec_dist.push_back(MATCHING);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "jaccard")) vec_dist.push_back(JACCARD);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "tanimoto")) vec_dist.push_back(TANIMOTO);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "dice")) vec_dist.push_back(DICE);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "antidice")) vec_dist.push_back(ANTIDICE);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "sneath")) vec_dist.push_back(SNEATH);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "hamman")) vec_dist.push_back(HAMMAN);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "phi")) vec_dist.push_back(PHI);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "anderberg")) vec_dist.push_back(ANDERBERG);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "gower")) vec_dist.push_back(GOWER);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "russel")) vec_dist.push_back(RUSSEL);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "yule")) vec_dist.push_back(YULE);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "ochiai")) vec_dist.push_back(OCHIAI);
		else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "kulczynski")) vec_dist.push_back(KULCZYNSKI);
		else
		{
			printf("[warning]: The distance measurement %s is unrecognized!\n ", vec_distStr[j].c_str());
			vec_distStr.erase(vec_distStr.begin() + j);
		}
	}
	if (!str_save_modelDir.empty() && !dir_exists(str_save_modelDir)) { std::string cmd = "mkdir " + str_save_modelDir; system(cmd.c_str()); }
	
	// validity check for Kmer length
	if (i_k <= 0)
	{
		printf("[error]: Kmer length should be positive! \n");
		print_usage_and_exit();
	}

	if (i_k <= 2 && containCvtree)
	{
		printf("[error]: Kmer length must exceed 2 in CVTree! \n");
		print_usage_and_exit();
	}

	// validity check for possible markov order
	if (vec_orderStr.empty())  vec_orderStr.push_back("0");
	if (vec_orderStr.size() != vec_fastaFiles.size())
	{
		if (vec_orderStr.size() > 1)
		{
			printf("[error]: The length of the order list should match the number of fasta files! \n");
			print_usage_and_exit();
		}
		else if (vec_orderStr.size() == 1)
		{
			while (vec_orderStr.size() < vec_fastaFiles.size()) vec_orderStr.push_back(vec_orderStr.at(0));
		}
	}
	for (int i = 0; i < vec_orderStr.size(); ++i)
	{
		int order = atoi(vec_orderStr[i].c_str());
		if (order >= i_k || order < -1)
		{
			printf("[error]: Markov Order should be non-negative and less than k or '-1' indicating auto inference! \n");
			std::cout << "k = " << i_k << "  order = " << order << std::endl;
			print_usage_and_exit();
		}
	}

	// validity check for input fasta files
	for (int i = 0; i < vec_fastaFiles.size(); ++i)
	{
		if (!file_exists(vec_fastaFiles[i]))
		{
			std::cout << "Input file not exist! Skip: " << vec_fastaFiles[i] << std::endl;
			vec_fastaFiles.erase(vec_fastaFiles.begin() + i);
			vec_orderStr.erase(vec_orderStr.begin() + i);
		}
	}
	printf("Finish parsing the arguments.\n");

	printf("Start pre-processing the hash of input fasta files...\n");
	// validity check and pre-processing for hashs of input fasta files
	for (int i = 0; i < vec_fastaFiles.size(); ++i)
	{
		int order = atoi(vec_orderStr[i].c_str()); 
		std::vector<int> cachedOrderVec;
		if (containChiSq) cachedOrderVec.push_back(i_k + 1);
		cachedOrderVec.push_back(i_k);
		if (containCvtree) { cachedOrderVec.push_back(i_k - 1); cachedOrderVec.push_back(i_k - 2); }

		if (0 == order) cachedOrderVec.push_back(1);
		else if (order > 0) { cachedOrderVec.push_back(order + 1); cachedOrderVec.push_back(order); }
		else { for (int minOrder = i_k-1; minOrder > 0 ; --minOrder) cachedOrderVec.push_back(minOrder); }

		std::string str_seqName = getFileName(vec_fastaFiles[i]);
		std::string str_saveDir = str_save_modelDir;
		if (!str_saveDir.empty() && !endsWith(str_saveDir, "/")) str_saveDir.append("/");
		str_saveDir.append("hash_").append(str_seqName).append("_L_").append(patch::to_string(i_lowerCnt)).append("_k_");

		for (int orderIdx = 0; orderIdx < cachedOrderVec.size(); ++orderIdx)
		{
			int currK = cachedOrderVec[orderIdx];
			std::string str_saveURL = str_saveDir + patch::to_string(currK); if (file_exists(str_saveURL)) continue;
			
			KmerModel* kmerModel = new KmerModel(currK, true);
			if (0 == orderIdx)
			{
				bool jellyfishSucceed = false;
				if (jellyfishValid)
				{
					std::string str_jfBinURL = str_saveURL; str_jfBinURL.append(".jf");
					std::string str_jfTabTxtURL = str_saveURL; str_jfTabTxtURL.append(".cnt");

					std::string lowerCntStr = " -L " + patch::to_string(i_lowerCnt);
					if (i_lowerCnt < 2) lowerCntStr = "";

					std::string cmd1 = str_jellyfishExeURL + " count -m " + patch::to_string(currK) + " -s 500M -t 20" + lowerCntStr + " -o " + str_jfBinURL + " " + vec_fastaFiles[i];
					system(cmd1.c_str());  std::cout << "Execute Command: " << cmd1 << std::endl;

					if (file_exists(str_jfBinURL))
					{
						std::cout << "Jellyfish succeed in count!" << std::endl;
						std::string cmd2 = str_jellyfishExeURL + " dump -t " + str_jfBinURL + lowerCntStr + " > " + str_jfTabTxtURL;
						system(cmd2.c_str());  std::cout << "Execute Command: " << cmd2 << std::endl;

						if (file_exists(str_jfTabTxtURL))
						{
							std::cout << "Jellyfish succeed in dump!" << std::endl;
							kmerModel->saveFromJellyFish(str_jfTabTxtURL, str_saveURL);
							jellyfishSucceed = true;
						}
					}
				}
				
				if (!jellyfishSucceed)
				{
                    // Remove jellyfish
					// std::cout << "Jellyfish not succeed! Now use slow counting!" << std::endl;
					kmerModel->saveFromFasta(currK, vec_fastaFiles[i], str_saveURL);
				}
			}
			else
			{
				int prevK = cachedOrderVec[orderIdx - 1];
				std::string str_prevOrderURL = str_saveDir + patch::to_string(prevK);
				kmerModel->saveFromLargerK(currK, prevK, str_prevOrderURL, str_saveURL);
			}
			delete kmerModel; std::cout << "Now save model to " << str_saveURL << std::endl;
		}
	}

	for (int i = 0; i < vec_fastaFiles.size(); ++i)
	{
		std::string str_seqName = getFileName(vec_fastaFiles[i]);

		std::string str_saveDir = str_save_modelDir;
		if (!str_saveDir.empty() && !endsWith(str_saveDir, "/"))
			str_saveDir.append("/");
		str_saveDir.append("hash_").append(str_seqName).append("_L_").append(patch::to_string(i_lowerCnt)).append("_k_");

		int order = atoi(vec_orderStr[i].c_str());
		if (order < 0) order = getEstMarkovOrder(i_k, str_saveDir, str_seqName); // need auto infer order

		vec_order.push_back(order);
		vec_saveURLlist.push_back(str_saveDir);
		vec_namelist.push_back(str_seqName);
	}

	for (int i = fasta1Idx; i < fasta1Idx + vec_fastaFiles1.size(); ++i) vec_namelist1.push_back(vec_namelist.at(i));
	for (int j = fasta2Idx; j < fasta2Idx + vec_fastaFiles2.size(); ++j) vec_namelist2.push_back(vec_namelist.at(j));
		

	printf("Finish pre-processing the hash of input fasta files.\n");

	endTime = clock();
	// std::cout << "Time Elapsed: " << ((float)endTime - (float)startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
	startTime = clock();

	printf("Start calculating the distance...\n");
	for (int distIdx = 0; distIdx < vec_dist.size(); ++distIdx)
	{
		dist currDist = vec_dist[distIdx];
		int hashK = i_k; if (CHISQ == currDist) hashK = i_k + 1;
		
		smat::Matrix<double>* simMat = new smat::Matrix<double>(vec_fastaFiles1.size(), vec_fastaFiles2.size(), 0);		
		KmerModel* src_kmerModel, *trgt_kmerModel;

		for (int i = fasta1Idx; i < fasta1Idx+vec_fastaFiles1.size(); ++i)
		{
			src_kmerModel = new KmerModel(hashK, singleStrain);
			src_kmerModel->load(hashK, vec_saveURLlist[i] + patch::to_string(hashK));

			for (int j = fasta2Idx; j < fasta2Idx + vec_fastaFiles2.size(); ++j)
			{
				trgt_kmerModel = new KmerModel(hashK, singleStrain);
				trgt_kmerModel->load(hashK, vec_saveURLlist[j] + patch::to_string(hashK));

				double distVal = 0;

				if (D2 == currDist || COSINE == currDist) distVal = DistFactory::getInstance()->getD2dist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (MATCHING == currDist) distVal = DistFactory::getInstance()->getMatchingdist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (JACCARD == currDist) distVal = DistFactory::getInstance()->getJaccarddist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (TANIMOTO == currDist) distVal = DistFactory::getInstance()->getTanimotodist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (DICE == currDist) distVal = DistFactory::getInstance()->getDicedist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (ANTIDICE == currDist) distVal = DistFactory::getInstance()->getAntidicedist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (SNEATH == currDist) distVal = DistFactory::getInstance()->getSneathdist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (HAMMAN == currDist) distVal = DistFactory::getInstance()->getHammandist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (PHI == currDist) distVal = DistFactory::getInstance()->getPhidist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (ANDERBERG == currDist) distVal = DistFactory::getInstance()->getAnderbergdist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (GOWER == currDist) distVal = DistFactory::getInstance()->getGowerdist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (RUSSEL == currDist) distVal = DistFactory::getInstance()->getRusseldist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (YULE == currDist) distVal = DistFactory::getInstance()->getYuledist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (OCHIAI == currDist) distVal = DistFactory::getInstance()->getOchiaidist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (KULCZYNSKI == currDist) distVal = DistFactory::getInstance()->getKulczynskidist(hashK, singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (D2STAR == currDist)
					distVal = DistFactory::getInstance()->getD2stardist(hashK, singleStrain, i_lowerCnt, vec_order.at(i), vec_order.at(j), vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				else if (D2SHEPP == currDist)
					distVal = DistFactory::getInstance()->getD2sheppdist(hashK, singleStrain, i_lowerCnt, vec_order.at(i), vec_order.at(j), vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				else if (CVtree == currDist)
					distVal = DistFactory::getInstance()->getHaodist(hashK, singleStrain, i_lowerCnt, vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				else if (Ma == currDist) distVal = DistFactory::getInstance()->getL1dist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (Eu == currDist) distVal = DistFactory::getInstance()->getL2dist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (Ch == currDist) distVal = DistFactory::getInstance()->getLInfdist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (FFP == currDist) distVal = DistFactory::getInstance()->getFFPdist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (Co_Phylog == currDist) distVal = DistFactory::getInstance()->getCoPhylogdist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (CHISQ == currDist) distVal = DistFactory::getInstance()->getChiSqdist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (PEARSON == currDist) distVal = DistFactory::getInstance()->getPearsondist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (CANBERRA == currDist) distVal = DistFactory::getInstance()->getCanberradist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (HAMMING == currDist) distVal = DistFactory::getInstance()->getHammingdist(hashK, singleStrain, src_kmerModel, trgt_kmerModel);
				else if (JS == currDist)
				{
					//if (vec_order.at(i) != vec_order.at(j)) throw std::runtime_error("JS expect same order! ");
					int maxOrder = vec_order.at(i); if(maxOrder < vec_order.at(j)) maxOrder = vec_order.at(j);
					distVal = DistFactory::getInstance()->getJensenShannondist(hashK, singleStrain, maxOrder, vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				}
				simMat->set(i - fasta1Idx, j - fasta2Idx, distVal);
				delete trgt_kmerModel;
			}

			delete src_kmerModel;
		}

		std::cout << "-------------------------------------------------" << std::endl;
		std::cout << "Dist: \t" << vec_distStr[distIdx] << std::endl;

		if (!str_outputFileURL.empty())
		{
			std::string postfixedOutputURL1 = str_outputFileURL;
			postfixedOutputURL1.append(".").append(vec_distStr[distIdx]).append(".plain");
			OutputWriter::getInstance()->writeToFile(PLAIN, simMat, &vec_namelist1, &vec_namelist2, postfixedOutputURL1);

		}
		//OutputWriter::getInstance()->writeToConsole(outputType, simMat, &vec_namelist);
		
		endTime = clock();
		// std::cout << "Time Elapsed: " << ((float)endTime - (float)startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
		startTime = clock();
	}
	printf("Finish calculating the distance.\n");
	printf("Done");
	//system("pause");
	return 0;
}
