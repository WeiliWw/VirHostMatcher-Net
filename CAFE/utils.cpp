#include "utils.h"

#include <limits>
#include <cmath>

#include <sys/stat.h>
#include <sys/types.h>

int nt2int(char char_arg_nt)
{
	int i_returnVal = -1;
	if (char_arg_nt == 'A' || char_arg_nt == 'a') i_returnVal = 0;	
	else if (char_arg_nt == 'C' || char_arg_nt == 'c') i_returnVal = 1;	
	else if (char_arg_nt == 'G' || char_arg_nt == 'g') i_returnVal = 2;	
	else if (char_arg_nt == 'T' || char_arg_nt == 't') i_returnVal = 3;
	return i_returnVal;
}

char nt2ComplementNt(char char_arg_nt)
{
	int c_returnChar = 'N';
	if (char_arg_nt == 'A' || char_arg_nt == 'a') c_returnChar = 'T';	
	else if (char_arg_nt == 'C' || char_arg_nt == 'c') c_returnChar = 'G';	
	else if (char_arg_nt == 'G' || char_arg_nt == 'g') c_returnChar = 'C';
	else if (char_arg_nt == 'T' || char_arg_nt == 't') c_returnChar = 'A';
	return c_returnChar;
}

unsigned long long nt2index(const std::string & str_arg_seq)
{
	//return updateIdx((unsigned long long)0, str_arg_seq);
	unsigned long long i_tmp_newIdx = 0;

	for (unsigned int i = 0; i < str_arg_seq.size(); ++i)
	{
		i_tmp_newIdx = (i_tmp_newIdx << 2) + nt2int(str_arg_seq.at(i));

		if (i_tmp_newIdx < 0)
			throw std::runtime_error(" left shifting cause overflow! ");
	}
	return i_tmp_newIdx;
}

std::string index2nt(unsigned long long i_arg_index, int i_arg_seqLength)
{
	if (i_arg_index < 0) return "";
		
	std::string str_tmp_dict = "ACGT";
	unsigned long long i_tmp_currIndex = i_arg_index;
	std::stringstream tmp_revStrStream;
	for (int idx = 0; idx < i_arg_seqLength; ++idx)
	{
		tmp_revStrStream << str_tmp_dict.at(i_tmp_currIndex % BASE);
		i_tmp_currIndex = (i_tmp_currIndex >> 2);
	}

	std::string str_revStr = tmp_revStrStream.str();
	std::stringstream tmp_strStream;
	for (std::string::reverse_iterator rit = str_revStr.rbegin(); rit != str_revStr.rend(); ++rit) tmp_strStream << *rit;
	return tmp_strStream.str();
}

unsigned long long index2revCompleIdx(unsigned long long i_arg_index, int i_arg_seqLength)
{
	if (i_arg_index < 0) return 0;

	std::string str_tmp_dict = "ACGT";
	unsigned long long i_tmp_currIndex = i_arg_index, revCompleIdx = 0;

	for (int idx = 0; idx < i_arg_seqLength; ++idx)
	{
		revCompleIdx = (revCompleIdx << 2) + (BASE - 1 - (i_tmp_currIndex % BASE));
		i_tmp_currIndex = (i_tmp_currIndex >> 2);
	}
	return revCompleIdx;
	//return nt2index(revComplementStr(index2nt(i_arg_index, i_arg_seqLength)));
}

std::string revComplementStr(const std::string & currStr)
{
	std::stringstream ss;
	for (int i = currStr.length() - 1; i >= 0; --i) ss << nt2ComplementNt(currStr[i]);
	return ss.str();
}

std::string complementStr(const std::string & currStr)
{
	std::stringstream ss;
	for (int i = 0; i < currStr.length(); ++i) ss << nt2ComplementNt(currStr[i]);
	return ss.str();
}

std::string revStr(const std::string & currStr)
{
	std::stringstream ss;
	for (int i = currStr.length() - 1; i >= 0; --i) ss << currStr[i];
	return ss.str();
}

std::string trim(std::string currStr)
{
	if (currStr.empty()) return currStr;

	currStr.erase(0, currStr.find_first_not_of(" "));
	currStr.erase(currStr.find_last_not_of(" ") + 1);
	return currStr;
}

void split(const std::string& currStr, std::string delim, std::vector<std::string > & ret)
{
	size_t last = 0;
	size_t index = currStr.find_first_of(delim, last);

	while (index != std::string::npos)
	{
		ret.push_back(trim(currStr.substr(last, index - last)));
		last = index + 1;
		index = currStr.find_first_of(delim, last);
	}
	if (index - last>0)  ret.push_back(trim(currStr.substr(last, index - last)));
}

std::string toLowerCase(std::string currStr)
{
	std::string str_tmp_out;
	std::transform(currStr.begin(), currStr.end(), std::back_inserter(str_tmp_out), ::tolower);
	return str_tmp_out;
}

bool endsWith(std::string const & currStr, std::string const & ending)
{
	if (ending.size() > currStr.size()) return false;
	return std::equal(ending.rbegin(), ending.rend(), currStr.rbegin());
}

bool file_exists(const std::string & str_arg_filename)
{
	std::ifstream f(str_arg_filename.c_str());
	if (f.good()) { f.close(); return true; }
	else { f.close(); return false; }
}

bool dir_exists(const std::string & str_arg_directory)
{
	struct stat buffer;
	return (stat(str_arg_directory.c_str(), &buffer) == 0);
}

std::string getFileName(const std::string & str_arg_path)
{
	std::string base = str_arg_path.substr(str_arg_path.find_last_of("/\\") + 1);
	const size_t period_idx = base.rfind('.');
	if (std::string::npos != period_idx) base.erase(period_idx);
	return base;
}


double log_sum(double log_a, double log_b)
{
	double v;
	if (log_a < log_b) v = log_b + log(1 + exp(log_a - log_b));
	else v = log_a + log(1 + exp(log_b - log_a));
	return v;
}


double log_normalize(double * array, int nlen)
{
	const double log_max = 100.0; // the log(maximum in double precision), make sure it is large enough.
	int argmax;
	double max_val = max(array, nlen, &argmax); //get the maximum value in the array to avoid overflow
	double log_shift = log_max - log(nlen + 1.0) - max_val;
	double sum = 0.0;
	for (int i = 0; i < nlen; i++) sum += exp(array[i] + log_shift); //shift it
		
	double log_norm = log(sum) - log_shift;
	for (int i = 0; i < nlen; i++) array[i] -= log_norm; //shift it back
	return log_norm;
}


double log_normalize(std::vector<double> & vec, int nlen)
{
	const double log_max = 100.0; // the log(maximum in double precision), make sure it is large enough.
	int argmax;
	double max_val = max_vec(vec, nlen, &argmax); //get the maximum value in the array to avoid overflow
	double log_shift = log_max - log(nlen + 1.0) - max_val;
	double sum = 0.0;
	for (int i = 0; i < nlen; i++) sum += exp(vec[i] + log_shift); //shift it
		
	double log_norm = log(sum) - log_shift;
	for (int i = 0; i < nlen; i++) vec[i] -= log_norm; //shift it back
	return log_norm;
}


double log_subtract(double log_a, double log_b)
{
	if (log_a < log_b) return -1000.0;
	double v;
	v = log_a + log(1 - exp(log_b - log_a));
	return v;
}

bool almostEquals(double a, double b)
{
	return std::fabs(a - b) <= std::numeric_limits<double>::epsilon();
}