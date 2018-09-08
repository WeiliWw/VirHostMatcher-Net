#include <iostream>
#include <unordered_map>
#include "ctpl_stl.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <thread>
#include <unistd.h>
#include <atomic>


std::atomic<int> X;
bool VALID = true;
std::unordered_map<char, int> nuc2num = {
    {'W',-1}, {'w',-1}, {'S',-1}, {'s',-1}, {'M',-1},
    {'m',-1}, {'K',-1}, {'k',-1}, {'R',-1}, {'r',-1},
    {'Y',-1}, {'y',-1}, {'B',-1}, {'b',-1}, {'D',-1},
    {'d',-1}, {'H',-1}, {'h',-1}, {'V',-1},
    {'v',-1}, {'N',-1}, {'n',-1},
    {'N', -1}, {'A', 0}, {'C', 1 }, {'G', 2}, {'T', 3},
    {'n', -1}, {'a', 0}, {'c', 1}, {'g', 2}, {'t', 3}
};

int revcomp(int num, int K){
    int nuc_rc = 0;
    int shift;
    for (int i=0; i<K; i++){
        shift = 2 * (K-i-1);
        nuc_rc += (3 - (num>>shift)&3) * pow(4, i);
    }
    return nuc_rc;
}


void count_one_read(int id, int K, std::string one_read, std::vector<std::atomic<int>> &count_array, bool Reverse) {
    int length = one_read.length();
    int mask = pow(2, (2*(K-1)))-1;
    int num = 0;
    int nuc_num = 0;
    std::unordered_map<char, int>::iterator search;
    int j = 0;
    char nuc;
    int rev = 0;
    for (int i=0;i<length;i++){ 
        nuc = one_read[i];
        search = nuc2num.find(nuc);
        if (search == nuc2num.end()){
            ::VALID = false;
            return;
        }
        else {
            nuc_num = search -> second;
        }
        if (nuc_num == -1){
            num = 0;
	    rev = 0;
            j = 0;
	    }
        else{
            if (j < (K-1)){
                num = num * 4 + nuc_num;
				if (Reverse) rev = (rev >> 2) + (3-nuc_num) * pow(4, K-1);
                j += 1;
			}
            else{
                num = (num&mask)<<2;
                num += nuc_num;
				count_array[num] ++;
				if (Reverse) {
					rev = (rev >> 2) + (3-nuc_num) * pow(4, K-1);
					count_array[rev] ++;
				}
			}
        }
    }
}

std::vector<std::atomic<int>> count(std::string filename, int K, int Num_Threads, bool Reverse) {
    ::VALID = true;
    const int SIZE = pow(4, K);
    std::vector<std::atomic<int>> count_array(SIZE);
    std::ifstream fs(filename);
    char temp_char[5000];
    std::string temp;
    std::string one_read;
    //int Num_Threads =  std::thread::hardware_concurrency();
    ctpl::thread_pool p(Num_Threads);
    while (::VALID) {
		fs.getline(temp_char, 5000, '\n');
		std::string temp(temp_char);
		if(temp[0] == '>'){
			one_read = one_read + "N";
		}
		else{
			one_read = one_read + temp;
		}
		if (one_read.length() >= 5000){
			p.push([one_read, K, Reverse, &count_array](int id){count_one_read(id, K, one_read, count_array, Reverse);});
			one_read = one_read.substr(one_read.length()-K+1);
        }
		if (fs.eof()) break;
		fs.clear();
	}
	p.push([one_read, K, Reverse, &count_array](int id){count_one_read(id, K, one_read, count_array, Reverse);});
	fs.close();
	p.stop(true);
	if (!::VALID)  count_array[0] = -1;
    return count_array;
}
