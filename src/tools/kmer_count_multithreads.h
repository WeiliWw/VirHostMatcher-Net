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
std::unordered_map<char, int> nuc2num = {
        {'A', 0}, {'C', 1 }, {'G', 2}, {'T', 3}, 
	{'a', 0}, {'c', 1}, {'g', 2}, {'t', 3}
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
    int j = 0;
    char nuc;
    int rev = 0;
    for (int i=0;i<length;i++){ 
        nuc = one_read[i];
        if (nuc == 'N'){
            num = 0;
            j = 0;
	}
        else{
            if (j < (K-1)){
                num = num * 4 + nuc2num[nuc];
                j += 1;
	    }
            else{
                num = (num&mask)<<2;
                num += nuc2num[nuc];
		count_array[num] ++;
		if (Reverse) {
		    rev = revcomp(num, K);
		    count_array[rev] ++;
		}
	   }
        }
    }
}

std::vector<std::atomic<int>> count(std::string filename, int K, bool Reverse, int Num_Threads) {
    //std::cout<<"Filename: "<<filename<<std::endl;
    const int SIZE = pow(4, K);
    std::vector<std::atomic<int>> count_array(SIZE);
    std::ifstream fs(filename);
    std::string temp;
    std::string one_read;
    ctpl::thread_pool p(Num_Threads);
    //std::cout << Num_Threads <<std::endl;
    while (getline(fs, temp)) {
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
    }
    p.push([one_read, K, Reverse, &count_array](int id){count_one_read(id, K, one_read, count_array, Reverse);});
    fs.close();
    return count_array;
    
}

