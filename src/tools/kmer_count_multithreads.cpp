#include "kmer_count_multithreads.h"

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

