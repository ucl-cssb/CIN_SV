#ifndef UTIL_HPP
#define UTIL_HPP


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <set>
#include <climits>
#include <cstdlib>

#include <unistd.h>

#include <boost/filesystem.hpp>

// For defining empirical_cumulative_distribution_function
// #include <iterator>
// #include <stdexcept>


using namespace std;


typedef map<pair<int, int>, int> pcn;   // copy number at a position
typedef map<pair<int, int>, double> dpcn;   // copy number at a position


const int NUM_CHR = 22;
const int NORM_PLOIDY = 2;
const int NUM_LOC = 5000;

// Number of chromsomes affected in a multipolar event
const int MULTI_NCHR = 16;

vector<int> CHR_LENGTHS;
vector<int> ARM_BOUNDS;
vector<int> CENT_STARTS;
vector<int> CENT_ENDS;
vector<int> TELO_ENDS1;
vector<int> TELO_ENDS2;


// equal probability of selecting each chromosome
vector<double> CHR_PROBS_vec(NUM_CHR, 1.0 / NUM_CHR);
// use pointer to be compatible with gsl
double* CHR_PROBS = &CHR_PROBS_vec[0];

const int NUM_SVTYPE = 6;

enum SStat_type{ALL};
enum Growth_type{ONLY_BIRTH, CHANGE_BIRTH, CHANGE_DEATH, CHANGE_BOTH};
enum Telo_type{NONTEL, PTEL, QTEL, COMPLETE};
enum Adj_type{INTERVAL, REF, VAR};   // 0: interval, 1: reference, 2: variant
enum SV_type{NONE, DEL, DUP, H2HINV, T2TINV, TRA};  // only for variant adjacency, TRA: intra-chromosomal
enum Junc_type{HEAD, TAIL};

gsl_rng * r;

unsigned long setup_rng(unsigned long set_seed){
  gsl_rng_env_setup();

  const gsl_rng_type* T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if(set_seed != 0){
    gsl_rng_set(r, set_seed);
    return(set_seed);
  }else{
    int t = time(NULL);
    int pid = getpid();
    long s = t*pid;
    //cout << "pid:" << "\t" << getpid() << endl;
    // cout << "seed:" << "\t" << t << "\t" << pid << "\t" << abs(s) << endl;
    std::cout << "seed:" << "\t" << abs(s) << std::endl;
    gsl_rng_set (r, abs(s));
    return(abs(s));
  }
}


// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
long unsigned myrng(long unsigned n){
  // This function returns a random integer from 0 to n-1 inclusive by scaling down and/or discarding samples from the generator r. All integers in the range [0,n-1] are produced with equal probability.
  return gsl_rng_uniform_int(r, n);
}


long unsigned (*fp)(long unsigned) = myrng;


// factorial for choose(n,k)
int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}


// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform(r);

  while(myrandom==0){
    myrandom = a + (b-a)*gsl_rng_uniform(r);
  }
  return myrandom;
}



// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const std::vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i=0; i<rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  std::vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i=0; i<rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  std::vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i=0; i<rates.size(); ++i){
    if( u < psum[i] ){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}


void set_outdir(string outdir, int verbose = 0){
    const char* path = outdir.c_str();
    boost::filesystem::path dir(path);
    if(boost::filesystem::create_directory(dir)){
        if(verbose > 0) cerr << "Directory Created: " << outdir <<endl;
    }
}


// Read strings separated by space into a vector.
// num: the number of element in the string. If not given, assume it is the first number in the string
template <typename T>
void get_vals_from_str(vector<T>& vals, string str_vals, int num = 0){
    assert(str_vals != "");
    stringstream ss(str_vals);
    if(num == 0)   ss >> num;
    for(int i = 0; i < num; i++){
        T d1;
        ss >> d1;
        vals.push_back(d1);
    }
    // cout << "vector from " << str_vals << ": ";
    // for(auto v : vals){
    //     cout << "\t" << v;
    // }
    // cout << endl;
}


// Read reference genome informaton (chr start end centromere_pos); there must be an empty line in the end
// chr1	0	2300000	p36.33	gneg
void read_genome_info(const string& filename, vector<int>& chr_lengths, vector<int>& arm_boundaries, vector<int>& centromere_starts, vector<int>& centromere_ends, vector<int>& telomere_ends1, vector<int>& telomere_ends2, int verbose = 0){
  verbose = 0;
  if(verbose) cout << "\nReading reference genome information from " << filename << endl;

  ifstream infile(filename.c_str());

  if(!infile.is_open()){
    std::cerr << "Error: open of input data unsuccessful: " << filename << std::endl;
    exit(1);
  }

  std::string line;

  getline(infile, line); // skip the first header line

  while(!getline(infile, line).eof()){
    if(verbose) cout << line << endl;
    if(line.empty()){
      continue;
    }

    std::vector<std::string> split;
    std::string buf;
    stringstream ss(line);
    while (ss >> buf) split.push_back(buf);
    assert(split.size() == 7);

    int len = atoi(split[1].c_str());
    chr_lengths.push_back(len);

    int pos = atoi(split[2].c_str());
    arm_boundaries.push_back(pos);

    pos = atoi(split[3].c_str());
    centromere_starts.push_back(pos);
    pos = atoi(split[4].c_str());
    centromere_ends.push_back(pos);

    pos = atoi(split[5].c_str());
    telomere_ends1.push_back(pos);
    pos = atoi(split[6].c_str());
    telomere_ends2.push_back(pos);
  }

  if(verbose){
    // cout << chr_lengths.size() << endl;
    int nchr = chr_lengths.size();
    for(int i = 0; i < nchr; i++){
      cout << i + 1 << "\t" << chr_lengths[i] << "\t" << arm_boundaries[i] << endl;
    }
  }
}

// read intervals, to facilitate overlap computation
void read_fragile_sites(const string& filename, vector<int>& fragile_sites){

}

// get the probability matrix for selecting chromosomes
// random, biased, fixed
// TODO: assign probability based on Dirichlet distribution? selected_chr specified from input?
void get_chr_prob(int type, vector<int>& selected_chr){
  if(type == 0){   // equal probability for each chromosome, or based on chromosome size?
    // for(int i = 0; i < NUM_CHR; i++){
    //   CHR_PROBS[i] = 1.0 / NUM_CHR;
    // }
    return;
  }

  int nb = 1;
  if(type == 1){ // biased
    nb = myrng(3);
    // randomly select nb chromosomes
    for(int i = 0; i < nb; i++){
      int csel = myrng(NUM_CHR);
      selected_chr.push_back(csel);
    }
  }else{ // fixed to one specific chromosome, chr1 by default,  if(type == 2)
    // # same probability bias as before for reproducibility
    // # choose g target chromosomes, -1 for index, advised to choose up to 3 targets max
    selected_chr.push_back(0);
  }

  double p0 = 2.0 / 3.0;
  double p1 = p0 / nb;      // p1 for selected chromosome
  double p2 = (1 - p0) / (NUM_CHR - nb);
  for(int i = 0; i < NUM_CHR; i++){
    CHR_PROBS[i] = p2;
  }
  for(int i = 0; i < selected_chr.size(); i++){
    CHR_PROBS[i] = p1;
  }
}


string get_side_string(int side){
  string side_chr = "+";  // HEAD
  if(side == TAIL){
    side_chr = "-";
  }
  return side_chr;
}


string get_adj_type_string(int type){
  string type_str = "";
  switch (type) {
    case 1: type_str = "REF"; break;
    case 2: type_str = "VAR"; break;
    default: type_str = "INTERVAL";
  }
  return type_str;
}


// 0: no telomere; 1: left telomere; 2: right telomere; 3: both telomeres
string get_telomere_type_string(int type){
  string type_str = "";
  switch (type) {
    case 1: type_str = "LEFT telomere"; break;
    case 2: type_str = "RIGHT telomere"; break;
    case 3: type_str = "BOTH telomeres"; break;
    default: type_str = "NO telomere";
  }
  return type_str;
}


string get_haplotype_string(int hap){
  if(hap == 0){
    return "A";
  }else{
    return "B";
  }
}


string get_sv_type_string(int type){
  string type_str = "";
  switch (type) {
    case DUP: type_str = "DUP"; break;
    case DEL: type_str = "DEL"; break;
    case TRA: type_str = "TRA"; break;
    case H2HINV: type_str = "H2HINV"; break;
    case T2TINV: type_str = "T2TINV"; break;
    default: type_str = "NO SV";
  }
  return type_str;
}



// Run length encoding a vector, find the size of the longest region with value x
int find_max_size(int x, const vector<int>& vec){
  int len = vec.size();
  vector<int> rle_vec_val;
  vector<int> rle_vec_len;

  for(int i = 0; i < len; i++){
    int count = 1;
    while(vec[i] == vec[i+1] && i < len - 1){
      count++;
      i++;
    }
    rle_vec_val.push_back(vec[i]);
    rle_vec_len.push_back(count);
  }

  bool zeros = std::all_of(vec.begin(), vec.end(), [](int i) { return i == 0; });

  if(zeros){
    return 0;
  }else{
    vector<int> x_lens;
    for(int i = 0; i < rle_vec_val.size(); i++){
      if(rle_vec_val[i] == x){
        x_lens.push_back(rle_vec_len[i]);
      }      
    }
    int max = *max_element(x_lens.begin(), x_lens.end());
    return max;
  }  
}

#endif
