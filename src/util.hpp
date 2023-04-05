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
// #include <numeric>   // for partial sum

#include <unistd.h>

#include <boost/filesystem.hpp>

// For defining empirical_cumulative_distribution_function
// #include <iterator>
// #include <stdexcept>


using namespace std;


typedef map<pair<int, int>, int> pcn;   // copy number at a position
typedef map<pair<int, int>, double> dpcn;   // copy number at a position


// A genomic position
struct pos_bin{
  int chr;
  int start;
  int end;  
};

struct pos_bp{
  int chr;
  int loc;
  int side;   

  bool operator<(pos_bp const &other) const {
    return (chr < other.chr ||
      (chr == other.chr && loc < other.loc));
  }
};


const int FAIL = 1;

const int NUM_CHR = 22;
const int NORM_PLOIDY = 2;
// const int NUM_LOC = 5000;
// const int BIN_SIZE = 500000; 

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

const int NUM_SVTYPE = 8;
// DEL: +/-, DUP: -/+, same orientation as DELREAL/DUPREAL, but copy number may not change
enum SV_type{NONE, DEL, DUP, H2HINV, T2TINV, BND, DELREAL, DUPREAL}; // only for variant adjacency, BND: intra-chromosomal

enum SStat_type{ALL};

enum Growth_type{ONLY_BIRTH, CHANGE_BIRTH, CHANGE_DEATH, CHANGE_BOTH};

enum Telo_type{NONTEL, PTEL, QTEL, COMPLETE};
enum Adj_type{INTERVAL, REF, VAR};   // 0: interval, 1: reference, 2: variant
enum Junc_type{HEAD, TAIL};


const double PROB_INTER = 1e-9;

const double SURVIVAL_D = 0.00039047;
const double SURVIVAL_C = -0.036132164;

// Chromosome and arm scores, taken from Davoli et al. Table S6A, opposite sign
const double CHR_SCORE[] = {-0.143640496, 0.638322635, 0.597508197, 0.106407616, -0.785208831, -0.664148445, 3.039521587, 1.650903175, 0.765873656, -1.23443224, 0.210103365,
                        1.720482377, -1.207617162, -0.712581034, -0.751608856, -1.277797927, -0.784673321, -1.428496154, 0.809097907, 1.780741874, 1.568732394, -1.576297101};

// use 0 for short p arms with few genes
// TSG-OG score, taken from Davoli et al. Table S6B (top genes only), opposite sign
const double ARM_SCORE[] = {-2.194424482, 1.226691224,	
-0.058765864,	0.722951857,	
-0.23241358,	2.77507441,	
0.343354369,	1.236932406,	
1.29642446,	-0.728377682,	
-0.841818493,	-0.783217918,	
5.195591398,	4.588576125,	
1.391378151,	1.26397449,	
-1.495356436,	1.979101476,	
-3.665105263,	-0.068322404,	
0.590530233,	-0.200811736,	
2.960242754,	1.661808926,	
0.,	-1.906547855,	  // 13
0.,	-0.784903448,	
0.,	-0.893062731,	
0.450231325,	-3.157515406,	
-2.803987461,	0.168961686,  
0.,	-2.868015464,	
1.632886207,	-0.908804749,	
0.602858757,	1.609560694,	
0.,	-0.990394366,	
0.,	-1.209528986};


// original TSG-OG-Ess score, to update
const double ARM_SCORE2[] = {0.181435341, 
-2.450861322, 
-2.217704595, 
-2.978829436, 
-2.485106996, 
-4.091920145, 
-1.348359223, 
-2.191723658, 
-3.234086331, 
-1.47302289, 
-1.811969178, 
-0.008782082, 
-7.470580645, 
-5.237576125, 
-2.092189076, 
-4.42784949, 
0.097757426, 
-3.55250738, 
2.751381579, 
-0.863211293, 
-2.884669767, 
-1.434382641, 
-4.166536232, 
-4.739354254, 
0.767049505, 
-1.931944828, 
-1.560444649, 
-3.274554217, 
0.53970028, 
0.649470219, 
-3.296263091, 
-0.403166667, 
1.991530928, 
-4.057684483, 
-0.415737467, 
-3.286858757, 
-3.772268786, 
0.249985915, 
0.389657005};


// https://www.sciencedirect.com/book/9780124046313/benign-and-pathological-chromosomal-imbalances
// The five human acrocentric chromosomes are numbered 13, 14, 15, 21, and 22. 
// They all have a cytogenetically similar short arm that is extremely gene-poor. 
// Their main contribution for the cell is that the acrocentric short arms are carriers of the nucleolus organizing regions (NOR) in subbands p12
// double arm_scores1[] = {-2.194424482, 1.226691224,	-0.058765864,	0.722951857,	-0.23241358,	2.77507441,	0.343354369,	1.236932406,	1.29642446,	-0.728377682,	-0.841818493,	-0.783217918,	5.195591398,	4.588576125,	1.391378151,	1.26397449,	-1.495356436,	1.979101476,	-3.665105263,	-0.068322404,	0.590530233,	-0.200811736,	2.960242754,	1.661808926,	-1.906547855,	-0.784903448,	-0.893062731,	0.450231325,	-3.157515406,	-2.803987461,	0.168961686,  0.,	-2.868015464,	1.632886207,	-0.908804749,	0.602858757,	1.609560694,	-0.990394366,	-1.209528986};
// Indexing of chromosome arms (the second component of each pair is 1 or 2 to refere to the p-arm or the q-arm, respectively)
// armlist:=[[1, 1], [1, 2], [2, 1], [2, 2], [3, 1], [3, 2], [4, 1], [4, 2], [5, 1], [5, 2], [6, 1], [6, 2], [7, 1], [7, 2], [8, 1], [8, 2], [9, 1], [9, 2], [10, 1], [10, 2], [11, 1], [11, 2], [12, 1], [12, 2],  [13, 2],  [14, 2],  [15, 2], [16, 1], [16, 2], [17, 1], [17, 2], [18, 1], [18, 2], [19, 1], [19, 2], [20, 1], [20, 2],  [21, 2],  [22, 2]]


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
  if(verbose) cout << "\nReading reference genome information from " << filename << endl;

  ifstream infile(filename);
  if(!infile.is_open()){
    std::cerr << "Error: open of input data unsuccessful: " << filename << std::endl;
    exit(FAIL);
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


// read bins, to be consistent with real data where copy numbers are often binned
void read_bins(const string& filename, vector<pos_bin>& bins, int verbose = 0){
  if(verbose) cout << "\nReading bins from " << filename << endl;

  ifstream infile(filename);
  if(!infile.is_open()){
    std::cerr << "Error: open of input data unsuccessful: " << filename << std::endl;
    exit(FAIL);
  }

  std::string line;
  getline(infile, line);  // skip header
  while(!getline(infile, line).eof()){
    // if(verbose) cout << line << endl;
    if(line.empty()){
      continue;
    }
    std::vector<std::string> split;
    std::string buf;
    stringstream ss(line);
    while (ss >> buf) split.push_back(buf);
    assert(split.size() == 3);

    int chr = atoi(split[0].c_str());
    if(chr > NUM_CHR) break;  
    int start = atoi(split[1].c_str());
    int end = atoi(split[2].c_str());
    pos_bin p{chr, start, end};
    bins.push_back(p);
  }
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

  int nb = 1;   // number of chromosomes with biased probability of DSBs
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

  double p0 = 2.0 / 3.0;  // ? 
  double p1 = p0 / nb;      // p1 for selected chromosome
  double p2 = (1 - p0) / (NUM_CHR - nb);
  for(int i = 0; i < NUM_CHR; i++){
    CHR_PROBS[i] = p2;
  }
  for(int i = 0; i < selected_chr.size(); i++){
    CHR_PROBS[i] = p1;
  }
}


void get_chr_prob_from_file(const string& filename, int verbose = 0){
  if(verbose) cout << "\nReading chromosome probability from " << filename << endl;

  ifstream infile(filename);
  if(!infile.is_open()){
    std::cerr << "Error: open of input data unsuccessful: " << filename << std::endl;
    exit(FAIL);
  }

  std::string line;
  int i = 0;
  while(!getline(infile, line).eof()){
    if(verbose) cout << line << endl;
    if(line.empty()){
      continue;
    }
    if(i >= NUM_CHR) break;   
    std::vector<std::string> split;
    std::string buf;
    stringstream ss(line);
    while (ss >> buf) split.push_back(buf);
    assert(split.size() == 1);

    double len = atof(split[0].c_str());
    CHR_PROBS[i++] = len;
  }
}


// file format: "chr1"    "pos1"    "strand1" "chr2"  "pos2"    "strand2"
// use set to ensure uniqueness
// vector<int>& intra_distance: distances of breakpoints on the same chromosome
// map<int, int>& inter_chrom: connection of different chromosomes
vector<pos_bp> get_bp_from_file(const string& filename, int verbose = 0){
  if(verbose) cout << "\nReading breakpoints from " << filename << endl;

  ifstream infile(filename);
  if(!infile.is_open()){
    std::cerr << "Error: open of input data unsuccessful: " << filename << std::endl;
    exit(FAIL);
  }

  std::string line;
  getline(infile, line);  // skip header

  set<pos_bp> bps0;
  while(!getline(infile, line).eof()){
    if(verbose) cout << line << endl;
    if(line.empty()){
      continue;
    } 
    std::vector<std::string> split;
    std::string buf;
    stringstream ss(line);
    while (ss >> buf) split.push_back(buf);
    assert(split.size() == 6);

    int chr = atoi(split[0].c_str());
    int pos = atoi(split[1].c_str());
    string strand = split[2].c_str();
    int side = 0;
    if(strand == "-") side = 1;    
    pos_bp bp{chr, pos, side};
    bps0.insert(bp);

    chr = atoi(split[3].c_str());
    pos = atoi(split[4].c_str());
    strand = split[5].c_str();
    side = 0;
    if(strand == "-") side = 1; 
    pos_bp bp2{chr, pos, side};
    bps0.insert(bp2);   
  }
  if(verbose > 1) cout << "Finish reading " << bps0.size() << " breakpoints " << endl;
  vector<pos_bp> bps(bps0.size());
  std::copy(bps0.begin(), bps0.end(), bps.begin());
  if(verbose > 1) cout << "Finish copying " << bps.size() << " breakpoints " << endl;
  
  return bps;
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
    case DUPREAL: type_str = "DUPREAL"; break;
    case DELREAL: type_str = "DELREAL"; break;    
    case BND: type_str = "BND"; break;
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


  // get survival probability of normal cell (compute once)
  double get_surv_prob_normal_chr(double selection_strength){
    double score = 0.0;
    for(int i = 0; i < NUM_CHR; i++){  
      score += CHR_SCORE[i] * 2;
    }
    double surv_prob = exp(SURVIVAL_D * score + SURVIVAL_C);
    surv_prob = pow(surv_prob, selection_strength);
    return surv_prob;
  }


  double get_surv_prob_normal_arm(double selection_strength){
    double score = 0.0;
    for(int i = 0; i < NUM_CHR * 2; i += 2){
        score += ARM_SCORE[i] * 2;
        score += ARM_SCORE[i + 1] * 2;
      }       

    double surv_prob = exp(SURVIVAL_D * score + SURVIVAL_C);
    surv_prob = pow(surv_prob, selection_strength);
    return surv_prob;
  }



#endif
