#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// #include "cell.hpp"
#include "clone.hpp"


using namespace std;


// Simulate structural variations during cell division

int main(int argc, char const *argv[]) {
    int n_cycle;
    double dsb_rate;
    int min_dsb, max_dsb;
    int n_unrepaired;
    string fchr;

    double leap_size;

    string num_cells; // number of samples to take at each side

    double birth_rate, death_rate;
    double mutation_rate;
    // parameters for mean of dup/del size distributions
    // int mean_gain_size, mean_loss_size;
    int loc_type;

    int model_ID;   // not making differences, fitness is used to introduce selection
    int use_alpha;  // when alpha is used, genotype differences are considered
    double fitness;
    int genotype_diff;
    int growth_type;
    int norm_by_bin;

    int stat_type;
    double frac_cutoff;
    double min_freq;
    double max_freq;
    double delta;
    int use_std;

    string outdir, suffix; // output
    int write_shatterseek, write_cicos, write_sumstats;

    unsigned long seed;

    int verbose;

    namespace po = boost::program_options;

    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ;

    // parameters that are essential to simulations
    po::options_description required("Required parameters");
    required.add_options()
      ("n_cycle,n", po::value<int>(&n_cycle)->default_value(2), "number of cell cycles")
      ("min_dsb", po::value<int>(&min_dsb)->default_value(0), "minimal number of double strand breaks")
      ("max_dsb", po::value<int>(&max_dsb)->default_value(40), "maximal number of double strand breaks")
      ("n_unrepaired", po::value<int>(&n_unrepaired)->default_value(5), "mean number of unrepaired double strand breaks")
      ("type_dsb", po::value<int>(&n_unrepaired)->default_value(0), "type of unrepaired double strand breaks. 0: random; 1: biased; 2: fixed")

      ("birth_rate,b", po::value<double>(&birth_rate)->default_value(1), "birth rate")
      ("death_rate,d", po::value<double>(&death_rate)->default_value(0), "death rate")

      ("dsb_rate,r", po::value<double>(&dsb_rate)->default_value(20), "rate of double strand break per division")

      ("odir,o", po::value<string>(&outdir)->required()->default_value("./"), "output directory")
       ;

    po::options_description optional("Optional parameters");
    optional.add_options()
      // input files which specifies #sampled cells and CNA informaton
      ("fchr", po::value<string>(&fchr)->default_value(""), "TSV file with chromosome size infomation")

      // options related to model of evolution
      ("model", po::value<int>(&model_ID)->default_value(0), "model of evolution. 0: neutral; 1: selection")
      ("use_alpha", po::value<int>(&use_alpha)->default_value(1), "whether or not to use alpha in selection model. 0: use selection coefficient; 1: use alpha")
      ("fitness,f", po::value<double>(&fitness)->default_value(0), "fitness values of mutatants")
      ("genotype_diff", po::value<int>(&genotype_diff)->default_value(3), "type of karyotype difference (L1 distance) in simulating selection. 0: no; 1: L1 distance of CN; 2: Hamming distance of CN; 3: number of new mutations")
      ("norm_by_bin", po::value<int>(&norm_by_bin)->default_value(0), "whether or not to normalize karyotype difference by number of bins (segments) in the genome. 0: no; 1: yes")
      ("growth_type,t", po::value<int>(&growth_type)->default_value(0), "Type of growth when adding selection. 0: only birth; 1: change birth rate; 2: change death rate; 3: change both birth or death rate")

      // options related to summary statistics
      ("stat_type,s", po::value<int>(&stat_type)->default_value(4), "type of summary statistics. 0: variance; 1: clone-pairwise differences; 2: average CNP; 3: complete CNP; 4: sample-pairwise differences")
      ("use_std", po::value<int>(&use_std)->default_value(0), "whether or not to use standard deviation of pairwise divergence. If not, the variance is output")
      // ("bp_cutoff", po::value<double>(&BP_CUTOFF)->default_value(BP_CUTOFF), "threshold to determine whether a breakpoint can be detected or not")
      // ("bin_cutoff", po::value<double>(&BIN_CUOFF)->default_value(BIN_CUOFF), "threshold to determine whether a bin can be detected as altered or not")
      // ("frac_cutoff", po::value<double>(&frac_cutoff)->default_value(0.5), "cutoff of counting breakpoints")

      // options related to output
      ("suffix", po::value<string>(&suffix)->default_value(""), "suffix of output file")
      ("write_cicos", po::value<int>(&write_cicos)->default_value(1), "whether or not to write files for circos plot")
      ("write_shatterseek", po::value<int>(&write_shatterseek)->default_value(1), "whether or not to write files for shatterseek")
      ("write_sumstats", po::value<int>(&write_sumstats)->default_value(0), "whether or not to write summary statistics")

      ("seed", po::value<unsigned long>(&seed)->default_value(0), "seed used for generating random numbers")
      ("verbose", po::value<int>(&verbose)->default_value(0), "verbose level (0: default, 1: print information of final cells; 2: print information of all cells)")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
        if(vm.count("help")){
            cout << cmdline_options << endl;
            return 1;
        }
        if(vm.count("version")){
            cout << "simcell [version 0.1], a program to simulate copy number variations along a stochastic branching tree" << endl;
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception& e) {
          cerr << e.what() << endl;
          return 1;
    }

    unsigned long rseed = setup_rng(seed);

    if(verbose > 0){
        cout << "Random seed: " << rseed << endl;
        cout << "Simulating copy numbers called from bulk WGS and WES samples of cells sampled from a patient crypt" << endl;
        cout << "Fitness:\t" << fitness << endl;
        cout << "Mutation rate:\t" << mutation_rate << endl;
    }

    // cout << "Using Boost "
    //   << BOOST_VERSION / 100000     << "."  // major version
    //   << BOOST_VERSION / 100 % 1000 << "."  // minor version
    //   << BOOST_VERSION % 100                // patch level
    //   << std::endl;

    read_genome_info(fchr, CHR_LENGTHS, ARM_BOUNDS, CENT_STARTS, CENT_ENDS, TELO_ENDS1, TELO_ENDS2);

    Model start_model(model_ID, genotype_diff, growth_type, fitness, use_alpha);
    Cell_ptr start_cell = new Cell(1, 0, birth_rate, death_rate, dsb_rate, 0);
    genome g(1);
    start_cell->g = g;
    Clone* s = new Clone(1, 0);
    s->grow_with_dsb(start_cell, start_model, n_cycle, n_unrepaired, verbose);

    for(auto cell : s->curr_cells){
      cell->g.calculate_segment_cn();
    }

    if(write_cicos){
      for(auto cell : s->curr_cells){
        string fname_sv = outdir +"/" + "sv_data_c" + to_string(cell->cell_ID) + suffix + ".tsv";
        cell->write_sv(fname_sv);
        string fname_cn = outdir +"/" + "cn_data_c" + to_string(cell->cell_ID) + suffix + ".tsv";
        cell->write_cnv(fname_cn);
      }
    }

    if(write_shatterseek){
      for(auto cell : s->curr_cells){
        string fname_sv_ss = outdir +"/" + "SVData_c" + to_string(cell->cell_ID) + suffix + ".tsv";
        string fname_cn_ss = outdir +"/" + "CNData_c" + to_string(cell->cell_ID) + suffix + ".tsv";
        cell->write_shatterseek(fname_sv_ss, fname_cn_ss);
      }
    }

    s->print_all_cells(s->curr_cells, 1);

    if(write_sumstats){

    }

    return 0;
}
