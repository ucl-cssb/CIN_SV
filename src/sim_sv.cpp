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
    // int n_cycle;
    int n_cell;
    // double dsb_rate;
    int min_dsb, max_dsb;
    int n_dsb;
    int n_unrepaired;
    int chr_prob;
    string fchr;

    double leap_size;

    string n_cells; // number of samples to take at each side

    double birth_rate, death_rate;

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
    int write_shatterseek, write_cicos, write_sumstats, write_genome;

    unsigned long seed;
    int track_all;
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
      // ("n_cycle", po::value<int>(&n_cycle)->default_value(2), "number of cell cycles") // only meaningful when tracking one child
      ("n_cell,n", po::value<int>(&n_cell)->default_value(2), "size of final population")

      // ("dsb_rate,r", po::value<double>(&dsb_rate)->default_value(20), "rate of double strand break per division")
      ("min_dsb", po::value<int>(&min_dsb)->default_value(0), "minimal number of double strand breaks")
      ("max_dsb", po::value<int>(&max_dsb)->default_value(40), "maximal number of double strand breaks")
      ("n_dsb", po::value<int>(&n_dsb)->default_value(20), "maximal number of double strand breaks")
      ("n_unrepaired", po::value<int>(&n_unrepaired)->default_value(5), "number of unrepaired double strand breaks")
      ("chr_prob", po::value<int>(&chr_prob)->default_value(0), "the probability of double strand breaks across chromosomes. 0: random; 1: biased; 2: fixed")

      ("birth_rate,b", po::value<double>(&birth_rate)->default_value(1), "birth rate")
      ("death_rate,d", po::value<double>(&death_rate)->default_value(0), "death rate")

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
      ("write_cicos", po::value<int>(&write_cicos)->default_value(0), "whether or not to write files for circos plot")
      ("write_shatterseek", po::value<int>(&write_shatterseek)->default_value(0), "whether or not to write files for shatterseek")
      ("write_sumstats", po::value<int>(&write_sumstats)->default_value(1), "whether or not to write summary statistics")
      ("write_genome", po::value<int>(&write_genome)->default_value(0), "whether or not to write derivative genome")

      ("track_all", po::value<int>(&track_all)->default_value(0), "whether or not to keep track of all cells")
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
            cout << "simsv [version 0.1], a program to simulate SVs by generating and repairing double strand breaks along a stochastic branching tree" << endl;
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception& e) {
          cerr << e.what() << endl;
          return 1;
    }

    unsigned long rseed = setup_rng(seed);

    cout << "Random seed: " << rseed << endl;

    // cout << "Using Boost "
    //   << BOOST_VERSION / 100000     << "."  // major version
    //   << BOOST_VERSION / 100 % 1000 << "."  // minor version
    //   << BOOST_VERSION % 100                // patch level
    //   << std::endl;

    read_genome_info(fchr, CHR_LENGTHS, ARM_BOUNDS, CENT_STARTS, CENT_ENDS, TELO_ENDS1, TELO_ENDS2, verbose);

    vector<int> selected_chr;
    get_chr_prob(chr_prob, selected_chr);
    // int diff = max_dsb - min_dsb;
    // int rdm = myrng(diff);
    // n_dsb = min_dsb + diff;
    // n_dsb = gsl_ran_poisson(r, dsb_rate);

    Model start_model(model_ID, genotype_diff, growth_type, fitness, use_alpha);
    Cell_ptr start_cell = new Cell(1, 0, birth_rate, death_rate, n_dsb, n_unrepaired, 0);
    genome g(1);
    start_cell->g = g;
    Clone* s = new Clone(1, 0);
    s->grow_with_dsb(start_cell, start_model, n_cell, track_all, verbose);

    vector<Cell_ptr> final_cells;
    if(track_all){
      cout << "Tracking all cells" << endl;
      final_cells = s->cells;
    }else{
      final_cells = s->curr_cells;
    }

    // merge segments with the same CN by haplotype
    for(auto cell : final_cells){
      cell->g.calculate_segment_cn();
    }

    string filetype = ".tsv";

    if(write_cicos){
      for(auto cell : final_cells){
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname_sv = outdir +"/" + "sv_data_c" + midfix + filetype;
        cell->write_sv(fname_sv);
        string fname_cn = outdir +"/" + "cn_data_c" + midfix + filetype;
        cell->write_cnv(fname_cn);
      }
    }

    // write CN and SV data to tsv - for ShatterSeek
    if(write_shatterseek){
      for(auto cell : final_cells){
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname_sv_ss = outdir +"/" + "SVData_c" + midfix + filetype;
        string fname_cn_ss = outdir +"/" + "CNData_c" + midfix + filetype;
        cell->write_shatterseek(fname_sv_ss, fname_cn_ss);
      }
    }

    s->print_all_cells(final_cells, verbose);

    // write summary statistics in the simulation
    string fname_stat_sim = outdir +"/" + "sumStats_sim" + filetype;
    ofstream fout(fname_stat_sim);
    fout << "nCell\tnDSB\tnUnrepair\tnComplex\tnMbreak\tnTelofusion\n";
    fout << final_cells.size() << "\t" << n_dsb << "\t" << n_unrepaired << "\t" + to_string(s->n_complex_path) + "\t" + to_string(s->n_path_break) + "\t" + to_string(s->n_telo_fusion) << endl;
    fout.close();

    if(write_sumstats){
      for(auto cell : final_cells){
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname_stat = outdir +"/" + "sumStats_total_c" + midfix + filetype;
        string fname_stat_chr = outdir +"/" + "sumStats_chrom_c" + midfix + filetype;
        cell->write_summary_stats(fname_stat, fname_stat_chr);
      }
    }

    if(write_genome){
      for(auto cell : final_cells){
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname = outdir +"/" + "genome_c" + midfix + filetype;
        cell->write_genome(fname);
      }
    }

    return 0;
}
