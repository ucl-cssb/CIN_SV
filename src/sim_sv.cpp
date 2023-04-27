#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "clone.hpp"

using namespace std;


// Simulate structural variations during cell division

void get_bin_number(string fbin, vector<pos_bin>& bins, vector<int>& bin_number, int verbose = 0){
    read_bins(fbin, bins, verbose);

    vector<int> bin_count(NUM_CHR, 0);
    if(verbose) cout << "There are " << bins.size() << " bins" << endl;
    for(auto bin : bins){
      // cout << bin.chr << " " << bin.start << " " << bin.end << endl;
      if(bin.chr > NUM_CHR) break;
      assert(bin.chr >= 1);
      bin_count[bin.chr - 1]++;    // real data starts from chr 1
    }

    // bin_number[0] = bin_count[0] when calling library function
    // std::partial_sum(bin_count.begin(), bin_count.end(), bin_number.begin());
    for(int i = 0; i < NUM_CHR - 1; i++){
      bin_number[i + 1] = bin_number[i] + bin_count[i];
    }

    if(verbose){
      cout << "number of bins for each chromosome" << endl;
      for(int i = 0; i < bin_count.size(); i++) {
        cout << i + 1 << "\t" << bin_count[i] << "\t" << bin_number[i] << endl;
      }        
    }   
}


void write_rck_output(vector<Cell_ptr>& final_cells, string outdir, int verbose = 0){
    for(auto cell : final_cells){
      map<int, set<int>> bps_by_chr;
      cell->g->get_bps_per_chr_orig(bps_by_chr, verbose);

      cell->g->calculate_segment_cn(bps_by_chr , verbose);

      //  and circos plot
      if(verbose > 0) cout << "\nWrite output for RCK graph" << endl;
      string subdir = outdir + "/" +  "c" + to_string(cell->cell_ID);
      boost::filesystem::path dir(subdir);
      boost::filesystem::create_directory(dir);
      string fname_cn = subdir + "/" + "rck.scnt.tsv";      
      string fname_sv = subdir + "/" + "rck.acnt.tsv";
      cell->write_rck(fname_cn, fname_sv, verbose);

      // same format as 
      if(verbose > 0) cout << "\nWrite output for chromosome graph plot" << endl;
      fname_cn = subdir + "/" + "copy_number.CN_opt.phased";      
      fname_sv = subdir + "/" + "SVs.CN_opt.phased";
      cell->write_plot(fname_cn, fname_sv, verbose);
    }  
}


void write_stat_bin_cn(Clone* s, vector<Cell_ptr>& final_cells, int num_loc, int verbose = 0){
    map<int, double> cell_ploidy;
    map<int, vector<double>> loc_cn;
    map<int, vector<double>> loc_cnA;
    map<int, vector<double>> loc_cnB;
    vector<int> ids;

    for(auto cell : final_cells){
        cell_ploidy[cell->cell_ID] = cell->ploidy;
        loc_cn[cell->cell_ID] = cell->g->bin_tcn;
        loc_cnA[cell->cell_ID] = cell->g->bin_cnA;
        loc_cnB[cell->cell_ID] = cell->g->bin_cnB;
        ids.push_back(cell->cell_ID);
    }
    if(verbose > 0) cout << "summary statistics for " << num_loc << " bins" << endl;

    double pga = s->get_pga(loc_cn, num_loc, cell_ploidy);
    double pgaA = s->get_pga(loc_cnA, num_loc, cell_ploidy, 1);
    double pgaB = s->get_pga(loc_cnB, num_loc, cell_ploidy, 1);
    pair<double, double> div = s->get_pairwise_divergence(ids, loc_cn, num_loc, cell_ploidy);
    pair<double, double> divA = s->get_pairwise_divergence(ids, loc_cnA, num_loc, cell_ploidy, 1);
    pair<double, double> divB = s->get_pairwise_divergence(ids, loc_cnB, num_loc, cell_ploidy, 1);
    cout << pga << "\t" << div.first << "\t" << div.second << "\t" << pgaA << "\t" << divA.first << "\t" << divA.second << "\t" << pgaB << "\t" << divB.first << "\t" << divB.second;
}


void write_stat_bp_freq(Clone* s, vector<Cell_ptr>& final_cells, int verbose = 0){
    // get number of new breakpoints locations of all current cells
    map<int, int> bp_freq;  // bp count: 1 to ncell, bp frequency
    s->get_bp_freq(final_cells, bp_freq);
    // based on https://stackoverflow.com/questions/9370945/finding-the-max-value-in-a-map
    // auto pr = std::max_element(std::begin(bp_freq), std::end(bp_freq), [](const pair<int, int>& p1, const pair<int, int>& p2){ return p1.second < p2.second;});
    int nbp = accumulate(std::begin(bp_freq), std::end(bp_freq), 0,
                                        [](const int previous, const std::pair<const int, int>& p)
                                        { return previous + p.second; });
    assert(bp_freq.size() == final_cells.size());
    for(auto bpq: bp_freq){
      if(verbose > 1) cout << "\t" << bpq.first << "\t" << bpq.second << endl;       
      // normalize count by the total number of counts to get similar scalce to PGA and DIV
      double freq = 0;
      if(nbp > 0){
        freq = (double) bpq.second / nbp;
      }
      cout << "\t" << freq;
    }  
}


int main(int argc, char const *argv[]){
    // int n_cycle;
    int n_cell;
    double dsb_rate, dsb_rate2;
    // int min_dsb, max_dsb;
    int n_dsb;
    double frac_unrepaired;
    double circular_prob;
    int chr_prob;
    // string target_chrs;
    string fchr_prob;
    string fchr, fbin, fbp;
    double n_local_frag;  // mean number of breakpoints introduced by local fragmentation during mitosis
    double frac_unrepaired_local;
    double prob_wgd;
    int pair_type; // type of joining pair of breakpoints
    double prob_correct_repaired;
    int div_break;   // ID of division when DSBs occurs
    int only_repair_new; // only repair new DSBs

    double birth_rate, death_rate;

    int model_ID;   
    int selection_type;
    int growth_type;
    double selection_strength;

    string outdir, suffix; // output
    int write_shatterseek, write_rck, write_sumstats, write_genome, write_bin, write_selection;
    int bin_level_sumstat; //

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
      ("div_break", po::value<int>(&div_break)->default_value(0), "maximum ID of cell division when double strand breaks occurs")
      ("dsb_rate,r", po::value<double>(&dsb_rate)->default_value(0), "mean number or rate of double strand breaks per cell division assuming gradutual evolution, which follows Poisson distribution")
      // ("dsb_rate2", po::value<double>(&dsb_rate2)->default_value(0), "mean constant rate of double strand breaks per division in gradutual evolution after a punctuated event")
      // ("min_dsb", po::value<int>(&min_dsb)->default_value(0), "minimal number of double strand breaks (lower bound of uniform distribution of the number of double strand breaks)")
      // ("max_dsb", po::value<int>(&max_dsb)->default_value(40), "maximal number of double strand breaks (upper bound of uniform distribution of the number of double strand breaks)")
      ("n_dsb", po::value<int>(&n_dsb)->default_value(20), "number of double strand breaks introduced in G1 (in catastrophic events, fixed at each cell cycle)")
      ("frac_unrepaired", po::value<double>(&frac_unrepaired)->default_value(0), "fraction of unrepaired double strand breaks in G1, default: all breaks are repaired.")
      ("prob_wgd", po::value<double>(&prob_wgd)->default_value(0), "probability of whole genome doubling")
      ("n_local_frag", po::value<double>(&n_local_frag)->default_value(0), "mean number of double strand breaks introduced by local fragmentation during mitosis")
      ("frac_unrepaired_local", po::value<double>(&frac_unrepaired_local)->default_value(1), "number of unrepaired double strand breaks in local fragmentation during mitosis, default: all breaks are not repaired")
      ("pair_type", po::value<int>(&pair_type)->default_value(0), "type of joining pair of breakpoints, default: randomly joined. 1: joining based on distance between breakpoints")  
      ("prob_correct_repaired", po::value<double>(&prob_correct_repaired)->default_value(0.5), "the probability of correctly repaired double strand breaks in G1")          
      ("chr_prob", po::value<int>(&chr_prob)->default_value(0), "the types of assigning probability of double strand breaks across chromosomes. 0: random; 1: biased; 2: fixed")
      ("only_repair_new", po::value<int>(&only_repair_new)->default_value(0), "whether or not to only repair new DSBs introduced in each cell cycle")
      ("circular_prob", po::value<double>(&circular_prob)->default_value(0), "the probability of a frament without centromere and telomeres forming circular DNA (ecDNA)")
      // ("target_chrs", po::value<string>(&target_chrs)->default_value(""), "biased chromosomes to introduce breaks, total number followed by ID of each chromosome")
      ("fchr_prob", po::value<string>(&fchr_prob)->default_value(""), "the file containing the probability of double strand breaks on each chromosome")

      ("birth_rate,b", po::value<double>(&birth_rate)->default_value(1), "birth rate")
      ("death_rate,d", po::value<double>(&death_rate)->default_value(0), "death rate")

      ("odir,o", po::value<string>(&outdir)->required()->default_value("./"), "output directory")
       ;
 
    po::options_description optional("Optional parameters");
    optional.add_options()
      // input files which specifies #sampled cells and CNA informaton
      ("fchr", po::value<string>(&fchr)->default_value(""), "TSV file with chromosome size infomation")
      ("fbin", po::value<string>(&fbin)->default_value(""), "TSV file with bin size infomation")
      ("fbp", po::value<string>(&fbp)->default_value(""), "TSV file with breakpoint infomation")

      // options related to model of evolution
      ("model", po::value<int>(&model_ID)->default_value(0), "model of evolution. 0: neutral; 1: selection")
      ("selection_type", po::value<int>(&selection_type)->default_value(0), "types used to define selection strength, 0: chr-level, 1: arm-level")
      ("growth_type,t", po::value<int>(&growth_type)->default_value(0), "type of growth when adding selection. 0: only birth; 1: change birth rate; 2: change death rate; 3: change both birth or death rate")
       ("selection_strength,d", po::value<double>(&selection_strength)->default_value(1), "strength of selection, the larger the value, the stronger the selection is")

      // options related to summary statistics

      // options related to output
      ("suffix", po::value<string>(&suffix)->default_value(""), "suffix of output file")
      ("bin_level_sumstat", po::value<int>(&bin_level_sumstat)->default_value(0), "whether or not to output summary statistics at bin level")
      ("write_rck", po::value<int>(&write_rck)->default_value(0), "whether or not to write files in RCK format")
      ("write_bin", po::value<int>(&write_bin)->default_value(0), "whether or not to write copy numbers in bin format")
      ("write_shatterseek", po::value<int>(&write_shatterseek)->default_value(0), "whether or not to write files for shatterseek")
      ("write_sumstats", po::value<int>(&write_sumstats)->default_value(1), "whether or not to write summary statistics")
      ("write_genome", po::value<int>(&write_genome)->default_value(0), "whether or not to write derivative genome")
      ("write_selection", po::value<int>(&write_selection)->default_value(0), "whether or not to write information on selection score for each cell")

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
    if(verbose > 0) cout << "Random seed: " << rseed << endl;

    // cout << "Using Boost "
    //   << BOOST_VERSION / 100000     << "."  // major version
    //   << BOOST_VERSION / 100 % 1000 << "."  // minor version
    //   << BOOST_VERSION % 100                // patch level
    //   << std::endl;

    read_genome_info(fchr, CHR_LENGTHS, ARM_BOUNDS, CENT_STARTS, CENT_ENDS, TELO_ENDS1, TELO_ENDS2, verbose);

    vector<int> selected_chr;
    // if(target_chrs == ""){
    if(fchr_prob == ""){
      get_chr_prob(chr_prob, selected_chr);
    }else{
      // get_vals_from_str(selected_chr, target_chrs);
      get_chr_prob_from_file(fchr_prob, verbose);
    }

    vector<pos_bin> bins;
    vector<int> bin_number(NUM_CHR, 0);  // store the cumulative number of bins for each chromosome   
    if(fbin != ""){
      get_bin_number(fbin, bins, bin_number, verbose);
    }

    vector<pos_bp> bps;
    vector<double> bp_fracs;
    if(fbp != ""){
      bps = get_bp_from_file(fbp, bp_fracs, verbose);
      // assert(bps.size() % 2 == 0);  // some breakpoints are duplicated
      // n_dsb = bps.size() / 2;
      if(verbose > 1) cout << "There are " << bps.size() << " known breakpoints " << endl; 
    }
      
    // int diff = max_dsb - min_dsb;
    // int rdm = myrng(diff);
    // n_dsb = min_dsb + diff;
    
    Model start_model(model_ID, selection_type, selection_strength, growth_type);
    int n_unrepaired = round(n_dsb * frac_unrepaired);
    Cell_ptr start_cell = new Cell(1, 0, birth_rate, death_rate, prob_wgd, dsb_rate, n_dsb, n_unrepaired, 0, div_break, only_repair_new);
    Clone* s = new Clone(1, 0, n_cell, bps, bp_fracs, frac_unrepaired, n_local_frag, frac_unrepaired_local, circular_prob, pair_type, prob_correct_repaired, track_all);
    
    if(verbose > 0){
      cout << "Start cell growth with " << n_unrepaired << " unrepaired DSBs" << endl;
      if(n_local_frag > 0) cout << " introducing local fragmentation randomly with mean number of breaks " << n_local_frag << endl;
    }
    s->grow_with_dsb(start_cell, start_model, verbose);
    if(verbose > 0) cout << "\nFinish cell growth " << endl;

    vector<Cell_ptr> final_cells;
    final_cells = s->curr_cells;
    assert(final_cells.size() == n_cell);
    // if(track_all){
    //   cout << "Tracking all cells" << endl;
    //   final_cells = s->cells;
    // }

    if(verbose > 0) cout << "\nComputing breakpoints on each chromosome across all cells for more intuitive comparison and visualization" << endl;
    // get breakpoints across cells first to align breakpoints for better visualization
    map<int, set<int>> bps_by_chr;
    for(auto cell : final_cells){
      if(verbose > 0){
        cout << "Cell " << cell->cell_ID << endl;
      }
      cell->g->get_bps_per_chr(bps_by_chr, verbose);
    }

    if(outdir != "./"){
      boost::filesystem::path dir(outdir);
      boost::filesystem::create_directories(dir);
    }

    string filetype = ".tsv";
    if(write_sumstats){
      // write summary statistics in the simulation
      string fname_stat_sim = outdir + "/" + "sumStats_sim" + filetype;
      ofstream fout(fname_stat_sim);
      fout << "nCell\tnDSB\tnUnrepair\tnComplex\tnMbreak\tnTelofusion\n";
      fout << final_cells.size() << "\t" << n_dsb << "\t" << n_unrepaired << "\t" + to_string(s->n_complex_path) + "\t" + to_string(s->n_path_break) + "\t" + to_string(s->n_telo_fusion) << endl;
      // fout.close();
    }

    string fname_sel = "";
    if(model_ID == 1 && write_selection){
        if(verbose > 0) cout << "\nWrite survival probability and selection coefficient of each cell" << endl;    
        fname_sel = outdir + "/" + "selection" + filetype;
        ofstream fout(fname_sel);
        fout << "cell\tprob_survival\tselection_coef\n";
        for(auto cell : final_cells){
          fout << cell->cell_ID  << "\t" << cell->surv_prob << "\t" << cell->fitness << endl;
        }
    }

    if(verbose > 0) cout << "\nComputing CN for each cell " << endl;


    for(auto cell : final_cells){
      // verbose = 1;
      if(verbose > 0){
        cout << "\nCell " << cell->cell_ID << endl;
      }
      // cell->print_bp_adj();
      // merge segments with the same CN by haplotype
      cell->g->calculate_segment_cn(bps_by_chr, verbose);

      if(bin_level_sumstat || write_bin){
        assert(fbin != "" && bins.size() > 0);
        cell->g->get_cn_bin(bins, bin_number, verbose);
      }

      // compute and print summary statistics, used for ABC
      cell->get_summary_stats(verbose);

      // write summary statistics for each cell
      if(write_sumstats){
        if(verbose > 0) cout << "\nWrite output for summary statistics" << endl;

        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname_stat = outdir +"/" + "sumStats_total_c" + midfix + filetype;
        string fname_stat_chr = outdir +"/" + "sumStats_chrom_c" + midfix + filetype;
        cell->write_summary_stats(fname_stat, fname_stat_chr, bin_level_sumstat, verbose);
      }

      // write CN and SV data to tsv - for ShatterSeek
      if(write_shatterseek){
        if(verbose > 0) cout << "\nWrite output for ShatterSeek" << endl;
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname_sv_ss = outdir + "/" + "SVData_c" + midfix + filetype;
        string fname_cn_ss = outdir + "/" + "CNData_c" + midfix + filetype;
        cell->write_shatterseek(fname_cn_ss, fname_sv_ss, verbose);
      }

      if(write_genome){
        if(verbose > 0) cout << "\nWrite output for the derivative genome" << endl;
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname = outdir +"/" + "genome_c" + midfix + filetype;
        cell->write_genome(fname);
      }

      if(write_bin){
        if(verbose > 0) cout << "\nWrite copy number in bins of fixed size" << endl;       
        string midfix = to_string(cell->cell_ID) + "_div" + to_string(cell->div_occur) + suffix;
        string fname_cn = outdir + "/" + "CNBin_c" + midfix + filetype;
        cell->write_cn_bin(fname_cn, bins, verbose);
      }
    }      
    
    // at the end to avoid conflict with merged copy numbers
    if(write_rck){
      if(verbose > 0) cout << "\nWrite output in RCK format" << endl;
      write_rck_output(final_cells, outdir, verbose);
    }

    if(track_all){
      string fname_tree = outdir + "/" + "cell_lineage.tsv";  
      ofstream fout_lineage(fname_tree);
      s->print_all_cells(s->cells, fout_lineage, verbose);
      // fout_lineage.close();
    }

    // output all the summary statistics in console for ABC 
    if(bin_level_sumstat){   
      int num_loc = bins.size();  
      write_stat_bin_cn(s, final_cells, num_loc, verbose);     
    }
    write_stat_bp_freq(s, final_cells, verbose);
    // for(auto cell : final_cells){  
    //   vector<pos_cn> dups;
    //   vector<pos_cn> dels;
    //   // seems not necessary as the breakpoints can be telled from copy number segments
    //   // cell->g->get_pseudo_adjacency(dups, dels, verbose);        
    //   cell->print_total_summary(dups, dels, verbose);
    // }
    cout << endl;

    delete s;
    gsl_rng_free(r);

    return 0;
}
