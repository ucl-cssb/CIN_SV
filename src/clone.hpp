#ifndef CLONE_HPP
#define CLONE_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>

#include <assert.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "util.hpp"
#include "cell.hpp"

using namespace std;


typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> VARIANCE;


// Node of a binary lineage tree, used for sampling certain number of nodes at each side
struct node {
    int data;
    node* left;
    node* right;
    int flag; // whether the node is available (dead) or not

    node(int data){
        this->data = data;
        this->flag = 0;
        left = NULL;
        right = NULL;
    }

    node(int data, int flag){
        this->data = data;
        this->flag = flag;
        left = NULL;
        right = NULL;
    }
};


node* search_node(int key, node* p){
    // cout << "search node " << key << endl;
    if(p == NULL) return NULL;

    if(key == p->data){
        return p;
    }
    else{
        node* t = search_node(key, p->left);
        if(t == NULL){
            return search_node(key, p->right);
        }
        else{
            return t;
        }
    }
}


// find all leaves below an internal node
void get_leaves_below(node* p, vector<int>& nodes){
    if(p == NULL) return;
    if(p->left == NULL && p->right == NULL){
        nodes.push_back(p->data);
    }else{
        get_leaves_below(p->left, nodes);
        get_leaves_below(p->right, nodes);
    }
}


// find all leaves below an internal node which are available
void get_leaves_below_avail(node* p, vector<int>& nodes){
    if(p == NULL) return;
    if(p->left == NULL && p->right == NULL && p->flag >= 0){
        nodes.push_back(p->data);
    }else{
        get_leaves_below_avail(p->left, nodes);
        get_leaves_below_avail(p->right, nodes);
    }
}


void destroy_tree(node* root){
    if(root != NULL)
    {
      destroy_tree(root->left);
      destroy_tree(root->right);
      delete root;
    }
}


// The model of CNA evolution
class Model{
public:
    int model_ID;
    int growth_type;   // how cells grow under selection
    int selection_type;
    double psurv_norm;

    Model(){
        model_ID = 0;
        selection_type = 0;
        growth_type = 0;
        psurv_norm = get_surv_prob_normal_chr();
    }

    Model(int model_ID, int selection_type, int growth_type){
        this->model_ID = model_ID;
        this->selection_type = selection_type;
        this->growth_type = growth_type;
        if(selection_type == 0){
            psurv_norm = get_surv_prob_normal_chr();
        }else{
            psurv_norm = get_surv_prob_normal_arm();
        }       
    }
};


/*
A Clone represents a population of cells
*/
class Clone{
public:
    int clone_ID;
    string name;

    // Used for tracking clone/deme expansion
    int parent_ID;

    vector<Cell_ptr> curr_cells;   // only available cells at present
    vector<Cell_ptr> cells;   // all cells at present
    vector<int> id_curr_cells;      // Only store the IDs of cell to save space

    // vector<pair<int, int>> lineage_cells;  // record lineage relationship of all cells
    // map<int, vector<double>> fits_by_time;

    // vector<string> bps; // breakpoints

    // double dsb_rate;       // assuming constant mutation rate for one clone
    // map<int, vector<int>> cnp_all[num_loc];     // store all CNVs in a big array
    // vector<int> id_by_loc;      // store cell ID sorted by location to facilitate sampling
    // map<int, pair<double, double>> grates;  // birth/death rates for cells with fitness (dis)advantages


    int ntot;   // total number of cells generated so far, used to increase cell ID
    // int n_cycle;   // number of cell cycles gone through
    int n_complex_path;  // number of paths with >1 centromeres, not useful when included for each cell as it is split into daughter cells
    int n_path_break;   // a complex path may have 2 or more centromeres and additional breaks are introduced to get single-centromere paths
    int n_telo_fusion;  // telomere removal and fusion as in BFB

    double time_end;    // ending time of the simulation

    int model;  // the model of evolution

    // ~Clone() = default;
    Clone(const Clone& other) = default;
    Clone& operator=(const Clone& other) = default;
    Clone& operator=(Clone&& other) = default;

    Clone(){
        clone_ID = 0;
        parent_ID = -1;
        name = "";
        time_end = 0;
        ntot = 0;
    }


    Clone(int cID, string cname, double time_start){
        clone_ID = cID;
        parent_ID = -1;
        name = cname;
        time_end = time_start;
        ntot = 0;
    }


    Clone(int cID, int pID){
        clone_ID = cID;
        parent_ID = pID;
        name = "";
        time_end = 0;
        ntot = 0;
    }


    ~Clone(){
        // cout << "Only delete available cells " << curr_cells.size() << endl;
        if(cells.size() > 0){
          // when all cells are recorded
          for(auto p : cells){
              // cout << p->cell_ID << "\t" << p->clone_ID << endl;
              delete p;
          }
        }else{
          for(auto p : curr_cells){
              // cout << p->cell_ID << "\t" << p->clone_ID << endl;
              delete p;
          }
        }
    }

    /*********************** functions related to simulating clonal growth **************************/
    void initialize(const Cell_ptr ncell, int verbose = 0){
        // this->cells.clear();
        this->id_curr_cells.clear();
        this->ntot = 1;

        this->id_curr_cells.push_back(ncell->cell_ID);
        this->time_end = ncell->time_occur;
    }


    /*
    Initialize the clone with clonal CNVs if exists
    */
    void initialize_with_dsb(const Cell_ptr ncell, const Model& model, int track_all = 0){
        this->curr_cells.clear();
        this->curr_cells.push_back(ncell);

        if(track_all){
          this->cells.clear();
          this->cells.push_back(ncell);
        }

        this->ntot = 1;
        this->model = model.model_ID;
        this->time_end = ncell->time_occur;
    }


    Cell_ptr get_cell_from_ID(vector<Cell_ptr>& cells, int cID){
        for(int i = 0; i < cells.size(); i++){
            Cell_ptr cell = cells[i];
            if(cell->cell_ID == cID) return cell;
        }
        return NULL;
    }


    int get_n_dsb(){
        int sum=0;
        for (auto cell: curr_cells){
            sum += cell->n_dsb;
        }
        return sum;
    }


    // Find  maximum birth rate (bmax) and maximum death rate (dmax) in the population
    double get_rmax(){
        double bmax = 0;
        double dmax = 0;
        for(unsigned int i = 0; i < curr_cells.size(); i++) {
            Cell_ptr ci = curr_cells[i];
            double bi = ci->birth_rate;
            double di = ci->death_rate;
            if(bi > bmax){
                bmax = bi;
            }
            if(di > dmax){
                dmax = di;
            }
        }
        // cout << "Maximum birth rate: " << bmax << endl;
        // cout << "Maximum death rate: " << dmax << endl
        return bmax + dmax;
    }


    /*
    * Key function when introducing selection into the stochastic branching process
    * Update the birth/death rate of a cell (only occur when there are new mutations in daughter cells) according to that of a reference cell
    * start_cell: the baseline for fitness changes (should be the cell with optimum karyotype)
    * Two models considered here: chr-level, arm_level
    * Use average total copy number for each chr or arm 
    */
    void update_cell_growth(Cell_ptr dcell, const Cell_ptr start_cell, const Model& model, int verbose = 0) {
        // Introduce selection to cells with CNAs using specified fitness
        double gdiff = 0;
        // assume start cell has optimum karyotype by default
        double b0 = start_cell->birth_rate;
        double d0 = start_cell->death_rate;

        if(verbose > 1){
            cout << "update birth-death rate for cell "<< dcell->cell_ID << " with " << dcell->birth_rate << "-" << dcell->death_rate << endl;
        }

        map<int, set<int>> bps_by_chr;
        dcell->g->get_bps_per_chr(bps_by_chr, verbose);
        dcell->g->calculate_segment_cn(bps_by_chr, verbose);

        dcell->get_surv_prob(model.selection_type, verbose);
        // the value may be 0
        if(dcell->surv_prob == 0){
            dcell->fitness = -0.999999;
        }else{
            dcell->fitness = log(dcell->surv_prob * (1 + dcell->surv_prob)) / log(model.psurv_norm * (1 + model.psurv_norm)) - 1;
        }
        
        if(verbose > 1){
            cout << "survival probability: " << dcell->surv_prob  << ", selection coefficient: " << dcell->fitness  << endl;
        }

        switch (model.growth_type) {
            case ONLY_BIRTH: {
                dcell->birth_rate = b0 * (1 + dcell->fitness);
                break;
            }
            case CHANGE_BIRTH: {
                double new_grate = (b0 - d0) * (1 + dcell->fitness);
                dcell->birth_rate = d0 + new_grate;
                dcell->death_rate = dcell->birth_rate - new_grate;
                break;
            }
            case CHANGE_DEATH: { 
                double new_grate = (b0 - d0) * (1 + dcell->fitness);
                double new_drate = (b0 - new_grate);
                dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                dcell->birth_rate = dcell->death_rate + new_grate;
                break;
            }
            case CHANGE_BOTH:{
                double new_grate = (b0 - d0) * (1 + dcell->fitness);
                double rn = runiform(r, 0, 1);
                if(rn < 0.5){
                    dcell->birth_rate = d0 + new_grate;
                    dcell->death_rate = dcell->birth_rate - new_grate;
                }else{               
                    double new_drate = (b0 - new_grate);
                    dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                    dcell->birth_rate = dcell->death_rate + new_grate;
                }
                break;
            }
            default: cout << "" << endl;  
        }  
    }


    /*
       This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm.
       intput:
        Nend -- the cell population size at the end
        ncell -- the starting cell. Given a cell in another clone, it can mimick migragation
        model -- 0: neutral, 1: gradual, 2: punctuated
        restart -- 1: start with a new cell; 0: continue with current state
       output:
        a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
     */
     void grow_with_dsb(const Cell_ptr ncell, const Model& model, int Nend, double frac_unrepaired, int mean_local_frag, double frac_unrepaired_local, double circular_prob, int track_all = 0, int verbose = 0, int restart = 1, double tend = DBL_MAX){
         // Initialize the simulation with one cell
         if(restart == 1) initialize_with_dsb(ncell, model, track_all);

         double t = 0;  // starting time, relative to the time of last end
         int nu = 0; // The number of new mutations, relative to the time of last end

         if(verbose > 0) cout << "\nSimulating tumour growth with CNAs under model " << model.model_ID << " at time " << ncell->time_occur + t << endl;

         while(this->curr_cells.size() < Nend) {
             if (this->curr_cells.size() == 0) {
                 t = 0;
                 nu = 0;
                 initialize_with_dsb(ncell, model, track_all);
                 continue;
             }

             // Choose a random cell from the current population
             int rindex = myrng(this->curr_cells.size());
             Cell_ptr rcell = this->curr_cells[rindex];
             double rbrate = rcell->birth_rate;
             double rdrate = rcell->death_rate;
             int rID = rcell->cell_ID;

             double rmax = get_rmax();
             // increase time
             double tau = -log(runiform(r, 0, 1));  // an exponentially distributed random variable
             double deltaT = tau/(rmax * this->curr_cells.size());
             t += deltaT;

             if(ncell->time_occur + t > tend){
                 t = tend - ncell->time_occur;
                 break;
             }

             // draw a random number
             double rb = runiform(r, 0, rmax);
             // cout << "random number " << rb << endl;
             if(rb < rbrate){
                 if(verbose > 0){
                    cout << "Select cell " << rID << " with " << rcell->n_dsb << " DSBs and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to divide" << endl;
                 }
                 // generate two daughter cells
                 int cID1 = this->ntot + 1;
                 Cell_ptr dcell1 = new Cell(cID1, rID, ncell->time_occur + t);
                 dcell1->copy_parent((*rcell));

                 int cID2 = this->ntot + 2;
                 Cell_ptr dcell2 = new Cell(cID2, rID, ncell->time_occur + t);
                 dcell2->copy_parent((*rcell));

                 // this->n_cycle++;
                 rcell->do_cell_cycle(dcell1, dcell2, frac_unrepaired, n_telo_fusion, n_complex_path, n_path_break, mean_local_frag, frac_unrepaired_local, circular_prob, verbose);

                 this->ntot = this->ntot + 2;

                 if(model.model_ID == 1){   // selection
                     if(verbose > 1){
                         cout << "Update the grow parameters of daughter cells " << endl;
                     }
                     update_cell_growth(dcell1, ncell, model, verbose);
                     update_cell_growth(dcell2, ncell, model, verbose);
                 }
                 // Remove the parent cell from the list of current cells
                 if(rID != 1 && !track_all){
                     delete (this->curr_cells[rindex]);
                     this->curr_cells[rindex] = NULL;
                 }
                 this->curr_cells.erase(this->curr_cells.begin() + rindex);
                 this->curr_cells.push_back(dcell1);
                 this->curr_cells.push_back(dcell2);

                 if(track_all){
                  this->cells.push_back(dcell1);
                  this->cells.push_back(dcell2);
                 }
             }else if(rb >= rbrate && rb < rbrate + rdrate){
                 // cout << " death event" << endl;
                 // if(verbose > 1){
                 //   cout << "Select cell " << rID << " with " << rcell->n_dsb << " DSBs and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to disappear" << endl;
                 // }
                 if(rID != 1 && !track_all){
                     delete (this->curr_cells[rindex]);
                     this->curr_cells[rindex] = NULL;
                 }
                 this->curr_cells.erase(this->curr_cells.begin() + rindex);
             }else{

             }
         }

         this->time_end += t;
         // this->n_novel_dsb  += nu;
     }


     /***************************************************************************************************************************/


    /*********************** functions related to output **************************/

    // treat single cell data as pseudo-bulk data and compute average ploidy 
    // distinguish gain/loss to account for different size distribution (not much help)
    // loc_cn stores the absolute copy numbers for each cell
    // frac_genome_alt = length(which(x!=y)) / length(x)
    double get_pga(const map<int, vector<int>>& loc_cn, int num_loc, int is_haplotype = 0){
        int nsample = loc_cn.size();
        // cout << "There are " << nsample << " samples" << endl;
        vector<int> alter_indicator_sep(num_loc, 0);
        double avg_nalter_sep = 0;
        int baseline = NORM_PLOIDY;
        if(is_haplotype){
            baseline = NORM_PLOIDY / 2;
        }

        for(auto s : loc_cn){
            for(int i = 0; i < num_loc; i++){
                if(s.second[i] != baseline){
                    alter_indicator_sep[i] += 1;
                }
            }
        }

        // exclude clonal regions
        for(int i = 0; i < num_loc; i++){
            if(alter_indicator_sep[i] > 0 && alter_indicator_sep[i] < nsample)
                avg_nalter_sep++;
        }
        // normalized by num_loc to account for different choices of bin size (in real data)
        avg_nalter_sep = (double) avg_nalter_sep / num_loc;

        return avg_nalter_sep;
    }


    // Divergence of CNAs across cells/clones was quantified by computing their pairwise divergence.
    // Specifically, this was the proportion of altered bins (copy number not equal to ploidy in either or both samples) that had different copy number in each sample.
    // ids: The ID of each clone (cell), sorted
    // loc_cn: the absolute copy numbers of each location (bin) along the genome
    pair<double, double> get_pairwise_divergence(const vector<int>& ids, map<int, vector<int>>& loc_cn, int num_loc, int is_haplotype = 0, int use_alter = 1, int verbose = 0){
        int ntotal = 0;
        vector<double> alters;
        double avg_alter = 0.0;
        VARIANCE ppalters;
        int baseline = NORM_PLOIDY;
        if(is_haplotype){
            baseline = NORM_PLOIDY / 2;
        }

        for(int i = 0; i < ids.size(); i++){
            vector<int> lchanges1 = loc_cn[ids[i]];
            for(int j = i + 1; j < ids.size(); j++){
                ntotal++;
                vector<int> lchanges2 = loc_cn[ids[j]];

                int num_alter = 0;
                int num_diff = 0;
                for(int k = 0; k < num_loc; k++){
                    if((lchanges1[k]) != baseline || (lchanges2[k]) != baseline){
                        num_alter++;
                        if((lchanges2[k]) != (lchanges1[k])){
                            num_diff++;
                        }
                    }
                }
                assert(num_alter >= num_diff);
                double prop_alter = 0.0;
                // when only altered regions are considered
                if(use_alter == 1 && num_alter > 0 && num_diff > 0){
                    prop_alter = (double) num_diff / num_alter;
                }else{  // use absolute bin size
                    prop_alter = num_diff;
                }
                alters.push_back(prop_alter);
                ppalters(prop_alter);
                // cout << num_diff << "\t" << num_alter << "\t" << prop_alter << endl;
            }
        }

        for(int i = 0; i < alters.size(); i++){
            avg_alter += alters[i];
            // cout << avg_alter << endl;
        }
        // alters.size() should equal to ntotal
        // cout << alters.size() << endl;

        // weird numbers when running on Mac
        avg_alter = avg_alter / alters.size();
        double var_alter = boost::accumulators::variance(ppalters);

        if(verbose > 1){
            cout << "total pairs of glands for computing pairwise divergence " << ntotal << endl;
        }

        pair<double, double> alt(avg_alter, var_alter);

        return alt;
    }


    /*
       This method prints out the copy numbers of each final cell in a clone
     */
    void print_all_cells(vector<Cell_ptr> cells, ofstream& fout, int verbose = 0){
        int Nend = cells.size();
        if(verbose > 0) cout << "Printing " << Nend << " cells" << endl;

        fout << "cell_ID\tparent_ID\n";
        for(unsigned int i = 0; i < Nend; i++) {
            Cell_ptr cell = cells[i];
            if(cell->parent_ID == 0) continue;
            fout << cell->cell_ID << "\t" << cell->parent_ID << endl;
        }
    }


    /*
    This function prints summary informaton about the simulated clone
     */
    void print_summary(string outfile) {
        ofstream out;
        out.setf(ios::fixed);
        out.setf(ios::showpoint);
        out.precision(9);
        out.open(outfile);

        double lambda = log(2);
        out << "Information for host population:" << endl;
        for(auto cell : curr_cells){
            if(cell->clone_ID == 0){
                cell->print_cell_info();
                if(model == 0) break;     // same rates under neutral evolution
            }
        }
        out << endl;
        // out << "\tNumber of new breaks: "<< n_novel_dsb << endl;
        out << "\tEnd time of simulation: " << time_end << endl;

        out.close();
    }


     // Print out all informaton for one clone
     void print_single_clone(Cell_ptr start_cell, Model model, string outdir, string suffix, int verbose = 1){
         double lambda = start_cell->birth_rate - start_cell->death_rate;
         double tend = log(curr_cells.size())/(lambda); // The latest time that a subclone occurs
         cout << "\nPrinting informaton for clone " << clone_ID << endl;

         cout << "Model of evolution: " << model.model_ID << endl;
         cout << "Initial Net growth rate: " << lambda << endl;
         cout << "Estimated simulation finish time (tumor doublings): " << tend << endl;

         string outfile = "";
         outfile = outdir + "summary" + suffix;
         cout << "Printing summary" << endl;
         print_summary(outfile);
     }


};


#endif
