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

#include "util.hpp"
#include "cell.hpp"

using namespace std;



// Node of a binary lineage tree, used for sampling certain number of nodes at each side
struct node {
    int data;
    node* left;
    node* right;
    int flag; // whether the node is available (dead) or not

    node(int data)
    {
        this->data = data;
        this->flag = 0;
        left = NULL;
        right = NULL;
    }

    node(int data, int flag)
    {
        this->data = data;
        this->flag = flag;
        left = NULL;
        right = NULL;
    }
};


node* search_node(int key, node* p)
{
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
    int genotype_diff;    // whether or not to consider genotype differences
    int growth_type;   // how cells grow under selection
    double fitness;    // the strength of selection
    int use_alpha;    // use the model based on alpha, in which fitness is alpha; Otherwise fitness is classic selection coefficient
    double min_diff = 0;  // The mininal genotype differences that will cause fitness changes
    int norm_by_bin = 1;  // normalized genotype differences by the number of bins

    Model(){
        model_ID = 0;
        genotype_diff = 0;
        growth_type = 0;
        use_alpha = 0;
        fitness = 0;
    }

    Model(int model_ID, int genotype_diff, int growth_type, double fitness){
        this->model_ID = model_ID;
        this->genotype_diff = genotype_diff;
        this->growth_type = growth_type;
        this->fitness = fitness;
    }

    Model(int model_ID, int genotype_diff, int growth_type, double fitness, int use_alpha, int norm_by_bin = 1){
        this->model_ID = model_ID;
        this->genotype_diff = genotype_diff;
        this->growth_type = growth_type;
        this->fitness = fitness;
        this->use_alpha = use_alpha;
        this->norm_by_bin = norm_by_bin;
    }
};


/*
A Clone represents a population of cells
*/
class Clone
{
public:
    int clone_ID;
    string name;

    // Used for tracking clone/deme expansion
    int parent_ID;

    vector<Cell_ptr> curr_cells;   // only available cells at present
    vector<Cell_ptr> cells;   // all cells at present
    vector<int> id_curr_cells;      // Only store the IDs of cell to save space

    // map<int, vector<double>> fits_by_time;

    // vector<string> bps; // breakpoints

    // double dsb_rate;       // assuming constant mutation rate for one clone
    // map<int, vector<int>> cnp_all[NUM_LOC];     // store all CNVs in a big array
    // vector<int> id_by_loc;      // store cell ID sorted by location to facilitate sampling
    // map<int, pair<double, double>> grates;  // birth/death rates for cells with fitness (dis)advantages


    int ntot;   // total number of cells generated so far, used to increase cell ID
    int n_cycle;   // number of cell cycles gone through

    double time_end;    // ending time of the simulation

    int model;  // the model of evolution

    Coord min_pos;
    Coord max_pos;

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
    void initialize_with_dsb(const Cell_ptr ncell, Model model, int verbose = 0){
        this->curr_cells.clear();
        this->curr_cells.push_back(ncell);
        this->ntot = 1;
        this->model = model.model_ID;
        this->time_end = ncell->time_occur;
    }


    /*
    Initialize the clone with clonal CNVs if exists, record all cells in the lineage
    nu: The number of new mutations
    */
    void initialize_with_dsb_cmpl(const Cell_ptr ncell, Model model, int verbose = 0){
        this->curr_cells.clear();
        this->cells.clear();
        this->curr_cells.push_back(ncell);
        this->cells.push_back(ncell);

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


    int get_num_dsb(){
        int sum=0;
        for (auto cell: curr_cells){
            sum += cell->num_dsb;
        }
        return sum;
    }


    // Find  maximum birth rate (bmax) and maximum death rate (dmax) in the population
    double get_rmax(){
        double bmax = 0;
        double dmax = 0;
        for(unsigned int i=0; i<curr_cells.size(); i++) {
            Cell_ptr ci = curr_cells[i];
            double bi = ci->birth_rate;
            double di = ci->death_rate;
            if (bi > bmax) {
                bmax = bi;
            }
            if(di > dmax) {
                dmax = di;
            }
        }
        // cout << "Maximum birth rate: " << bmax << endl;
        // cout << "Maximum death rate: " << dmax << endl
        return bmax + dmax;
    }


    // Update the location of new-born daughter cells
    Coord get_neighbor_position(Coord p_pos, double sd){
        Coord d_pos;

        gsl_vector *mu = gsl_vector_alloc(3);
        gsl_vector *res = gsl_vector_alloc(3);
        gsl_matrix *lower_triangle = gsl_matrix_alloc(3, 3);

        gsl_vector_set(mu, 0, p_pos.x);
        gsl_vector_set(mu, 1, p_pos.y);
        gsl_vector_set(mu, 2, p_pos.z);

        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                if(i == j){
                    gsl_matrix_set(lower_triangle, i, j, sd * sd);
                }else{
                    gsl_matrix_set(lower_triangle, i, j, 0);
                }
            }
        }

        gsl_ran_multivariate_gaussian(r, mu, lower_triangle, res);

        d_pos.x = gsl_vector_get(res, 0);
        d_pos.y = gsl_vector_get(res, 1);
        d_pos.z = gsl_vector_get(res, 2);

        return d_pos;
    }


    /*
    * Key function when introducing selection into the stochastic branching process
    * Update the birth/death rate of a cell (only occur when there are new mutations in daughter cells) according to that of a reference cell
    * TODO: Different fitness values may apply at different locations of the genome
    * start_cell: the baseline for fitness changes (should be the cell with optimum karyotype)
    * Two models considered here:
    * When genotype_diff > 0, selection model is based on alpha (negative selection -- positive alpha)
    * Otherwise, selection model is based on classic selection coefficient (negative selection -- negative s)
    * Genotype differences are considered to impose different selective advantages to different cells
    */
    void update_cell_growth(Cell_ptr dcell, const Cell_ptr start_cell, const Model& model, int loc_type, int verbose = 0) {
        // Introduce selection to cells with CNAs using specified fitness
        double gdiff = 0;
        // assume start cell has optimum karyotype by default
        double b0 = start_cell->birth_rate;
        double d0 = start_cell->death_rate;


        if(verbose > 1){
            cout << "new birth-death rate for cell "<< dcell->cell_ID << " is " << dcell->birth_rate << "-" << dcell->death_rate << " with fitness " << model.fitness << " and genotype difference " << gdiff << endl;
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
     void grow_with_dsb(const Cell_ptr ncell, const Model& model, int n_cycle, int n_unrepaired, double leap_size=0, int verbose = 0, int restart = 1, double tend = DBL_MAX){
         // Initialize the simulation with one cell
         if(restart == 1) initialize_with_dsb(ncell, model, verbose);

         double t = 0;  // starting time, relative to the time of last end
         int mut_ID = ncell->num_dsb;  // used to distinguish different mutations
         int nu = 0; // The number of new mutations, relative to the time of last end
         int num_dsb_event = 0; // count the number of times a CNA event is introduced

         if(verbose > 0) cout << "\nSimulating tumour growth with CNAs under model " << model.model_ID << " at time " << ncell->time_occur + t << endl;

         while(this->n_cycle < n_cycle) {
             if (this->curr_cells.size() == 0) {
                 t = 0;
                 mut_ID = ncell->num_dsb;
                 nu = 0;
                 initialize_with_dsb(ncell, model, verbose);
                 continue;
             }
             // print_all_cells(this->curr_cells, verbose);

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

             if(ncell->time_occur + t > tend && this->curr_cells.size() >= MIN_NCELL){
                 t = tend - ncell->time_occur;
                 break;
             }

             // draw a random number
             double rb = runiform(r, 0, rmax);
             // cout << "random number " << rb << endl;
             if(rb < rbrate){
                 // if(verbose > 1){
                      cout << "Select cell " << rID << " with " << rcell->num_dsb << " mutations and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to divide" << endl;
                 // }
                 // generate two daughter cells
                 int cID1 = this->ntot + 1;
                 Cell_ptr dcell1 = new Cell(cID1, rID, ncell->time_occur + t);
                 dcell1->copy_parent((*rcell));

                 int cID2 = this->ntot + 2;
                 Cell_ptr dcell2 = new Cell(cID2, rID, ncell->time_occur + t);
                 dcell2->copy_parent((*rcell));
                 // dcell2->pos = get_neighbor_position(rcell->pos, POS_SIGMA);

                this->n_cycle++;
                rcell->do_cell_cycle(n_unrepaired, dcell1, dcell2);


                this->ntot = this->ntot + 2;

                 // daughter cells aquire nu new mutations, where nu ~ Poisson(dsb_rate)
                 if(model.fitness != 0)
                 {
                     if(verbose > 1){
                         cout << "Update the grow parameters of daughter cells " << endl;
                     }
                     update_cell_growth(dcell1, ncell, model, verbose);
                     update_cell_growth(dcell2, ncell, model, verbose);
                 }
                 // Remove the parent cell from the list of current cells
                 if(rID != 1){
                     delete (this->curr_cells[rindex]);
                     this->curr_cells[rindex] = NULL;
                 }

                 this->curr_cells.erase(this->curr_cells.begin() + rindex);
                 this->curr_cells.push_back(dcell1);
                 this->curr_cells.push_back(dcell2);
             }
             // death event if b<r<b+d
             else if(rb >= rbrate && rb < rbrate + rdrate) {
                 // cout << " death event" << endl;
                 // if(verbose > 1){
                 //   cout << "Select cell " << rID << " with " << rcell->num_dsb << " mutations and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to disappear" << endl;
                 // }
                 if(rID != 1){
                     delete (this->curr_cells[rindex]);
                     this->curr_cells[rindex] = NULL;
                 }
                 this->curr_cells.erase(this->curr_cells.begin() + rindex);
             }else{

             }
         }

         this->time_end += t;
         // this->num_novel_dsb  += nu;

         if(verbose > 0){
             // cout << "Generated " << nu << " mutations during time " << t << endl;
             if(verbose > 1) print_all_cells(this->curr_cells, verbose);
         }
     }



     /***************************************************************************************************************************/


    /*********************** functions related to output **************************/

    /*
       This method prints out the copy numbers of each final cell in a clone
     */
    void print_all_cells(vector<Cell_ptr> cells, int verbose = 0){
        int num_cell = cells.size();
        if(verbose > 0) cout << "Printing " << num_cell << " cells" << endl;

        for(unsigned int i = 0; i < num_cell; i++) {
                Cell_ptr cell = cells[i];
                cout << cell->cell_ID << "\t" << cell->parent_ID << endl;
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
                lambda = cell->birth_rate - cell->death_rate;
                out << "\tCell ID: " << cell->cell_ID << endl;
                out << "\tMutation rate: " << cell->dsb_rate << endl;

                out << "\tBirth rate: " << cell->birth_rate << endl;
                out << "\tDeath rate: " << cell->death_rate << endl;

                out << "\tEffective break rate (μ/β): " << cell->dsb_rate / ((cell->birth_rate-cell->death_rate)/cell->birth_rate) << endl;
                out << endl;
                if(model == 0) break;     // same rates under neutral evolution
            }
        }
        out << endl;
        // out << "\tNumber of new breaks: "<< num_novel_dsb << endl;
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
         cout << "Mutation rate: " << start_cell->dsb_rate << endl;
         cout << "Estimated simulation finish time (tumor doublings): " << tend << endl;

         string outfile = "";
         outfile = outdir + "summary" + suffix;
         cout << "Printing summary" << endl;
         print_summary(outfile);
     }


};


#endif
