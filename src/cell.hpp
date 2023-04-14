#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>
#include <unordered_map>

#include <assert.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

// #include "util.hpp"
#include "genome.hpp"


using namespace std;

struct bp_interval{
  int bp;
  int chr;
  int haplotype;
  int left_jid;
  int right_jid;
  bool is_inverted;
};


class Mutation{
public:
    int mut_ID;

    int type;   // Mutation type. 0: SNV 1: CNV
    double time_occur;

    int cell_ID;
    int chr;  // chromosome on which the mutation occurs
    int arm;    // 0: no informaton, 1: p; 2: q
    int reciprocal;

    int start;  // start position
    int end;    // end position
    int size;

    int number;  // copy of this mutation

    ~Mutation() = default;
    Mutation(const Mutation& other) = default;
    Mutation(Mutation&& other) = default;
    Mutation& operator=(const Mutation& other) = default;
    Mutation& operator=(Mutation&& other) = default;

    Mutation(){
        mut_ID = 0;
        time_occur = 0;
        number = 1;
    }


    Mutation(int mut_ID, double time_occur){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;
        this->number = 1;
    }

    Mutation(int mut_ID, double time_occur, int chr, int arm, int type, int reciprocal){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;

        this->chr = chr;
        this->arm = arm;
        this->type = type;
        this->reciprocal = reciprocal;

        this->number = 1;
    }

    Mutation(int mut_ID, double time_occur, int chr, int start, int end, int arm, int type, int reciprocal){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;

        this->chr = chr;
        this->start = start;
        this->end = end;
        this->arm = arm;
        this->type = type;
        this->reciprocal = reciprocal;

        this->number = 1;
    }
};

class Cell;

// typedef shared_ptr<Cell> Cell_ptr;
typedef Cell* Cell_ptr;


class Cell{
public:
    int cell_ID;
    int parent_ID;
    int clone_ID;

    // int flag;   // whether the cell is alive or not. 0:new, 1:divided, -1: death
    double time_occur;
    int div_occur;  // occur at which division
    int div_break;
    int only_repair_new;

    // parameters related to cell growth
    double birth_rate;
    double death_rate;

    // chr-level and arm-level CNs for computing fitness values
    vector<double> chr_tcns;
    vector<double> arm_tcns;
    double surv_prob;   // survival probability
    double fitness;   // selection coefficient as in the classic definition

    // parameters related to DSB generation
    int dsb_rate;  // rate of double strand breaks during cell division
    int n_dsb;    // number of double strand breaks during cell division, use exact number for external control over each cell
    int n_unrepaired; // number of unrepaired double strand breaks

    int n_bp;  // each breakpoint may have at most two connections
    int n_del;
    int n_dup;
    int n_h2h;
    int n_t2t;
    int n_tra;
    int n_circle;
    int n_ecdna;  // circular chromosomes without centromeres
    int n_nocentro;  // linear chromosomes without centromeres

    // count events by chromosome
    map<int, set<vector<int>>> chr_bp_unique;  // use position to be consistent with real data
    map<int, vector<int>> chr_type_num;
    map<int, int> chr_n_osc2;
    map<int, int> chr_n_osc3;
    set<string> bp_unique; // the location of unique breakpoints
    set<vector<int>> dbs_unique; // may be fewer than nDSB + nMbreak due to uneven distribution into daughter cells

    genome* g; // each cell has a genome

    Cell(const Cell& other) = default;
    Cell(Cell&& other) = default;
    Cell& operator=(const Cell& other) = default;
    Cell& operator=(Cell&& other) = default;

    ~Cell(){
      delete g;
    }

    Cell(){
        cell_ID = 0;
        parent_ID = 0;
        clone_ID = 0;

        time_occur = 0;
        div_occur = 0; // whether the cell is divided or not

        birth_rate = log(2);
        death_rate = 0;

        dsb_rate = 0;
        n_dsb = 0;
        n_unrepaired = 0;

        n_bp = 0;  // each breakpoint may have at most two connections
        n_del = 0;
        n_dup = 0;
        n_h2h = 0;
        n_t2t = 0;
        n_tra = 0;
        n_circle = 0;  
        n_ecdna = 0;
        n_nocentro = 0;
    }         


    // used when creating daughter cells
    Cell(int cell_ID, int parent_ID, double time_occur){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;

        this->time_occur = time_occur;

        // build an empty genome to be filled based on parent genome
        this->g = new genome();
        this->g->cell_ID = cell_ID;
        // this->div_occur = 0;
        
        // this->birth_rate = birth_rate;
        // this->death_rate = death_rate;
        //
        // this->dsb_rate = 0;
        // this->n_dsb  = 0;
        this->n_bp = 0;  // each breakpoint may have at most two connections
        this->n_del = 0;
        this->n_dup = 0;
        this->n_h2h = 0;
        this->n_t2t = 0;
        this->n_tra = 0;
        this->n_circle = 0; 
        this->n_ecdna = 0;     
        this->n_nocentro = 0;   
    }


    Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, int dsb_rate, int n_dsb, int n_unrepaired, double time_occur, int div_break, int only_repair_new){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;
        this->clone_ID = 0;

        this->time_occur = time_occur;
        this->div_occur = 0;

        this->birth_rate = birth_rate;
        this->death_rate = death_rate;

        this->dsb_rate = dsb_rate;
        this->n_dsb  = n_dsb;
        this->n_unrepaired = n_unrepaired;
        this->div_break = div_break;
        this->only_repair_new = 0;

        this->g = new genome(cell_ID);

        this->n_bp = 0;  // each breakpoint may have at most two connections
        this->n_del = 0;
        this->n_dup = 0;
        this->n_h2h = 0;
        this->n_t2t = 0;
        this->n_tra = 0;
        this->n_circle = 0; 
        this->n_ecdna = 0; 
        this->n_nocentro = 0;              
    }


    Cell_ptr get_parent(vector<Cell_ptr> cells){
        for(int i = 0; i < cells.size(); i++){
            Cell_ptr cell = cells[i];
            if(cell->cell_ID == parent_ID) return cell;
        }
        return NULL;
    }


    // initialization after daughter cells are generated
    void copy_parent(const Cell& ncell){
        this->clone_ID = ncell.clone_ID;

        // this->time_occur = ncell.time_occur;
        this->div_occur = ncell.div_occur + 1;
        this->div_break = ncell.div_break;

        this->birth_rate = ncell.birth_rate;
        this->death_rate = ncell.death_rate;

        this->dsb_rate = ncell.dsb_rate;
        this->n_dsb = ncell.n_dsb;
        this->n_unrepaired = ncell.n_unrepaired;

        // genome information will be recomputed later
    }


    /*********************** functions related to mutation generations **************************/
    // pass a path and related objects to both daughter cells, assuming correct implicit copy, decrpetated after explicit copy
    void inherit_path_both(path* p, Cell_ptr dcell1, Cell_ptr dcell2){
      path* p1 = new path(*p);
      p1->cell_ID = dcell1->cell_ID;
      dcell1->g->paths[p1->id] = p1;

      path* p2 = new path(*p);
      p2->cell_ID = dcell2->cell_ID;
      dcell2->g->paths[p2->id] = p2;

      // copy related breakpoints
      for(auto jid : p->nodes){
        breakpoint* j = g->breakpoints[jid];

        breakpoint* j1 = new breakpoint(*j);
        j1->cell_ID = dcell1->cell_ID;
        dcell1->g->breakpoints[jid] = j1;

        breakpoint* j2 = new breakpoint(*j);
        j2->cell_ID = dcell2->cell_ID;
        dcell2->g->breakpoints[jid] = j2;
      }

      // copy related adjacencies
      for(auto aid : p->edges){
        adjacency* a = g->adjacencies[aid];

        adjacency* a1 = new adjacency(*a);
        a1->cell_ID = dcell1->cell_ID;
        dcell1->g->adjacencies[aid] = a1;

        adjacency* a2 = new adjacency(*a);
        a2->cell_ID = dcell2->cell_ID;
        dcell2->g->adjacencies[aid] = a2;
      }
    }


    /*
    * Copy a path and associated breakpoints and adjacencies to a daughter cell
    * Paths in the parent cell will be released so need to create a new path for the daughter cell
    * Assume path p is within g->paths
    */
    void inherit_path_one(path* p, Cell_ptr dcell1, int verbose = 0){
      // keep track of which daughter cell the path is distributed to get balanced distribution
      // https://stackoverflow.com/questions/63934567/vector-of-pointers-why-does-changing-the-pointer-externally-not-change-the-vect
      p->child_cell_ID = dcell1->cell_ID;
      g->paths[p->id] = p; 
      if(verbose > 1){
        cout << "inherit path " << p->id + 1 << " to cell " << dcell1->cell_ID << "\t" << g->paths[p->id]->child_cell_ID << endl;
        p->print();
      }

      // keep same path ID in the daughter cell
      path* p1 = new path(*p);  
      p1->child_cell_ID = -1;   
      p1->cell_ID = dcell1->cell_ID;
      dcell1->g->paths[p1->id] = p1;

      // copy related breakpoints
      for(auto jid : p->nodes){
        breakpoint* j = g->breakpoints[jid];
        breakpoint* j1 = new breakpoint(*j);
        j1->cell_ID = dcell1->cell_ID;
        dcell1->g->check_duplicated_breakpoint(j1);
        dcell1->g->breakpoints[j1->id] = j1;

        if(verbose > 2){
          cout << "copied breakpoint " << j1->id << endl;
          j1->print();
        }        
      }

      // copy related adjacencies
      for(auto aid : p->edges){
        adjacency* a = g->adjacencies[aid];
        adjacency* a1 = new adjacency(*a);
        a1->cell_ID = dcell1->cell_ID;
        dcell1->g->check_duplicated_adjacency(a1);
        dcell1->g->adjacencies[a1->id] = a1;

        if(verbose > 2){
          cout << "copied adjacency " << a1->id << endl;
          a1->print();
        }
      }

      if(verbose > 1){
        cout << "path " << p->id + 1 << " inherited to cell " << dcell1->cell_ID << "\t" << g->paths[p->id]->child_cell_ID << endl;
        p->print();
        p1->print();
      }
    }


    // distribute a path to either daughter cell randomly
    void distribute_path_random(path* p, Cell_ptr dcell1, Cell_ptr dcell2, int verbose = 0){
        string desc = "";
        if(p->n_centromere == 0){
          desc = " (without centromere) ";
        }else{
          desc = " (with centromere) ";
        }

        double u = runiform(r, 0, 1);
        if(u < 0.5){
          if(verbose > 0){
            cout << "Distributing copied path " << p->id + 1 << desc << "to daughter cell " << dcell1->cell_ID << endl;
            g->write_path(p, cout);
          }
          inherit_path_one(p, dcell1, verbose);
        }else{
          if(verbose > 0){
            cout << "Distributing copied path " << p->id + 1 << desc << "to daughter cell " << dcell2->cell_ID << endl;
            g->write_path(p, cout);
          }
          inherit_path_one(p, dcell2, verbose);
        }
        if(verbose > 1) cout << "\nFinish random path distributing\n" << endl;
    }


    void break_circle(path* p2split, int verbose = 0){
        int nnode = p2split->nodes.size();
        breakpoint* js = g->breakpoints[p2split->nodes[0]];
        breakpoint* je = g->breakpoints[p2split->nodes[nnode - 1]];
        int aid = p2split->edges.back();
        // the novel connection causing a circle
        adjacency* al = g->adjacencies[aid];
        int js_left = js->left_jid;
        int js_right = js->right_jid;
        int al1 = al->junc_id1;
        int al2 = al->junc_id2;

        if(js_left == al1 || js_left == al2){
          js->left_jid = -1;
          if(je->right_jid == js->id){
            je->right_jid = -1;
          }else{
            assert(je->left_jid == js->id);
            je->left_jid = -1;
          }         
        }else{
          assert(js_right == al1 || js_right == al2);
          js->right_jid = -1;
          if(je->right_jid == js->id){
            je->right_jid = -1;
          }else{
            assert(je->left_jid == js->id);
            je->left_jid = -1;
          }
        }

        g->adjacencies.erase(aid);
        // p2split->nodes.pop_back();
        p2split->edges.pop_back();

        if(verbose > 1){
          cout << "remove last connection " << aid << " with left breakpoint " << al1 << " and right breakpoint " << al2 << " as the path is circular" << endl;
          p2split->print();
          js->print();
          je->print();
        }
    }


    // Find two intervals containing each centromere
    void find_centromeric_intervals(path* p2split, vector<pair<int, int>>& adjID_pair_withcent, int verbose = 0){
      int ncent = 0;   // count #centromeres traversed
      int prev_adj = p2split->edges[0];
      adjacency* adj = g->adjacencies[prev_adj];
      if(adj->is_centromeric){
        ncent++;
      }

      for(int i = 1; i < p2split->edges.size(); i++){
        int curr_adj = p2split->edges[i];
        adj = g->adjacencies[curr_adj];
        if(adj->is_centromeric){  // can be true only for intervals
          // cout << "is_centromeric" << endl;
          ncent++;
          if(ncent == 2){
            adjID_pair_withcent.push_back(pair<int, int>(prev_adj, curr_adj));
            ncent = 1;  // restart from current interval
          }
          prev_adj = curr_adj;
        }
      }

      if(verbose > 1){
        cout << "Adjacent intervals with " << adjID_pair_withcent.size() + 1 << " centromeres: " << endl;
        for(auto pair : adjID_pair_withcent){
          cout << pair.first << "\t" << pair.second << endl;
          g->adjacencies[pair.first]->print();
          g->adjacencies[pair.second]->print();
        }
      }
    }


    void get_left_interval(adjacency* adj, int& ps, int& pe, int& chr, int verbose = 0){
        int jid1 = adj->junc_id1;
        int jid2 = adj->junc_id2;
        assert(g->breakpoints[jid1]->chr == g->breakpoints[jid2]->chr);

        if(verbose > 1){
          cout << "add break on left side of second interval " << adj->id << endl;
          adj->print();
          g->breakpoints[jid1]->print();
          g->breakpoints[jid2]->print();
        }

        chr = g->breakpoints[jid1]->chr;
        int pos1 = g->breakpoints[jid1]->pos;
        int pos2 = g->breakpoints[jid2]->pos;
        if(adj->is_inverted){
          if(pos1 < pos2){
            ps = CENT_ENDS[chr];
            pe = pos2;
          }else{
            ps = pos2;
            pe = CENT_STARTS[chr];
          }
        }else{
          if(pos1 < pos2){
            ps = pos1;
            pe = CENT_STARTS[chr];
          }else{
            ps = CENT_ENDS[chr];
            pe = pos1;
          }
        }
    }


    void get_right_interval(adjacency* adj, int& ps, int& pe, int& chr, int verbose = 0){
        int jid1 = adj->junc_id1;
        int jid2 = adj->junc_id2;
        assert(g->breakpoints[jid1]->chr == g->breakpoints[jid2]->chr);

        if(verbose > 1){
          cout << "add break on right side of first interval " << adj->id << endl;
          adj->print();
          g->breakpoints[jid1]->print();
          g->breakpoints[jid2]->print();
        }

        chr = g->breakpoints[jid1]->chr;
        int pos1 = g->breakpoints[jid1]->pos;
        int pos2 = g->breakpoints[jid2]->pos;
        if(adj->is_inverted){
          if(pos1 < pos2){
            ps = pos1;
            pe = CENT_STARTS[chr];
          }else{
            ps = CENT_ENDS[chr];
            pe = pos1;
          }
        }else{
          if(pos1 < pos2){
            ps = CENT_ENDS[chr];
            pe = pos2;
          }else{    // inverted region
            ps = pos2;
            pe = CENT_STARTS[chr];
          }
        }      
    }

    // randombly choose an interval to introduce the break
    // when breaking paths with multiple centromeres, a breakpoint must be to the other side of another breakpoint
    // For each adjacent interval, breakpoint is either at right of 1st interval or left of 2nd interval when finding all breakpoints before splitting intervals
    // Here, left and right should be based on the direction of the interval on the path, so need to check inversion
    void find_breakpoint_in_interval(const vector<pair<int, int>>& adjID_pair_withcent, vector<bp_interval>& bps_in_interval, int verbose = 0){
      int chr = -1;
      int haplotype = -1;
      int left_jid = -1;
      int right_jid = -1;

      // vector<int> old_aids;
      for(auto pair : adjID_pair_withcent){
        int sel = myrng(2);  // 0 or 1
        adjacency* adj = g->adjacencies[pair.first];
        assert(adj->is_centromeric);
        assert(g->breakpoints[adj->junc_id1]->chr == g->breakpoints[adj->junc_id2]->chr);
        assert(g->breakpoints[adj->junc_id1]->haplotype == g->breakpoints[adj->junc_id2]->haplotype);
        int bp, ps, pe;        

        if(sel == 0){  // break on left side of second interval
          adj = g->adjacencies[pair.second];   // should not be ptel
          get_left_interval(adj, ps, pe, chr, verbose);
        }else{ // break on right side of first interval
          // should not be qtel
          get_right_interval(adj, ps, pe, chr, verbose);
        }

        if(verbose > 1) cout << ps << "\t" << pe << endl;
        // assert(ps < pe);     // it may not be possible to introduce a breakpoint without breaking centromere
        if(ps >= pe){
          if(sel == 0){
            adj = g->adjacencies[pair.first]; 
            get_right_interval(adj, ps, pe, chr, verbose);
          }else{
            adj = g->adjacencies[pair.second];   // should not be ptel
            get_left_interval(adj, ps, pe, chr, verbose);
          }
        }

        bp = (int)runiform(r, ps, pe);
        left_jid = adj->junc_id1;
        right_jid = adj->junc_id2;
        int haplotype = g->breakpoints[adj->junc_id1]->haplotype;
        bp_interval bi{bp, chr, haplotype, left_jid, right_jid, adj->is_inverted};
        bps_in_interval.push_back(bi);
        // old_aids.push_back(adj->id);
      }
    }


    void set_path_breakpoints(const vector<bp_interval>& bps_in_interval, vector<int>& junc_starts, int verbose){
      int jid = g->breakpoints.rbegin()->first + 1;
      int aid = g->adjacencies.rbegin()->first + 1;
      // some interval may be chosen twice and be broken at 1st choice
      int prev_left_jid = -1;
      // int prev_right_jid = -1;
      int i = 0;

      for(auto bi : bps_in_interval){
        int bp = bi.bp;
        int chr = bi.chr;
        int haplotype = bi.haplotype;
        // neighbours from bi are before any new breaks and may be changed if a break occur in this interval
        int left_jid = bi.left_jid;
        int right_jid = bi.right_jid;
        // cout << "cell ID before adding breakpoints " << this->cell_ID << endl;
        // cout << g->breakpoints.size() << " breakpoints " << endl;
        // update left_jid if an interval is chosen twice
        if(prev_left_jid == left_jid){
          if(bi.is_inverted){
            right_jid = jid - 2;
          }else{
            left_jid = jid - 1;
          }
        }

        if(verbose > 0){
          cout << "\tpath break " << i + 1 << " at position " << bp << " chr " << chr + 1 << " haplotype " << get_haplotype_string (haplotype) << endl;
          if(verbose > 1) cout << "\twith left breakpoint " << left_jid << " and right breakpoint " << right_jid << endl;
        }

        i++;
        g->add_new_breakpoint(jid, aid, chr, bp, haplotype, left_jid, right_jid, false, verbose);
        // one breakpoint is end of a path, the other is the start of another path, depending on the orientations
        // jid has increased by 2 after 
        if(bi.is_inverted){
            if(verbose > 0){
              cout << "\tInverted breakpoint interval " << endl;
            }
            junc_starts.push_back(jid - 2);
        }else{
            junc_starts.push_back(jid - 1); // only add 2nd breakpoint to avoid repeated tranverse
        }

        // update number of DSBs per chrom per haplotype
        // int idx = chr + haplotype * NUM_CHR;
        // g->vec_n_mitosis_break[idx] += 1;
        prev_left_jid = left_jid;
        // prev_right_jid = right_jid;
      }
    }


    // construct a new path starting from a specific breakpoint
    void get_path_from_bp(path* p, breakpoint* js, double frac_unrepaired, double circular_prob, int verbose = 0) {
        js->path_ID = p->id;
        p->nodes.push_back(js->id);
        if(verbose > 1){
          cout << "\nstart connecting path from " << js->id << endl;
          js->print();
        }
        // there will be no cycles
        g->get_connected_breakpoint(js, p, g->adjacencies, false, verbose);

        if(verbose > 1) {
          cout << "setting circular path type at probability " << circular_prob << endl;
        }          
        double u2 = runiform(r, 0, 1);
        if(u2 < circular_prob){
            g->set_path_type(p, true, verbose);
        }else{
            g->set_path_type(p, false, verbose);
        } 

        assert(p->n_centromere <= 1);

        if(verbose > 1){
          cout << "\nconstructed a new path " << p->id + 1 << " at breakpoint " << js->id << endl;
          js->print();
        }
    }

 
    // new_bps: new breakpoints introduced during fragmentation after repairing
    // for simplification, assume all breakpoints are not repaired in mitosis
    bool path_local_fragmentation(path* p, int n_local_frag, vector<breakpoint*>& new_bps, int verbose = 0){
        int nbreak = gsl_ran_poisson(r, n_local_frag);   // expected number of breaks, more reasonable with different number of breaks
        // int nbreak = n_local_frag;
        int nbreak_succ = 0; // real number of breaks during fragmentation, as some breaks may be infeasible due to constraints
        if(verbose > 1){
          cout << "Fragmentize path " << p->id + 1 << " in cell " << cell_ID << endl;
          p->print();
        }
        // generate DSBs on the selected path, form nbreak+1 new paths
        int bcount = 0;  
        int ntrial = 0;
        int nfail = 0;
        vector<breakpoint*> local_bps;
        int jid = g->breakpoints.rbegin()->first + 1;
        int aid = g->adjacencies.rbegin()->first + 1;
        while(bcount != nbreak){
          // randomly pick an edge
          int eid = myrng(p->edges.size());
          adjacency* e = g->adjacencies[p->edges[eid]];
          bool is_inverted = e->is_inverted;
          if(verbose > 1){
            cout << "selected edge " << eid << endl;
            e->print();
          }

          int p1 = g->breakpoints[e->junc_id1]->pos;
          int p2 = g->breakpoints[e->junc_id2]->pos;
          int ps = p1;
          int pe = p2;
          if(p1 > p2){
            ps = p2;
            pe = p1;
          }
          if(e->type == INTERVAL){
            // double u2 = runiform(r, 0, 1);
            // if(u2 < 0.5){  // break this edge
                if(verbose > 1){
                  cout << "break on edge " << e->id << endl;
                  g->breakpoints[e->junc_id1]->print();
                  g->breakpoints[e->junc_id2]->print();
                  e->print();
                }

                bcount++;

                int chr = g->breakpoints[e->junc_id1]->chr;
                int haplotype = g->breakpoints[e->junc_id1]->haplotype;
                assert(g->breakpoints[e->junc_id1]->chr == g->breakpoints[e->junc_id2]->chr);
                assert(g->breakpoints[e->junc_id1]->haplotype == g->breakpoints[e->junc_id2]->haplotype);
                int bp = 0;
                int left_jid = -1;
                int right_jid = -1;

                do{
                  ntrial++;
                  bp = (int)runiform(r, ps, pe);
                  left_jid = e->junc_id1;
                  right_jid = e->junc_id2;
                }while((bp <= TELO_ENDS1[chr] || bp >= TELO_ENDS2[chr] || (bp >= CENT_STARTS[chr] && bp <= CENT_ENDS[chr])) && ntrial <= nbreak);

                // it may be impossible to split the region without breaking telomore or centromere
                if(ntrial > nbreak){
                  nfail++;
                  continue;
                }

                nbreak_succ++;

                // remove old edges, add new edges
                g->add_new_breakpoint(jid, aid, chr, bp, haplotype, left_jid, right_jid, false, verbose);

                // record both breakpoints for repairing
                local_bps.push_back(g->breakpoints[jid - 2]);
                local_bps.push_back(g->breakpoints[jid - 1]);

                // always add the right (based on the direction of the edge) most breakpoint to avoid duplicated path
                if(is_inverted){
                  breakpoint* bp1 = g->breakpoints[jid - 2];
                  assert(bp1->id == jid - 2);
                  new_bps.push_back(bp1);
                  if(verbose > 1){
                    cout << "Inverted breakpoint interval " << endl;
                    cout << "storing new breakpoints with current jid " << jid << endl;
                    bp1->print();
                  }
                }else{
                  breakpoint* bp2 = g->breakpoints[jid - 1];
                  assert(bp2->id == jid - 1);
                  new_bps.push_back(bp2);
                  if(verbose > 1){
                    cout << "storing new breakpoints with current jid " << jid << endl;
                    bp2->print();
                  }
                }

                if(verbose > 1){
                  cout << "Updating path with new edges " << endl;
                }

                // no need to add new breakpoints, no need to keep edge order
                // what matters is the new segment where new breaks can be introduced
                p->edges.erase(p->edges.begin() + eid);
                p->edges.push_back(aid - 2);
                p->edges.push_back(aid - 1);
            // }
          }
        }
        if(nfail == nbreak){  // all trials to add breakpoints failed
          return false;
        }else{
          return true;    
        }      
      }


    void fragment_path(path* p, int& pid, int n_local_frag, double frac_unrepaired, double circular_prob, Cell_ptr dcell1, Cell_ptr dcell2, int verbose = 0){
        // p will be broken into new paths
        vector<breakpoint*> new_bps;
        bool splited = path_local_fragmentation(p, n_local_frag, new_bps, verbose);
        if(!splited){
          g->validate_path(p);
          g->paths[p->id] = p;
          distribute_path_random(p, dcell1, dcell2, verbose);
          return;
        }
        new_bps.push_back(g->breakpoints[p->nodes[0]]);

        int i = 0;
        for(auto njs: new_bps){ 
          if(verbose > 1){
            cout << "\nGenerating " << ++i << "th new path out of local fragmentation from breakpoint " << njs->id << endl;
            njs->print();
          }  

          path* p2 = new path(++pid, this->cell_ID, COMPLETE);              
          get_path_from_bp(p2, njs, frac_unrepaired, circular_prob, verbose);

          if(verbose > 1){
            cout << "Path " << p2->id + 1 << " from breakpoint " << njs->id << endl;               
            p2->print();
          }
          g->validate_path(p2);
          g->paths[p2->id] = p2;
          distribute_path_random(p2, dcell1, dcell2, verbose);
        }
    }


    // called when there are >2 centromeres
    void split_path_from_bp(int& pid, breakpoint* js, int n_local_frag, double frac_unrepaired, double circular_prob, Cell_ptr dcell1, Cell_ptr dcell2, int verbose = 0){
        //  each path will have a unique ID
        path* p = new path(++pid, this->cell_ID, COMPLETE);
        get_path_from_bp(p, js, frac_unrepaired, circular_prob, verbose);

        if(n_local_frag == 0){
          if(verbose > 1){
            cout << "\nThe new path has one centromere after simple break" << endl;
          }
          g->validate_path(p);
          g->paths[p->id] = p;
          distribute_path_random(p, dcell1, dcell2, verbose);
        }else{
          // randomly fragmentatize a path
          double u = runiform(r, 0, 1);
          if(u < 0.5){
              if(verbose > 0){
                cout << "\nThe new path has local fragmentation" << endl;
              }
              fragment_path(p, pid, n_local_frag, frac_unrepaired, circular_prob, dcell1, dcell2, verbose);
          }else{
              g->validate_path(p);
              g->paths[p->id] = p;
              distribute_path_random(p, dcell1, dcell2, verbose);
          }
        }
    }

    // split the paths to two daughter cells respectively when there are two paths/centromeres 
    void split_new_path_2centro(path* p, int& cell_ID1, int& cell_ID2, Cell_ptr dcell1, Cell_ptr dcell2, int verbose = 0){
        if(cell_ID1 == -1){       
          double u = runiform(r, 0, 1); 
          if(u < 0.5){     
            inherit_path_one(p, dcell1, verbose);             
            cell_ID1 = dcell1->cell_ID;
          }else{
            inherit_path_one(p, dcell2, verbose);
            cell_ID1 = dcell2->cell_ID;
          }
          if(verbose > 0) cout << "path " << p->id + 1 << " ditributed to daughter cell " << cell_ID1 << endl;
        }else{
          if(cell_ID1 == dcell1->cell_ID){     
            inherit_path_one(p, dcell2, verbose);
            cell_ID2 = dcell2->cell_ID;
          }else{
            inherit_path_one(p, dcell1, verbose);
            cell_ID2 = dcell1->cell_ID;
          }
          if(verbose > 0) cout << "path " << p->id + 1 << " ditributed to daughter cell " << cell_ID2 << endl;              
        }
    }


    // called when number of centromeres in a path exceeds 1 (may be caused by random joining in G1)
    // choose breakpoint location so that each path contains one centromere, introduce new breakpoints and form connections
    // assume the path has 2 telomeres at the ends
    // p2split will be split into n_centromere paths, which are distributed randomly
    void segregate_polycentric(path* p2split, Cell_ptr dcell1, Cell_ptr dcell2, int& pid, int n_local_frag, double frac_unrepaired, double circular_prob, int verbose = 0){
      // verbose = 2;
      if(verbose > 0){
        cout << "\nbreaking a complex path with >1 centromeres" << endl;
        g->write_path(p2split, cout);
      }
      if(verbose > 1){
        p2split->print();
      }

      // remove last connection if the path is circular
      if(p2split->is_circle){
        if(verbose > 1) cout << "\nbreaking a circular path" << endl;
        break_circle(p2split, verbose);
      }

      int nbreak = p2split->n_centromere - 1;

      vector<pair<int, int>> adjID_pair_withcent;
      if(verbose > 1) cout << "\nfinding intervals with centromeres" << endl;
      find_centromeric_intervals(p2split, adjID_pair_withcent, verbose);
      assert(adjID_pair_withcent.size() == nbreak);

      vector<bp_interval> bps_in_interval;   // breakpoints at each interval
      if(verbose > 1) cout << "\nselecting breakpoints in the centromeric intervals" << endl;
      find_breakpoint_in_interval(adjID_pair_withcent, bps_in_interval, verbose);
      assert(bps_in_interval.size() == nbreak);

      vector<int> junc_starts;   
      if(verbose > 1) cout << "\nsetting start breakpoints for new paths" << endl;
      set_path_breakpoints(bps_in_interval, junc_starts, verbose);
      // The path starting from both path breaks will be distributed twice when adding both sides of a break
      // assert(junc_starts.size() == nbreak * 2);
      assert(junc_starts.size() == nbreak);  // one junc_start lead to a path

      // add path start and end to junc_starts to avoid single traversation of terminal intervals when the original path is not dupliated
      // breakpoint* ps = g->breakpoints[p2split->nodes[0]];
      // breakpoint* pe = g->breakpoints[p2split->nodes[p2split->nodes.size() - 1]];
      junc_starts.push_back(p2split->nodes[0]);
      // junc_starts.push_back(pe);
      // split paths at once to save cost
      // no need to rejoin paths after splitting, just distribute to different daughter cells
      if(verbose > 1){
        cout << "\nStart breakpoints of split paths: " << endl;
        for(auto j : junc_starts){
          g->breakpoints[j]->print();
        }
        cout << "\n#paths before splitting: " << g->paths.size() << endl;
      }

      if(verbose > 0) cout << "\nSegregating into at least " << junc_starts.size() << " subpaths for path " << p2split->id + 1 << endl;

      // recompute path from each breakpoint to avoid complexity of tracking each broken adjacency
      // each js leads a new path with only one centromere
      if(junc_starts.size() == 2){
        int cell_ID1 = -1; // record the cell ID inheriting 1st new path
        int cell_ID2 = -1;            
        if(n_local_frag == 0){
          if(verbose > 0) cout << "No local fragmentation\n";
          for(auto js : junc_starts){
            path* p = new path(++pid, this->cell_ID, COMPLETE);
            get_path_from_bp(p, g->breakpoints[js], frac_unrepaired, circular_prob, verbose);
            g->validate_path(p);
            g->paths[p->id] = p;          
            split_new_path_2centro(p, cell_ID1, cell_ID2, dcell1, dcell2, verbose);   
          }   
        }else{      
          // randomly fragmentatize a path
          for(auto js : junc_starts){
            path* p = new path(++pid, this->cell_ID, COMPLETE);
            get_path_from_bp(p, g->breakpoints[js], frac_unrepaired, circular_prob, verbose);

            double u = runiform(r, 0, 1);
            if(u < 0.5){
              if(verbose > 0){
                cout << "\nThe new path has local fragmentation with probability 0.5" << endl;
              }
              fragment_path(p, pid, n_local_frag, frac_unrepaired, circular_prob, dcell1, dcell2, verbose);
            }else{
              if(verbose > 0){
                cout << "\nThe new path has no local fragmentation with probability 0.5" << endl;
              }             
              g->validate_path(p);
              g->paths[p->id] = p;          
              split_new_path_2centro(p, cell_ID1, cell_ID2, dcell1, dcell2, verbose); 
            }  
          }        
        }
      }else{
        for(auto jid : junc_starts){
          if(verbose > 0){
            cout << "\nspliting path starting from breakpoint " << jid << endl;
            g->breakpoints[jid]->print();
          }
          split_path_from_bp(pid, g->breakpoints[jid], n_local_frag, frac_unrepaired, circular_prob, dcell1, dcell2, verbose);
        }
      }
    }


    // introduce DSB and repair breaks (adding variant adjacencies)
    void g1(vector<pos_bp>& bps, vector<double>& bp_fracs, double frac_unrepaired, double circular_prob, int pair_type = 0, double prob_correct_repaired = 0, int verbose = 0){
      // verbose = 1;
      if(verbose > 0) cout << "\nIntroducing " << n_dsb << " DSBs" << endl;
      vector<breakpoint*> junc2repair;

      g->generate_dsb(n_dsb, bps, bp_fracs, junc2repair, verbose);
      
      // connect segments randomly to get variant adjacencies (repair DSBs)
      // only breakpoints with missing connections need to be repaired
      if(!only_repair_new){
        junc2repair.clear();
        for(auto jm : g->breakpoints){
          breakpoint* j = jm.second;
          if(!j->is_repaired){   // breakpoints from previous cycles may not be repaired
            junc2repair.push_back(j);
            if(verbose > 0) j->print();
          }
        }
      }

      if(verbose > 0) cout << junc2repair.size() << " breakpoints to repair" << endl;

      // assume two junction points are related to one DSB
      // each repairing connects two breakpoints, so not all breakpoints will be repaired if there is an odd number of breakpoints
      n_unrepaired = round(junc2repair.size() / 2 * frac_unrepaired);
      int n_torepair = junc2repair.size() / 2 - n_unrepaired;
      if(verbose > 0) cout << "\nRepairing " << n_torepair << " DSBs" << endl;      
      if(n_torepair > 0){
        g->repair_dsb(n_unrepaired, junc2repair, pair_type, prob_correct_repaired, verbose);
      }

      if(verbose > 1){
        cout << "\nRemaing breakpoints to repair: " << endl;
        for(auto j: junc2repair){
            j->print();
        }
        cout << "\nConnecting breakpoints to paths" << endl;
      }

      if(verbose > 0) cout << "\nGetting the derivative genome" << endl;
      g->get_derivative_genome(circular_prob, verbose);

      if(verbose > 1){
        // all regions should be included with CN 2, may have different orders
        cout << "\nAll paths in the derivative genome" << endl;
        for(auto p : g->paths){
          p.second->print();
          g->write_path(p.second, cout);
        }
        // TODO: add CN sanity
        cout << "\nEnd of G1 phase" << endl;
      }
    }


    // S phase and G2 implemented together
    void sphase_g2(int& n_telo_fusion, int verbose = 0){
      // verbose = 1;
      int max_pID = g->find_max_pathID();
      map<int, path*> orig_paths = this->g->paths;  // new paths added when duplicating paths
      for(auto pp : orig_paths){
        path* p = pp.second;
        // replicates incomplete paths by adding new breakpoints (S) and path fusion
        if(p->type != COMPLETE && !p->is_circle){
          // if(p->type == NONTEL && p->n_centromere > 0){
            // biologically feasible to have such event ?
            // continue;   // skip these paths to avoid circular DNA with centromere
          // }
          if(verbose > 0){
            cout << "Duplicating incomplete path " << p->id + 1 << endl;
            g->write_path(p, cout);
          }

          // TODO: distinguish fusion at which cell cycle to see its linkage with chromothripsis?
          // an complete path must lost at least one telomere
          // if(p->type == PTEL || p->type == QTEL){
            n_telo_fusion += 1;
          // }

          if(verbose > 0){
            cout << "fusion of path " << p->id + 1 << endl;
            p->print();
            g->write_path(p, cout);
          }

          // keep the same path ID, as it is extended with its copy
          g->duplicate_path_fusion(*p, verbose);

          if(verbose > 1){
            cout << "after path fusion" << endl;
            p->print();
          }
        }else{ // make a copy to avoid mixing with paths split from incomplete fusions in mitosis (a complete path may have >1 centromeres)
            vector<int> edges_copy;
            vector<int> nodes_copy;
            g->duplicate_path(*p, edges_copy, nodes_copy, false, verbose);

            // ID of new paths starting from maximum value of previous paths
            path* p2 = new path(++max_pID, this->cell_ID, p->type);
            for(auto n : nodes_copy){
              g->breakpoints[n]->path_ID = p2->id;
            }
            for(auto e : edges_copy){
              g->adjacencies[e]->path_ID = p2->id;
            }
            p2->edges = edges_copy;
            p2->nodes = nodes_copy;
            p2->n_centromere = p->n_centromere;
            p2->is_circle = p->is_circle;
            p2->sibling = p->id;
            p->sibling = p2->id;

            g->validate_path(p2);
            g->paths[p2->id] = p2;

            if(verbose > 1){
              cout << "copy path " << p->id + 1 << " to path " << p2->id + 1 << endl;
              p->print();
              g->write_path(p, cout);
              p2->print();
              g->write_path(p2, cout);
            }
        }
      }
    }


    // randomly distribute between daughter cells, depend on #centromeres
    void mitosis(Cell_ptr dcell1, Cell_ptr dcell2, int& n_complex_path, int& n_path_break, int n_local_frag, double frac_unrepaired, double circular_prob, int verbose = 0){
      // verbose = 2;
      if(verbose > 1){
        cout << g->paths.size() << " paths before split (may include duplicated path distributed to daughter cells)" << endl;
      }

      // get the max Path ID to define new IDs for random distribution
      int max_pID = g->find_max_pathID();
      map<int, path*> orig_paths = this->g->paths;  // new paths added when splitting complex paths
      for(auto pp : orig_paths){
        path* p = pp.second;
        assert(pp.first == p->id);
        if(verbose > 1){
          cout << "\ndistributing path " << p->id + 1 << endl;
          g->write_path(p, cout);
        }
        if(verbose > 1) p->print();
        // a circular path may have even number of nodes during path replacation
        // assert(p->nodes.size() % 2 == 0);
        if((p->nodes.size() % 2 != 0)){
          cout << "path node size is not even!" << endl;
          p->print();
          exit(FAIL);
        }
        if(p->n_centromere == 1){  // balanced distribution, may be circular
          // keep the sibling path, distribution to different cells
          assert(p->sibling != -1); 
          if(verbose > 0){
            cout << "\nbalanced distribution of path " << p->id + 1 << " with one centromere in cell " << cell_ID << endl;       
            cout << "distributing the same copy to different cells" << endl;
            // cout << "sibling path " << p->sibling + 1 << " ditributed to daughter cell " << g->paths[p->sibling]->child_cell_ID << endl;
          }

          if(g->paths[p->sibling]->child_cell_ID == dcell1->cell_ID){             
            inherit_path_one(p, dcell2, verbose);
            if(verbose > 0) cout << "path " << p->id + 1 << " ditributed to daughter cell " << dcell2->cell_ID << endl;
          }else{
            inherit_path_one(p, dcell1, verbose);
            if(verbose > 0) cout << "path " << p->id + 1 << " ditributed to daughter cell " << dcell1->cell_ID << endl;              
          }
        }else if(p->n_centromere == 0){  // may be lost finally in some cells
           if(verbose > 0) cout << "\ndistribution of path " << p->id + 1 << " without centromere in cell " << cell_ID << endl;
          //  assume all such paths have been duplicated correctly
           distribute_path_random(p, dcell1, dcell2, verbose);          
        }else{
          // if this path has >1 centromeres, maybe circular, duplicated and random breaking for each copy
          if(verbose > 0) cout << "\ncomplex distribution of path " << p->id + 1 << " with >1 centromeres in cell " << cell_ID << endl;
          n_complex_path += 1;
          n_path_break += p->n_centromere - 1; 
          segregate_polycentric(p, dcell1, dcell2, max_pID, n_local_frag, frac_unrepaired, circular_prob, verbose);
          if(verbose > 0) cout << "\nFinish complex splitting\n" << endl;
        }
      }

      if(verbose > 1){
        // path ID follows those from parent cell
        cout << "\n#path after splitting in 1st daughter cell: " << dcell1->g->paths.size() << endl;
        dcell1->print_paths();
        cout << "\n#path after splitting in 2nd daughter cell: " << dcell2->g->paths.size() << endl;
        dcell2->print_paths();
      }
    }



    // assume segments have been generated
    // count number of oscillating CNs by chromosome
    void get_oscillating_CN(int verbose = 0){  
      if(verbose > 1) cout << "count number of oscillating CNs by chromosome\n";     
      for(int chr = 0; chr < NUM_CHR; chr++){
        vector<int> tcns;
        vector<int> tdiff;
        int nseg = g->chr_segments[chr].size();
        if(verbose > 1) cout << chr << "\t" << nseg << endl;
        if(nseg == 0){
          chr_n_osc2[chr] = 0;
          chr_n_osc3[chr] = 0;
          continue;
        }
        if(nseg == 1){
          chr_n_osc2[chr] = 1;
          chr_n_osc3[chr] = 1;
          continue;
        }
        vector<int> v3_states(nseg, 0);
        vector<int> v2_states(nseg, 0);

        for(auto s : g->chr_segments[chr]){
          // s->print();
          int tcn = s->cnA + s->cnB;
          tcns.push_back(tcn);
        }

        for(int i = 0; i < nseg - 2; i++){
          int dcn = tcns[i] - tcns[i + 2];
          if(dcn == 0){
            v3_states[i] = 1;
            v2_states[i] = 1;
          }else if(abs(dcn) == 1){
            v3_states[i] = 1;
          }else{

          }
          if(verbose > 1) cout << dcn << "\t" << v2_states[i] << "\t" << v3_states[i] << endl;
        }

        int n_osc2 = find_max_size(1, v2_states) + 2;
        int n_osc3 = find_max_size(1, v3_states) + 2;
        chr_n_osc2[chr] = n_osc2;
        chr_n_osc3[chr] = n_osc3;
      }

    }


   bool is_copy_changed(int chr, int start, int end, const vector<pos_cn>& dups, int verbose = 0){
    for(auto d : dups){
      if(verbose > 1){
        cout << d.chr + 1 << "\t" << d.start << "\t" << d.end << "\t" << d.cnA << "\t" << d.cnB << endl;
      }
      if(d.chr != chr){
        continue;
      }else{  // on the same chromosome
        if(start >= d.start && end <= d.end){
            return true;
        } 
      }
    }

    return false;
   }


    // get average total CN for each chromosome, scaled by segment size
    void get_chr_tcn(int verbose = 0){
        chr_tcns.clear();
        map<int, vector<double>> chr_cns;
        // intialize to avoid chromosome loss 
        for(int i = 0; i < NUM_CHR; i++){
          vector<double> cns;
          chr_cns[i] = cns;
        }
        for(auto sg : g->chr_segments){
          for(auto s : sg.second){ 
            // string extra = "cn={'c1': {'A': " + to_string(s->cnA) + ", 'B': " + to_string(s->cnB) + "}}";
            // [) interval to avoid manual overlapping
            // string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + extra + "\n";
            // cout << line;               
            int tcn = s->cnA + s->cnB;
            double weighted_tcn = (double) tcn * (s->end - s->start + 1) / CHR_LENGTHS[s->chr];
            chr_cns[s->chr].push_back(weighted_tcn);
            if(verbose > 1) cout << s->chr << "\t" << s->start << "\t" << s->end << "\t" << tcn << "\t" << weighted_tcn << endl;      
          }
        }
        // cout << chr_cns.size() << endl;
        assert(chr_cns.size() == NUM_CHR);

        for(auto ct : chr_cns){
          vector<double> cns = ct.second;
          double acn = 0.0;
          if(cns.size() > 0) acn = accumulate(cns.begin(), cns.end(), 0.0);
          if(verbose > 1){
            cout << ct.first << "\t";
            for(auto cn : cns){
              cout << "\t" << cn;
            }   
            cout << "\t" << acn << endl;
          }           
          chr_tcns.push_back((acn));
        }
        // cout << chr_tcns.size() << endl;
        assert(chr_tcns.size() == NUM_CHR);
    }


    // need to determine arm boundary
    void get_arm_tcn(int verbose = 0){
        arm_tcns.clear();
        map<pair<int, int>, vector<double>> arm_cns;
        for(int i = 0; i < NUM_CHR; i++){
          vector<double> cns;
          pair<int, int> key(i, 0);
          arm_cns[key] = cns;
          pair<int, int> key1(i, 1);
          arm_cns[key1] = cns;          
        }       
        for(auto sg : g->chr_segments){
          for(auto s : sg.second){  
            // string extra = "cn={'c1': {'A': " + to_string(s->cnA) + ", 'B': " + to_string(s->cnB) + "}}";
            // // [) interval to avoid manual overlapping
            // string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + extra + "\n";
            // cout << line;             
            int tcn = s->cnA + s->cnB;
            pair<int, int> key(s->chr, 0);
            pair<int, int> key1(s->chr, 1);
            double weighted_tcn = 0.0;
            double weighted_tcn1 = 0.0;
            int seg_size = s->end - s->start + 1;
            int parm_size = ARM_BOUNDS[s->chr];
            int qarm_size = CHR_LENGTHS[s->chr] - ARM_BOUNDS[s->chr] + 1;
            if(s->start < ARM_BOUNDS[s->chr] && s->end < ARM_BOUNDS[s->chr]){
              weighted_tcn = (double) tcn * seg_size / parm_size;
              arm_cns[key].push_back(weighted_tcn);
            }else if(s->start > ARM_BOUNDS[s->chr] && s->end > ARM_BOUNDS[s->chr]){
              weighted_tcn1 = (double) tcn * seg_size / qarm_size;
              arm_cns[key1].push_back(weighted_tcn1);             
            }else{
              int sizep = ARM_BOUNDS[s->chr] - s->start + 1;
              weighted_tcn = (double) tcn * sizep / parm_size;
              arm_cns[key].push_back(weighted_tcn);

              int sizeq = s->end - ARM_BOUNDS[s->chr] + 1;
              weighted_tcn1 = (double) tcn * sizeq / qarm_size;
              arm_cns[key1].push_back(weighted_tcn1);                 
            }   
            if(verbose > 1) cout << s->chr << "\t" << s->start << "\t" << s->end << "\t" << ARM_BOUNDS[s->chr] << "\t" << tcn << "\t" << weighted_tcn << "\t" << weighted_tcn1 << endl;        
          }
        }
        assert(arm_cns.size() == NUM_CHR * 2);

        for(auto ct : arm_cns){
          vector<double> cns = ct.second;      
          double acn = 0.0;
          if(cns.size() > 0) acn = accumulate(cns.begin(), cns.end(), 0.0);
          if(verbose > 1){
            cout << ct.first.first << "\t" << ct.first.second;
            for(auto cn : cns){
              cout << "\t" << cn;
            }   
            cout << "\t" << acn << endl;
          }          
          arm_tcns.push_back((acn));
        }
        // cout << arm_tcns.size() << endl;
        assert(arm_tcns.size() == NUM_CHR * 2);      
    }


    // according to the formula in Laughney et al, 2015
    void get_surv_prob(int selection_type, double selection_strength, int verbose){     
      double score = 0.0;
      if(selection_type == 0){
        get_chr_tcn(verbose);
        for(int i = 0; i < NUM_CHR; i++){
          score += CHR_SCORE[i] * chr_tcns[i];
          if(verbose > 1) cout << "chr" << "\t" << CHR_SCORE[i] << "\t" << chr_tcns[i] << endl;
        }
        surv_prob = exp(SURVIVAL_D * score + SURVIVAL_C);  
      }else if(selection_type == 1){
        get_arm_tcn(verbose);
        for(int i = 0; i < NUM_CHR * 2; i++){
          score += ARM_SCORE[i] * arm_tcns[i];
          if(verbose > 1) cout << "chr" << i / 2 << "\tarm" << i % 2 << "\t" << ARM_SCORE[i]  << "\t" << arm_tcns[i] << endl;
        }   
        surv_prob = exp(SURVIVAL_D * score + SURVIVAL_C);    
      }else{  // selection based on individual SVs, fitness proportional to SV size
        cout << "not supported selection type!" << endl;
        exit(FAIL);
      }
      surv_prob = pow(surv_prob,  selection_strength);
      if(verbose > 1) cout << "computing survival probability of cell " << cell_ID << " with selection type " << selection_type << " and original score " << score << endl;
    }

    // summary for the whole genome of each cell, computed from the set of adjacencies
    void get_summary_stats(int verbose = 0){
      if(verbose > 1) cout << "computing summary stats\n";
      get_oscillating_CN(verbose);

      for(int i  = 0; i < NUM_CHR; i++){
        vector<int> n_sv(NUM_SVTYPE, 0);
        chr_type_num[i] = n_sv;
      }

      for(auto adjm : g->adjacencies){
        adjacency* adj = adjm.second;
        breakpoint* bp1 = g->breakpoints[adj->junc_id1];
        breakpoint* bp2 = g->breakpoints[adj->junc_id2];
        int type = adj->sv_type;
        if(adj->type != VAR) continue;   // only need to consider variant adjacencies
        int chr1 = bp1->chr;
        int chr2 = bp2->chr;

        int pos1 = bp1->pos;
        int pos2 = bp2->pos;

        // exclude chromosome ends when counting breakpoints (chr ends based on position should only appear in intervals)
        // each breakpoint connects with just one variant adjacency
        if(pos1 != 1 && pos1 != CHR_LENGTHS[chr1]){
          n_bp += 1;
          bp_unique.insert(to_string(chr1) + "_" + to_string(pos1));
         
          int haplotype1 = bp1->haplotype;
          vector<int> dbs1{chr1, pos1, haplotype1};
          dbs_unique.insert(dbs1);
          chr_bp_unique[chr1].insert(dbs1);
        }

        if(pos2 != 1 && pos2 != CHR_LENGTHS[chr2]){
          n_bp += 1;
          bp_unique.insert(to_string(chr2) + "_" + to_string(pos2));
         
          int haplotype2 = bp2->haplotype;
          vector<int> dbs2{chr2, pos2, haplotype2};
          dbs_unique.insert(dbs2);
          chr_bp_unique[chr2].insert(dbs2);  // should be unique when counting directly
        }

        if(chr1 != chr2){
          if(type != BND){
              cout << "Weird type of adjacency!" << endl;
              adj->print();
              bp1->print();
              bp2->print();
              exit(FAIL);
          }
          chr_type_num[chr1][type] += 1;
          chr_type_num[chr2][type] += 1;
        }else{
          chr_type_num[chr1][type] += 1;
        }

        // int start = min(pos1, pos2);
        // int end = max(pos1, pos2);
        // string pos = to_string(chr1) + "_" + to_string(start) + "_" + to_string(end);
        switch(type){
          case DUP: {
            // if(verbose > 1) cout << "candidate DUP: " << chr1 + 1 << "\t" << start << "\t" << end << endl;
            // if(pos_dup_by_adj.count(pos)){  // really duplicated
            // if(is_copy_changed(chr1, start, end, dups, verbose)){  // really duplicated
            n_dup += 1; 
              // pos_dup_by_adj.insert(pos);
            // }else{              
            //   chr_type_num[chr1][type] -= 1;
            //   adj->sv_type = OTHER;
            // }
            break;
          }
          case DEL:{
            // if(verbose > 1) cout << "candidate DEL: " << chr1 + 1 << "\t" << start << "\t" << end << endl;
            // // if(!pos_del_by_adj.count(pos)){
            // if(is_copy_changed(chr1, start, end, dels, verbose)){  // really deleted
            n_del += 1; 
              // pos_del_by_adj.insert(pos);
            // }else{
            //   chr_type_num[chr1][type] -= 1;
            //   adj->sv_type = OTHER;
            // }
            break;
          }
          case H2HINV: n_h2h += 1; break;
          case T2TINV: n_t2t += 1; break;
          case BND: n_tra += 1; break;
          default: ;
        }
      }        

      for(auto pp : g->paths){
        path* p = pp.second;

        if(p->n_centromere > 1){
          cout << "Path " << p->id + 1 << " has " << p->n_centromere << " centromeres!" << endl;
          p->print();
          exit(FAIL);
        }

        if(p->is_circle){
          n_circle += 1;
          if(p->n_centromere == 0) n_ecdna += 1;
        }else{
          if(p->n_centromere == 0) n_nocentro += 1;
        }
      }
    }

    // print summary information for each cell whose values are used for inference
    void print_total_summary(vector<pos_cn>& dups, vector<pos_cn> dels, int verbose = 0){
      // set<string> pos_del_by_adj;
      // set<string> pos_dup_by_adj;
      // compute summary statistics based on adjacency
      // get_summary_stats(pos_dup_by_adj, pos_del_by_adj, verbose);
      // get_summary_stats(verbose);

      if(verbose > 1) cout << "appending real DUP/DEL\n";   
      // set<string>::iterator piter;
      vector<pos_cn>::iterator citer = dups.begin();
      while(citer != dups.end()){
        pos_cn d = *citer;
        // string pos = to_string(d.chr) + "_" + to_string(d.start) + "_" + to_string(d.end);
        // piter = find(pos_dup_by_adj.begin(), pos_dup_by_adj.end(), pos);
        // if(piter == pos_dup_by_adj.end()){
          chr_type_num[d.chr][DUPREAL] += 1;
          n_dup += 1; 
          citer++;       
        // }else{
        //   citer = dups.erase(citer);  // remove duplicated ones
        // }
      }

      citer = dels.begin();
      while(citer != dels.end()){
        pos_cn d = *citer;
        // string pos = to_string(d.chr) + "_" + to_string(d.start) + "_" + to_string(d.end);
        // piter = find(pos_del_by_adj.begin(), pos_del_by_adj.end(), pos);
        // if(piter == pos_del_by_adj.end()){
          chr_type_num[d.chr][DELREAL] += 1;
          n_del += 1;  
          citer++;      
        // }else{
        //   citer = dels.erase(citer);  // remove duplicated ones
        // }
      }
 

      // int nSV = n_dup + n_del + n_h2h + n_t2t;
      cout << div_occur << "\t" << cell_ID;
      //  << "\t" << 0 << "\t" << bp_unique.size() << "\t" << nSV << endl;
      for(int i = 0; i < NUM_CHR; i++){
        int nSV_chr = chr_type_num[i][DEL] + chr_type_num[i][DUP] + chr_type_num[i][H2HINV] + chr_type_num[i][T2TINV];
        // cout << div_occur << "\t" << cell_ID << "\t" << i+1 << "\t" << chr_n_bp[i] << "\t" << nSV_chr << endl;
        cout << "\t" << chr_bp_unique[i].size() << "\t" << chr_n_osc2[i] << "\t" << chr_n_osc3[i]  << "\t" << nSV_chr;
      }
      cout << endl;      
    }


    // a whole cycle of cell division
    // duplication and repair of its genome
    // Path IDs are reencoded in the next cell cylce
    void do_cell_cycle(Cell_ptr dcell1, Cell_ptr dcell2, vector<pos_bp>& bps, vector<double>& bp_fracs, double frac_unrepaired, int& n_telo_fusion, int& n_complex_path, int& n_path_break, int n_local_frag, double frac_unrepaired_local, double circular_prob, int pair_type = 0, double prob_correct_repaired = 0, int verbose = 0){
      // verbose = 2;
      // introduces DSBs into the genome prior to G1 repairs
      if(verbose > 1){
        cout << "Original genome: " << endl;
        g->print();
      }

      if(verbose > 0) cout << "\nG1 phase" << endl;

      if(dsb_rate > 0){  // different for each cycle
        n_dsb = gsl_ran_poisson(r, dsb_rate);
      }  
      // DSB rate also applies for puncpunctuated events 
      if(div_occur > div_break){
        if(verbose > 1) cout << "No new breaks in division " << div_occur << " in cell " << cell_ID << endl;
        n_dsb = 0;
        // n_local_frag = 0;
      }

      g1(bps, bp_fracs, frac_unrepaired, circular_prob, pair_type, prob_correct_repaired, verbose);

      if(verbose > 0) cout << "\nSphase and G2" << endl;
      sphase_g2(n_telo_fusion, verbose);

      if(verbose > 0) cout << "\nMitosis phase" << endl;
      mitosis(dcell1, dcell2, n_complex_path, n_path_break, n_local_frag, frac_unrepaired_local, circular_prob, verbose);

      if(verbose > 2){
        cout << "\nnumber of DSBs in G1: " << n_dsb << endl;
        this->print_bp_adj();
        dcell1->print_bp_adj();
        dcell2->print_bp_adj();
        cout << g->paths.size() << "\t" << dcell1->g->paths.size() << "\t" << dcell2->g->paths.size() << endl;
        print_paths();
        dcell1->print_paths();
        dcell2->print_paths();        
      }

      if(verbose > 0){
        cout << "\nOne cell cycle finish!\n";
        cout << "\nUnrepaired DSBs in this cycle\n";
        for(auto jm : g->breakpoints){
          breakpoint* j = jm.second;
          if(!j->is_repaired){   // breakpoints from previous cycles may not be repaired
            j->print();
          }
        }
      }
    }


    /***********************************  below are functions related to output  ***********************************/

    // Follow format of RCK output
    // CN: "chr", "start", "end", "extra"
    // SV: "chr1", "coord1", "strand1", "chr2", "coord2", "strand2",	"extra"
    void write_rck(string fname_cn, string fname_sv, int verbose = 0){
      ofstream fout_cn(fname_cn);
      string header = "chr\tstart\tend\textra\n";
      fout_cn << header;

      for(auto sg : g->chr_segments){
        for(auto s : sg.second){
          // s->print();
          // use cn1 to be consistent with patterns in Jabba
          string extra = "cn={'c1': {'A': " + to_string(s->cnA) + ", 'B': " + to_string(s->cnB) + "}}";
          // [) interval to avoid manual overlapping
          string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + extra + "\n";
          fout_cn << line;
        }
      }  
      fout_cn.close();

      ofstream fout(fname_sv);
      header = "aid\tchr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\textra\n";
      fout << header;
      g->get_adjacency_CN();
      int id_alt = 0;
      int id_ref = 0;
      string idx = "";
      for(auto adjm : g->adjacency_CNs){
        adj_pos ap = adjm.first;
        adj_cn ac = adjm.second;
        
        if(ap.type == "R"){
          id_ref++;
          idx = ap.type + to_string(id_ref);
        }else{
          id_alt++;
          idx = to_string(id_alt);
        }

        int cn_AA = ac.cnAA;
        int cn_AB = ac.cnAB;
        int cn_BA = ac.cnBA;
        int cn_BB = ac.cnBB;
        int tcn = cn_AA + cn_AB + cn_BA + cn_BB;
        // if(tcn == 0) continue;
        string extra = "aid=" + (idx) + ";cn={'c1': {'AA': " + to_string(cn_AA) + ", 'AB': "+ to_string(cn_AB) + ", 'BA': " + to_string(cn_BA) + ", 'BB': " + to_string(cn_BB) + "}};at=" + ap.type;
        string line = (idx) + "\t" + to_string(ap.chr1 + 1) + "\t" + to_string(ap.pos1) + "\t" + get_side_string(ap.strand1) + "\t" + to_string(ap.chr2 + 1) + "\t" + to_string(ap.pos2) + "\t" + get_side_string(ap.strand2) + "\t" + extra + "\n";
        fout << line;
      }

      // if(verbose > 1) cout << "recomputing segment copy number for each cell to get consecutive regions" << endl;
      // map<int, set<int>> bps_by_chr;
      // g->get_bps_per_chr(bps_by_chr, verbose);
      // g->calculate_segment_cn(bps_by_chr, verbose);
      // for(int chr = 0; chr < NUM_CHR; chr++){
      //   for(auto s: g->chr_segments[chr]){
      //     if(s->cnA > 1 || s->cnB > 1){
      //       string extra = "sv_type=DUP";
      //       string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + "-" + "\t" + to_string(s->chr + 1) + "\t" + to_string(s->end) + "\t" + "+" + "\t" + extra + "\n";
      //       fout << line;
      //     }
      //   }
      // }

      // for(auto d: dups){
      //     string extra = "sv_type=DUP";
      //     string line = to_string(d.chr + 1) + "\t" + to_string(d.start) + "\t" + "-" + "\t" + to_string(d.chr + 1) + "\t" + to_string(d.end) + "\t" + "+" + "\t" + extra + "\n";
      //     fout << line;       
      // }

      // for(auto d: dels){
      //     string extra = "sv_type=DEL";
      //     string line = to_string(d.chr + 1) + "\t" + to_string(d.start) + "\t" + "+" + "\t" + to_string(d.chr + 1) + "\t" + to_string(d.end) + "\t" + "-" + "\t" + extra + "\n";
      //     fout << line;       
      // }

      fout.close();
    }

   // Follow format of RCK output
    // CN: Chromosome  Start.bp    End.bp allele_1 allele_2
    // SV: <SV>	2	213454835	2	213584799	3to5	0	1	0	1
    void write_plot(string fname_cn, string fname_sv, int verbose = 0){
      ofstream fout_cn(fname_cn);
      string header = "Chromosome\tStart.bp\tEnd.bp\tallele_1\tallele_2\n";
      fout_cn << header;
      for(auto sg : g->chr_segments){
        for(auto s : sg.second){
          // [) interval to avoid manual overlapping
          string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + to_string(s->cnA)  + "\t" + to_string(s->cnB) + "\n";
          fout_cn << line;
        }
      }  
      fout_cn.close();

      ofstream fout(fname_sv);
      g->get_adjacency_CN();
      int id_alt = 0;
      int id_ref = 0;
      string idx = "";
      for(auto adjm : g->adjacency_CNs){
        adj_pos ap = adjm.first;
        adj_cn ac = adjm.second;
        
        int cn_AA = ac.cnAA;
        int cn_AB = ac.cnAB;
        int cn_BA = ac.cnBA;
        int cn_BB = ac.cnBB;

        string direction = "";
        if(ap.strand1 == TAIL && ap.strand2 == TAIL){
          direction = "3to3";
        }else if(ap.strand1 == HEAD && ap.strand2 == HEAD){
          direction = "5to5";
        }else if(ap.strand1 == TAIL && ap.strand2 == HEAD){
          direction = "3to5";
        }else{
          assert(ap.strand1 == HEAD && ap.strand2 == TAIL);
          direction = "5to3";
        }
        string line = "<SV>\t" + to_string(ap.chr1 + 1) + "\t" + to_string(ap.pos1) + "\t" + to_string(ap.chr2 + 1) + "\t" + to_string(ap.pos2) + "\t" + direction + "\t" + to_string(cn_AA) + "\t" + to_string(cn_AB) + "\t" + to_string(cn_BA) + "\t" + to_string(cn_BB) + "\n";
        fout << line;
      }

      fout.close();
    }


    // output in a way to facilitate understanding of CN and SV data
    void write_genome(string fname){
        ofstream fout(fname);
        // fout << g->paths.size() << endl;
        // weird path ID, starting from 1
        fout << "ID\tshape\ttype\tNcentromere\tnodes\n";
        for(auto p: g->paths){
          g->write_path(p.second, fout);
        }
    }


    // Assume the summary statistics have been computed
    void write_summary_stats(string fname, string fname_chr){     
      assert(dbs_unique.size() >= 0);
      assert(bp_unique.size() >= 0);
      assert(chr_bp_unique.size() >= 0);
      assert(chr_type_num.size() >= 0);

      // writing summary statistics for all chromosomes
      ofstream fout(fname);
      // nBP includes chromosome ends
      fout << "cycleID\tcellID\tnDSB\tnUnrepair\tnDBS_unique\tnBP_unique\tnBP\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\tnCircle\tnECDNA\tnNoCentromere\n";
      // writing summary statistics for each chromosome
      ofstream fout_chr(fname_chr);
      fout_chr << "cycleID\tcellID\tchr\tnBP\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\n";

      // for(auto dbs : dbs_unique){
      //   cout << dbs[0] << "\t" << dbs[1] << "\t" << dbs[2] << "\n";
      // }

      fout << to_string(div_occur) + "\t" + to_string(cell_ID) + "\t" + to_string(n_dsb) + "\t" + to_string(n_unrepaired)  + "\t" + to_string(dbs_unique.size()) + "\t" + to_string(bp_unique.size()) + "\t" + to_string(n_bp) + "\t" + to_string(n_del) + "\t" + to_string(n_dup) + "\t" + to_string(n_h2h) + "\t" + to_string(n_t2t)  + "\t" + to_string(n_tra) + "\t" + to_string(n_circle) + "\t" + to_string(n_ecdna) + "\t" + to_string(n_nocentro) +"\n";
      fout.close();

      for(int i = 0; i < NUM_CHR; i++){
        fout_chr << to_string(div_occur) + "\t" + to_string(cell_ID) + "\t" + to_string(i + 1) + "\t" << chr_bp_unique[i].size() << "\t" << chr_type_num[i][DEL] << "\t" << chr_type_num[i][DUP] << "\t" << chr_type_num[i][H2HINV] << "\t" << chr_type_num[i][T2TINV] << "\t" << chr_type_num[i][BND] << "\n";
      }
      fout_chr.close();
    }


    // // summary for each chromosome
    // void write_summary_stats(string fname, string fname_chr){
    //   ofstream fout(fname);
    //   fout << "cycleID\tnDSB\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\tnComplex\n";
    //   int n_del = 0;
    //   int n_dup = 0;
    //   int n_h2h = 0;
    //   int n_t2t = 0;
    //   int n_tra = 0;

    //   ofstream fout_chr(fname_chr);
    //   fout_chr << "chr\tnDSB\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\n";
    //   for(int i  = 0; i < NUM_CHR; i++){
    //     n_del += chr_type_num[i][DEL];
    //     n_dup += chr_type_num[i][DUP];
    //     n_h2h += chr_type_num[i][H2HINV];
    //     n_t2t += chr_type_num[i][T2TINV];
    //     n_tra += chr_type_num[i][BND];

    //     fout_chr << i + 1 << "\t" << g->vec_n_dsb[i] << "\t" << chr_type_num[i][DEL] << "\t" << chr_type_num[i][DUP] << "\t" << chr_type_num[i][H2HINV] << "\t" << chr_type_num[i][T2TINV] << "\t" << chr_type_num[i][BND] << endl;
    //   }
    //   fout_chr.close();

    //   string line = to_strins->end - 1g(div_occur) + "\t" + to_string(n_dsb) + "\t" + to_string(n_del) + "\t" + to_string(n_dup) + "\t" + to_string(n_h2h) + "\t" + to_string(n_t2t)  + "\t" + to_string(n_tra) + "\t" + to_string(n_complex_path) + "\n";
    //   fout << line;
    //   fout.close();
    // }


    void write_cn_bin(string fname_cn, const vector<pos_bin>& bins, int verbose = 0){
      ofstream fout_cn(fname_cn);
      string header = "chromosome\tstart\tend\ttotal_cn\tcnA\tcnB\n";
      fout_cn << header;  

      for(int i = 0; i < bins.size(); i++){
        pos_bin bin = bins[i];
        int tcn = round(g->bin_tcn[i]);
        int cnA = round(g->bin_cnA[i]);
        int cnB = round(g->bin_cnB[i]);
        string line = to_string(bin.chr) + "\t" + to_string(bin.start) + "\t" + to_string(bin.end) + "\t" + to_string(tcn) + "\t" + to_string(cnA) + "\t" + to_string(cnB)+ "\n";
        fout_cn << line;
      } 

      fout_cn.close(); 
    }


    // SVtype (character): type of SV, encoded as: DEL (deletion-like; +/-), DUP (duplication-like; -/+), h2hINV (head-to-head inversion; +/+), and t2tINV (tail-to-tail inversion; -/-).
    void write_shatterseek(string fname_cn, string fname_sv, const vector<pos_cn>& dups, const vector<pos_cn>& dels, int verbose = 0){
      // "chromosome", "start", "end", "total_cn"
      ofstream fout_cn(fname_cn);
      string header = "chromosome\tstart\tend\ttotal_cn\tcnA\tcnB\n";
      fout_cn << header;
      for(auto sg : g->chr_segments){
        for(auto s: sg.second){
          // s->print();
          int tcn = s->cnA + s->cnB;
          string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + to_string(tcn) + "\t" + to_string(s->cnA) + "\t" + to_string(s->cnB)+ "\n";
          fout_cn << line;
        }
      }
      fout_cn.close();

      // "chrom1", "start1", "strand1", "chrom2", "end2", "strand2", "svclass"
      ofstream fout(fname_sv);
      header = "chrom1\tstart1\tstrand1\tchrom2\tend2\tstrand2\tsvclass\n";
      fout << header;
      for(auto adjm : g->adjacencies){
        adjacency* adj = adjm.second;
        if(adj->sv_type == NONE) continue;
        string extra = get_sv_type_string(adj->sv_type);
        string line = g->get_sv_string(adj, extra);
        fout << line;
      }

      for(auto d: dups){
          string extra = "DUPREAL";
          string line = to_string(d.chr + 1) + "\t" + to_string(d.start) + "\t" + "-" + "\t" + to_string(d.chr + 1) + "\t" + to_string(d.end) + "\t" + "+" + "\t" + extra + "\n";
          fout << line;       
      }

      for(auto d: dels){
          string extra = "DELREAL";
          string line = to_string(d.chr + 1) + "\t" + to_string(d.start) + "\t" + "+" + "\t" + to_string(d.chr + 1) + "\t" + to_string(d.end) + "\t" + "-" + "\t" + extra + "\n";
          fout << line;       
      }

      fout.close();

    }


    void print_bp_adj(){
      cout << "breakpoints in cell " << cell_ID << endl;
      for(auto jm : g->breakpoints){
        breakpoint* j = jm.second;
        j->print();
      }

      cout << "adjacencies in cell " << cell_ID << endl;
      for(auto am : g->adjacencies){
        adjacency* a = am.second;
        if(a->type == INTERVAL){
          a->print_interval(g->breakpoints);
          a->print();
          g->breakpoints[a->junc_id1]->print();
          g->breakpoints[a->junc_id2]->print();
        }
      }
    }


    void print_paths(){
      cout << "\nAll paths in cell " << cell_ID << endl;
      vector<int> pids;
      for(auto p : g->paths){
        pids.push_back(p.first);
        p.second->print();
      }

      unordered_map<int, int> freq;
      for(auto i: pids) {
          freq[i]++;
      }
 
      for(auto f : freq){
        if(f.second > 1){
          cout << "duplicated path: " << f.first << "\t" << f.second << endl;
          exit(FAIL);
        }       
      }
    }

    void print_cell_info(){
        cout << "Cell " << this->cell_ID << endl;
        cout << "\t parent " << this->parent_ID << endl;
        cout << "\t in clone " << this->clone_ID << endl;

        cout << "\t Occur at division " << this->div_occur << endl;

        cout << "\t birth_rate " << this->birth_rate << endl;
        cout << "\t death_rate " << this->death_rate << endl;

        cout << "\t rate of double strand breaks " << this->dsb_rate << endl;
        cout << "\t number of double strand breaks " << this->n_dsb << endl;
        cout << "\t number of unrepaired breaks " << this->n_unrepaired << endl;
        cout << "\t number of unique breakpoints " << this->bp_unique.size() << endl;
    }

};


#endif
