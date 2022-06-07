#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>

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

    // parameters related to cell growth
    double birth_rate;
    double death_rate;

    // parameters related to DSB generation
    // double dsb_rate;  // use exact number for external control over each cell
    int n_dsb;    // number of double strand breaks during cell division
    int n_unrepaired; // number of unrepaired double strand breaks


    int n_bp = 0;  // each breakpoint may have at most two connections
    int n_del = 0;
    int n_dup = 0;
    int n_h2h = 0;
    int n_t2t = 0;
    int n_tra = 0;
    int n_circle = 0;

    // count events by chromosome
    map<int, set<vector<int>>> chr_bp_unique;  // use position to be consistent with real data
    map<int, vector<int>> chr_type_num;
    map<int, int> chr_n_osc2;
    map<int, int> chr_n_osc3;
    set<int> bp_unique; // a DBS may introduce two breakpoints with different IDs
    set<vector<int>> dbs_unique; // may be fewer than nDSB + nMbreak due to uneven distribution into daughter cells

    genome g; // each cell has a genome

    ~Cell() = default;
    Cell(const Cell& other) = default;
    Cell(Cell&& other) = default;
    Cell& operator=(const Cell& other) = default;
    Cell& operator=(Cell&& other) = default;


    Cell() {
        cell_ID = 0;
        parent_ID = 0;
        clone_ID = 0;

        time_occur = 0;
        div_occur = 0; // whether the cell is sampled or not

        birth_rate = log(2);
        death_rate = 0;

        // dsb_rate = 0;
        n_dsb  = 0;
        n_unrepaired = 0;
    }


    Cell(int cell_ID, int parent_ID, double time_occur) {
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;

        this->time_occur = time_occur;

        genome g(cell_ID);
        this->g = g;
        // this->div_occur = 0;
        //
        // this->birth_rate = birth_rate;
        // this->death_rate = death_rate;
        //
        // this->dsb_rate = 0;
        // this->n_dsb  = 0;
    }


    Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, int n_dsb, int n_unrepaired, double time_occur){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;
        this->clone_ID = 0;

        this->time_occur = time_occur;
        this->div_occur = 0;

        this->birth_rate = birth_rate;
        this->death_rate = death_rate;

        // this->dsb_rate = dsb_rate;
        this->n_dsb  = n_dsb;
        this->n_unrepaired = n_unrepaired;
        // this->n_complex_path = 0;
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

        this->birth_rate = ncell.birth_rate;
        this->death_rate = ncell.death_rate;

        // this->dsb_rate = ncell.dsb_rate;
        this->n_dsb = ncell.n_dsb;
        this->n_unrepaired = ncell.n_unrepaired;
        // this->n_complex_path = ncell.n_complex_path;
        
        this->g.paths.clear();
        this->g.breakpoints.clear();
        this->g.adjacencies.clear();
        this->g.cell_ID = this->cell_ID;
    }


    /*********************** functions related to mutation generations **************************/
    // pass a path and related objects to a daughter
    void inherit_path_both(path* p, Cell_ptr dcell1, Cell_ptr dcell2){
      path* p1 = new path(*p);
      p1->cell_ID = dcell1->cell_ID;
      dcell1->g.paths.push_back(p1);

      path* p2 = new path(*p);
      p2->cell_ID = dcell2->cell_ID;
      dcell2->g.paths.push_back(p2);

      // copy related breakpoints
      for(auto jid : p->nodes){
        breakpoint* j = g.breakpoints[jid];

        breakpoint* j1 = new breakpoint(*j);
        j1->cell_ID = dcell1->cell_ID;
        dcell1->g.breakpoints[jid] = j1;

        breakpoint* j2 = new breakpoint(*j);
        j2->cell_ID = dcell2->cell_ID;
        dcell2->g.breakpoints[jid] = j2;
      }

      // copy related adjacencies
      for(auto aid : p->edges){
        adjacency* a = g.adjacencies[aid];

        adjacency* a1 = new adjacency(*a);
        a1->cell_ID = dcell1->cell_ID;
        dcell1->g.adjacencies[aid] = a1;

        adjacency* a2 = new adjacency(*a);
        a2->cell_ID = dcell2->cell_ID;
        dcell2->g.adjacencies[aid] = a2;
      }
    }


    int find_max_bpID(){
      vector<int> bp_IDs(g.breakpoints.size(), 0);
      int i = 0;

      for(auto bp : g.breakpoints){
        bp_IDs[i++] = bp.first;
      }

      int max_bpID = *max_element(bp_IDs.begin(), bp_IDs.end());

      return max_bpID;
    }


    int find_max_adjID(){
      vector<int> adj_IDs(g.adjacencies.size(), 0);
      int i = 0;

      for(auto adj : g.adjacencies){
        adj_IDs[i++] = adj.first;
      }

      int max_adjID = *max_element(adj_IDs.begin(), adj_IDs.end());

      return max_adjID;
    }


    int find_max_pathID(){
      vector<int> path_IDs(g.paths.size(), 0);
      int i = 0;

      for(auto p : g.paths){
        path_IDs[i++] = p->id;
      }

      int max_pID = *max_element(path_IDs.begin(), path_IDs.end());

      return max_pID;
    }


    // There may be two copies of the same path assigned to the same daughter cell 
    // need to assign different IDs for duplicated path (indicated by n1) and not affect the path to process later
    void inherit_path_one(path* p, Cell_ptr dcell1, int& max_pID, int ncopy = 1, int verbose = 0){
      if(verbose > 1){
        cout << "inherit path " << p->id + 1 << " to cell " << dcell1->cell_ID << endl;
        p->print();
      }
      
      path* p1 = new path(*p);
      p1->cell_ID = dcell1->cell_ID;
      dcell1->g.paths.push_back(p1);

      // copy related breakpoints
      int max_bpID = -1;  // clear nodes in case new IDs are created
      if(ncopy > 1){
        max_bpID = find_max_bpID();
        p1->nodes.clear();     
      }  
      map<int, int> jid_map;  // mapping of copied ID  
      jid_map[-1] = -1;
      for(auto jid : p->nodes){
        breakpoint* j = g.breakpoints[jid];
        if(verbose > 1){
          cout << "copy breakpoint " << j->id << endl;
          j->print();
        }
        breakpoint* j1 = new breakpoint(*j);
        j1->cell_ID = dcell1->cell_ID;

        int njid = jid;
        if(ncopy > 1){
          njid = ++max_bpID;   // plus one first 
          jid_map[j->id] = njid;
          j1->id = njid;        
          p1->nodes.push_back(njid);
          // add the copied breakpoint to parent cell to keep track all and avoid ID reusing
          // point to old breakpoints to avoid early release of the copied breakpoint
          breakpoint* j2 = new breakpoint(*j);
          g.breakpoints[njid] = j2;
        }

        if(verbose > 1){  
          cout << "copied breakpoint " << j1->id << endl;
          j1->print();
        }       
        dcell1->g.breakpoints[njid] = j1;  
      }

      // update neighbors after all breakpoints are copied
      if(ncopy > 1){
        for(auto njid : p1->nodes){
          breakpoint* j1 = dcell1->g.breakpoints[njid];
          j1->left_jid = jid_map[j1->left_jid];
          j1->right_jid = jid_map[j1->right_jid];
        }
      }

      // copy related adjacencies
      int max_adjID = -1;
      if(ncopy > 1){
        max_adjID = find_max_adjID();
        p1->edges.clear(); 
      }
      for(auto aid : p->edges){
        adjacency* a = g.adjacencies[aid];
        if(verbose > 1){
          cout << "copy adjacency " << a->id << endl;
          a->print();
        }
        adjacency* a1 = new adjacency(*a);
        a1->cell_ID = dcell1->cell_ID;

        int naid = aid;
        if(ncopy > 1){
          naid = ++max_adjID;  
          a1->id = naid;  
          // update junction IDs of copied adjacency
          a1->junc_id1 = jid_map[a->junc_id1];
          a1->junc_id2 = jid_map[a->junc_id2];
          p1->edges.push_back(naid);
          adjacency* a2 = new adjacency(*a);
          g.adjacencies[naid] = a2;
        }

        if(verbose > 1){
          cout << "copied adjacency " << a1->id << endl;
          a1->print();
        }         
        dcell1->g.adjacencies[naid] = a1;
      }

      // update path ID
      if(ncopy > 1){ 
         assert(max_pID > 0);     
         p1->id = ++max_pID;
         if(verbose > 1){
           cout << "new path ID " << p1->id + 1 << endl;
           p1->print();
         }
      }
    }


    // distribute a path to either daughter cell
    void distribute_path(path* p, Cell_ptr dcell1, Cell_ptr dcell2, int& max_pID, int verbose = 0){
        double u = runiform(r, 0, 1);
        assert(max_pID < 0);  // no need to update path ID
        if(u < 0.5){
          if(verbose > 0){
            cout << "Distributing split path " << p->id + 1 << " to daughter cell " << dcell1->cell_ID << endl;
            write_path(p, cout);
          }
          inherit_path_one(p, dcell1, max_pID, 1, verbose);
        }else{
          if(verbose > 0){
            cout << "Distributing split path " << p->id + 1 << " to daughter cell " << dcell2->cell_ID << endl;
            write_path(p, cout);
          }           
          inherit_path_one(p, dcell2, max_pID, 1, verbose);
        }
        if(verbose > 1) cout << "Finish path distributing" << endl;
    }


    // called when number of centromeres in a path exceeds 1 (may be caused by random joining in G1)
    // choose breakpoint location so that each path contains one centromere, introduce new breakpoints and form connections
    // assume the path has 2 telomeres at the ends
    // p2split will be split into n_centromere paths, which are distributed randomly
    void segregate_polycentric(path* p2split, Cell_ptr dcell1, Cell_ptr dcell2, int& pid, vector<path*>& new_paths, int verbose = 0){
      // verbose = 2;      
      if(verbose > 0) write_path(p2split, cout);
      if(verbose > 1){      
        p2split->print();
      }

      // remove last connection if the path is circular 
      if(p2split->is_circle){
        breakpoint* js = g.breakpoints[p2split->nodes[0]];       
        int aid = p2split->edges.back();
        adjacency* al = g.adjacencies[aid];
        int js_left = js->left_jid;
        int js_right = js->right_jid;
        int al1 = al->junc_id1;
        int al2 = al->junc_id2;
        if(js_left == al1 || js_left == al2){
          js->left_jid = -1;  
        }else{
          assert(js_right == al1 || js_right == al2);
          js->right_jid = -1;  
        }       
        g.adjacencies.erase(aid);
        p2split->nodes.pop_back();
        p2split->edges.pop_back();
        if(verbose > 1){
          cout << "remove last connection " << aid << " with left breakpoint " << al1 << " and right breakpoint " << al2 << " as the path is circular" << endl;
          p2split->print();
        }
      }
      // Find two intervals containing each centromere
      int nbreak = p2split->n_centromere - 1;     
      vector<pair<int, int>> adjID_pair_withcent;
      int ncent = 0;   // count #centromeres traversed
      int prev_adj = p2split->edges[0];
      adjacency* adj = g.adjacencies[prev_adj];
      if(adj->is_centromeric){
        ncent++;
      }

      for(int i = 1; i < p2split->edges.size(); i++){
        int curr_adj = p2split->edges[i];
        adj = g.adjacencies[curr_adj];
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
          g.adjacencies[pair.first]->print();
          g.adjacencies[pair.second]->print();
        }
      }
    
      assert(adjID_pair_withcent.size() == nbreak);

      // randombly choose an interval to introduce the break
      // when breaking paths with multiple centromeres, a breakpoint must be to the other side of another breakpoint 
      // For each adjacent interval, breakpoint is either at right of 1st interval or left of 2nd interval when finding all breakpoints before splitting intervals
      // Here, left and right should be based on the direction of the interval on the path, so need to check inversion
      int chr = -1;
      int haplotype = -1;
      int left_jid = -1;
      int right_jid = -1;
      vector<bp_interval> bps_in_interval;   // breakpoints at each interval
      // vector<int> old_aids;
      for(auto pair : adjID_pair_withcent){
        int sel = myrng(2);  // 0 or 1
        adjacency* adj = g.adjacencies[pair.first];
        assert(adj->is_centromeric);
        assert(g.breakpoints[adj->junc_id1]->chr == g.breakpoints[adj->junc_id2]->chr);
        assert(g.breakpoints[adj->junc_id1]->haplotype == g.breakpoints[adj->junc_id2]->haplotype);     
        int bp, ps, pe;

        if(sel == 0){  // break on left side of second interval
          adj = g.adjacencies[pair.second];   // should not be ptel
          int jid1 = adj->junc_id1;
          int jid2 = adj->junc_id2;
          assert(g.breakpoints[jid1]->chr == g.breakpoints[jid2]->chr);         
          if(verbose > 1){
            cout << "add break on left side of second interval " << adj->id << endl;
            adj->print();
            g.breakpoints[jid1]->print();
            g.breakpoints[jid2]->print();
          }
          chr = g.breakpoints[jid1]->chr;
          int pos1 = g.breakpoints[jid1]->pos;
          int pos2 = g.breakpoints[jid2]->pos;          
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
        }else{ // break on right side of first interval
          // should not be qtel
          int jid1 = adj->junc_id1;
          int jid2 = adj->junc_id2;
          assert(g.breakpoints[jid1]->chr == g.breakpoints[jid2]->chr);
          if(verbose > 1){
            cout << "add break on right side of first interval " << adj->id << endl;
            adj->print();
            g.breakpoints[jid1]->print();
            g.breakpoints[jid2]->print();            
          }         
          chr = g.breakpoints[jid1]->chr; 
          int pos1 = g.breakpoints[jid1]->pos;
          int pos2 = g.breakpoints[jid2]->pos;
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

        if(verbose > 1) cout << ps << "\t" << pe << endl;
        assert(ps < pe);
        bp = (int)runiform(r, ps, pe);
        left_jid = adj->junc_id1;
        right_jid = adj->junc_id2;
        int haplotype = g.breakpoints[adj->junc_id1]->haplotype;
        bp_interval bi{bp, chr, haplotype, left_jid, right_jid, adj->is_inverted};
        bps_in_interval.push_back(bi);
        // old_aids.push_back(adj->id);
      }

      assert(bps_in_interval.size() == nbreak);

      int jid = g.breakpoints.rbegin()->first + 1;
      int aid = g.adjacencies.rbegin()->first + 1;
      // some interval may be chosen twice and be broken at 1st choice
      int prev_left_jid = -1;
      // int prev_right_jid = -1;
      int i = 0;
      vector<breakpoint*> junc_starts;  // start breakpoints for new paths

      for(auto bi : bps_in_interval){
        int bp = bi.bp;
        int chr = bi.chr;
        int haplotype = bi.haplotype;
        // neighbours from bi are before any new breaks and may be changed if a break occur in this interval
        int left_jid = bi.left_jid;
        int right_jid = bi.right_jid;
        // cout << "cell ID before adding breakpoints " << this->cell_ID << endl;
        // cout << g.breakpoints.size() << " breakpoints " << endl;
        // update left_jid if an interval is chosen twice
        if(prev_left_jid == left_jid){  
          if(bi.is_inverted){
            right_jid = jid - 2;
          }else{
            left_jid = jid - 1;
          }
        }
 
        if(verbose > 0){
          cout << "   path break " << i + 1 << " at position " << bp << " chr " << chr + 1 << " haplotype " << get_haplotype_string (haplotype) << endl;
          if(verbose > 1) cout << "   with left breakpoint " << left_jid << " and right breakpoint " << right_jid << endl;
        }
       
        i++;
        g.add_new_breakpoint(jid, aid, chr, bp, haplotype, left_jid, right_jid, false, verbose);
        // one breakpoint is end of a path, the other is the start of another path, depending on the orientations
        // jid has increased by 2 after add_new_breakpoint
        junc_starts.push_back(g.breakpoints[jid - 2]);
        junc_starts.push_back(g.breakpoints[jid - 1]);

        // update number of DSBs per chrom per haplotype
        // int idx = chr + haplotype * NUM_CHR;
        // g.vec_n_mitosis_break[idx] += 1;
        prev_left_jid = left_jid;
        // prev_right_jid = right_jid;
      }

      assert(junc_starts.size() == nbreak * 2);

      // split paths at once to save cost
      // no need to rejoin paths after splitting, just distribute to different daughter cells
      if(verbose > 1){
        cout << "\nStart breakpoints of split paths: " << endl;
        for(auto j : junc_starts){
          cout << j->id << " at path " << j->path_ID + 1 << endl;
        }
        cout << "#paths before: " << g.paths.size() << endl;
      }

      int nnode = 0;
      // recompute path from each breakpoint to avoid complexity of tracking each broken adjacency
      // The path starting from both path breaks will be distributed twice
      if(verbose > 0) cout << "Segregating into " << junc_starts.size() << " subpaths for path " << p2split->id + 1 << endl;
      for(auto js : junc_starts){
        path* p = new path(++pid, this->cell_ID, COMPLETE);
        new_paths.push_back(p);  

        if(verbose > 1){
          cout << "new path " << p->id + 1 << " at breakpoint " << js->id << endl;
          js->print();
        }

        js->path_ID = p->id;
        p->nodes.push_back(js->id);
        // there will be no cycles
        g.get_connected_breakpoint(js, p, g.adjacencies, false);
        g.set_path_type(p, verbose);
        // p->print();

        int dcell_sel = myrng(2);  // 0 or 1
        // each path will have a unique ID
        int max_pID = -1;
        distribute_path(p, dcell1, dcell2, max_pID, verbose);
        // path p will not be in current cell
        nnode += p->nodes.size();
      }

      assert(nnode = p2split->nodes.size() + 2 * nbreak);     
    }


    // introduce DSB and repair breaks (adding variant adjacencies)
    void g1(int n_dsb, int n_unrepaired, int verbose = 0){
      // verbose = 1;
       if(verbose > 0) cout << "\nIntroducing " << n_dsb << " DSBs" << endl;
      g.generate_dsb(n_dsb, verbose);

      int n_torepair = n_dsb - n_unrepaired;
      if(verbose > 0) cout << "\nRepairing " << n_torepair << " DSBs" << endl;
      vector<breakpoint*> junc2repair;
      if(n_torepair > 0){      
        g.repair_dsb(n_unrepaired, junc2repair, verbose);
      }

      if(verbose > 1){
        cout << "\nRemaing breakpoints to repair: " << endl;
        for(auto j: junc2repair){
            j->print();
        }
        cout << "\nConnecting breakpoints to paths" << endl;
      }

      if(verbose > 0) cout << "\nGetting the derivative genome" << endl;
      g.get_derivative_genome(verbose);

      if(verbose > 0){
        cout << "\nAll paths in the derivative genome" << endl;
        for(auto p: g.paths){
          p->print();
        }
        cout << "\nEnd of G1 phase" << endl;
      }
    }


    // S phase and G2 implemented together
    void sphase_g2(int& n_telo_fusion, int verbose = 0){
      // verbose = 1;     
      for(auto p : g.paths){
        if(!p->is_circle){        
          assert(p->nodes.size() % 2 == 0);
        }
        assert(p->nodes.size() == p->edges.size() + 1);

        // replicates incomplete paths by adding new breakpoints (S) and path fusion
        if(p->type != COMPLETE && !p->is_circle){  
          if(p->type == NONTEL && p->n_centromere > 0) continue;   // skip these paths to avoid circular DNA with centromere
          if(verbose > 0){
            cout << "Duplicating incomplete path " << p->id + 1 << endl;
            write_path(p, cout);
          }

          // TODO: distinguish fusion at which cell cycle to see its linkage with chromothripsis?
          if(p->type == PTEL || p->type == QTEL){
            n_telo_fusion += 1;
          }

          if(verbose > 1){ 
            cout << "before" << endl;         
            p->print();
          }

          g.duplicate_path(*p, verbose);

          if(verbose > 1){
            cout << "after" << endl; 
            p->print();
          }
        }
      }
    }


    void random_split(path* p, Cell_ptr dcell1, Cell_ptr dcell2, int& max_pID, int verbose = 0){
      // verbose = 1;
      int n1 = 0;  // times this path is distributed to cell 1
      int n2 = 0;  // times this path is distributed to cell 2

      if(verbose > 0) cout << "\nRandom distribution of path " << p->id + 1 << endl;

      for(int i = 0; i < 2; i++){  // imitate duplication       
        string desc = "";
        if(p->n_centromere == 0){
          desc = " (without centromere) ";
        }else{
          desc = " (circular) ";
        }

        double u = runiform(r, 0, 1);
        // need to change path ID when the same path gets distributed to the same cell twice
        if(u < 0.5){
          n1++;
          if(verbose > 0) cout << "Distributing path " << p->id + 1 << desc << " to daughter cell " << dcell1->cell_ID << endl;
          inherit_path_one(p, dcell1, max_pID, n1, verbose);
        }else{
          n2++;
          if(verbose > 0) cout << "Distributing path " << p->id + 1 << desc << " to daughter cell " << dcell2->cell_ID << endl;
          inherit_path_one(p, dcell2, max_pID, n2, verbose);
        }
      }
    }


    // randomly distribute segments between daughter cells, depend on #centromeres
    void mitosis(Cell_ptr dcell1, Cell_ptr dcell2, int&  n_complex_path, int& n_path_break, int verbose = 0){
      // verbose = 2;
      if(verbose > 1){
        cout << g.paths.size() << " paths before split (may include duplicated path distributed to daughter cells)" << endl;
        for(auto p : g.paths){
          p->print();
        }
      }

      // get the max Path ID to define new IDs for random distribution
      int max_pID = find_max_pathID();
      vector<path*> new_paths;
      for(auto p : this->g.paths){
        if(verbose > 1) p->print();
        if(!p->is_circle){   // a circular path may have even number of nodes during path replacation
          assert(p->nodes.size() % 2 == 0);
        }

        if(p->n_centromere == 1){  // balanced distribution
          if(verbose > 1) cout << "balanced distribution of path " << p->id + 1 << endl;
          inherit_path_both(p, dcell1, dcell2);
        }else if(p->n_centromere == 0){  // may be lost finally in some cells
          // delay later to get all the breakpoints
        }else{
          // if this path has >1 centromeres, maybe circular, duplicated and random breaking for each copy
          if(verbose > 1) cout << "complex distribution of path " << p->id + 1 << endl;
          n_complex_path += 1;
          n_path_break += p->n_centromere - 1;
          // for(int i = 0; i < 2; i++){ // do it twice to iminate duplication ?
          // }         
          segregate_polycentric(p, dcell1, dcell2, max_pID, new_paths, verbose);
          if(verbose > 1) cout << "\nFinish complex splitting" << endl;
        }
      }

      g.paths.insert(g.paths.begin(), new_paths.begin(), new_paths.end());
      max_pID = find_max_pathID();  // new paths may be added before
      for(auto p : this->g.paths){
         if(p->n_centromere == 0){
           random_split(p, dcell1, dcell2, max_pID, verbose);
           if(verbose > 0) write_path(p, cout);
         }
      }

      if(verbose > 1){
        // path ID follows those from parent cell
        cout << "\npath after split in 1st daughter cell" << endl;
        cout << "#path after splitting: " << dcell1->g.paths.size() << endl;
        for(auto p : dcell1->g.paths){
          p->print();
        }
        cout << "\npath after split in 2nd daughter cell" << endl;
        cout << "#path after splitting: " << dcell2->g.paths.size() << endl;
        for(auto p : dcell2->g.paths){
          p->print();
        }
      }
    }


    // Follow format of RCK output
    // "chr1", "coord1", "strand1", "chr2", "coord2", "strand2",	"extra"
    void write_sv(string fname){
      // string fname = "sv_data.tsv"
      ofstream fout(fname);
      string header = "chr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\textra\n";
      fout << header;

      for(auto adjm : g.adjacencies){
        adjacency* adj = adjm.second;
        if(adj->sv_type == NONE) continue;

        breakpoint* j1 = g.breakpoints[adj->junc_id1];
        breakpoint* j2 = g.breakpoints[adj->junc_id2];
        string type = "sv_type=" + get_sv_type_string(adj->sv_type);
        string extra = type;
        // int cn_AA = 0;
        // int cn_AB = 0;
        // int cn_BA = 0;
        // int cn_BB = 0;
        // string cn = "cn={'c1':{'AA': " + to_string(cn_AA) + "'AB': "+ to_string(cn_AB) + ", 'BA':" + to_string(cn_BA) + ", 'BB':" + to_string(cn_BB) + "}";
        // type = type + ";" + cn;
        string line = to_string(j1->chr + 1) + "\t" + to_string(j1->pos) + "\t" + get_side_string(j1->side) + "\t" + to_string(j2->chr + 1) + "\t" + to_string(j2->pos) + "\t" + get_side_string(j2->side) + "\t" + extra + "\n";
        fout << line;
      }

      fout.close();
    }


    // "chr", "start", "end", "extra"
    void write_cnv(string fname){
      // string fname = "cn_data.tsv"
      ofstream fout(fname);
      string header = "chr\tstart\tend\textra\n";
      fout << header;

      for(auto s : g.segments){
        // s->print();
        string extra = "cn={'1': {'A': " + to_string(s->cnA) + ", 'B': " + to_string(s->cnB) + "}}";
        string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + extra + "\n";
        fout << line;
      }

      fout.close();
    }


    void write_path(path* p, ostream& fout){
        // cout << "writing path" << endl;
        string shape = "linear";
        if(p->is_circle){
          shape = "circular";
        }    
             
        fout << p->id + 1 << "\t" << shape << "\t" << get_telomere_type_string(p->type) << "\t" << p->n_centromere << "\t";

        for(auto e : p->edges){
          // cout << "\nedge " << e << endl;
          adjacency *adj = g.adjacencies[e];
          // adj->print();
          breakpoint *j1 = g.breakpoints[adj->junc_id1];
          breakpoint *j2 = g.breakpoints[adj->junc_id2];           
          if(adj->is_inverted){
            j1 = g.breakpoints[adj->junc_id2];  
            j2 = g.breakpoints[adj->junc_id1];
          }

          fout << (j1->chr % NUM_CHR) + 1 << ":" << get_haplotype_string(j1->haplotype) << ":" << j1->pos << "-" << (j2->chr % NUM_CHR) + 1 << ":" << get_haplotype_string(j2->haplotype) << ":" << j2->pos << ",";
        }   
        fout << endl;       
    }

    
    // output in a way to facilitate understanding of CN and SV data
    void write_genome(string fname){
        // print_bp_adj();
        // g.get_derivative_genome();  // to check consistencies of paths, breakpoints, adjacencies
        
        ofstream fout(fname);
        // fout << g.paths.size() << endl;
        // weird path ID, starting from 1
        fout << "ID\tshape\ttype\tNcentromere\tnodes\n";
        for(auto p: g.paths){
          write_path(p, fout);
        }
    }


    

    // count number of oscillating CNs by chromosome
    void get_oscillating_CN(int verbose = 0){
      // group segments by chr 
      map<int, vector<segment*>> chr_segments;
      for(auto s : g.segments){
        chr_segments[s->chr].push_back(s);
      }

      for(int chr = 0; chr < NUM_CHR; chr++){
        vector<int> tcns;
        vector<int> tdiff;        
        int nseg = chr_segments[chr].size();
        if(verbose > 0) cout << chr << "\t" << nseg << endl;
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

        for(auto s : chr_segments[chr]){
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
          if(verbose > 0) cout << dcn << "\t" << v2_states[i] << "\t" << v3_states[i] << endl;
        }


        int n_osc2 = find_max_size(1, v2_states) + 2;
        int n_osc3 = find_max_size(1, v3_states) + 2;
        chr_n_osc2[chr] = n_osc2;
        chr_n_osc3[chr] = n_osc3;
      }

    }


    // summary for the whole genome of each cell, computed from the set of adjacencies
    void get_summary_stats(int verbose = 0){
      get_oscillating_CN(verbose);

      for(int i  = 0; i < NUM_CHR; i++){
        vector<int> n_sv(NUM_SVTYPE, 0);
        chr_type_num[i] = n_sv;
      }

      for(auto adjm : g.adjacencies){
        adjacency* adj = adjm.second;
        int type = adj->sv_type;
        if(adj->type != VAR) continue;   // only need to consider variant adjacencies 
        int chr1 = g.breakpoints[adj->junc_id1]->chr;
        int chr2 = g.breakpoints[adj->junc_id2]->chr;

        // exclude chromosome end when counting breakpoints (chr_ends should only appear in intervals)
        // each breakpoint connects with just one variant adjacency
        if(find(g.chr_ends.begin(), g.chr_ends.end(), adj->junc_id1) == g.chr_ends.end()){
          n_bp += 1;
          bp_unique.insert(adj->junc_id1);
         

          int pos1 = g.breakpoints[adj->junc_id1]->pos;
          int haplotype1 = g.breakpoints[adj->junc_id1]->haplotype;
          vector<int> dbs1{chr1, pos1, haplotype1};
          dbs_unique.insert(dbs1);
          chr_bp_unique[chr1].insert(dbs1);
        }

        if(find(g.chr_ends.begin(), g.chr_ends.end(), adj->junc_id2) == g.chr_ends.end()){
          n_bp += 1;
          bp_unique.insert(adj->junc_id2);

          int pos2 = g.breakpoints[adj->junc_id2]->pos;
          int haplotype2 = g.breakpoints[adj->junc_id2]->haplotype;
          vector<int> dbs2{chr2, pos2, haplotype2};
          dbs_unique.insert(dbs2);
          chr_bp_unique[chr2].insert(dbs2);  // should be unique when counting directly
        }

        if(chr1 != chr2){
          assert(type == TRA);
          chr_type_num[chr1][type] += 1;
          chr_type_num[chr2][type] += 1;
        }else{
          chr_type_num[chr1][type] += 1;
        }

        switch(type){
          case DUP: n_dup += 1; break;
          case DEL: n_del += 1; break;
          case H2HINV: n_h2h += 1; break;
          case T2TINV: n_t2t += 1; break;
          case TRA: n_tra += 1; break;
          default: ;
        }
      }

      for(auto p : g.paths){
        if(p->n_centromere > 1){
          cout << "Path " << p->id + 1 << " has " << p->n_centromere << " centromeres!" << endl;
          p->print();
          exit(-1);
        }
        if(p->is_circle) n_circle += 1;
      }
    }

    
    // Assume the summary statistics have been computed
    void write_summary_stats(string fname, string fname_chr){
      // writing summary statistics for all chromosomes
      ofstream fout(fname);
      // nBP includes chromosome ends
      fout << "cycleID\tcellID\tnDSB\tnUnrepair\tnDBS_unique\tnBP_unique\tnBP\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\tnCircle\n";
      // writing summary statistics for each chromosome
      ofstream fout_chr(fname_chr);
      fout_chr << "cycleID\tcellID\tchr\tnBP\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\n";
   
      // for(auto dbs : dbs_unique){
      //   cout << dbs[0] << "\t" << dbs[1] << "\t" << dbs[2] << "\n";
      // }

      fout << to_string(div_occur) + "\t" + to_string(cell_ID) + "\t" + to_string(n_dsb) + "\t" + to_string(n_unrepaired)  + "\t" + to_string(dbs_unique.size()) + "\t" + to_string(bp_unique.size()) + "\t" + to_string(n_bp) + "\t" + to_string(n_del) + "\t" + to_string(n_dup) + "\t" + to_string(n_h2h) + "\t" + to_string(n_t2t)  + "\t" + to_string(n_tra) + "\t" + to_string(n_circle) + "\n";
      fout.close();

      for(int i = 0; i < NUM_CHR; i++){
        fout_chr << to_string(div_occur) + "\t" + to_string(cell_ID) + "\t" + to_string(i + 1) + "\t" << chr_bp_unique[i].size() << "\t" << chr_type_num[i][DEL] << "\t" << chr_type_num[i][DUP] << "\t" << chr_type_num[i][H2HINV] << "\t" << chr_type_num[i][T2TINV] << "\t" << chr_type_num[i][TRA] << "\n";
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
    //     n_tra += chr_type_num[i][TRA];

    //     fout_chr << i + 1 << "\t" << g.vec_n_dsb[i] << "\t" << chr_type_num[i][DEL] << "\t" << chr_type_num[i][DUP] << "\t" << chr_type_num[i][H2HINV] << "\t" << chr_type_num[i][T2TINV] << "\t" << chr_type_num[i][TRA] << endl;
    //   }
    //   fout_chr.close();

    //   string line = to_string(div_occur) + "\t" + to_string(n_dsb) + "\t" + to_string(n_del) + "\t" + to_string(n_dup) + "\t" + to_string(n_h2h) + "\t" + to_string(n_t2t)  + "\t" + to_string(n_tra) + "\t" + to_string(n_complex_path) + "\n";
    //   fout << line;
    //   fout.close();
    // }


    // SVtype (character): type of SV, encoded as: DEL (deletion-like; +/-), DUP (duplication-like; -/+), h2hINV (head-to-head inversion; +/+), and t2tINV (tail-to-tail inversion; -/-).
    void write_shatterseek(string fname_sv, string fname_cn){
      // "chrom1", "start1", "strand1", "chrom2", "end2", "strand2", "svclass"
      ofstream fout(fname_sv);
      string header = "chrom1\tstart1\tstrand1\tchrom2\tend2\tstrand2\tsvclass\n";
      fout << header;
      for(auto adjm : g.adjacencies){
        adjacency* adj = adjm.second;
        int type = adj->sv_type;
        if(type == NONE) continue;
        breakpoint* j1 = g.breakpoints[adj->junc_id1];
        breakpoint* j2 = g.breakpoints[adj->junc_id2];

        string line = to_string(j1->chr + 1) + "\t" + to_string(j1->pos) + "\t" + get_side_string(j1->side) + "\t" + to_string(j2->chr + 1) + "\t" + to_string(j2->pos) + "\t" + get_side_string(j2->side) + "\t" + get_sv_type_string(type) + "\n";
        fout << line;
      }
      fout.close();

      // "chromosome", "start", "end", "total_cn"
      ofstream fout_cn(fname_cn);
      header = "chromosome\tstart\tend\ttotal_cn\n";
      fout_cn << header;
      for(auto s : g.segments){
        // s->print();
        int tcn = s->cnA + s->cnB;
        string line = to_string(s->chr + 1) + "\t" + to_string(s->start) + "\t" + to_string(s->end) + "\t" + to_string(tcn) + "\n";
        fout_cn << line;
      }
      fout_cn.close();
    }

  
    void print_bp_adj(){
      cout << "breakpoints in daughter cell " << cell_ID << endl;
      for(auto jm : g.breakpoints){
        breakpoint* j = jm.second;
        j->print();
      }

      cout << "adjacencies in cell " << cell_ID << endl;
      for(auto am : g.adjacencies){
        adjacency* a = am.second;
        if(a->type == INTERVAL){
          a->print_interval(g.breakpoints);
          a->print();
          g.breakpoints[a->junc_id1]->print();
          g.breakpoints[a->junc_id2]->print();
        }
      }   
    }

    // a whole cycle of cell division
    // duplication and repair of its genome
    void do_cell_cycle(Cell_ptr dcell1, Cell_ptr dcell2, int& n_telo_fusion, int& n_complex_path, int& n_path_break, int verbose = 0){
      // verbose = 2;
      // introduces DSBs into the genome prior to G1 repairs
      if(verbose > 1){
        cout << "Original genome: " << endl;
        g.print();
      }

      if(verbose > 0) cout << "\nG1 phase" << endl;
      g1(n_dsb, n_unrepaired, verbose);

      if(verbose > 0) cout << "\nSphase and G2" << endl;
      sphase_g2(n_telo_fusion, verbose);

      if(verbose > 0) cout << "\nMitosis phase" << endl;
      mitosis(dcell1, dcell2, n_complex_path, n_path_break, verbose);

      if(verbose > 2){
        cout << "\nnumber of DSBs: " << n_dsb << endl;
        this->print_bp_adj();
        dcell1->print_bp_adj();
        dcell2->print_bp_adj();
      }

      // if(verbose > 1) cout << "computing segments in current cell" << endl;
      // g.calculate_segment_cn(verbose);
      // if(verbose > 1) cout << "\ncomputing segments in daughter cell " << dcell1->cell_ID << endl;
      // dcell1->g.calculate_segment_cn(verbose);
      // if(verbose > 1) cout << "\ncomputing segments in daughter cell " << dcell2->cell_ID << endl;
      // dcell2->g.calculate_segment_cn(verbose);

      if(verbose > 0) cout << "\nOne cell cycle finish!\n";
    }


/*********************************************************************************/

    /*********************** function related to output generations **************************/

    void print_cell_info(){
        cout << "Cell " << this->cell_ID << endl;
        cout << "\t parent " << this->parent_ID << endl;
        cout << "\t in clone " << this->clone_ID << endl;

        cout << "\t Occur at division " << this->div_occur << endl;

        cout << "\t birth_rate " << this->birth_rate << endl;
        cout << "\t death_rate " << this->death_rate << endl;

        cout << "\t number of double strand breaks " << this->n_dsb << endl;
        cout << "\t number of unrepaired breaks " << this->n_unrepaired << endl;
    }
};


#endif
