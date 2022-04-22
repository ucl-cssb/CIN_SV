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
        // n_complex_path = 0;
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


    void inherit_path_one(path* p, Cell_ptr dcell1, int verbose = 0){
      if(verbose) cout << "inherit path " << p->id << " to cell " << dcell1->cell_ID << endl;
      path* p1 = new path(*p);
      p1->cell_ID = dcell1->cell_ID;
      dcell1->g.paths.push_back(p1);

      // copy related breakpoints
      for(auto jid : p->nodes){
        breakpoint* j = g.breakpoints[jid];
        if(verbose) cout << "copy breakpoint " << j->id << endl;
        breakpoint* j1 = new breakpoint(*j);
        j1->cell_ID = dcell1->cell_ID;
        dcell1->g.breakpoints[jid] = j1;
      }

      // copy related adjacencies
      for(auto aid : p->edges){
        adjacency* a = g.adjacencies[aid];
        if(verbose) cout << "copy adjacency " << a->id << endl;
        adjacency* a1 = new adjacency(*a);
        a1->cell_ID = dcell1->cell_ID;
        dcell1->g.adjacencies[aid] = a1;
      }

    }


    // called when number of centromeres in a path exceeds 1
    // function chooses breakpoint location so that each path contains one centromere, introduce new breakpoints and form connections
    // assume the path has 2 telomeres at the ends
    void segregate_polycentric(path* p2split, Cell_ptr dcell1, Cell_ptr dcell2, int& pid, int verbose = 0){
      // verbose = 1;
      // p will be split into n_centromere paths, which are distributed randomly
      int nbreak = p2split->n_centromere - 1;
      if(verbose){
        cout << "complicate segregation into " << p2split->n_centromere << " subpaths for path " << p2split->id << endl;
        p2split->print();
      }

      // Find two intervals containing each centromere
      vector<pair<int, int>> adjID_pair_withcent;
      vector<breakpoint*> junc_starts;  // start breakpoints for new paths
      int ncent = 0;   // count #centromeres traversed

      int prev_adj = p2split->edges[0];
      adjacency* adj = g.adjacencies[prev_adj];
      // adj->print();
      if(adj->is_centromeric){
        ncent++;
      }
      // avoid using original ends
      // if(g.breakpoints[adj->junc_id1]->is_end){
      //   junc_starts.push_back(g.breakpoints[adj->junc_id1]);
      // }else{
      //   assert(g.breakpoints[adj->junc_id2]->is_end);
      //   junc_starts.push_back(g.breakpoints[adj->junc_id2]);
      // }


      for(int i = 1; i < p2split->edges.size(); i++){
        int curr_adj = p2split->edges[i];
        adjacency* adj = g.adjacencies[curr_adj];
        // adj->print();
        if(adj->is_centromeric){  // can be true only for intervals
          // cout << "is_centromeric" << endl;
          ncent++;
          if(ncent == 2){
            adjID_pair_withcent.push_back(pair<int, int>(prev_adj, curr_adj));
            // restart from current interval
            ncent = 1;
          }
          prev_adj = curr_adj;
        }
      }

      if(verbose){
        cout << "Adjacent intervals with " << adjID_pair_withcent.size() + 1 << " centromeres: " << endl;
        for(auto pair : adjID_pair_withcent){
          cout << pair.first << "\t" << pair.second << endl;
          g.adjacencies[pair.first]->print();
          g.adjacencies[pair.second]->print();
        }
      }

      assert(adjID_pair_withcent.size() == nbreak);

      // randombly choose an interval to introduce the break
      // when breaking paths with multiple centromeres, the breakpoints must be to the other side of a breakpoint at the end
      // For each adjacent interval, breakpoint is either at right of 1st interval or left of 2nd interval when finding all breakpoints before splitting intervals
      // Here, left and right should be based on the direction
      int chr = -1;
      int haplotype = -1;
      int left_jid = -1;
      int right_jid = -1;
      vector<bp_interval> bps_in_interval;   // breakpoints at each interval
      for(auto pair : adjID_pair_withcent){
        int sel = myrng(2);  // 0 or 1
        adjacency* adj = g.adjacencies[pair.first];
        assert(adj->is_centromeric);
        // assert(g.breakpoints[adj->junc_id1]->chr == g.breakpoints[adj->junc_id2]->chr);
        // assert(g.breakpoints[adj->junc_id1]->haplotype == g.breakpoints[adj->junc_id2]->haplotype);
       
        int bp, ps, pe;
        if(sel == 0){  // break on left side of second interval
          adj = g.adjacencies[pair.second];   // should not be ptel
          // cout << "add break on left side of second interval " << adj->id << endl;
          // adj->print();
          int jid = adj->junc_id1;
          if(adj->is_inverted){
            jid = adj->junc_id2;
          }
          chr = g.breakpoints[jid]->chr;
          pe = CENT_STARTS[chr];
          ps = g.breakpoints[jid]->pos;
          bp = (int)runiform(r, ps, pe);
        }else{ // break on right side of first interval
          // should not be qtel
          // cout << "add break on right side of first interval " << adj->id << endl;
          // adj->print();
          int jid = adj->junc_id2;
          if(adj->is_inverted){
            jid = adj->junc_id1;
          }
          chr = g.breakpoints[jid]->chr;
          ps = CENT_ENDS[chr];
          pe = g.breakpoints[jid]->pos;
          bp = (int)runiform(r, ps, pe);
        }
        left_jid = adj->junc_id1;
        right_jid = adj->junc_id2;
        int haplotype = g.breakpoints[adj->junc_id1]->haplotype;
        bp_interval bi{bp, chr, haplotype, left_jid, right_jid};
        bps_in_interval.push_back(bi);
      }

      assert(bps_in_interval.size() == nbreak);

      int jid = g.breakpoints.rbegin()->first + 1;
      int aid = g.adjacencies.rbegin()->first + 1;
      // some interval may be chosen twice and be broken at 1st choice
      int prev_left_jid = -1;
      for(auto bi : bps_in_interval){
        int bp = bi.bp;
        int chr = bi.chr;
        int haplotype = bi.haplotype;
        int left_jid = bi.left_jid;
        int right_jid = bi.right_jid;
        // cout << "cell ID before adding breakpoints " << this->cell_ID << endl;
        // cout << g.breakpoints.size() << " breakpoints " << endl;
        // update left_jid if an interval is chosen twice
        if(prev_left_jid == left_jid){
          left_jid = jid - 1;
        }
        g.add_new_breakpoint(this->cell_ID, jid, aid, chr, bp, haplotype, left_jid, right_jid, verbose);
        // g.breakpoints[jid - 1]->print();
        // g.breakpoints[jid - 2]->print();
        // one breakpoint is end of a path, the other is the start of another path, depending on the orientations
        // jid has increased by 2 after add_new_breakpoint
        junc_starts.push_back(g.breakpoints[jid - 1]);
        junc_starts.push_back(g.breakpoints[jid - 2]);

        // update number of DSBs per chrom per haplotype
        // int idx = chr + haplotype * NUM_CHR;
        // g.vec_n_mitosis_break[idx] += 1;
        prev_left_jid = left_jid;
      }

      assert(junc_starts.size() == p2split->n_centromere);

      // split paths at once to save cost
      // no need to rejoin paths after splitting, just distribute to different daughter cells
      if(verbose){
        cout << "Start breakpoints of split paths: " << endl;
        for(auto j : junc_starts){
          cout << j->id << " at path " << j->path_ID << endl;
        }
        cout << "#path before: " << g.paths.size() << endl;
      }

      // use ID order in parent cell for convenience
      // reset path ID for p2split
      // for(auto jid : p2split->nodes){
      //   g.breakpoints[jid]->print();
      //   g.breakpoints[jid]->path_ID = -1;
      // }
      // for(auto aid : p2split->edges){
      //   g.adjacencies[aid]->print();
      //   g.adjacencies[aid]->path_ID = -1;
      // }

      int nnode = 0;
      for(auto js : junc_starts){
        path* p = new path(pid++, this->cell_ID, COMPLETE);

        if(verbose){
          cout << "new path " << p->id << " at breakpoint " << js->id << endl;
          js->print();
        }

        js->path_ID = p->id;
        p->nodes.push_back(js->id);
        // there will be no cycles
        g.get_connected_breakpoint(js, p, g.adjacencies, false);

        // end nodes should be telomere
        int size = p->nodes.size();
        if(g.breakpoints[p->nodes[0]]->id == g.breakpoints[p->nodes[size-1]]->id){
          p->is_circle = true;
        }
        if(js->is_end){  // pTel
          p->type = PTEL;
        }else if(g.breakpoints[ p->nodes[size-1]]->is_end){  // qTel
          p->type = QTEL;
        }else{     // nonTel
          p->type = NONTEL;
        }

        // p->print();

        int dcell_sel = myrng(2);  // 0 or 1
        double u = runiform(r, 0, 1);
        if(u < 0.5){
          if(verbose) cout << "distribute split path " << p->id << " to daughter cell 1" << endl;
          inherit_path_one(p, dcell1, verbose);
        }else{
          if(verbose) cout << "distribute split path " << p->id << " to daughter cell 2" << endl;
          inherit_path_one(p, dcell2, verbose);
        }
        // path p will not be in current cell
        nnode += p->nodes.size();
      }

      assert(nnode = p2split->nodes.size() + 2 * nbreak);
    }


    // introduce DSB and repair breaks (adding variant adjacencies)
    void g1(int n_dsb, int n_unrepaired, int verbose = 0){
      cout << "Introducing " << n_dsb << " DSBs" << endl;
      g.generate_dsb(n_dsb, verbose);

      vector<breakpoint*> junc2repair;
      g.repair_dsb(n_unrepaired, junc2repair);

      if(verbose){
        cout << "Remaing breakpoints to repair: " << endl;
        for(auto j: junc2repair){
            j->print();
        }
      }

      cout << "Connecting breakpoints to paths" << endl;
      g.get_derivative_genome();
      if(verbose){
        for(auto p: g.paths){
          p->print();
        }
      }

      cout << "End of G1 phase" << endl;
    }


    // S phase and G2 implemented together
    void sphase_g2(int verbose = 0){
      // replicates incomplete paths by adding new breakpoints (S) and path fusion
      for(auto p : g.paths){
        if(p->is_circle){
          assert(p->nodes.size() % 2 != 0);
        }else{
          assert(p->nodes.size() % 2 == 0);
        }
        assert(p->nodes.size() == p->edges.size() + 1);

        if(p->type != COMPLETE && !p->is_circle){
          if(verbose){
            cout << "\nduplicate path " << p->id << endl;
            p->print();
          }
          g.duplicate_path(*p);

          if(verbose) p->print();
        }
      }
    }

    // randomly distribute segments between daughter cells, depend on #centromeres
    void mitosis(Cell_ptr dcell1, Cell_ptr dcell2, int&  n_complex_path, int& n_path_break, int verbose = 0){
      dcell1->g.paths.clear();
      dcell2->g.paths.clear();
      dcell1->g.breakpoints.clear();
      dcell2->g.breakpoints.clear();
      dcell1->g.adjacencies.clear();
      dcell2->g.adjacencies.clear();

      dcell1->g.cell_ID = dcell1->cell_ID;
      dcell2->g.cell_ID = dcell2->cell_ID;

      if(verbose){
        cout << g.paths.size() << " paths before split (may include duplicated path distributed to daughter cells)" << endl;
        for(auto p : g.paths){
          p->print();
        }
      }

      int pid = g.paths.size();
      for(auto p : this->g.paths){
        if(p->is_circle){
          assert(p->nodes.size() % 2 != 0);
        }else{
          assert(p->nodes.size() % 2 == 0);
        }
        assert(p->nodes.size() == p->edges.size() + 1);

        if(p->n_centromere == 1){  // balanced distribution
          if(verbose) cout << "balanced distribution of path " << p->id << endl;
          inherit_path_both(p, dcell1, dcell2);
        }else if(p->n_centromere == 0 || p->is_circle){  // may be lost finally
          if(verbose) cout << "random distribution of path " << p->id << endl;
          double u = runiform(r, 0, 1);
          if(u < 0.5){
            inherit_path_one(p, dcell1);
          }else{
            inherit_path_one(p, dcell2);
          }
        }else{
          // if this path is complete, it needs to be duplicated or not?
          if(verbose) cout << "complex distribution of path " << p->id << endl;
          n_complex_path += 1;
          n_path_break += p->n_centromere - 1;
          segregate_polycentric(p, dcell1, dcell2, pid, verbose);
        }
      }

      if(verbose){
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


    // summary for the whole genome of each cell, computed from the set of adjacencies
    // TO exclude chromosome end when counting breakpoints
    void write_summary_stats(string fname, string fname_chr){
      ofstream fout(fname);
      fout << "cycleID\tcellID\tnDSB\tnUnrepair\tnDBS_unique\tnBP_unique\tnBP\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\tnCircle\n";

      ofstream fout_chr(fname_chr);
      fout_chr << "cycleID\tcellID\tchr\tnBP\tnDel\tnDup\tnH2HInv\tnT2TInv\tnTra\n";

      int n_bp = 0;  // each breakpoint may have at most two connections
      int n_del = 0;
      int n_dup = 0;
      int n_h2h = 0;
      int n_t2t = 0;
      int n_tra = 0;
      int n_circle = 0;
      // int n_complex_path = 0; // should be 0 in the final cells

      // count events by chromosome
      map<int, int> chr_n_bp;
      map<int, vector<int>> chr_type_num;
      for(int i  = 0; i < NUM_CHR; i++){
        chr_n_bp[i] = 0;
        vector<int> n_sv(NUM_SVTYPE, 0);
        chr_type_num[i] = n_sv;
      }

      set<int> bp_unique; // a DBS may introduce two breakpoints with different IDs
      set<vector<int>> dbs_unique; // may be fewer than nDSB + nMbreak due to uneven distribution into daughter cells
      for(auto adjm : g.adjacencies){
        adjacency* adj = adjm.second;
        int type = adj->sv_type;
        int chr1 = g.breakpoints[adj->junc_id1]->chr;
        int chr2 = g.breakpoints[adj->junc_id2]->chr;

        if(find(g.chr_ends.begin(), g.chr_ends.end(), adj->junc_id1) == g.chr_ends.end()){
          n_bp += 1;
          bp_unique.insert(adj->junc_id1);

          int pos1 = g.breakpoints[adj->junc_id1]->pos;
          int haplotype1 = g.breakpoints[adj->junc_id1]->haplotype;
          vector<int> dbs1{chr1, pos1, haplotype1};
          dbs_unique.insert(dbs1);
        }

        if(find(g.chr_ends.begin(), g.chr_ends.end(), adj->junc_id2) == g.chr_ends.end()){
          n_bp += 1;
          bp_unique.insert(adj->junc_id2);

          int pos2 = g.breakpoints[adj->junc_id2]->pos;
          int haplotype2 = g.breakpoints[adj->junc_id2]->haplotype;
          vector<int> dbs2{chr2, pos2, haplotype2};
          dbs_unique.insert(dbs2);
        }

        if(chr1 != chr2){
          assert(type == TRA);
          chr_type_num[chr1][type] += 1;
          chr_type_num[chr2][type] += 1;
          chr_n_bp[chr1] += 1;
          chr_n_bp[chr2] += 1;
        }else{
          chr_type_num[chr1][type] += 1;
          chr_n_bp[chr1] += 2;
        }

        switch(type){
          case DUP: n_dup += 1;
          case DEL: n_del += 1;
          case H2HINV: n_h2h += 1;
          case T2TINV: n_t2t += 1;
          case TRA: n_tra += 1;
          default: ;
        }
      }

      for(auto p : g.paths){
        if(p->n_centromere > 1){
          cout << "Path " << p->id << " has " << p->n_centromere << " centromeres!" << endl;
          p->print();
        }
        if(p->is_circle) n_circle += 1;
      }

      // for(auto dbs : dbs_unique){
      //   cout << dbs[0] << "\t" << dbs[1] << "\t" << dbs[2] << "\n";
      // }

      fout << to_string(div_occur) + "\t" + to_string(cell_ID) + "\t" + to_string(n_dsb) + "\t" + to_string(n_unrepaired)  + "\t" + to_string(dbs_unique.size()) + "\t" + to_string(bp_unique.size()) + "\t" + to_string(n_bp) + "\t" + to_string(n_del) + "\t" + to_string(n_dup) + "\t" + to_string(n_h2h) + "\t" + to_string(n_t2t)  + "\t" + to_string(n_tra) + "\t" + to_string(n_circle) + "\n";
      fout.close();

      for(int i  = 0; i < NUM_CHR; i++){
        fout_chr << to_string(div_occur) + "\t" + to_string(cell_ID) + "\t" + to_string(i + 1) + "\t" << chr_n_bp[i] << "\t" << chr_type_num[i][DEL] << "\t" << chr_type_num[i][DUP] << "\t" << chr_type_num[i][H2HINV] << "\t" << chr_type_num[i][T2TINV] << "\t" << chr_type_num[i][TRA] << "\n";
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


    // a whole cycle of cell division
    // duplication and repair of its genome
    void do_cell_cycle(Cell_ptr dcell1, Cell_ptr dcell2, int&  n_complex_path, int& n_path_break, int verbose = 0){
      // introduces DSBs into the genome prior to G1 repairs
      if(verbose){
        cout << "Original genome: " << endl;
        g.print();
      }

      cout << "\nG1 phase" << endl;
      g1(n_dsb, n_unrepaired, verbose);

      cout << "\nSphase and G2" << endl;
      sphase_g2(verbose);

      cout << "\nmitosis phase" << endl;
      mitosis(dcell1, dcell2, n_complex_path, n_path_break, verbose);

      if(verbose){
        cout << "\nnumber of DSBs: " << n_dsb << endl;

        cout << "adjacencies in current cell" << endl;
        for(auto am : g.adjacencies){
          adjacency* a = am.second;
          a->print();
        }

        cout << "adjacencies in daughter cell 1" << endl;
        for(auto am : dcell1->g.adjacencies){
          adjacency* a = am.second;
          if(a->type == INTERVAL){
            a->print_interval(dcell1->g.breakpoints);
            // a->print();
            // dcell1->g.breakpoints[a->junc_id1]->print();
            // dcell1->g.breakpoints[a->junc_id2]->print();
          }
        }

        cout << "adjacencies in daughter cell 2" << endl;
        for(auto am : dcell2->g.adjacencies){
          adjacency* a = am.second;
          if(a->type == INTERVAL){
            a->print_interval(dcell2->g.breakpoints);
            // a->print();
            // dcell2->g.breakpoints[a->junc_id1]->print();
            // dcell2->g.breakpoints[a->junc_id2]->print();
          }
        }

        cout << "breakpoints in current cell" << endl;
        for(auto jm : g.breakpoints){
          breakpoint* j = jm.second;
          j->print();
        }

        cout << "breakpoints in daughter cell 1" << endl;
        for(auto jm : dcell1->g.breakpoints){
          breakpoint* j = jm.second;
          j->print();
        }

        cout << "breakpoints in daughter cell 2" << endl;
        for(auto jm : dcell2->g.breakpoints){
          breakpoint* j = jm.second;
          j->print();
        }
      }

      cout << "computing segments in current cell" << endl;
      g.calculate_segment_cn(verbose);

      cout << "computing segments in daughter cell 1" << endl;
      dcell1->g.calculate_segment_cn(verbose);

      cout << "computing segments in daughter cell 2" << endl;
      dcell2->g.calculate_segment_cn(verbose);
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
