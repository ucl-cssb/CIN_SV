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


struct Coord{
    double x;
    double y;
    double z;

    Coord() {}
    Coord(double x, double y, double z) : x(x), y(y), z(z) {}

    // sort on all axis
    // bool operator<(const Coord &o) const {
    //     if (x != o.x) {
    //      return x < o.x;
    //     }
    //     if (y != o.y) {
    //      return y < o.y;
    //     }
    //     return z < o.z;
    // }

    bool operator<(const Coord &o) const {
        return x < o.x;
    }
};


class Mutation
{
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


class Cell
{
public:
    int cell_ID;
    int parent_ID;
    int clone_ID;

    // int flag;   // whether the cell is alive or not. 0:new, 1:divided, -1: death
    double time_occur;

    int is_sampled;

    // parameters related to cell growth
    double birth_rate;
    double death_rate;

    Coord pos{0, 0, 0};

    // parameters related to mutation generations
    double dsb_rate;

    // parameters related to storing mutations
    int num_dsb;    // number of double strand breaks during cell division

    // record CNs of each segment
    // record SVs
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
        is_sampled = 0; // whether the cell is sampled or not

        birth_rate = log(2);
        death_rate = 0;

        dsb_rate = 0;
        num_dsb  = 0;
    }


    Cell(int cell_ID, int parent_ID, double time_occur) {
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;

        this->time_occur = time_occur;

        genome g(cell_ID);
        this->g = g;
        // this->is_sampled = 0;
        //
        // this->birth_rate = birth_rate;
        // this->death_rate = death_rate;
        //
        // this->dsb_rate = 0;
        // this->num_dsb  = 0;
    }


    Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, double dsb_rate, double time_occur){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;
        this->clone_ID = 0;

        this->time_occur = time_occur;
        this->is_sampled = 0;

        this->birth_rate = birth_rate;
        this->death_rate = death_rate;

        this->dsb_rate = dsb_rate;
        this->num_dsb  = 0;
    }


    bool operator<(const Cell &o) const {
        return pos < o.pos;
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
        this->is_sampled = ncell.is_sampled;
        this->pos = ncell.pos;

        this->birth_rate = ncell.birth_rate;
        this->death_rate = ncell.death_rate;

        this->dsb_rate = ncell.dsb_rate;
        this->num_dsb = ncell.num_dsb;

    }


    void update_loc(double offset_x, double offset_y, double offset_z) {
        pos.x += offset_x;
        pos.y += offset_y;
        pos.z += offset_z;
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

      // copy related junctions
      for(auto jid : p->nodes){
        junction* j = g.junctions[jid];

        junction* j1 = new junction(*j);
        j1->cell_ID = dcell1->cell_ID;
        dcell1->g.junctions[jid] = j1;

        junction* j2 = new junction(*j);
        j2->cell_ID = dcell2->cell_ID;
        dcell2->g.junctions[jid] = j2;
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


    void inherit_path_one(path* p, Cell_ptr dcell1){
      cout << "inherit path " << p->id << " to cell " << dcell1->cell_ID << endl;
      path* p1 = new path(*p);
      p1->cell_ID = dcell1->cell_ID;
      dcell1->g.paths.push_back(p1);

      // copy related junctions
      for(auto jid : p->nodes){
        junction* j = g.junctions[jid];
        cout << "copy junction " << j->id << endl;
        junction* j1 = new junction(*j);
        j1->cell_ID = dcell1->cell_ID;
        dcell1->g.junctions[jid] = j1;
      }

      // copy related adjacencies
      for(auto aid : p->edges){
        adjacency* a = g.adjacencies[aid];
        cout << "copy adjacency " << a->id << endl;
        adjacency* a1 = new adjacency(*a);
        a1->cell_ID = dcell1->cell_ID;
        dcell1->g.adjacencies[aid] = a1;
      }

    }


    // called when number of centromeres in a path exceeds 1
    // function chooses breakpoint location, introduce new breakpoints and form connections
    // assume the path has 2 telomeres at the ends
    void segregate_polycentric(path* p2split, Cell_ptr dcell1, Cell_ptr dcell2){
      bool debug = true;
      // p will be split into nbreak paths, which are distributed randomly
      int nbreak = p2split->num_centromere - 1;
      if(debug){
        cout << "complicate segregation into " << nbreak << " subpaths for path " << p2split->id << endl;
        p2split->print();
      }

      // Find two intervals containing each centromere
      vector<pair<int, int>> adjID_pair_withcent;
      vector<junction*> junc_starts;  // start junctions for new paths
      int ncent = 0;   // count #centromeres traversed
      int prev_adj = p2split->edges[0];
      adjacency* adj = g.adjacencies[prev_adj];
      // adj->print();
      if(adj->is_centromeric){
        ncent++;
      }
      if(g.junctions[adj->junc_id1]->is_end){
        junc_starts.push_back(g.junctions[adj->junc_id1]);
      }else{
        assert(g.junctions[adj->junc_id2]->is_end);
        junc_starts.push_back(g.junctions[adj->junc_id2]);
      }


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

      if(debug){
        cout << "Adjacent intervals with " << adjID_pair_withcent.size() + 1 << " centromeres: " << endl;
        for(auto pair : adjID_pair_withcent){
          cout << pair.first << "\t" << pair.second << endl;
          g.adjacencies[pair.first]->print();
          g.adjacencies[pair.second]->print();
        }
      }

      assert(adjID_pair_withcent.size() == nbreak);

      // randombly choose an interval to introduce the break
      // for each adjacent interval, breakpoint is either at right of 1st interval or left of 2nd interval
      // find all breakpoints before splitting intervals
      int chr = -1;
      int haplotype = -1;
      int left_jid = -1;
      int right_jid = -1;
      vector<bp_interval> bps_in_interval;   // breakpoints at each interval
      for(auto pair : adjID_pair_withcent){
        int sel = myrng(2);  // 0 or 1
        adjacency* adj = g.adjacencies[pair.first];
        assert(adj->is_centromeric);
        // assert(g.junctions[adj->junc_id1]->chr == g.junctions[adj->junc_id2]->chr);
        // assert(g.junctions[adj->junc_id1]->haplotype == g.junctions[adj->junc_id2]->haplotype);

        int bp, ps, pe;
        if(sel == 0){  // break on left side of second interval
          adj = g.adjacencies[pair.second];   // should not be ptel
          // cout << "add break on left side of second interval " << adj->id << endl;
          // adj->print();
          chr = g.junctions[adj->junc_id1]->chr;
          pe = CENT_STARTS[chr];
          ps = g.junctions[adj->junc_id1]->pos;
          bp = (int)runiform(r, ps, pe);
        }else{ // break on right side of first interval
          // should not be qtel
          // cout << "add break on right side of first interval " << adj->id << endl;
          // adj->print();
          chr = g.junctions[adj->junc_id1]->chr;
          ps = CENT_ENDS[chr];
          pe = g.junctions[adj->junc_id2]->pos;
          bp = (int)runiform(r, ps, pe);
        }
        left_jid = adj->junc_id1;
        right_jid = adj->junc_id2;
        int haplotype = g.junctions[adj->junc_id1]->haplotype;
        bp_interval bi{bp, chr, haplotype, left_jid, right_jid};
        bps_in_interval.push_back(bi);
      }

      assert(bps_in_interval.size() == nbreak);

      int jid = g.junctions.rbegin()->first + 1;
      int aid = g.adjacencies.rbegin()->first + 1;
      // some interval may be chosen twice and be broken at 1st choice
      int prev_left_jid = -1;
      for(auto bi : bps_in_interval){
        int bp = bi.bp;
        int chr = bi.chr;
        int haplotype = bi.haplotype;
        int left_jid = bi.left_jid;
        int right_jid = bi.right_jid;
        // cout << "cell ID before adding junctions " << this->cell_ID << endl;
        // cout << g.junctions.size() << " junctions " << endl;
        // update left_jid if an interval is chosen twice
        if(prev_left_jid == left_jid){
          left_jid = jid - 1;
        }
        g.add_new_junction(this->cell_ID, jid, aid, chr, bp, haplotype, left_jid, right_jid);
        // g.junctions[jid - 1]->print();
        // g.junctions[jid - 2]->print();
        // 1st junction is end of a path, 2nd is the start of another path
        // jid has increased by 2 after add_new_junction
        junc_starts.push_back(g.junctions[jid - 1]);

        // update number of DSBs per chrom per haplotype
        int idx = 2 * chr + haplotype;
        g.vec_num_mitosis_break[idx] += 1;
        prev_left_jid = left_jid;
      }

      assert(junc_starts.size() == p2split->num_centromere);

      // split paths at once to save cost
      // no need to rejoin paths after splitting, just distribute to different daughter cells
      if(debug){
        cout << "Start junctions of split paths: " << endl;
        for(auto j : junc_starts){
          cout << j->id << " at path " << j->path_ID << endl;
        }
        cout << "#path before: " << g.paths.size() << endl;
      }

      // use ID order in parent cell for convenience
      // reset path ID for p2split
      // for(auto jid : p2split->nodes){
      //   g.junctions[jid]->print();
      //   g.junctions[jid]->path_ID = -1;
      // }
      // for(auto aid : p2split->edges){
      //   g.adjacencies[aid]->print();
      //   g.adjacencies[aid]->path_ID = -1;
      // }

      int pid = g.paths.size();
      int nnode = 0;
      for(auto js : junc_starts){
        path* p = new path(pid++, this->cell_ID, COMPLETE);
        cout << "new path " << p->id << " at junction " << js->id << endl;
        js->print();
        js->path_ID = p->id;
        p->nodes.push_back(js->id);
        // there will be no cycles
        g.get_connected_junction(js, p, g.adjacencies, false);

        // end nodes should be telomere
        int size = p->nodes.size();
        if(g.junctions[p->nodes[0]]->id == g.junctions[p->nodes[size-1]]->id){
          p->is_circle = true;
        }
        if(js->is_end){  // pTel
          p->type = PTEL;
        }else if(g.junctions[ p->nodes[size-1]]->is_end){  // qTel
          p->type = QTEL;
        }else{     // nonTel
          p->type = NONTEL;
        }

        p->print();

        int dcell_sel = myrng(2);  // 0 or 1
        double u = runiform(r, 0, 1);
        if(u < 0.5){
          cout << "distribute split path " << p->id << " to daughter cell 1" << endl;
          inherit_path_one(p, dcell1);
        }else{
          cout << "distribute split path " << p->id << " to daughter cell 2" << endl;
          inherit_path_one(p, dcell2);
        }
        // path p will not be in current cell
        nnode += p->nodes.size();
      }

      assert(nnode = p2split->nodes.size() + 2 * nbreak);
    }


    // introduce DSB and repair breaks (adding variant adjacencies)
    void g1(int n_dsb, int n_unrepaired){
      bool debug = false;

      cout << "Introducing " << n_dsb << " DSBs" << endl;
      g.generate_dsb(n_dsb);

      vector<junction*> junc2repair;
      g.repair_dsb(n_unrepaired, junc2repair);
      cout << "Remaing junctions to repair: " << endl;
      for(auto j: junc2repair){
          j->print();
      }

      cout << "Connecting junctions to paths" << endl;
      g.get_derivative_genome();
      for(auto p: g.paths){
        p->print();
      }

      cout << "End of G1 phase" << endl;
    }


    // S phase and G2 implemented together
    void sphase_g2(){
      // replicates incomplete paths by adding new junctions (S) and path fusion
      for(auto p : g.paths){
        if(p->is_circle){
          assert(p->nodes.size() % 2 != 0);
        }else{
          assert(p->nodes.size() % 2 == 0);
        }
        assert(p->nodes.size() == p->edges.size() + 1);

        if(p->type != COMPLETE && !p->is_circle){
          cout << "\nduplicate path " << p->id << endl;
          p->print();
          g.duplicate_path(*p);
          p->print();
        }
      }
    }

    // randomly distribute segments between daughter cells
    void mitosis(Cell_ptr dcell1, Cell_ptr dcell2){
      bool debug = true;

      // depend on #centromeres
      dcell1->g.paths.clear();
      dcell2->g.paths.clear();
      dcell1->g.junctions.clear();
      dcell2->g.junctions.clear();
      dcell1->g.adjacencies.clear();
      dcell2->g.adjacencies.clear();

      dcell1->g.cell_ID = dcell1->cell_ID;
      dcell2->g.cell_ID = dcell2->cell_ID;

      if(debug){
        cout << g.paths.size() << " paths before split (may include duplicated path distributed to daughter cells)" << endl;
        for(auto p : g.paths){
          p->print();
        }
      }

      for(auto p : this->g.paths){
        if(p->is_circle){
          assert(p->nodes.size() % 2 != 0);
        }else{
          assert(p->nodes.size() % 2 == 0);
        }
        assert(p->nodes.size() == p->edges.size() + 1);

        if(p->num_centromere == 1){  // balanced distribution
          cout << "balanced distribution of path " << p->id << endl;

          inherit_path_both(p, dcell1, dcell2);

        }else if(p->num_centromere == 0 || p->is_circle){  // may be lost finally
          cout << "random distribution of path " << p->id << endl;
          double u = runiform(r, 0, 1);
          if(u < 0.5){
            inherit_path_one(p, dcell1);
          }else{
            inherit_path_one(p, dcell2);
          }
        }else{
          cout << "complex distribution of path " << p->id << endl;
          segregate_polycentric(p, dcell1, dcell2);
        }
      }

      if(debug){
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
        junction* j1 = g.junctions[adj->junc_id1];
        junction* j2 = g.junctions[adj->junc_id2];
        if(adj->sv_type == NONE) continue;
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


    void write_summary_stats(string dir){
      // whole genome
      ofstream fout("sumStats_total.tsv");
      string header = "nDSB\tnDel\tnInv\tnIns\tnDup\tcycleID\tnBiasedChroms\tnMitosisBreaks\n";
      fout << header;
      fout.close();

      // each chromosome
      ofstream fout_chr("sumStats_chrom.tsv");
      header = "chr\tnDSB\tnOsc\tnDel\tnIns\tnInv\tnDup";
      fout_chr << header;
      for(int i  = 0; i < NUM_CHR; i++){
        string line = "";
        fout_chr << line << endl;
      }
      fout_chr.close();
    }


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
        junction* j1 = g.junctions[adj->junc_id1];
        junction* j2 = g.junctions[adj->junc_id2];

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
    void do_cell_cycle(int n_unrepaired, Cell_ptr dcell1, Cell_ptr dcell2){
      bool debug = true;

      // introduces DSBs into the genome prior to G1 repairs
      // int diff = max_dsb - min_dsb;
      // int rdm = myrng(diff);
      // num_dsb = min_dsb + diff;
      num_dsb = gsl_ran_poisson(r, dsb_rate);

      cout << "Original genome: " << endl;
      g.print();

      cout << "\nG1 phase" << endl;

      g1(num_dsb, n_unrepaired);

      cout << "\nSphase and G2" << endl;
      sphase_g2();

      cout << "\nmitosis phase" << endl;
      mitosis(dcell1, dcell2);

      if(debug){
        cout << "\nnumber of DSBs: " << num_dsb << endl;

        cout << "adjacencies in current cell" << endl;
        for(auto am : g.adjacencies){
          adjacency* a = am.second;
          a->print();
        }

        cout << "adjacencies in daughter cell 1" << endl;
        for(auto am : dcell1->g.adjacencies){
          adjacency* a = am.second;
          if(a->type == INTERVAL){
            a->print_interval(dcell1->g.junctions);
            // a->print();
            // dcell1->g.junctions[a->junc_id1]->print();
            // dcell1->g.junctions[a->junc_id2]->print();
          }
        }

        cout << "adjacencies in daughter cell 2" << endl;
        for(auto am : dcell2->g.adjacencies){
          adjacency* a = am.second;
          if(a->type == INTERVAL){
            a->print_interval(dcell2->g.junctions);
            // a->print();
            // dcell2->g.junctions[a->junc_id1]->print();
            // dcell2->g.junctions[a->junc_id2]->print();
          }
        }

        cout << "junctions in current cell" << endl;
        for(auto jm : g.junctions){
          junction* j = jm.second;
          j->print();
        }

        cout << "junctions in daughter cell 1" << endl;
        for(auto jm : dcell1->g.junctions){
          junction* j = jm.second;
          j->print();
        }

        cout << "junctions in daughter cell 2" << endl;
        for(auto jm : dcell2->g.junctions){
          junction* j = jm.second;
          j->print();
        }
      }

      cout << "segments in current cell" << endl;
      g.calculate_segment_cn();

      cout << "segments in daughter cell 1" << endl;
      dcell1->g.calculate_segment_cn();

      cout << "segments in daughter cell 2" << endl;
      dcell2->g.calculate_segment_cn();
    }


/*********************************************************************************/

    /*********************** function related to output generations **************************/

    void print_cell_info(){
        cout << "All infomation of cell " << this->cell_ID << endl;
        cout << "\t parent " << this->parent_ID << endl;
        cout << "\t in clone " << this->clone_ID << endl;

        // cout << "\t is_sampled " << this->is_sampled << endl;

        // cout << "\t position " << this->pos.x << endl;

        cout << "\t birth_rate " << this->birth_rate << endl;
        cout << "\t death_rate " << this->death_rate << endl;

        cout << "\t dsb_rate " << this->dsb_rate << endl;
    }
};


#endif
