// genome.hpp

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>


using namespace std;

#include "util.hpp"

#include "breakpoint.hpp"
#include "adjacency.hpp"
#include "segment.hpp"
#include "path.hpp"


class genome {
public:
  int cell_ID;

  // breakpoint (breakpoint) and adjacency define the genome graph
  // use map to access each element by its unique ID and avoid duplications
  map<int, breakpoint*> breakpoints;
  map<int, adjacency*> adjacencies;
  map<int, path*> paths;

  // map<pair<int, int>, vector<breakpoint*>> junc_map;  // group breakpoints by haplotype and chr to facilitate neighbor searching
  map<pair<int, int>, vector<interval>>  cn_by_chr_hap;
  map<pair<int, int>, vector<interval>>  cn_by_chr_hap_merged;
  map<int, vector<segment*>> chr_segments;     // group segments by chr
  vector<double> bin_tcn;    // use double to avoid data loss on bins spanning segments
  vector<double> bin_cnA;
  vector<double> bin_cnB;

  map<adj_pos, adj_cn> adjacency_CNs;   // group adjacencies by location

  // number of DSBs on each chr each haplotype
  // vector<int> vec_n_dsb;   // hard to maintain across cell divisions due to random path distributions
  // breaks introduced due to multiple centromeres during mitosis
  // vector<int> vec_n_mitosis_break;  // seem not biologically meaningful
  // map<int, vector<int>> chr_type_num;  // chr, num for each SV type indexed by TYPE

  ~genome(){
    // cout << "start releasing pointers" << endl;

    for(auto p : paths){
      delete p.second;
    }
    for(auto am : adjacencies){
      delete am.second;
    }
    for(auto jm : breakpoints){
      delete jm.second;
    }
    for(auto sg : chr_segments){
      for(auto s : sg.second){
        delete s;
      }
    }

    // cout << "finish releasing pointers" << endl;
  }

  genome(){

  }

  // only breakpoint at ends and interval adjacencies
  genome(int cell_ID){
    this->cell_ID = cell_ID;
    // create breakpoints on chromosome ends, for each haplotype
    int jid = 0;
    int aid = 0;
    int pid = 0;
    for(int chr = 0; chr < NUM_CHR; chr++){  //  start from 0 for access convenience
      // vec_n_dsb.push_back(0);
      // vector<int> n_sv(NUM_SVTYPE, 0);
      // chr_type_num[chr] = n_sv;

      bool is_end = true;
      bool is_repaired = true;

      // two breakpoints are generated
      int haplotype = 0;
      // int cell_ID, int id, int chr, int pos, int side, int haplotype, bool is_end, bool is_repaired
      breakpoint* j1 = new breakpoint(cell_ID, jid++, chr, 1, TAIL, haplotype, is_end, is_repaired);
      breakpoint* j2 = new breakpoint(cell_ID, jid++, chr, CHR_LENGTHS[chr], HEAD, haplotype, is_end, is_repaired);
      j1->right_jid = j2->id;
      j2->left_jid = j1->id;
      j1->path_ID = pid;
      j2->path_ID = pid;
      breakpoints[j1->id] = j1;
      breakpoints[j2->id] = j2;
      // interval
      adjacency* adj = new adjacency(cell_ID, aid++, pid, j1->id, j2->id, INTERVAL, COMPLETE, NONE);
      adjacencies[adj->id] = adj;
      // path
      path* p = new path(pid++, cell_ID, COMPLETE);
      p->nodes.push_back(j1->id);
      p->nodes.push_back(j2->id);
      p->edges.push_back(adj->id);
      // validate_path(p);
      paths[p->id] = p;

      haplotype = 1;
      breakpoint* j3 = new breakpoint(cell_ID, jid++, chr, 1, TAIL, haplotype, is_end, is_repaired);
      breakpoint* j4 = new breakpoint(cell_ID, jid++, chr, CHR_LENGTHS[chr], HEAD, haplotype, is_end, is_repaired);
      j3->right_jid = j4->id;
      j4->left_jid = j3->id;
      j3->path_ID = pid;
      j4->path_ID = pid;
      breakpoints[j3->id] = j3;
      breakpoints[j4->id] = j4;
      // interval
      adjacency* adj2 = new adjacency(cell_ID, aid++, pid, j3->id, j4->id, INTERVAL, COMPLETE, NONE);
      adjacencies[adj2->id] = adj2;
      // path
      path* p2 = new path(pid++, cell_ID, COMPLETE);
      p2->nodes.push_back(j3->id);
      p2->nodes.push_back(j4->id);
      p2->edges.push_back(adj2->id);
      // validate_path(p2);
      paths[p2->id] = p2;
    }

    // get_breakpoint_map();
  }


  genome(const genome& _g2) {

  }

  /*********************** functions related to initialization **************************/ 
   
  // only connect unrepaired DSBs later in the repair stage
  void set_initial_bp_status(breakpoint* j1, breakpoint* j2, int& aid, int chr, int haplotype, int bp, int side, map<pair<int, int>, vector<int>>& chr_bps, map<pos_hap, int>& bp_jids, int verbose = 0){
    pair<int, int> chp(chr, haplotype);
    insert_sorted_vec(chr_bps[chp], bp);
    // find the postions of the two breakpoints to the left and right of bp in chr_bps
    pair<int, int> adj_bps = find_adjacent_values(chr_bps[chp], bp);
    // assert(adj_bps.first < bp && adj_bps.second > bp);
    if(adj_bps.first >= bp || adj_bps.second <= bp){
      cout << "wrong neighbors for breakpoint " << bp << " on chr " << chr + 1 << " haplotype " << haplotype << ": " << adj_bps.first << ", " << adj_bps.second << endl;
      print_bp_id(chr_bps, bp_jids);
      exit(FAIL);
    }

    pos_hap a1 = {chr, haplotype, adj_bps.first, TAIL}; 
    assert(bp_jids.find(a1) != bp_jids.end());   
    int j1l_id = bp_jids[a1];
    breakpoint* j1l = breakpoints[j1l_id];     
    if(verbose > 1){
      cout << "breakpoint at the left of " << bp << " is " << chr + 1 << " " << haplotype << " " << adj_bps.first << " " << TAIL << " with ID " << j1l_id << endl;
      j1l->print();
    }
    
    pos_hap a2 = {chr, haplotype, adj_bps.second, HEAD};  
    assert(bp_jids.find(a2) != bp_jids.end());   
    int j2r_id = bp_jids[a2];
    breakpoint* j2r = breakpoints[j2r_id];     
    if(verbose > 1){
      cout << "breakpoint at the right of " << bp << " is "  << chr + 1 << " " << haplotype << " " << adj_bps.second << " " << HEAD << " with ID " << j2r_id << endl;
      j2r->print();
    } 
      
    // add adjacency later
    if(side == HEAD){
        j1->is_repaired = true;
    }else{
        j2->is_repaired = true;
        j1->pos--;
        j2->pos--;
    }  

    breakpoints[j1->id] = j1;
    breakpoints[j2->id] = j2;    

    // add the bp pair later to avoid confusing neighbour finding
    insert_sorted_vec(chr_bps[chp], bp + 1);

    pos_hap p1 = {chr, haplotype, bp, HEAD};
    bp_jids[p1] = j1->id;
    pos_hap p2 = {chr, haplotype, bp + 1, TAIL};
    bp_jids[p2] = j2->id;

    // add the new interval adjacency     
    update_adjacency(j1, j2, j1l, j2r, aid, false, verbose);
  }


  // generate DSBs and repair according to known events
  // need to track all breakpoints on a chromosome to link them
  void intialize_with_svs(const vector<pos_sv>& svs_common, int verbose = 0){ 
    map<pair<int, int>, vector<int>> chr_bps;   
    map<pos_hap, int> bp_jids;    
    intialize_chr_bp_jid(chr_bps, bp_jids);   

    int jid = breakpoints.rbegin()->first + 1;
    int aid = adjacencies.rbegin()->first + 1;

    if(verbose > 1){
      cout << jid << " breakpoints before introducing " << svs_common.size() << " initial DSBs" << endl;
    }

    bool is_end = false;
    bool is_repaired = false;

    for(auto pos_sv : svs_common){
      if(verbose > 1){
        cout << pos_sv.chr1 << "\t" << pos_sv.loc1 << "\t" << pos_sv.side1 << "\t" << pos_sv.chr2 << "\t" << pos_sv.loc2 << "\t" << pos_sv.side2 << endl;
      }

      // add first DSB
      int bp = pos_sv.loc1;
      int chr = pos_sv.chr1;
      int side = pos_sv.side1; 
      int haplotype = myrng(2);
      breakpoint* j1 = new breakpoint(cell_ID, jid++, chr, bp, HEAD, haplotype, is_end, is_repaired);
      breakpoint* j2 = new breakpoint(cell_ID, jid++, chr, bp + 1, TAIL, haplotype, is_end, is_repaired);
      set_initial_bp_status(j1, j2, aid, chr, haplotype, bp, side, chr_bps, bp_jids, verbose);

      // add second DSB
      bp = pos_sv.loc2;
      chr = pos_sv.chr2;
      side = pos_sv.side2;
      haplotype = myrng(2);
      breakpoint* j3 = new breakpoint(cell_ID, jid++, chr, bp, HEAD, haplotype, is_end, is_repaired);
      breakpoint* j4 = new breakpoint(cell_ID, jid++, chr, bp + 1, TAIL, haplotype, is_end, is_repaired);
      set_initial_bp_status(j3, j4, aid, chr, haplotype, bp, side, chr_bps, bp_jids, verbose);

      if(verbose > 1) printf("adding the known adjacency\n");
      // add adjacencies based on the sides
      if(pos_sv.side1 == HEAD && pos_sv.side2 == HEAD){
        add_new_adjacency(j1, j3, aid, verbose);
        assert(j1->is_repaired == true && j3->is_repaired == true);
      }else if(pos_sv.side1 == HEAD && pos_sv.side2 == TAIL){
        add_new_adjacency(j1, j4, aid, verbose);
        assert(j1->is_repaired == true && j4->is_repaired == true);
      }else if(pos_sv.side1 == TAIL && pos_sv.side2 == TAIL){
        add_new_adjacency(j2, j4, aid, verbose);
        assert(j2->is_repaired == true && j4->is_repaired == true);
      }else{
        assert(pos_sv.side1 == TAIL && pos_sv.side2 == HEAD);
        add_new_adjacency(j2, j3, aid, verbose);
        assert(j2->is_repaired == true && j3->is_repaired == true);
      }      
    }
  }

  /*********************** functions related to output **************************/ 
  void print(){
    cout << breakpoints.size() << " breakpoints in the genome: " << endl;
    for(auto j : breakpoints){
      (j.second)->print();
    }

    // cout << "breakpoints in the genome (grouped by haplotype and chr): " << endl;
    // for(auto jm : junc_map){
    //   cout << "haplotype " << jm.first.first << " chr " << jm.first.second << ": ";
    //   for(auto j : jm.second){
    //     j->print();
    //   }
    // }

    cout << adjacencies.size() << " adjacencies in the genome: " << endl;
    for(auto a : adjacencies){
      (a.second)->print();
    }

    cout << paths.size() << " paths in the genome: " << endl;
    for(auto p : paths){
      p.second->print();
    }
  }


  void print_full_path(path* p){
    p->print();
    write_path(p, cout);
    for(auto n : p->nodes){
      breakpoints[n]->print();
    }
    for(auto e : p->edges){
      adjacencies[e]->print();
    }      
  }


  void write_path(path* p, ostream& fout){
    // cout << "writing path" << endl;
    string shape = "linear";
    if(p->is_circle){
      shape = "circular";
    }

    fout << p->id + 1 << "\t" << shape << "\t" << get_telomere_type_string(p->type) << "\t" << p->n_centromere << "\t";

    int gsize = 0;
    for(auto e : p->edges){
      // cout << "\nedge " << e << endl;
      adjacency *adj = adjacencies[e];
      // adj->print();
      breakpoint *j1 = breakpoints[adj->junc_id1];
      breakpoint *j2 = breakpoints[adj->junc_id2];
      if(adj->is_inverted){
        j1 = breakpoints[adj->junc_id2];
        j2 = breakpoints[adj->junc_id1];
      }
      if(adj->type == 0){  // only consider intervals
          gsize += abs(j2->pos - j1->pos) + 1;
      }     
      fout << (j1->chr % NUM_CHR) + 1 << ":" << get_haplotype_string(j1->haplotype) << ":" << j1->pos << "-" << (j2->chr % NUM_CHR) + 1 << ":" << get_haplotype_string(j2->haplotype) << ":" << j2->pos << ",";
    }
    fout << "\t" << gsize << endl;
  }


   /*********************** functions related to finding IDs **************************/ 

  // Find the maximum values of breakpoint IDs
  int find_max_bpID(){
    vector<int> bp_IDs(breakpoints.size(), 0);
    int i = 0;

    for(auto bp : breakpoints){
      bp_IDs[i++] = bp.first;
    }

    int max_bpID = *max_element(bp_IDs.begin(), bp_IDs.end());

    return max_bpID;
  }

  // Find the maximum values of adjacency IDs
  int find_max_adjID(){
    vector<int> adj_IDs(adjacencies.size(), 0);
    int i = 0;

    for(auto adj : adjacencies){
      adj_IDs[i++] = adj.first;
    }

    int max_adjID = *max_element(adj_IDs.begin(), adj_IDs.end());

    return max_adjID;
  }

  // Find the maximum values of path IDs
  int find_max_pathID(){
    vector<int> path_IDs(paths.size(), 0);
    int i = 0;

    for(auto p : paths){
      path_IDs[i++] = p.first;
    }

    int max_pID = *max_element(path_IDs.begin(), path_IDs.end());

    return max_pID;
  }


  // // find breakpoints by haplotype and chr
  // void get_breakpoint_map(){
  //   for(auto jm : breakpoints){
  //     breakpoint* j = jm.second;
  //     int haplotype = j->haplotype;
  //     int chr = j->chr;
  //     pair<int, int> key(chr, haplotype);
  //     junc_map[key].push_back(j);
  //   }
  // }


/********************* functions related to DSB introduction and repair *******************/


// This function merges the rearranged genome to get all the unique genomic intervals
// used to find feasible regions for introducing breakpoints and computing copy numbers
// assume adjacencies contain all edges
// when intervals are at circular DNAs, it may cause breaks on ecDNA or ecDNA structure change. 
void get_unique_interval(int verbose = 0){
  cn_by_chr_hap.clear();

  if(verbose > 1) cout << "Collecting all the unique genomic intervals\n";
  // get CN for each interval
  // chr, start, end, haplotype : cn
  map<haplotype_pos, int> cn_by_pos;  // default value is 0
  for(auto am : adjacencies){
    adjacency* a = am.second;
    assert(a != NULL);
    if(a->type == INTERVAL){
      breakpoint* j1 = breakpoints[a->junc_id1];
      breakpoint* j2 = breakpoints[a->junc_id2];
      if(verbose > 2){
        a->print();
        j1->print();
        j2->print();
      }

      int chr = j1->chr;
      int haplotype = j1->haplotype;
      assert(j1->chr == j2->chr && j1->haplotype == j2->haplotype);
      int start = j1->pos;
      int end = j2->pos;
      int jid_start = j1->id;
      int jid_end = j2->id;
      haplotype_pos key{chr, haplotype, start, end, jid_start, jid_end};
      cn_by_pos[key]++;
    }
  }

  if(verbose > 1){
    cout << "\nOriginal interval and copy numbers after cell cycle" << endl;
    cout << "cell\tchr\thaplotype\tstart\tend\tJID_start\tJID_end\tCN\n";
    for(auto cnp : cn_by_pos){
      haplotype_pos key = cnp.first;
      int cn = cnp.second;
      int chr = key.chr;
      int haplotype = key.haplotype;
      int start = key.start;
      int end = key.end;
      int jid_start = key.jid_start;
      int jid_end = key.jid_end;
      cout << cell_ID << "\t" << chr + 1  << "\t" << get_haplotype_string(haplotype) << "\t" << start << "\t" << end << "\t" << jid_start << "\t" << jid_end << "\t" << cn << "\n";
    }
  }

  for(auto cnp : cn_by_pos){
    haplotype_pos key = cnp.first;
    int cn = cnp.second;
    int chr = key.chr;
    int haplotype = key.haplotype;
    int start = key.start;
    int end = key.end;
    int jid_start = key.jid_start;
    int jid_end = key.jid_end;
    pair<int, int> key1(chr, haplotype);
    interval intl{start, end, cn, jid_start, jid_end};
    cn_by_chr_hap[key1].push_back(intl);
  }

  if(verbose > 1){
    cout << "\nUnique intervals by chr and haplotype" << endl;
    cout << "cell\tchr\thaplotype\tstart\tend\tJID_start\tJID_end\tCN\n";
    for(auto cnp : cn_by_chr_hap){
      int chr = cnp.first.first;
      int haplotype = cnp.first.second;
      for(auto intl : cnp.second){
        int start = intl.start;
        int end = intl.end;
        int cn = intl.cn;
        int jid_start = intl.jid_start;
        int jid_end = intl.jid_end;
        cout << cell_ID << "\t" << chr + 1 << "\t" << get_haplotype_string(haplotype) << "\t" << start << "\t" << end << "\t" << jid_start << "\t" << jid_end << "\t" << cn << "\n";
      }
    }
  }
}


  // assume all the breakpoints are known, sample without replacement when bins is not empty
  void get_random_bp(vector<pos_bp>& bps, vector<double>& bp_fracs, int& chr, int& bp, int& haplotype, int& left_jid, int& right_jid, int verbose = 0){
      bool is_insertable;
      // keep centromere and telomere intact for simiplicity
      do{
        is_insertable = true;
        pos_bp bp_sel;
        if(bps.size() > 0){
          if(verbose > 1){
            cout << "sampling from known breakpoints (by their frequency)" << endl;
          }
          // random_shuffle(bps.begin(), bps.end(), myrng);
          // bp_sel = bps.back();
          gsl_ran_discrete_t* dis_bp = gsl_ran_discrete_preproc(bp_fracs.size(), &bp_fracs[0]);
          int idx = gsl_ran_discrete(r, dis_bp);
          bp_sel = bps[idx];
          chr = bp_sel.chr;
          bps.pop_back();
        }else{
          if(verbose > 1){
            cout << "random sampling on the genome" << endl;
          }
          gsl_ran_discrete_t* dis_loc = gsl_ran_discrete_preproc(NUM_CHR, CHR_PROBS);
          chr = gsl_ran_discrete(r, dis_loc);
        }

        // each chr has two haplotypes, which haplotype the breakpoint is on
        haplotype = myrng(2);
        // cout << " select chr " << chr << " haplotype " << haplotype << endl;
        pair<int, int> key(chr, haplotype);

        if(cn_by_chr_hap.find(key) == cn_by_chr_hap.end()){
          is_insertable = false;
          // cout << "not in genome" << endl;
          continue;
        }

        vector<interval> intls = cn_by_chr_hap[key];
        int nintl = intls.size();
        if(verbose > 1) cout << nintl << " intervals for chr " << chr << " haplotype " << haplotype << endl;
        if(nintl == 0){
          is_insertable = false;
          continue;
        }

        // randomly select an interval
        int sel_intl = (int)runiform(r, 0, nintl);
        interval intl = intls[sel_intl];
        if(bps.size() > 0){
          bp = bp_sel.loc;
          if(bp < intl.start || bp > intl.end){
            is_insertable = false;
            continue;
          }
        }else{
          bp = (int)runiform(r, intl.start, intl.end);
          // cout << "Not enough breakpoints to sample!\n";
          // exit(FAIL);
        }

        left_jid = intl.jid_start;
        right_jid = intl.jid_end;
        // cout << "select interval "  << sel_intl << ": " << intl.start << ", " << intl.end  << " at bp " << bp << " with left breakpoint " << left_jid << " and right breakpoint " << right_jid << endl;

      }while(!is_insertable || bp <= TELO_ENDS1[chr] || bp >= TELO_ENDS2[chr] || (bp >= CENT_STARTS[chr] && bp <= CENT_ENDS[chr]));
  }

  // This function introduces a new breakpoint at position bp and bp+1 and updates the relevant adjacencies
  // jid and aid are passed by reference to automatically get the new IDs for next new breakpoints
  void add_new_breakpoint(int& jid, int &aid, int chr, int bp, int haplotype, int j1l_id = -1, int j2r_id = -1, bool reset_pathID = false, int verbose = 0){
    bool is_end = false;
    // a new breakpoint will generate two new segment j1, j2
    // (j1, j2) still unconnected, need to be repaired
    bool is_repaired = false;
    // two breakpoints at the same? position are generated, with different sides (orientation)
    breakpoint* j1 = new breakpoint(cell_ID, jid++, chr, bp, HEAD, haplotype, is_end, is_repaired);
    breakpoint* j2 = new breakpoint(cell_ID, jid++, chr, bp + 1, TAIL, haplotype, is_end, is_repaired);
    // only connect DSBs later in the repair stage
    breakpoints[j1->id] = j1;
    breakpoints[j2->id] = j2;

    // add new smaller intervals (j1l, j1), (j2, j2r)
    assert(j1l_id != -1);
    assert(j2r_id != -1);
    breakpoint* j1l = breakpoints[j1l_id];
    breakpoint* j2r = breakpoints[j2r_id];
    assert(j1l != NULL);
    assert(j2r != NULL);

    if(verbose > 1){
      cout << "\nadding a new breakpoint " << j1->id << "\t" << j2->id << " at chr " << chr + 1 << " position " << bp << " haplotype " << get_haplotype_string(haplotype) << " between breakpoint " << j1l_id << " and " << j2r_id << endl;
    }

    // the other side of j1l and j2r are unaffected
    // record intervals
    // pair<int, int> key(chr, haplotype);
    // junc_map[key].push_back(j1);
    // junc_map[key].push_back(j2);
    //
    // vector<breakpoint*> juncs = junc_map[key];
    // sort(juncs.begin(), juncs.end(), compare_bp_pos);
    update_adjacency(j1, j2, j1l, j2r, aid, false, verbose);

    // remove from paths
    // cout << "adjacencies after removing j1l to j2r" << endl;
    // for(auto a : adjacencies){
    //   a->print();
    // if(update_path){  // remove and add relevant adjacencies in paths
    //   assert(aid > 0);
    //   paths.edges.push_back(adj1)
    // }
    // }
  }


  // This function imitates DSB by introducing breakpoints (breakpoints) across the whole genome, following infinite sites assumption
  // a breakpoint will not disappear once exist
  // n_dsb: number of breakpoints for each chromosome
  // TODO: add interval DSB probability based on overlapping with fragile sites
  void generate_dsb(int n_dsb, vector<pos_bp>& bps, vector<double>& bp_fracs, vector<breakpoint*>& junc2repair, int verbose = 0){
    // group breakpoints by haplotype and chr for sorting
    // get_breakpoint_map();
    if(n_dsb <= 0){
      if(verbose > 1){
        cout << "NO DSBs in this cycle!" << endl;
      }
      return;
    }

    // if(bps.size() > 0 && bps.size() < n_dsb){
    //   cout << n_dsb << " DSBs to introduce, but only " << bps.size() << " DSBs in the provided breakpoint file!" << endl;
    //   exit(FAIL);
    // }

    int jid = breakpoints.rbegin()->first + 1;
    int aid = adjacencies.rbegin()->first + 1;
    if(verbose > 1){
      cout << jid << " breakpoints before introducing DSB" << endl;
    }

    for(int i = 0; i < n_dsb; i++){
      // a breakpoint may be at centromere or telomere
      // a chromosome may lost some regions, only feasible on available segments
      int bp = 0;
      int chr = -1;
      int haplotype = -1;
      int left_jid = -1;
      int right_jid = -1;

       // available intervals in the haplotype after a new dsb
      get_unique_interval(verbose);

      get_random_bp(bps, bp_fracs, chr, bp, haplotype, left_jid, right_jid, verbose);

      // update number of DSBs per chrom per haplotype
      // int idx = chr + haplotype * NUM_CHR;
      // vec_n_dsb[chr] += 1;

      if(verbose > 0){
        cout << "   break " << i + 1 << " at position " << bp << " chr " << chr + 1 << " haplotype " << get_haplotype_string(haplotype) << endl;
        if(verbose > 1) cout << "   with left breakpoint " << left_jid << " and right breakpoint " << right_jid << endl;
      }

      add_new_breakpoint(jid, aid, chr, bp, haplotype, left_jid, right_jid, false, verbose);
      junc2repair.push_back(breakpoints[jid - 2]);
      junc2repair.push_back(breakpoints[jid - 1]);
    }
  }


  // An auxiliary function called by repair_dsb
  // Find another breakpoint to repair one breakpoint based on the type of rejoining: random or based on distance
  // The distance is based on the positions of the breakpoints in the reference genome for simplicity, as all the coordinates are based on the reference.
  // TODO: consider the distance on the rearranged chromosome.
  breakpoint* get_dsb_pair(breakpoint* j1, vector<breakpoint*>& junc2repair, int pair_type = 0, double prob_correct_repaired = 0, int verbose = 0){
    if(verbose > 1){
      cout << "start ";
      j1->print();
    }
    breakpoint* j2 = NULL;
    if(pair_type == 0){ // randomly join two breakpoints
      int j2id = myrng(junc2repair.size());
      j2 = junc2repair[j2id];
    }else{
      // probability is proportional to 1 / distance
      int nstate = junc2repair.size();
      if(nstate == 1){
        return junc2repair[nstate - 1];
      }
      // manual normalization to avoid constant correct repairing
      double *probs = new double[nstate];
      memset(probs, 0.0, nstate);
      long double sum_prob = 0.0;    // sum of probabilities for non-adjacent breakpoints
      int ndist1 = 0;
      for(int i = 0; i < nstate; i++){  // a breakpoint may have multiple copies
        breakpoint* j = junc2repair[i];
        if(verbose > 1) j->print();
        int distance = 0;
        if(j->chr == j1->chr){
          distance = abs(j->pos - j1->pos);
          if(distance == 0){  // duplicates of a breakpoint
            probs[i] = PROB_SELF;
            sum_prob += probs[i];
          }else if(distance == 1){   // adjacent breakpoints
            probs[i] = -1;  // use -1 to avoid float comparison of very small values to 0
            ndist1++;
          }else{
            probs[i] = (double) 1 / (distance);
            sum_prob += probs[i];
          }
        }else{
          probs[i] = PROB_INTER;
          sum_prob += probs[i];
        }
        if(verbose > 1) cout << distance << "\t" << probs[i] << endl;
      }

      if(sum_prob == 0){  // only breakpoints with the adjacent breakpoint, probably from some duplicates
          for(int i = 0; i < nstate; i++){
            probs[i] =  (double) 1 / nstate;
            if(verbose > 1) cout << probs[i] << endl;
          }
      }else{
        for(int i = 0; i < nstate; i++){
          if(probs[i] == -1){
            assert(ndist1 >= 1);
            probs[i] = prob_correct_repaired / ndist1;
          }else{
            if(ndist1 == 0){
              assert(sum_prob > 0);
              probs[i] = probs[i] / sum_prob;
            }else{
              probs[i] = (1 - prob_correct_repaired) * probs[i] / sum_prob;
            }
          }
          if(verbose > 1) cout << probs[i] << endl;
        }
      }
      // if(verbose > 0) cout << sum_prob << "\t" << sum_prob2 << endl;
      // assert(fabs(sum_prob2 - 1) < PROB_INTER);    # the difference may be > 0.01, not matter too much, as probs will be normalized to sum 1 by gsl_ran_discrete_preproc

      // sample based on distance to j1
      gsl_ran_discrete_t* dis = gsl_ran_discrete_preproc(nstate, probs);
      int sel = gsl_ran_discrete(r, dis);
      j2 = junc2repair[sel];
    }
    assert(j2 != NULL);

    return j2;
  }


  // This function repairs DSBs according to the user-specified fraction of unrepaired DSBs
  // n_unrepaired: number of unrepaired DSBs, each DSB introduces two breakpoints
  // junc2repair: all the unrepaired DSBs
  void repair_dsb(int n_unrepaired, vector<breakpoint*>& junc2repair, int pair_type = 0, double prob_correct_repaired = 0, int verbose = 0){
    // verbose = 1;
    assert(n_unrepaired >= 0);
    // store segments without telomeres
    int aid = adjacencies.rbegin()->first + 1;

    if(verbose > 1){
      cout << "#breakpoints to repair: " << junc2repair.size() << endl;
      for(auto j : junc2repair){
        j->print();
      }

      cout << adjacencies.size() << " adjacencies before repairing dsb" << endl;
      for(auto a : adjacencies){
        (a.second)->print();
      }
    }

    // checks how many remaining breaks there are in G1
    // if less than a pre-specified value then G1 will stop
    // each breakpoint has at most one missing connection
    // breakpoints should not be at genome end, at least two breakpoints are needed to form a connection
    while(junc2repair.size() > n_unrepaired * 2 && junc2repair.size() >= 2){
      int j1id = myrng(junc2repair.size());
      breakpoint* j1 = junc2repair[j1id];
      junc2repair.erase(std::remove(junc2repair.begin(), junc2repair.end(), j1), junc2repair.end());
      // cout << "#breakpoints remaining after 1st repair: " << junc2repair.size() << endl;
      j1->is_repaired = true;

      // int j2id = myrng(junc2repair.size());
      // breakpoint* j2 = junc2repair[j2id];
      breakpoint* j2 = get_dsb_pair(j1, junc2repair, pair_type, prob_correct_repaired, verbose);
      junc2repair.erase(std::remove(junc2repair.begin(), junc2repair.end(), j2), junc2repair.end());
      // cout << "#breakpoints remaining after 2nd repair: " << junc2repair.size() << endl;
      j2->is_repaired = true;

      if(verbose > 1){
        cout << "repair breakpoint " << j1->id << "\t" << j2->id << endl;
        j1->print();
        j2->print();
      }

      if(verbose > 0){
        cout << "   connecting breakpoint " << j1->chr + 1 << ":" << get_haplotype_string(j1->haplotype) << ":" << j1->pos << get_side_string(j1->side) << " with " << j2->chr + 1 << ":" << get_haplotype_string(j2->haplotype) << ":" << j2->pos << get_side_string(j2->side) << endl;
      }

      add_new_adjacency(j1, j2, aid, verbose);  

    } // while

    if(verbose > 1){
      cout << breakpoints.size() << " breakpoints after repairing dsb" << endl;
      for(auto j : breakpoints){
        (j.second)->print();
      }
      cout << adjacencies.size() << " adjacencies after repairing dsb" << endl;
      for(auto a : adjacencies){
        (a.second)->print_interval(breakpoints);
      }
      cout << junc2repair.size() << " breakpoints unrepaired" << endl;
      for(auto j : junc2repair){
        j->print();
      }
    }
  }

  /********************* functions related to adjacency operation *******************/

  // add novel variant adjacency between two new breakpoints
  void add_new_adjacency(breakpoint* j1, breakpoint* j2, int& aid, int verbose = 0){
      // each breakpoint node should has degree 2
      adjacency* adj = NULL;

      // j1 and j2 should be in creasing order of coordinates when repairing
      breakpoint* jm = NULL;
      if(j1->chr == j2->chr && j1->pos > j2->pos){
        if(verbose > 1){
          cout << "switch j1 and j2 to have increasing coordinates" << endl;
          j1->print();
          j2->print();
        }

        jm = j1;
        j1 = j2;
        j2 = jm;

        if(verbose > 1){
          cout << "switched j1 and j2 with increasing coordinates" << endl;
          j1->print();
          j2->print();
        }
      }

      SV_type sv_type = NONE;
      // determine direction by breakpoint side
      if(j1->chr == j2->chr && abs(j1->pos - j2->pos) == 1 && j1->haplotype == j2->haplotype && j1->side != j2->side){
        // reference adjacency (should be rare)
        adj = new adjacency(cell_ID, aid++, -1, j1->id, j2->id, REF, NONTEL, NONE);
        adj->set_telomeric_type(breakpoints);
        adjacencies[adj->id] = adj;
        // adj->print();
        if(j1->left_jid == -1){
          j1->left_jid = j2->id;
        }else{
          j1->right_jid = j2->id;
        }
        if(j2->left_jid == -1){
          j2->left_jid = j1->id;
        }else{
          j2->right_jid = j1->id;
        }
        // cout << "before updating sv type count " << chr_type_num[j1->chr][sv_type] << endl;
        // chr_type_num[j1->chr][sv_type] += 1;
        // cout << "after updating sv type count " << chr_type_num[j1->chr][sv_type] << endl;
      }else{
        // variant adjacencies
        sv_type = set_var_type(j1, j2, verbose);
        adj = new adjacency(cell_ID, aid++, -1, j1->id, j2->id, VAR, NONTEL, sv_type);
        adj->set_telomeric_type(breakpoints);
        adjacencies[adj->id] = adj;
      } // else
      if(verbose > 1){
          cout << "adding new adjacency " << aid << endl;
          adj->print();
          j1->print();
          j2->print();
      } 
  }


  // remove adjacency connected by j1 and j2
  // update the two resultant adjacencies due to the breakage of the adjacency
  int remove_adjacency(int j1, int j2, adjacency* adj1, adjacency* adj2, int verbose = 0){
    for(std::map<int, adjacency*>::iterator it = adjacencies.begin(); it != adjacencies.end(); ++it){
      adjacency* adj = it->second;
      if((adj->junc_id1 == j1 && adj->junc_id2 == j2) || (adj->junc_id2 == j1 && adj->junc_id1 == j2)){
        int aid = adj->id;
        assert(it->first == aid);

        if(verbose > 1){
          cout << "remove adjacency " << aid << ": " << j1 << "-" << j2 << endl;
          adj->print();
        }
        adj1->is_inverted = adj->is_inverted;
        adj2->is_inverted = adj->is_inverted;

        delete adj;
        adjacencies.erase(it);

        return aid;
      }
    }
    return -1;
  }


  void remove_adjacency_by_ID(int aid, int verbose = 0){
    adjacency* adj = adjacencies[aid];
    if(verbose > 1){
      cout << "remove adjacency " << aid << endl;
      adj->print();
    }
    delete adj;
    adjacencies.erase(aid);
  }


  void update_adjacency(breakpoint* j1, breakpoint* j2, breakpoint* j1l, breakpoint* j2r, int& aid, bool reset_pathID = false, int verbose = 0){
    if(reset_pathID){
      j1->path_ID = -1;
      j2->path_ID = -1;
    }else{
      j1->path_ID = j1l->path_ID;
      j2->path_ID = j2r->path_ID;
    }

    adjacency* adj1 = new adjacency(cell_ID, aid++, j1->path_ID, j1l->id, j1->id, INTERVAL, NONTEL, NONE);
    adj1->set_telomeric_type(breakpoints);
    adjacencies[adj1->id] = adj1;
    j1l->right_jid = j1->id;
    j1->left_jid = j1l->id;
    // j1->right_jid = -1;
    
    adjacency* adj2 = new adjacency(cell_ID, aid++, j2->path_ID, j2->id, j2r->id, INTERVAL, NONTEL, NONE);
    adj2->set_telomeric_type(breakpoints);
    adjacencies[adj2->id] = adj2;
    j2->right_jid = j2r->id;
    j2->left_jid = -1;
    j2r->left_jid = j2->id;

    // remove old larger intervals
    // cout << "adjacencies before removing j1l to j2r" << endl;
    // for(auto a : adjacencies){
    //   a->print();
    // }
    remove_adjacency(j1l->id, j2r->id, adj1, adj2, verbose);  

    if(verbose > 1){
      cout << "\nupdated four breakpoints affected by a new DSB (breakpoint at the left of j1, j1, j2, breakpoint at the right of j2) " << endl;
      j1l->print();
      j1->print();      
      j2->print();
      j2r->print();
      cout << "\nadded two adjacencies introduced by break at j1, j2" << endl;
      adj1->print();
      adj2->print();
      cout << endl;
    }
  }


  /********************* functions related to path reconstruction *******************/
  
  // used when constructing a new path starting from a specific breakpoint
  void set_path_type(path* p, bool set_circle = true, int verbose = 0){
    // verbose = 0;
    int size = p->nodes.size();
    int start = p->nodes[0];
    int end = p->nodes[size - 1];

     // end nodes should be telomere
    if(breakpoints[start]->is_end && breakpoints[end]->is_end){
      p->type = COMPLETE;
    }else if(breakpoints[start]->is_end){  // pTel
      p->type = PTEL;
    }else if(breakpoints[end]->is_end){  // qTel
      p->type = QTEL;
    }else{     // nonTel
      p->type = NONTEL;
      assert(!p->is_circle);
      assert(p->nodes.size() > 1);
      // set to circular if there is no centromere
      if(p->n_centromere == 0 && set_circle){
        set_circle_by_edge(p, start, end, verbose);
      }
    }
    if(verbose > 1) p->print();
  }


  // a connected path of breakpoints starting from js
  // is_forward: true -- left to right (5' to 3', p-tel to q-tel), 1 -- right to left
  // check left breakpoint if the right breakpoint is the same as jp
  // check_interval: check whether the 1st adjacency is an interval
  // circular path may have the start node duplicated, which will be removed later
  void get_connected_breakpoint(breakpoint* js, path* p, map<int, adjacency*>& adjacencies, bool check_interval = false, int verbose = 0){
    // verbose = 1;
    int nei = js->right_jid;
    if(verbose > 1) cout << js->id << " at path " << p->id + 1 << ", traversing from right breakpoint " << nei << endl;

    if(nei == -1){
      nei = js->left_jid;
      if(verbose > 1) cout << js->id << " at path " << p->id + 1 << ", traversing from left breakpoint (right is end) " << nei << endl;
    }

    bool not_interval = false;
    if(check_interval){
      if(verbose > 2){
        cout << "checking whether path from " << js->id << " to " << nei << " is an interval\n";
        for(auto adj : adjacencies){
          adj.second->print();
        }
      }
      vector<int> aids;
      get_adjacency_ID(aids, js->id, nei, adjacencies);
      if(verbose > 1){
        cout << aids.size() << " adjacencies inbetween\n";
        adjacencies[aids[0]]->print();
      }
      if(aids.size() == 1 && adjacencies[aids[0]]->type != INTERVAL){
        if(verbose > 1) cout << "not an interval!\n";
        not_interval = true;
      }
    }

    if(not_interval && nei == js->right_jid){
      nei = js->left_jid;
      if(verbose > 1) cout << js->id << " at path " << p->id + 1 << ", traversing from left breakpoint (right is not an interval)" << nei << endl;
    }

    breakpoint* jp = js;
    int prev_atype = -1;  // record adjacency types to ensure they are alternative
    while(nei != -1){
      if(verbose > 1) cout << "current breakpoint " << nei << endl;
      breakpoint* jn = breakpoints[nei];     
      bool succ = update_path_by_adj(p, jp, jn, adjacencies, prev_atype, verbose);
      if(verbose > 1){
        js->print();
        jn->print();
      }

      if(!succ)  break;   // stop when the path cannot be extended any more

      nei = jn->right_jid;
      if(verbose > 1) cout << jn->id << ", traversing from right breakpoint " << nei << endl;
      vector<int> aids;
      get_adjacency_ID(aids, jn->id, nei, adjacencies);
      if(nei == jp->id && aids.size() == 1){  // right breakpoint has been visited
        nei = jn->left_jid;
        if(verbose > 1) cout << jn->id << ", traversing from left breakpoint (right breakpoint has been visited) " << nei << endl;
      }
      jp = jn;
    }
  }


  // set a path to be a cycle after fragmentation by connecting start and end node directly 
  // only used when the input parameter set_circle is true so that most circles are formed naturally
  void set_circle_by_edge(path* p, int start, int end, int verbose = 0){
    // int start_edge = p->edges[0];
    // int end_edge = p->edges[p->edges.size() - 1];
    // if(adjacencies[start_edge]->junc_id1 == adjacencies[end_edge]->junc_id1 ||
    //   adjacencies[start_edge]->junc_id1 == adjacencies[end_edge]->junc_id2 ||
    //   adjacencies[start_edge]->junc_id2 == adjacencies[end_edge]->junc_id1 ||
    //   adjacencies[start_edge]->junc_id2 == adjacencies[end_edge]->junc_id2){
    //     if(verbose > 0) cout << "set circle obtained by connecting two different nodes \n";
    //     p->is_circle = true;
    // }
    // verbose = 1;
    if(verbose > 0){
      cout << "\nset circular type for path " << p->id + 1 << endl;
    }

    if((!(breakpoints[start]->left_jid == -1 || breakpoints[start]->right_jid == -1))||
    (!(breakpoints[end]->left_jid == -1 || breakpoints[end]->right_jid == -1))){
        cout << "The end node must be disconnected at one side !" << endl;
        print_full_path(p);
        exit(FAIL);
    }
    // assert(!breakpoints[start]->is_repaired);   // is_repaired may not be set properly
    if(breakpoints[start]->left_jid == -1){
      breakpoints[start]->left_jid = end;
    }else{
      breakpoints[start]->right_jid = end;
    }
    // connect end node to start node
    // assert(!breakpoints[end]->is_repaired);
    if(breakpoints[end]->left_jid == -1){
      breakpoints[end]->left_jid = start;
    }else{
      breakpoints[end]->right_jid = start;
    }
    
    // add another adjacency from end to start
    breakpoint* j1 = breakpoints[end];
    breakpoint* j2 = breakpoints[start];
    breakpoint* jm = NULL;
    // ensure the coordinates of j1 is smaller than j2 to set SV type
    if(j1->chr == j2->chr && j1->pos > j2->pos){
      // if(verbose > 1){
      //   cout << "to switch j1 and j2 to have increasing coordinates" << endl;
      //   j1->print();
      //   j2->print();
      // }
      jm = j1;
      j1 = j2;
      j2 = jm;
      // if(verbose > 1){
      //   cout << "switched j1 and j2 with increasing coordinates" << endl;
      //   j1->print();
      //   j2->print();
      // }
    }

    SV_type sv_type = set_var_type(j1, j2, verbose);
    int aid = adjacencies.rbegin()->first + 1;   // id of the new edge    
    adjacency* adj = new adjacency(cell_ID, aid, p->id, end, start, VAR, NONTEL, sv_type);
    bool is_inverted = false;
    if(p->edges.size() == 1){
      is_inverted = !adjacencies[p->edges[0]]->is_inverted;
    }
    adj->is_inverted = is_inverted;  // from end to start
    assert(adjacencies.find(aid) == adjacencies.end());
    adjacencies[aid] = adj;

    // p->nodes.push_back(start);
    p->edges.push_back(aid);
    assert(p->edges.size() == p->nodes.size());

    breakpoints[end]->is_repaired = true;
    breakpoints[start]->is_repaired = true;
    p->is_circle = true;

    if(verbose > 0){
      cout << "set circle manually by adding an adjacency" << endl;
      print_full_path(p);
    }
  }


  // # defines paths that start and end on a telomere prior to S phase
  // # i.e. these paths are fully connected upon completion of G1
  // circular_prob used for artificial circles as the path will be recomputed in each cycle whatever it is derived
  void get_derivative_genome(double circular_prob, int verbose = 0){
    // verbose = 1;
    // remove previous path connections
    paths.clear();
    for(auto jm : breakpoints){
      (jm.second)->path_ID = -1;
    }
    for(auto am : adjacencies){
      (am.second)->path_ID = -1;
    }
    int pid = 0;

    // Find all breakpoints at the left telomere (pTel)
    vector<breakpoint*> junc_starts;
    for(auto jm : breakpoints){
      breakpoint* j = jm.second;
      if(j->pos == 0 && j->is_end){
        junc_starts.push_back(j);
      }
    }

    if(verbose > 1){
      cout << "\nbreakpoints in the genome:" << endl;
      for(auto jm: breakpoints){
        breakpoint* j = jm.second;
        j->print();
      }
      cout << "\nbreakpoints at the left telomere:";
      for(auto js : junc_starts){
        cout << " " << js->id;
      }
      cout << endl;
    }

    // a path must alternate between interval and variant/reference adjacency
    for(auto js : junc_starts){  // 1st region must be interval when starting from telomere
      if(js->path_ID >= 0) continue;   // some nodes may be already in a path

      path* p = new path(pid++, cell_ID, NONTEL);
      js->path_ID = p->id;
      p->nodes.push_back(js->id);
      if(verbose > 1) cout << "\nadding a new path " << p->id + 1 << endl;

      get_connected_breakpoint(js, p, adjacencies, false, verbose);

      int size = p->nodes.size();
      if(breakpoints[p->nodes[0]]->id == breakpoints[p->nodes[size-1]]->id){
        if(verbose > 0) cout << "set circle to itself\n";
        p->nodes.pop_back();
        p->is_circle = true;
      }
      // no need as this is assumed to be formed during fragmentation and the path will remain complete in the subsequent cycles
      // if(!p->is_circle && circular_prob > 0 && p->edges.size() > 1){ 
      //   set_circle_by_edge(p, verbose);
      // }
      if(breakpoints[p->nodes[size-1]]->is_end){
        p->type = COMPLETE; // end nodes should also be telomere
      }else{ 
        p->type = PTEL; // pTel
      }

      if(!validate_path(p)){
        cout << "Wrong path from junc_starts breakpoint " << js->id << endl;
        exit(FAIL);
      }

      if(verbose > 1){
        cout << "path after validation\n";
        p->print();
      }
      paths[p->id] = p;
    }

    if(verbose > 1){
      cout << "\npath starting at the left telomere " << endl;
      for(auto p : paths){
        p.second->print();
      }
    }

    // store unconnected paths for parsing to S phase; should be qTel
    vector<breakpoint*> junc_unconnected_lends;
    for(auto jm : breakpoints){
      breakpoint* j = jm.second;
      if(j->path_ID < 0 && j->pos > 0 && j->is_end){
        junc_unconnected_lends.push_back(j);
      }
    }

    if(verbose > 1){
      cout << "\nqTel breakpoints not fully connected:";
      for(auto js : junc_unconnected_lends){
        cout << " " << js->id;
      }
      cout << endl;
    }

    for(auto js : junc_unconnected_lends){
      if(js->path_ID >= 0) continue;   // some nodes may be already in a path

      path* p = new path(pid++, cell_ID, NONTEL);
      js->path_ID = p->id;
      p->nodes.push_back(js->id);
      if(verbose > 1) cout << "\nadding a new path " << p->id + 1 << endl;

      get_connected_breakpoint(js, p, adjacencies, false, verbose);

      // end nodes should be telomere
      int size = p->nodes.size();
      if(breakpoints[p->nodes[0]]->id == breakpoints[p->nodes[size-1]]->id){
        if(verbose > 0) cout << "set circle to itself\n";
        p->nodes.pop_back();
        p->is_circle = true;
      }
      // if(!p->is_circle && circular_prob > 0 && p->edges.size() > 1){ 
      //   set_circle_by_edge(p, verbose);
      // }      
      if(breakpoints[p->nodes[size-1]]->is_end){
        p->type = COMPLETE;
      }else{
        p->type = QTEL;
      }

      if(!validate_path(p)){
        cout << "Wrong path from junc_unconnected_lends breakpoint " << js->id << endl;
        exit(FAIL);
      }
      if(verbose > 1){
        cout << "path after validation\n";
        p->print();
      }
      paths[p->id] = p;
    }

    if(verbose > 1){
      cout << "\npath starting at the right telomere " << endl;
      for(auto p : paths){
        p.second->print();
      }
    }

    vector<breakpoint*> junc_unconnected;
    // may be repaired when inherited from parent
    for(auto jm : breakpoints){
      breakpoint* j = jm.second;
      if(j->path_ID < 0 && (j->left_jid == -1 || j->right_jid == -1)){
        assert(!j->is_end);
        junc_unconnected.push_back(j);
      }
    }

    if(verbose > 1){
      cout << "\nnon-end breakpoints not fully connected:";
      for(auto js : junc_unconnected){
        cout << " " << js->id;
      }
      cout << endl;
      for(auto js : junc_unconnected){
        breakpoints[js->id]->print();
      }
    }

    // an interval may be self-connected after repairing, forming a circle
    for(auto js : junc_unconnected){  // need additional check to ensure 1st region is interval
      if(js->path_ID >= 0) continue;   // some nodes may be already in a path

      path* p = new path(pid++, cell_ID, NONTEL);
      js->path_ID = p->id;
      p->nodes.push_back(js->id);
      if(verbose > 1) cout << "\nadding a new path " << p->id + 1 << endl;

      // check whether the 1st region is interval or not
      get_connected_breakpoint(js, p, adjacencies, true, verbose);

      // end nodes should be telomere
      int size = p->nodes.size();
      if(breakpoints[p->nodes[0]]->id == breakpoints[p->nodes[size-1]]->id){
        if(verbose > 0) cout << "set circle to itself\n";
        p->nodes.pop_back();
        p->is_circle = true;
      }
      // if(!p->is_circle && circular_prob > 0 && p->edges.size() > 1){ 
      //   set_circle_by_edge(p, verbose);
      // }      
      assert(!breakpoints[p->nodes[size-1]]->is_end);

      if(!validate_path(p)){
        cout << "Wrong path from junc_unconnected breakpoint " << js->id << endl;
        exit(FAIL);
      }
      if(verbose > 1){
        cout << "path after validation\n";
        p->print();
      }
      paths[p->id] = p;
    }

    vector<breakpoint*> junc_remained;
    for(auto jm : breakpoints){
      breakpoint* j = jm.second;
      if(j->path_ID < 0){
        assert(!j->is_end);
        junc_remained.push_back(j);
      }
    }

    if(verbose > 1){
      cout << "\nnon-end breakpoints not in a path:";
      for(auto js : junc_remained){
        cout << " " << js->id;
      }
      cout << endl;
      for(auto js : junc_remained){
        breakpoints[js->id]->print();
      }
    }

    for(auto js : junc_remained){  // need additional check to ensure 1st region is interval
      if(js->path_ID >= 0) continue;   // some nodes may be already in a path

      path* p = new path(pid++, cell_ID, NONTEL);
      js->path_ID = p->id;
      p->nodes.push_back(js->id);
      if(verbose > 1) cout << "\nadding a new path " << p->id + 1 << endl;

      // check whether the 1st region is interval or not
      get_connected_breakpoint(js, p, adjacencies, true, verbose);

      // end nodes should be telomere
      int size = p->nodes.size();
      if(breakpoints[p->nodes[0]]->id == breakpoints[p->nodes[size-1]]->id){
        if(verbose > 0) cout << "set circle to itself\n";
        p->nodes.pop_back();
        p->is_circle = true;
      }
      // if(!p->is_circle && circular_prob > 0 && p->edges.size() > 1){ 
      //   set_circle_by_edge(p, verbose);
      // }     
      assert(!breakpoints[p->nodes[size-1]]->is_end);

      if(!validate_path(p)){
        cout << "Wrong path from junc_remained breakpoint " << js->id << endl;
        exit(FAIL);
      }
      if(verbose > 1){
        cout << "path after validation\n";
        p->print();
      }
      paths[p->id] = p;
    }

    check_derived_genome();
  }


// update the other neighbor of last breakpoint in the duplicated path
void update_end_junc(int last_jid_orig, int last_jid) {
  if(!(breakpoints[last_jid_orig]->right_jid == breakpoints[last_jid]->right_jid ||
  breakpoints[last_jid_orig]->left_jid == breakpoints[last_jid]->left_jid)){
    cout << "Weird breakpoint neighbours!" << endl;
    breakpoints[last_jid_orig]->print();
    breakpoints[last_jid]->print();
    exit(FAIL);
  }

  assert(breakpoints[last_jid_orig]->right_jid == breakpoints[last_jid]->right_jid ||
  breakpoints[last_jid_orig]->left_jid == breakpoints[last_jid]->left_jid);

  if(breakpoints[last_jid_orig]->right_jid == breakpoints[last_jid]->right_jid){
    assert(breakpoints[last_jid_orig]->right_jid == -1);
    breakpoints[last_jid_orig]->right_jid = last_jid;
    breakpoints[last_jid]->right_jid = last_jid_orig;
  }else{
    assert(breakpoints[last_jid_orig]->left_jid == -1);
    breakpoints[last_jid_orig]->left_jid = last_jid;
    breakpoints[last_jid]->left_jid = last_jid_orig;
  }
}


// when invert_adj is true, it is used for fusion path and hence copied adjacency has direction inverted from original
// return the last copied breakpoint to form path fusion
breakpoint* duplicate_path(path& p, vector<int>& edges_copy, vector<int>& nodes_copy, bool invert_adj = false, int verbose = 0){
  // assert(p.nodes.size() == p.edges.size() + 1);
  if(verbose > 1){
    cout << "\nduplicate path " << p.id + 1 << endl;
    p.print();
  }
  assert(p.nodes.size() % 2 == 0);

  int curr_eid = 0;
  adjacency* adj_copy_prev = NULL;
  adjacency* adj_prev = NULL;
  breakpoint* junc_copy_prev = NULL;
  breakpoint* junc_prev = NULL;
  map<int, int> junc_id_map;  // a map of breakpoint ids to link copied breakpoints together
  // duplicate by edge (adjacency) to get node orders
  // nodes in forward order, adjacency of interval may be inverted
  // duplicated adjacency will have direction inverted from original direction
  for(int i = 0; i < p.edges.size(); i++){
    curr_eid = p.edges[i];

    adjacency* adj = adjacencies[curr_eid];
    adjacency* adj_copy = new adjacency(*adj);
    adj_copy->id = adjacencies.rbegin()->first + 1;
    adjacencies[adj_copy->id] = adj_copy;

    // duplicate node (breakpoint) and update adj_copy
    assert(p.nodes[i] == adj->junc_id1 || p.nodes[i] == adj->junc_id2);
    breakpoint* bp = breakpoints[p.nodes[i]];
    breakpoint* junc_copy = new breakpoint(*bp);
    junc_copy->id = breakpoints.rbegin()->first + 1;
    breakpoints[junc_copy->id] = junc_copy;
    junc_id_map[bp->id] = junc_copy->id;

    if(verbose > 1){
      cout << "copy adjacency " << adj->id << endl;
      adj->print();
      adj_copy->print();
      cout << "copy breakpoint " << bp->id << endl;
      bp->print();
      junc_copy->print();
    }

    if(junc_copy_prev != NULL){
      update_adj_junc(adj_copy_prev, junc_copy_prev, junc_copy, verbose);
      if(verbose > 1){
        cout << "adjacency " << adj_copy_prev->id << " after updating junctions" << endl;
        adj_copy_prev->print();
      }
      // set direction of copied adjacency after junction updating based on original direction
      if(invert_adj){
        if(verbose > 1){
          cout << "Invert copied path for fusion" << endl;
        }
        adj_copy_prev->is_inverted = !(adj_prev->is_inverted);
      }
    }

    edges_copy.push_back(adj_copy->id);
    nodes_copy.push_back(junc_copy->id);

    junc_copy_prev = junc_copy;
    adj_copy_prev = adj_copy;
    junc_prev = bp;
    adj_prev = adj;

  }

  assert(curr_eid >= 0);
  // duplicate last node (breakpoint) when the path is not circular, since #node = #edge + 1
  breakpoint* junc_copy = NULL;
  if(!p.is_circle){
    breakpoint* bp = breakpoints[p.nodes[p.nodes.size()-1]];
    junc_copy = new breakpoint(*bp);
    junc_copy->id = breakpoints.rbegin()->first + 1;
    breakpoints[junc_copy->id] = junc_copy;
    junc_id_map[bp->id] = junc_copy->id;
    nodes_copy.push_back(junc_copy->id);

    if(verbose > 1){
      cout << "copy breakpoint " << bp->id << endl;
      bp->print();
      junc_copy->print();
    }
  }else{
    assert(nodes_copy.size() == edges_copy.size());
    junc_copy = breakpoints[nodes_copy[0]];
  }

  // the neighbor updated depends on direction
  update_adj_junc(adj_copy_prev, junc_copy_prev, junc_copy, verbose);
  if(verbose > 1){
    cout << "adjacency " << adj_copy_prev->id << " after updating junctions" << endl;
    adj_copy_prev->print();
  }
  if(invert_adj){
    adj_copy_prev->is_inverted = !(adj_prev->is_inverted);
  }

  int nbp = nodes_copy.size();
  for(int i = 0; i < nbp; i++){
    int jid = nodes_copy[i];
    if(verbose > 1){
      cout << "update neighbor of node " << jid << endl;
    }
    if(breakpoints[jid]->left_jid >= 0){
      breakpoints[jid]->left_jid = junc_id_map[breakpoints[jid]->left_jid];
    }
    if(breakpoints[jid]->right_jid >= 0){
      breakpoints[jid]->right_jid = junc_id_map[breakpoints[jid]->right_jid];
    }
  }

  if(verbose > 1){
    cout << "all copied nodes: " << endl;
    for(int i = 0; i < edges_copy.size(); i++){
      int aid = edges_copy[i];
      adjacency* a = adjacencies[aid];
      a->print();
    }
    for(int i = 0; i < nodes_copy.size(); i++){
      int jid = nodes_copy[i];
      breakpoint* j= breakpoints[jid];
      j->print();
    }
  }

  return junc_copy;
}

// assume p->type != COMPLETE
// duplicate a non-circular incomplete path p and related objects (breakpoints and adjacencies)
void duplicate_path_fusion(path& p, int verbose = 0){
  vector<int> edges_copy;
  vector<int> nodes_copy;
  breakpoint* junc_copy = duplicate_path(p, edges_copy, nodes_copy, true, verbose);

  // merge duplicated path with original path
  p.n_centromere = p.n_centromere * 2;
  if(p.type == PTEL || p.type == NONTEL){  // pTel -- connect last breakpoint
    int last_jid_orig = p.nodes[p.nodes.size() - 1];
    int last_jid = junc_copy->id;   // breakpoint copied as forward

    // add connection between original path and its copy at the end
    int aid = adjacencies.rbegin()->first + 1;
    SV_type sv_type = T2TINV;
    if(breakpoints[last_jid_orig]->side == HEAD && breakpoints[last_jid]->side == HEAD){
      sv_type = H2HINV;
    }
    adjacency* adj_var = new adjacency(cell_ID, aid, p.id, last_jid_orig, last_jid, VAR, NONTEL, sv_type);
    adjacencies[adj_var->id] = adj_var;

    // one neighbor of the copied breakpoint must have pointed to its neighbour copy
    // only the neighbor not updated will be equal to original node's neighbor
    update_end_junc(last_jid_orig, last_jid);

    if((breakpoints[last_jid_orig]->left_jid < 0 || breakpoints[last_jid_orig]->right_jid < 0)){
      cout << "Weird breakpoint neighbours after updating!" << endl;
      breakpoints[last_jid_orig]->print();
      breakpoints[last_jid]->print();
      exit(FAIL);
    }

    assert(breakpoints[last_jid_orig]->left_jid >= 0 && breakpoints[last_jid_orig]->right_jid >= 0);
    assert(breakpoints[last_jid]->left_jid >= 0 && breakpoints[last_jid]->right_jid >= 0);

    breakpoints[last_jid_orig]->is_repaired = true;
    breakpoints[last_jid]->is_repaired = true;

    p.edges.push_back(adj_var->id);
    p.edges.insert(p.edges.end(), edges_copy.rbegin(), edges_copy.rend());
    p.nodes.insert(p.nodes.end(), nodes_copy.rbegin(), nodes_copy.rend());

    if(p.type == NONTEL){  // add one more connection to form a ring
      int first_jid_orig = p.nodes[0];
      int first_jid = nodes_copy[0];

      int aid = adjacencies.rbegin()->first + 1;
      SV_type sv_type = T2TINV;
      if(breakpoints[first_jid_orig]->side == HEAD && breakpoints[first_jid]->side == HEAD){
        sv_type = H2HINV;
      }
      adjacency* adj_var = new adjacency(cell_ID, aid, p.id, first_jid_orig, first_jid, VAR, NONTEL, sv_type);
      adjacencies[adj_var->id] = adj_var;

      update_end_junc(first_jid_orig, first_jid);

      assert(breakpoints[first_jid_orig]->left_jid >= 0 && breakpoints[first_jid_orig]->right_jid >= 0);
      assert(breakpoints[first_jid]->left_jid >= 0 && breakpoints[first_jid]->right_jid >= 0);

      breakpoints[first_jid_orig]->is_repaired = true;
      breakpoints[first_jid]->is_repaired = true;

      p.edges.push_back(adj_var->id);
      if(verbose > 0) cout << "set circle from fusion" << endl;
      p.is_circle = true;
    }
  }else{
    assert(p.type == QTEL);  // qTel -- connect first breakpoint, path start from qTel
    int last_jid_orig = p.nodes[p.nodes.size() - 1];
    int last_jid = junc_copy->id;

    int aid = adjacencies.rbegin()->first + 1;
    SV_type sv_type = T2TINV;
    if(breakpoints[last_jid_orig]->side == HEAD && breakpoints[last_jid]->side == HEAD){
      sv_type = H2HINV;
    }
    adjacency* adj_var = new adjacency(cell_ID, aid, p.id, last_jid_orig, last_jid, VAR, NONTEL, sv_type);
    adjacencies[adj_var->id] = adj_var;

    update_end_junc(last_jid_orig, last_jid);

    assert(breakpoints[last_jid_orig]->left_jid >= 0 && breakpoints[last_jid_orig]->right_jid >= 0);
    assert(breakpoints[last_jid]->left_jid >= 0 && breakpoints[last_jid]->right_jid >= 0);

    breakpoints[last_jid_orig]->is_repaired = true;
    breakpoints[last_jid]->is_repaired = true;

    p.edges.push_back(adj_var->id);
    p.edges.insert(p.edges.end(), edges_copy.rbegin(), edges_copy.rend());
    p.nodes.insert(p.nodes.end(), nodes_copy.rbegin(), nodes_copy.rend());
  }

  if(p.type == PTEL || p.type == QTEL){
    p.type = COMPLETE;
  }

  if(verbose > 1){
    cout << "duplicated path with fusion " << endl;
    p.print();
    for(auto jid : p.nodes){
      breakpoints[jid]->print();
    }
    for(auto aid : p.edges){
      adjacencies[aid]->print();
    }
  }
}



void remove_path(int path_ID, int verbose = 0){
  path* p = paths[path_ID];
  if(verbose > 1){
    cout << "remove path " << path_ID << endl;
    p->print();
  }
  delete p;
  paths.erase(path_ID);
}


/********************************** functions related to summarizing ouput  ***********************************/


// where there are several intervals starting from the same position, some of them may be adjacent to the next one
void get_merged_interval(int verbose = 0){
    cn_by_chr_hap_merged.clear();

    get_unique_interval(verbose);

    if(verbose > 1) cout << "Merging continous intervals on same chr and haplotype" << endl;

    for(auto cnp : cn_by_chr_hap){
      vector<interval> intls = cnp.second;
      sort(intls.begin(), intls.end());

      interval intl_prev = intls[0];
      cn_by_chr_hap_merged[cnp.first].push_back(intl_prev);
      if(verbose > 1) cout << "interval " << 1 << "\t" << intl_prev.start << "\t" << intl_prev.end << "\t" << intl_prev.cn << endl;

      bool merged = false;
      // cout << intls.size() << endl;
      for(int i = 1; i < intls.size(); i++){
        interval intl_curr = intls[i];
        if(verbose > 1) cout << "interval " << i + 1 << "\t" << intl_curr.start << "\t" << intl_curr.end << "\t" << intl_curr.cn << endl;

        // no intervals with the same position
        if(intl_prev.end + 1 == intl_curr.start && intl_prev.cn == intl_curr.cn){
          interval intl_merged{intl_prev.start, intl_curr.end, intl_prev.cn};
          if(verbose > 1) cout << "interval after merging with previous one\t" << intl_merged.start << "\t" << intl_merged.end << "\t" << intl_merged.cn << endl;
          cn_by_chr_hap_merged[cnp.first].pop_back();
          cn_by_chr_hap_merged[cnp.first].push_back(intl_merged);
          intl_prev = intl_merged;
        }else{
          if(verbose > 1) cout << "\t" << intl_curr.start << "\t" << intl_curr.end << "\t" << intl_curr.cn << endl;
          cn_by_chr_hap_merged[cnp.first].push_back(intl_curr);
          intl_prev = intl_curr;
        }
      }
    }

    if(verbose > 1){
      cout << "\nMerged intervals" << endl;
      cout << "cell\tchr\thaplotype\tstart\tend\tJID_start\tJID_end\tCN\n";
      for(auto cnp : cn_by_chr_hap_merged){
        int chr = cnp.first.first;
        int haplotype = cnp.first.second;
        // cout << cell_ID << "\t" << chr << "\t" << haplotype << endl;
        vector<interval> intls = cnp.second;
        for(int i = 0; i < intls.size(); i++){
          interval intl = intls[i];
          cout << cell_ID << "\t" << chr + 1 << "\t" << get_haplotype_string(haplotype) << "\t" << intl.start << "\t" << intl.end  << "\t" << intl.jid_start << "\t" << intl.jid_end << "\t"  << intl.cn << "\n";
        }
      }
    }
  }

// get breakpoints after merging intervals with the same copy number
// TODO: To avoid segments completely lost in one cell not shown, add chr end when needed
void get_bps_per_chr(map<int, set<int>>& bps_by_chr, int verbose = 0){
  get_merged_interval(verbose);

  // split regions on same chr to get total CN (for shatterseek input)
  for(auto cnp : cn_by_chr_hap_merged){
    int chr = cnp.first.first;
    int haplotype = cnp.first.second;
    vector<interval> intls = cnp.second;
    // find all breakpoints
    set<int> bps;
    for(auto intl : cnp.second){
      int start = intl.start;
      int end = intl.end;
      bps_by_chr[chr].insert(start);
      bps_by_chr[chr].insert(end);
    }
  }

  if(verbose > 1){
    cout << "\nbreakpoints for each chr after merging regions with the same CN" << endl;
    for(auto bpc : bps_by_chr){
      int chr = bpc.first;
      cout << chr + 1 << ":";
      for(auto bp : bpc.second){
        cout << " " << bp;
      }
      cout << endl;
    }
  }
}


// a correctly repaired breakpoint will also be reported
void get_bps_per_chr_orig(map<int, set<int>>& bps_by_chr, int verbose = 0){
  // split regions on same chr to get total CN (for shatterseek input)
  for(auto cnp : cn_by_chr_hap){
    int chr = cnp.first.first;
    int haplotype = cnp.first.second;
    vector<interval> intls = cnp.second;
    // find all breakpoints
    set<int> bps;
    for(auto intl : cnp.second){
      int start = intl.start;
      int end = intl.end;
      bps_by_chr[chr].insert(start);
      bps_by_chr[chr].insert(end);
    }
  }

  // add end points of each chromosome to avoid missing segments with copy number 0
  for(int chr = 0; chr < NUM_CHR; chr++){
    bps_by_chr[chr].insert(1);
    bps_by_chr[chr].insert(CHR_LENGTHS[chr]);
  }

  if(verbose > 1){
    cout << "\nbreakpoints for each chr after merging regions with the same CN" << endl;
    for(auto bpc : bps_by_chr){
      int chr = bpc.first;
      cout << chr + 1 << ":";
      for(auto bp : bpc.second){
        cout << " " << bp;
      }
      cout << endl;
    }
  }
}


// Compute allele-specific and total copy number for non-overlapping segments
void calculate_segment_cn(map<int, set<int>>& bps_by_chr, int verbose = 0){
  // verbose = 1;
  // assert(this->chr_segments.size() == 0); // only compute once for a cell
  if(this->chr_segments.size() > 0){
    for(auto sg : chr_segments){
      for(auto s : sg.second){
        delete s;
      }
    }
    chr_segments.clear();
  }
  int sid = 0;

  // split the regions according to the breakpoints
  // one record for each interval across haplotypes
  for(auto bpc : bps_by_chr){
    int chr = bpc.first;
    set<int> bps = bpc.second;  // should be sorted
    vector<pair<int, int>> intls;

    set<int>::iterator it = bps.begin();
    int pos1 = *it;
    it++;
    for(; it != bps.end(); it++){
      int pos2 = *it;
      if(abs(pos1 - pos2) == 1){  // consecutive breakpoints
        pos1 = pos2;
        continue;
      }
      pair<int, int> intl(pos1, pos2);
      intls.push_back(intl);
      pos1 = pos2;
    }

    if(verbose > 1){
      cout << "\nintervals for chr" << chr + 1 << ":";
      for(auto intl : intls){
        cout << " (" << intl.first << "," << intl.second << ")";
      }
      cout << endl;

      // find cns for each interval in each haplotype
      cout << "CNs for haplotype A" << endl;
    }

    map<pair<int, int>, int> cn_A;
    pair<int, int> key(chr, 0);
    vector<interval> cnsA = cn_by_chr_hap_merged[key];
    for(auto intl : intls){
      for(auto cnp : cnsA){
        int start = cnp.start;
        int end = cnp.end;
        int cn = cnp.cn;
        if(intl.first >= start && intl.second <= end){
          cn_A[intl] += cn;
        }
      }
    }

    if(verbose > 1) cout << "CNs for haplotype B" << endl;
    map<pair<int, int>, int> cn_B;
    pair<int, int> key1(chr, 1);
    vector<interval> cnsB = cn_by_chr_hap_merged[key1];
    for(auto intl : intls){
      for(auto cnp : cnsB){
        int start = cnp.start;
        int end = cnp.end;
        int cn = cnp.cn;
        if(intl.first >= start && intl.second <= end){
          cn_B[intl] += cn;
        }
      }
    }

    if(verbose > 1) cout << "CNs for segment" << endl;
    for(auto intl : intls){
      segment* s = new segment(sid++, cell_ID, chr, intl.first, intl.second, cn_A[intl], cn_B[intl]);
      this->chr_segments[chr].push_back(s);
      if(verbose > 1) s->print();
    }
  }
}


string get_sv_string(adjacency* adj, string extra){
    assert(adj->sv_type != NONE);

    breakpoint* j1 = breakpoints[adj->junc_id1];
    breakpoint* j2 = breakpoints[adj->junc_id2];

    breakpoint* jm = NULL;
    if(j1->chr == j2->chr && j1->pos > j2->pos){
      jm = j1;
      j1 = j2;
      j2 = jm;
    }

    string line = to_string(j1->chr + 1) + "\t" + to_string(j1->pos) + "\t" + get_side_string(j1->side) + "\t" + to_string(j2->chr + 1) + "\t" + to_string(j2->pos) + "\t" + get_side_string(j2->side) + "\t" + extra + "\n";

    return line;
}


// merge consecutive positions with the same copy number at both haplotypes
void get_merged_svs(const vector<pos_cn>& orig_pos, vector<pos_cn>& merged_pos, int verbose = 0){
    pos_cn intl_prev = orig_pos[0];
    merged_pos.push_back(intl_prev);
    if(verbose > 1) cout << "interval " << 1 << "\t" << intl_prev.chr << "\t" << intl_prev.start << "\t" << intl_prev.end << "\t" << intl_prev.cnA << "\t" << intl_prev.cnB << endl;

    bool merged = false;
    // cout << intls.size() << endl;
    for(int i = 1; i < orig_pos.size(); i++){
      pos_cn intl_curr = orig_pos[i];
      if(verbose > 1) cout << "interval " << i + 1 << "\t" << intl_curr.chr << "\t" << intl_curr.start << "\t" << intl_curr.end << "\t" << intl_curr.cnA << "\t" << intl_curr.cnB << endl;

      // if(intl_prev.chr == intl_curr.chr && intl_prev.end == intl_curr.start && intl_prev.cnA == intl_curr.cnA  && intl_prev.cnB == intl_curr.cnB)
      // only check total CN
      if(intl_prev.chr == intl_curr.chr && intl_prev.end == intl_curr.start && intl_prev.cnA + intl_prev.cnB == intl_curr.cnA + intl_curr.cnB){
        pos_cn intl_merged{intl_prev.chr, intl_prev.start, intl_curr.end, intl_prev.cnA, intl_prev.cnB};
        if(verbose > 1) cout << "interval after merging with previous one\t" << intl_merged.chr <<intl_merged.start << "\t" << intl_merged.end << "\t" << intl_merged.cnA << "\t" << intl_merged.cnB << endl;
        merged_pos.pop_back();
        merged_pos.push_back(intl_merged);
        intl_prev = intl_merged;
      }else{
        if(verbose > 1) cout << "interval without merging with previous one\t" << intl_curr.chr << "\t" << intl_curr.start << "\t" << intl_curr.end << "\t" << intl_curr.cnA << "\t" << "\t" << intl_curr.cnB << endl;
        merged_pos.push_back(intl_curr);
        intl_prev = intl_curr;
      }
    }
}

// Get duplicated or deleted region to set DUP/DEL-like patterns based on copy number, which may be caused during mitosis as a result of imbalanced distribution
void get_pseudo_adjacency(vector<pos_cn>& merged_dups_pos, vector<pos_cn>& merged_dels_pos, int verbose = 0){
    assert(chr_segments.size() > 0);
    assert(merged_dups_pos.size() == 0);
    assert(merged_dels_pos.size() == 0);

    if(verbose > 1) cout << "Set additional DUP/DEL-like patterns based on copy number" << endl;
    vector<pos_cn> dups_pos;
    vector<pos_cn> dels_pos;
    // original sets of positions may be more fragmented due to global breakpoints across cells
    for(auto segs : chr_segments){
      for(auto s: segs.second){
        if(verbose > 1) s->print();
        pos_cn pos{s->chr, s->start, s->end, s->cnA, s->cnB};
        int tcn = s->cnA + s->cnB;
        if(tcn > 2){
          dups_pos.push_back(pos);
        }else if(tcn < 2){
          dels_pos.push_back(pos);
        }else{

        }
      }
    }
    // cout << dups_pos.size() << endl;
    if(dups_pos.size() > 0){
      get_merged_svs(dups_pos, merged_dups_pos, verbose);
    }
    if(dels_pos.size() > 0){
      get_merged_svs(dels_pos, merged_dels_pos, verbose);
    }
  }


void update_adj_cn(int hap1, int hap2, adj_cn& ac){
    if(hap1 == 0 && hap2 == 0){
      ac.cnAA++;
    }else if(hap1 == 0 && hap2 == 1){
      ac.cnAB++;
    }else if(hap1 == 1 && hap2 == 0){
      ac.cnBA++;
    }else{
      assert(hap1 == 1 && hap2 == 1);
      ac.cnBB++;
    }
}


// Get adjacency haplotype-specific CNs for AA, AB, BA, BB
void get_adjacency_CN(){
  adjacency_CNs.clear();

  string at = "";
  for(auto p: paths){
    for(auto e : p.second->edges){
      // cout << "\nedge " << e << endl;
      adjacency *adj = adjacencies[e];
      // adj->print();

      if(adj->type == 0){ // 0: interval, 1: reference, 2: variant
        continue;
      }

      // https://github.com/aganezov/RCK/blob/master/docs/Adjacencies.md
      if(adj->type == 1){
        at = "R";  // reference
      }else{
        at = "N"; // novel
      }

      breakpoint *j1 = breakpoints[adj->junc_id1];
      breakpoint *j2 = breakpoints[adj->junc_id2];
      // make sure adjacencies connecting the same two positions are grouped together
      // hard to compare positions by both chr and pos, assume is_inverted is correct
      if(adj->is_inverted){
        j1 = breakpoints[adj->junc_id2];
        j2 = breakpoints[adj->junc_id1];
      }

      // 0: A, 1: B
      int hap1 = j1->haplotype;
      int hap2 = j2->haplotype;

      adj_pos ap{j1->chr, j1->pos, j1->side, j2->chr, j2->pos, j2->side, at};

      if(adjacency_CNs.count(ap) > 0){
        adj_cn ac = adjacency_CNs[ap];
        // update CN
        update_adj_cn(hap1, hap2, ac);
        adjacency_CNs[ap] = ac;
      }else{
        int cnAA = 0;
        int cnAB = 0;
        int cnBA = 0;
        int cnBB = 0;
        adj_cn ac{cnAA, cnAB, cnBA, cnBB};
        update_adj_cn(hap1, hap2, ac);
        adjacency_CNs[ap] = ac;
      }
    }
  }
}


void get_cn_bin(const vector<pos_bin>& bins, const vector<int>& bin_number, int verbose = 0){
  vector<int> bin_count(bins.size(), 0);   // used to count number of segments falling into the same bin
  for(int i = 0; i < bins.size(); i++){
    bin_tcn.push_back(0.0);
    bin_cnA.push_back(0.0);
    bin_cnB.push_back(0.0);
  }

  int bin_size = bins[0].end - bins[0].start + 1; // bin size used in copy number calling

  for(auto sg : chr_segments){
    for(auto s: sg.second){
      // compute bin number
      int start_bin = floor(s->start / bin_size) + bin_number[s->chr];
      int end_bin = floor(s->end / bin_size) + bin_number[s->chr];

      //check bins at the boundaries to compute fraction of coverage
      int cov_start = bins[start_bin].end - s->start + 1;
      bin_cnA[start_bin] += (double) s->cnA * cov_start / bin_size;
      bin_cnB[start_bin] += (double) s->cnB * cov_start / bin_size;
      bin_count[start_bin]++;

      if(end_bin > start_bin){
        // check bin at chromosome boundaries
        int cov_end = s->end - bins[end_bin].start + 1;
        int cov_size = bin_size;
        if(CHR_LENGTHS[s->chr] > bins[end_bin].start && CHR_LENGTHS[s->chr] < bins[end_bin].end){
          cov_size = CHR_LENGTHS[s->chr] - bins[end_bin].start + 1;
          if(verbose > 1) cout << "last bin at chr " << s->chr + 1 << " with length " << CHR_LENGTHS[s->chr] << endl;
        }
        double cnA = (double) s->cnA * cov_end / cov_size;
        double cnB = (double) s->cnB * cov_end / cov_size;
        bin_cnA[end_bin] += cnA;
        bin_cnB[end_bin] += cnB;
        bin_count[end_bin]++;

        for(int j = start_bin + 1; j < end_bin; j++){
          bin_cnA[j] = s->cnA;
          bin_cnB[j] = s->cnB;
          bin_count[j]++;
          if(bin_count[j] > 1){
            cout << "multiple segments in a bin: " << cov_start << " " << cov_end << " " << j << endl;
            exit(1);
          }
        }

          if(verbose > 1){
          cout << s->chr + 1 << "\t" << s->start << "\t" << s->end << "\t" << s->cnA << "\t" << s->cnB << "\t" <<  bins[start_bin].start << "\t" << bins[start_bin].end << "\t"  << bins[end_bin].start << "\t"  << bins[end_bin].end << "\t" << cov_start << "\t" << cov_end << "\t" << cnA << "\t" << cnB << endl;
          for(int j = start_bin; j <= end_bin; j++){
            cout << bin_cnA[j] << "\t" << bin_cnB[j] << "\t" << bin_count[j] << endl;
          }
        }
      }
    }
  }
  // all segments should have unique positions, so bin_count should be equal to 1
  // compute average CN for each bin
  for(int i = 0; i < bins.size(); i++){
    if(bin_count[i] == 0 && bin_tcn[i] ==0 && bin_cnA[i] == 0 && bin_cnB[i] == 0){
      bin_cnA[i] = 1;
      bin_cnB[i] = 1;
    }
    bin_tcn[i] = bin_cnA[i] + bin_cnB[i];
  }
}


/********************* functions related to sanity checking ******************/

// check whether two intervals defined by j1 and j2 are overlapping
bool is_interval_overlap(breakpoint* j1, breakpoint* j2){
  int i1_start = j1->pos;
  int i1_end = breakpoints[j1->right_jid]->pos;
  int i2_start = breakpoints[j2->left_jid]->pos;
  int i2_end = j2->pos;

  int im = 0;
  if(i1_start > i1_end){
    im = i1_start;
    i1_start = i1_end;
    i1_end = im;
  }

  if(i2_start > i2_end){
    im = i2_start;
    i2_start = i2_end;
    i2_end = im;
  }

  if((i1_start < i2_end && i1_start > i2_start) || (i1_end < i2_end && i1_end > i2_start)){
    return true;
  }else{
    return false;
  }
}

void check_duplicated_breakpoint(breakpoint* bp){
  // check path is not in current paths
  vector<int> ids;
  for(auto bp : breakpoints){
    ids.push_back(bp.first);
  }
  if(find(ids.begin(), ids.end(), bp->id) != ids.end()){
    cout << "breakpoint already included!" << endl;
    bp->print();
    exit(FAIL);
  }
}


void check_duplicated_adjacency(adjacency* adj){
  // check path is not in current paths
  vector<int> ids;
  for(auto adj : adjacencies){
    ids.push_back(adj.first);
  }
  if(find(ids.begin(), ids.end(), adj->id) != ids.end()){
    cout << "adjacency already included!" << endl;
    adj->print();
    exit(FAIL);
  }
}


void check_duplicated_path(path* p){
  // check path is not in current paths
  vector<int> pids;
  for(auto p : paths){
    pids.push_back(p.first);
  }
  if(find(pids.begin(), pids.end(), p->id) != pids.end()){
    cout << "path already included!" << endl;
    p->print();
    exit(FAIL);
  }
}


// check the consistency of nodes and edges in a path of the genome 
// after calling get_derivative_genome, get_path_from_bp
bool validate_path(path* p){
  check_duplicated_path(p);

  // remove duplicated start node  
  if(p->nodes[0] == p->nodes[p->nodes.size() - 1]){
    // if(!p->is_circle){  // function may be called when circle tag is not set
    //   cout << "circlurar path should have the right tag!" << endl;
    //   p->print();
    //   exit(FAIL);
    // }  
    if(p->nodes.size() % 2 == 0){
      cout << "Incorrect number of nodes in circular path " << p->id + 1 << endl;
      print_full_path(p);
    }
    // p->nodes.pop_back();   // assume duplicated node has been removed
  }

  if(p->nodes.size() % 2 != 0){   // must be satisfied due to alternative path types
    cout << "Incorrect number of nodes in path " << p->id + 1 << endl;
    print_full_path(p);
    // exit(FAIL);
    return false;
  }

  if(!p->is_circle){
    if(p->nodes.size() != p->edges.size() + 1){
      // set_circle_by_edge(p, 0);  // correct before reporting error
      cout << "Incorrect number of edges in linear path " << p->id + 1 << endl;
      print_full_path(p);
      // exit(FAIL);
      return false;
    }
  }else{
      if(p->nodes.size() != p->edges.size()){
      cout << "Incorrect number of edges in circular path " << p->id + 1 << endl;
      print_full_path(p);
      // exit(FAIL);
      return false;
    }
  }
  return true;
}


// make sure all breakpoints and adjacencies are in some paths
void check_derived_genome(){
    for(auto bp : breakpoints){
      breakpoint* b = bp.second;
      if(b->path_ID < 0){
        cout << "\nwrong path connection with missing breakpoint!" << endl;
        b->print();
        exit(FAIL);
      }
    }

    for(auto am : adjacencies){
      adjacency* a = am.second;
      if(a->path_ID < 0){
        cout << "\nwrong path connection with missing adjacency!" << endl;
        a->print();
        exit(FAIL);
      }
    }
}


};
