// genome.hpp

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>


using namespace std;

#include "util.hpp"


struct interval{
  int start;
  int end;
  int cn;
  int jid_start;  // keep breakpoint ID for tracking
  int jid_end;

  bool operator<(interval const &other) const {
    return (start < other.start ||
      (start == other.start && end < other.end));
  }
};


struct haplotype_pos{
  int chr;
  int haplotype;
  int start;
  int end;
  int jid_start;  // keep breakpoint ID for tracking
  int jid_end;

  bool operator<(haplotype_pos const &other) const {
    return (chr < other.chr ||
    (chr == other.chr && haplotype < other.haplotype) ||
    (chr == other.chr && haplotype == other.haplotype && start < other.start) ||
    (chr == other.chr && haplotype == other.haplotype && start == other.start && end < other.end));
  }
};


// an interval on a chromosome, used in summarizing CNs across genome
class segment{
public:
  int id;
  int cell_ID;
  int chr;
  int start;
  int end;
  // int haplotype;
  int cnA;
  int cnB;
  // int type;  // telomere or not. 0: normal, 1: telomere, 2: centromere
  // bool is_inverted;

  segment(int id, int cell_ID, int chr, int start, int end, int cnA, int cnB){
    this->id = id;
    this->cell_ID = cell_ID;
    this->chr = chr;
    this->start = start;
    this->end = end;
    this->cnA = cnA;
    this->cnB = cnB;
  }

  segment(const segment& _s2){
    chr    = _s2.chr;
    id = _s2.id;
  }


  void print(){
    cout << "segment " << id + 1 << " in cell " << cell_ID << " at chr " << chr + 1 << ": " << start << " - " << end << " with copy number " << cnA << " and " << cnB << endl;
  }

  // this method checks whether a node is telomeric
  // the terminal segments
  bool is_telomere(segment seg){
      int chr = seg.chr;
      // int chr_size = CHR_LENGTHS[chr];
      int telo_end1 = TELO_ENDS1[chr];
      int telo_end2 = TELO_ENDS2[chr];

      if(seg.start < telo_end1 || seg.end > telo_end2){
        return true;
      }
      return false;
  }

  // iterates along list of breakpoints, assigns centromeric = True property
  // based on breakpoint position relating to true genomic centromere position
  bool is_centromere(segment seg){
    int chr = seg.chr;
    // int chr_size = CHR_LENGTHS[chr];
    int cent_start = CENT_STARTS[chr];
    int cent_end = CENT_ENDS[chr];

    if((seg.start < cent_start && seg.end > cent_start)
      || (seg.start > cent_start && seg.end < cent_end)
      || (seg.start < cent_end && seg.end > cent_end)){
      return true;
    }
    return false;
  }

  // checks whether a segment (accessed via the breakpoint) is inverted; as inferred via the list of paths
  // if segment is inverted, node property "inv" = True
  bool is_inverted(){
    return false;
  }

};


// a breakpoint in the chromosome, introduced by a DSB, haplotype-specific,
// node in the genome graph
// one instance for each copy
// TODO: rename to be more accurate breakpoint
class breakpoint{
public:
  int id;   // keep increasing when new breakpoints are added
  int cell_ID;
  int path_ID = -1;  // a breakpoint can only belong to one path
  int chr;
  int pos;
  int side;   // 0 left or 1 right, strand in RCK, orientation
  int haplotype;
  // int cn;  // there may be multiple copies due to replication??
  bool is_repaired; // end breakpoint has only one connection (not necessarily telomere)
  // bool is_visited;  // whether in a path or not, used to determine circular paths
  int left_jid = -1;
  int right_jid = -1;
  bool is_end;  // whether it is at the end of genome or not, correspond to telomere when DSBs are not introduced to telomere regions
  // bool is_centromeric; // whether it overlaps with centromere, correspond to centromere when DSBs are not introduced to centromere regions

  breakpoint(){

  }


  breakpoint(int cell_ID, int id, int chr, int pos, int side, int haplotype, bool is_end, bool is_repaired){
    this->cell_ID = cell_ID;
    this->id = id;
    this->chr = chr;
    this->pos = pos;
    this->side = side;
    this->haplotype = haplotype;
    this->is_end = is_end;
    this->is_repaired = is_repaired;
    // this->is_visited = is_visited;
  }


  breakpoint(const breakpoint& bp){
    cell_ID = bp.cell_ID;
    id = bp.id;
    path_ID = bp.path_ID;
    chr = bp.chr;
    pos = bp.pos;
    side = bp.side;
    haplotype = bp.haplotype;
    is_end = bp.is_end;
    is_repaired = bp.is_repaired;
    left_jid = bp.left_jid;
    right_jid = bp.right_jid;
    is_end = bp.is_end;
  }


  bool operator==(const breakpoint &j) const {
      return id == j.id;
  }


  void print() {
    string conn = "not repaired";
    if(is_repaired){
      conn = "repaired";
    }
    string end = "not at genome end";
    if(is_end){
      end = "at genome end";
    }
    cout << "breakpoint " << id << " in cell " << cell_ID << " at path " << path_ID + 1 << " at chr " << chr + 1 << " position " << pos << ", side " << get_side_string(side) << ", haplotype " << get_haplotype_string(haplotype) << ", " << conn << ", " << end << ", left breakpoint " << left_jid << ", right breakpoint " << right_jid << endl;
  }


  void find_adjacent_breakpoint(){
    // left sided breakpoints

    // right sided breakpoints

    // closest breakpoint

    // furtherest breakpoint

    // telomere
  }

};


bool compare_bp_pos(breakpoint* j1, breakpoint* j2){
  assert(j1->haplotype == j2->haplotype);
  assert(j1->chr == j2->chr);
  return j1->pos < j2->pos;
}


bool compare_bp_pos_by_chr(breakpoint* j1, breakpoint* j2){
  assert(j1->chr == j2->chr);
  return (j1->pos < j2->pos) || (j1->pos == j2->pos && j1->haplotype < j2->haplotype) ;
}


bool compare_bp_pos_by_haplotype(breakpoint* j1, breakpoint* j2){
  assert(j1->chr == j2->chr && j1->pos == j2->pos);
  return (j1->haplotype < j2->haplotype);
}


// used to output SVs
class adjacency{
public:
  int id;
  int cell_ID;
  int path_ID;
  int junc_id1;
  int junc_id2;
  // int cn;
  int type;   // 0: interval, 1: reference, 2: variant
  bool is_centromeric; // for interval only
  int telomeric_type;  // for interval only; 0: no telomere; 1: left telomere; 2: right telomere; 3: both telomeres
  bool is_inverted;    // indicative for interval only, as it is hard to define direction of variant adjacency
  int sv_type;   // assigned later based on the two breakpoints

  adjacency(int cell_ID, int id, int path_ID, int junc_id1, int junc_id2, int type, int telomeric_type, int sv_type = NONE){
    this->id = id;
    this->cell_ID = cell_ID;
    this->path_ID = path_ID;
    this->junc_id1 = junc_id1;
    this->junc_id2 = junc_id2;
    this->type = type;
    this->telomeric_type = telomeric_type;
    this->sv_type = sv_type;
    this->is_inverted = false;
    this->is_centromeric = false;
  }

  adjacency(const adjacency& adj){
    id = adj.id;
    cell_ID = adj.cell_ID;
    path_ID = adj.path_ID;
    junc_id1 = adj.junc_id1;
    junc_id2 = adj.junc_id2;
    type = adj.type;
    is_centromeric = adj.is_centromeric;
    telomeric_type = adj.telomeric_type;
    is_inverted = adj.is_inverted;
    sv_type = adj.sv_type;
  }

  void print() {
    string cent = "non-centromere";
    if(is_centromeric){
      cent = "centromere";
    }
    string invt = "forward";
    if(is_inverted){
      invt = "inverted";
    }
    // path 0 means the adjacency is not in any path
    cout << "Adjacency " << id << " in cell " << cell_ID << " at path " << path_ID + 1 << " with left breakpoint " << junc_id1 << " and right breakpoint " << junc_id2 << ", " << cent << ", " << invt << ", " << get_adj_type_string(type) << ", " << get_telomere_type_string(telomeric_type) << ", " << get_sv_type_string(sv_type) << endl;
  }


  // print interval with position
  void print_interval(map<int, breakpoint*>& breakpoints){
    string cent = "non-centromere";
    if(is_centromeric){
      cent = "centromere";
    }
    string invt = "forward";
    if(is_inverted){
      invt = "inverted";
    }
    cout << "Adjacency " << id << " in cell " << cell_ID << " at path " << path_ID + 1 << " with left breakpoint " << junc_id1 << " at " << breakpoints[junc_id1]->chr + 1 << "-" << breakpoints[junc_id1]->pos << " and right breakpoint " << junc_id2 << " at " << breakpoints[junc_id2]->chr + 1 << "-" << breakpoints[junc_id2]->pos << ", " << cent << ", " << invt << ", " << get_telomere_type_string(telomeric_type) << ", " << get_sv_type_string(sv_type) << endl;
  }

  // this method checks whether an interval adjacency is telomeric
  // the terminal segments
  void set_telomeric_type(map<int, breakpoint*>& breakpoints){
    if(type != INTERVAL){
      telomeric_type = NONTEL;
      return;
    }

    int chr = breakpoints[junc_id1]->chr;
    int chr2 = breakpoints[junc_id2]->chr;
    assert(chr == chr2);

    // int chr_size = CHR_LENGTHS[chr];
    int telo_end1 = TELO_ENDS1[chr];
    int telo_end2 = TELO_ENDS2[chr];
    int start = breakpoints[junc_id1]->pos;
    int end = breakpoints[junc_id2]->pos;

    if(start <= telo_end1 && end >= telo_end2){
      telomeric_type = COMPLETE;
    }else if(start > telo_end1 && end >= telo_end2){
      telomeric_type = QTEL;
    }else if(start <= telo_end1 && end < telo_end2){
      telomeric_type = PTEL;
    }else{  // start > telo_end1 && end < telo_end2
      telomeric_type = NONTEL;
    }
  }

  // iterates along list of breakpoints, assigns centromeric = True property
  // based on breakpoint position relating to true genomic centromere position
  void set_centromeric_status(map<int, breakpoint*>& breakpoints){
    if(type != 0){
      is_centromeric = false;
      return;
    }

    int chr = breakpoints[junc_id1]->chr;
    int chr2 = breakpoints[junc_id2]->chr;
    assert(chr == chr2);

    // int chr_size = CHR_LENGTHS[chr];
    int cent_start = CENT_STARTS[chr];
    int cent_end = CENT_ENDS[chr];
    int start = breakpoints[junc_id1]->pos;
    int end = breakpoints[junc_id2]->pos;

    if((start < cent_start && end > cent_start)
      || (start > cent_start && end < cent_end)
      || (start < cent_end && end > cent_end)){
      is_centromeric = true;
    }
    is_centromeric = false;
  }
};


// a path is a consecutive genomic segment (linear chromosome)
// a connected component in the interval-adjacency graph
class path{
public:
  int id;   // + 1 when printing to facilitate counting
  int cell_ID;
  // int start;
  // int end;
  vector<int> nodes;  // breakpoint IDs
  vector<int> edges;  // adjacency IDs following the derivative chromosome
  // int size;  // number of nodes in the path
  int type;  // 0: nonTel; 1: pTel; 2: qTel; 3: ptel to tel (complete)
  int n_centromere;   // counts how many centromeres there are in a defined path
  bool is_circle;
  int sibling;   // Its copy in S phase, used in keep tracking of balanced distribution
  int child_cell_ID;  // The cell ID of its copy in next generation, may be more than one for circular paths, used in keep tracking of balanced distribution

  path(int id, int cell_ID, int type){
    this->id = id;
    this->cell_ID = cell_ID;
    this->type = type;
    this->n_centromere = 0;
    this->is_circle = false;
    this->sibling = -1;
    this->child_cell_ID = -1;
  }

  void print(){
    string shape = "linear";
    if(is_circle){
      shape = "circular";
    }
    string sibling_str = "";
    if(sibling != -1){
        sibling_str = ", sibling path ID " + to_string(sibling + 1);
    }
    cout << "Path " << id + 1 << " in cell " << cell_ID << ", " << shape << sibling_str << ", with " << n_centromere << " centromere and " << get_telomere_type_string(type) << "; " << nodes.size() << " nodes: ";
    for(auto n : nodes){
      cout << n << ",";
    }
    cout << "; " << edges.size() << " edges: ";
    for(auto e : edges){
      cout << e << ",";
    }
    cout << endl;
  }


  // split the path at (left_jid, right_jid), not used for now
  void split_path(int j, path* p1, path* p2, int left_jid, int right_jid, const vector<int>& old_aids, bool is_inverted, const vector<breakpoint*>& junc_new, const vector<adjacency*>& adj_new, int verbose = 0){

    int j1 = junc_new[0 + j * 2]->id;
    int j2 = junc_new[1 + j * 2]->id;
    int adj1 = adj_new[0 + j * 2]->id;
    int adj2 = adj_new[1 + j * 2]->id;
    int j1l = left_jid;
    int j2r = right_jid;
    int aid = old_aids[0 + j * 2];

    if(is_inverted){
      if(verbose > 1) cout << "inverted" << endl;
      j1 = junc_new[1 + j * 2]->id;
      j2 = junc_new[0 + j * 2]->id;
      adj1 = adj_new[1 + j * 2]->id;
      adj2 = adj_new[0 + j * 2]->id;
      j1l = right_jid;
      j2r = left_jid;
    }

    if(verbose > 1){
      cout << "splitting path into two via old adjacency " << aid << "\t" << j1 << "\t" << j2 << "\t" << adj1 << "\t" << adj2 << "\t" << j1l << "\t" << j2r << "\t" << endl;
    }

    // copy nodes until left_jid to p1, the remaining part from j2 belong to p2
    if(verbose > 1) cout << "copy nodes" << endl;
    bool p1_end = false;
    p2->nodes.push_back(j2);
    for(auto node : nodes){
      if(node == j1l){
        p1_end = true;
        continue;
      }
      if(!p1_end){
        p1->nodes.push_back(node);
      }else{
        p2->nodes.push_back(node);
      }
    }
    p1->nodes.push_back(j1l);
    p1->nodes.push_back(j1);

    if(verbose > 1) cout << "copy edges" << endl;
    p1_end = false;
    p2->edges.push_back(adj2);
    for(auto edge : edges){
      if(edge == aid){
        p1_end = true;
        continue;
      }
      if(!p1_end){
        p1->edges.push_back(edge);
      }else{
        p2->edges.push_back(edge);
      }
    }
    p1->edges.push_back(adj1);

    p1->n_centromere = 1;
    p2->n_centromere = 1;
  }


  // copy a subpath of p from breakpoint j1 to j2 (not used)
  void copy_subpath(path p, int j1, int j2, const vector<adjacency>& adjacencies){
    nodes.clear();
    edges.clear();

    int curr_eid = -1;
    int start_eid = -1;
    // find the id of start breakpoint
    for(int i = 0; i < p.edges.size(); i++){
      curr_eid = p.edges[i];
      if(adjacencies[curr_eid].junc_id1 == j1){
        start_eid = i;
        break;
      }
    }

    for(int i = start_eid; i < p.edges.size(); i++){
      curr_eid = p.edges[i];
      if(adjacencies[curr_eid].junc_id1 == j2){
        break;
      }else{
        edges.push_back(curr_eid);
        nodes.push_back(adjacencies[curr_eid].junc_id1);
      }
    }

    assert(curr_eid > 0);
    // last node
    nodes.push_back(adjacencies[curr_eid].junc_id2);

    // type based on property of j1 and j2
  }


  // counts how many centromeres there are in a defined path, not used for now
  int get_ncent(){
    int ncent = 0;
    return ncent;
  }

};




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
  vector<segment*> segments;
  // int n_adj;   // total number of adjacencies introduced in history, used to get unique adjacency ID;

  // number of DSBs on each chr each haplotype
  // vector<int> vec_n_dsb;   // hard to maintain across cell divisions due to random path distributions
  // breaks introduced due to multiple centromeres during mitosis
  // vector<int> vec_n_mitosis_break;  // seem not biologically meaningful
  // map<int, vector<int>> chr_type_num;  // chr, num for each SV type indexed by TYPE

  vector<int> chr_ends;

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
    for(auto s : segments){
      delete s;
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
    for(int chr = 0; chr < NUM_CHR; chr++){
      // vec_n_dsb.push_back(0);
      // vector<int> n_sv(NUM_SVTYPE, 0);
      // chr_type_num[chr] = n_sv;

      bool is_end = true;
      bool is_repaired = true;

      // two breakpoints are generated
      int haplotype = 0;
      // int cell_ID, int id, int chr, int pos, int side, int haplotype, bool is_end, bool is_repaired
      breakpoint* j1 = new breakpoint(cell_ID, jid++, chr, 0, TAIL, haplotype, is_end, is_repaired);
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
      validate_path(p);
      paths[p->id] = p;

      haplotype = 1;
      breakpoint* j3 = new breakpoint(cell_ID, jid++, chr, 0, TAIL, haplotype, is_end, is_repaired);
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
      validate_path(p2);
      paths[p2->id] = p2;

      chr_ends.insert(chr_ends.end(), {j1->id, j2->id, j3->id, j4->id});
    }
    // n_adj = adjacencies.size();

    // get_breakpoint_map();
  };


  genome(const genome& _g2) {

  }


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


  int find_max_adjID(){
    vector<int> adj_IDs(adjacencies.size(), 0);
    int i = 0;

    for(auto adj : adjacencies){
      adj_IDs[i++] = adj.first;
    }

    int max_adjID = *max_element(adj_IDs.begin(), adj_IDs.end());

    return max_adjID;
  }


  int find_max_pathID(){
    vector<int> path_IDs(paths.size(), 0);
    int i = 0;

    for(auto p : paths){
      path_IDs[i++] = p.first;
    }

    int max_pID = *max_element(path_IDs.begin(), path_IDs.end());

    return max_pID;
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
      adjacency *adj = adjacencies[e];
      // adj->print();
      breakpoint *j1 = breakpoints[adj->junc_id1];
      breakpoint *j2 = breakpoints[adj->junc_id2];
      if(adj->is_inverted){
        j1 = breakpoints[adj->junc_id2];
        j2 = breakpoints[adj->junc_id1];
      }

      fout << (j1->chr % NUM_CHR) + 1 << ":" << get_haplotype_string(j1->haplotype) << ":" << j1->pos << "-" << (j2->chr % NUM_CHR) + 1 << ":" << get_haplotype_string(j2->haplotype) << ":" << j2->pos << ",";
    }
    fout << endl;
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


  // find right neighbor of j1 if exists
  // assume juncs stores sorted breakpoints on the same haplotype and chr
  breakpoint* get_left_breakpoint(breakpoint* j1, const vector<breakpoint*>& juncs){
    for(int i = 0; i < juncs.size(); i++){
      if(juncs[i]->id == j1->id && i - 1 >= 0){
        // cout << "left jid " << juncs[i - 1]->id << endl;
        return juncs[i - 1];
      }
    }
    return NULL;
  }


  // find right neighbor of j1 if exists
  // assume juncs stores sorted breakpoints on the same haplotype and chr
  breakpoint* get_right_breakpoint(breakpoint* j1, const vector<breakpoint*>& juncs){
    for(int i = 0; i < juncs.size(); i++){
      if(juncs[i]->id == j1->id && i + 1 < juncs.size()){
        // cout << "right jid " << juncs[i + 1]->id << endl;
        return juncs[i + 1];
      }
    }
    return NULL;
  }


  // return the index of the adjacency specified by two breakpoints
  // an adjacency may appear inverted in a path
  // if one interval forms a circle, there will be two adjacencies between the breakpoints (one interval, one variant)
  void get_adjacency_ID(vector<int>& aids, int j1, int j2, map<int, adjacency*>& adjacencies){
    // cout << "get adjacency id for breakpoint " << j1 << " and " << j2 << endl;
    for(auto adjm : adjacencies){
      adjacency* adj = adjm.second;
      // adj->print();
      if(adj->junc_id1 == j1 && adj->junc_id2 == j2){
        // cout << "forward" << endl;
        aids.push_back(adj->id);
      }else if(adj->junc_id2 == j1 && adj->junc_id1 == j2){  // inverted
        // cout << "inverted" << endl;
        adj->is_inverted = true;
        aids.push_back(adj->id);
      }
    }
  }


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


  void remove_path(int path_ID, int verbose = 0){
    // for(std::vector<path*>::iterator it = paths.begin() ; it != paths.end(); ++it){
    //   if((*it)->id == path_ID){
    //     delete *it;
    //     paths.erase((it));
    //   }
    // }
    path* p = paths[path_ID];
    if(verbose > 1){
      cout << "remove path " << path_ID << endl;
      p->print();
    }
    delete p;
    paths.erase(path_ID);
  }


  // jid and aid are passed by reference to automatically get the new IDs for next new breakpoints
  // , bool update_path = false
  void add_new_breakpoint(int& jid, int &aid, int chr, int bp, int haplotype, int j1l_id = -1, int j2r_id = -1, bool reset_pathID = false, int verbose = 0){
    bool is_end = false;
    // a new breakpoint will generate two new segments j1, j2
    // (j1, j2) still unconnected, need to be repaired
    bool is_repaired = false;

    // two breakpoints at the same position are generated, with different sides (orientation)
    breakpoint* j1 = new breakpoint(cell_ID, jid++, chr, bp, HEAD, haplotype, is_end, is_repaired);
    breakpoint* j2 = new breakpoint(cell_ID, jid++, chr, bp, TAIL, haplotype, is_end, is_repaired);
    // only connect DSBs later in the repair stage
    j1->right_jid = -1;
    j2->left_jid = -1;
    breakpoints[j1->id] = j1;
    breakpoints[j2->id] = j2;

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

    // add new smaller intervals (j1l, j1), (j2, j2r)
    assert(j1l_id != -1);
    assert(j2r_id != -1);
    breakpoint* j1l = breakpoints[j1l_id];
    breakpoint* j2r = breakpoints[j2r_id];
    assert(j1l != NULL);
    assert(j2r != NULL);

    if(reset_pathID){
      j1->path_ID = -1;
    }else{
      j1->path_ID = j1l->path_ID;
    }
    adjacency* adj1 = new adjacency(cell_ID, aid++, j1->path_ID, j1l->id, j1->id, INTERVAL, NONTEL, NONE);
    adj1->set_telomeric_type(breakpoints);
    adjacencies[adj1->id] = adj1;
    // n_adj += 1;
    j1l->right_jid = j1->id;
    j1->left_jid = j1l->id;

    if(reset_pathID){
      j2->path_ID = -1;
    }else{
      j2->path_ID = j2r->path_ID;
    }
    adjacency* adj2 = new adjacency(cell_ID, aid++, j2->path_ID, j2->id, j2r->id, INTERVAL, NONTEL, NONE);
    adj2->set_telomeric_type(breakpoints);
    adjacencies[adj2->id] = adj2;
    j2->right_jid = j2r->id;
    j2r->left_jid = j2->id;
    // n_adj += 1;

    if(verbose > 1){
      cout << "\nadded two adjacencies introduced by break " << endl;
      adj1->print();
      adj2->print();
      cout << "\nupdated four breakpoints affected by break " << endl;
      j1l->print();
      cout << "breakpoint at the left of j1: ";
      j1l->print();
      j1->print();
      j2->print();
      cout << "breakpoint at the right of j2: ";
      j2r->print();
    }

    // remove old larger intervals
    // cout << "adjacencies before removing j1l to j2r" << endl;
    // for(auto a : adjacencies){
    //   a->print();
    // }
    remove_adjacency(j1l->id, j2r->id, adj1, adj2, verbose);
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

  // iminate DSB by introducing breakpoints (breakpoints) across the whole genome, following infinite sites assumption
  // a breakpoint will not disappear once exist
  // n_dsb: number of breakpoints for each chromosome
  // TODO: add interval DSB probability based on overlapping with FSs
  void generate_dsb(int n_dsb, int verbose = 0){
    // group breakpoints by haplotype and chr for sorting
    // get_breakpoint_map();
    gsl_ran_discrete_t* dis_loc = gsl_ran_discrete_preproc(NUM_CHR, CHR_PROBS);
    int jid = breakpoints.rbegin()->first + 1;
    int aid = adjacencies.rbegin()->first + 1;
    if(verbose > 1){
      cout << jid << " breakpoints before introducing DSB" << endl;
    }

    for(int i = 0; i < n_dsb; i = i + 1){
      // a breakpoint may be at centromere or telomere
      // a chromosome may lost some regions, only feasible on available segments
      int bp = 0;
      int chr = -1;
      int haplotype = -1;
      int left_jid = -1;
      int right_jid = -1;

      // available intervals in the haplotype after a new dsb
      get_unique_interval(verbose);

      bool is_insertable;
      // keep centromere and telomere intact for simiplicity
      do{
        is_insertable = true;
        chr = gsl_ran_discrete(r, dis_loc);
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
        bp = (int)runiform(r, intl.start, intl.end);
        left_jid = intl.jid_start;
        right_jid = intl.jid_end;
        // cout << "select interval "  << sel_intl << ": " << intl.start << ", " << intl.end  << " at bp " << bp << " with left breakpoint " << left_jid << " and right breakpoint " << right_jid << endl;

      }while(!is_insertable || bp <= TELO_ENDS1[chr] || bp >= TELO_ENDS2[chr] || (bp >= CENT_STARTS[chr] && bp <= CENT_ENDS[chr]));

      // update number of DSBs per chrom per haplotype
      // int idx = chr + haplotype * NUM_CHR;
      // vec_n_dsb[chr] += 1;

      if(verbose > 0){
        cout << "   break " << i + 1 << " at position " << bp << " chr " << chr + 1 << " haplotype " << get_haplotype_string (haplotype) << endl;
        if(verbose > 1) cout << "   with left breakpoint " << left_jid << " and right breakpoint " << right_jid << endl;
      }

      add_new_breakpoint(jid, aid, chr, bp, haplotype, left_jid, right_jid, false, verbose);
    }
  }


  SV_type set_var_type(breakpoint* j1, breakpoint* j2){
        SV_type sv_type = NONE;
        if(j1->side == HEAD && j2->side == TAIL){
          assert(j1->left_jid >= 0);
          assert(j2->right_jid >= 0);
          j1->right_jid = j2->id;
          j2->left_jid = j1->id;
          // DEL (deletion-like; +/-)
          sv_type = DEL;
          // interval in the reference genome is lost
          if(j1->chr != j2->chr){
            sv_type = TRA;
            // chr_type_num[j1->chr][sv_type] += 1;
            // chr_type_num[j2->chr][sv_type] += 1;
          }else{
            // chr_type_num[j1->chr][sv_type] += 1;
          }
        }else if(j1->side == HEAD && j2->side == HEAD){ //h2hINV
          assert(j1->left_jid >= 0);
          assert(j2->left_jid >= 0);
          j1->right_jid = j2->id;
          j2->right_jid = j1->id;
          // h2hINV (head-to-head inversion; +/+)
          sv_type = H2HINV;
          if(j1->chr != j2->chr){
            sv_type = TRA;
            // chr_type_num[j1->chr][sv_type] += 1;
            // chr_type_num[j2->chr][sv_type] += 1;
          }else{
            // chr_type_num[j1->chr][sv_type] += 1;
          }
        }else if(j1->side == TAIL && j2->side == HEAD){
          assert(j1->right_jid >= 0);
          assert(j2->left_jid >= 0);
          j1->left_jid = j2->id;
          j2->right_jid = j1->id;
          // DUP (duplication-like; -/+)
          sv_type = DUP;
          if(j1->chr != j2->chr){
            sv_type = TRA;
            // chr_type_num[j1->chr][sv_type] += 1;
            // chr_type_num[j2->chr][sv_type] += 1;
          }else{
            // chr_type_num[j1->chr][sv_type] += 1;
          }
        }else{  //if(j1->side == TAIL && j2->side == TAIL)
          assert(j1->right_jid >= 0);
          assert(j2->right_jid >= 0);
          j1->left_jid = j2->id;
          j2->left_jid = j1->id;
          // t2tINV (tail-to-tail inversion; -/-)
          sv_type = T2TINV;
          if(j1->chr != j2->chr){
            sv_type = TRA;
            // chr_type_num[j1->chr][sv_type] += 1;
            // chr_type_num[j2->chr][sv_type] += 1;
          }else{
            // chr_type_num[j1->chr][sv_type] += 1;
          }
        }

        return sv_type;
  }

  // n_unrepaired: number of unrepaired DSBs, each DSB introduces two breakpoints
  void repair_dsb(int n_unrepaired, vector<breakpoint*>& junc2repair, int verbose = 0){
    // verbose = 1;
    assert(n_unrepaired >= 0);
    // store segments without telomeres
    int aid = adjacencies.rbegin()->first + 1;

    // connect segments randomly to get variant adjacencies (repair DSBs)
    // only breakpoints with missing connections need to be repaired
    // vector<breakpoint> junc2repair(breakpoints.begin(), breakpoints.end());
    for(auto jm : breakpoints){
      breakpoint* j = jm.second;
      if(!j->is_repaired){
        junc2repair.push_back(j);
      }
    }

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

      int j2id = myrng(junc2repair.size());
      breakpoint* j2 = junc2repair[j2id];
      junc2repair.erase(std::remove(junc2repair.begin(), junc2repair.end(), j2), junc2repair.end());
      // cout << "#breakpoints remaining after 2nd repair: " << junc2repair.size() << endl;
      j2->is_repaired = true;

      if(verbose > 1){
        cout << "repair breakpoint " << j1->id << "\t" << j2->id << endl;
        j1->print();
        j2->print();
      }

      if(verbose > 0){
        cout << "   connecting " << j1->chr + 1 << ":" << get_haplotype_string(j1->haplotype) << ":" << j1->pos << get_side_string(j1->side) << " with " << j2->chr + 1 << ":" << get_haplotype_string(j2->haplotype) << ":" << j2->pos << get_side_string(j2->side) << endl;
      }

      // each breakpoint node should has degree 2
      adjacency* adj = NULL;

      SV_type sv_type = NONE;
      // determine direction by breakpoint side
      if(j1->chr == j2->chr && j1->pos == j2->pos && j1->haplotype == j2->haplotype && j1->side != j2->side){
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
        sv_type = set_var_type(j1, j2);
        adj = new adjacency(cell_ID, aid++, -1, j1->id, j2->id, VAR, NONTEL, sv_type);
        adj->set_telomeric_type(breakpoints);
        adjacencies[adj->id] = adj;
        // adj->print();
      } // else

      // this->n_adj++;
    } // while

    if(verbose > 1){
      cout << breakpoints.size() << " breakpoints after repairing dsb" << endl;
      for(auto j : breakpoints){
        (j.second)->print();
      }
      cout << adjacencies.size() << " adjacencies after repairing dsb" << endl;
      for(auto a : adjacencies){
        (a.second)->print();
      }
    }
  }


  void set_path_type(path* p, bool set_circle = true, int verbose = 0){
    if(verbose > 1){
      cout << "set type for path " << p->id + 1 << endl;
      p->print();
      write_path(p, cout);
    }

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
      // set to circular if there is no centromere
      if(p->n_centromere == 0 && set_circle && start != end){
        if((breakpoints[start]->left_jid == -1 || breakpoints[start]->right_jid == -1)){
          if(breakpoints[start]->left_jid == -1){
            breakpoints[start]->left_jid = end;
          }else{
            breakpoints[start]->right_jid = end;
          }
        }

        if((breakpoints[end]->left_jid == -1 || breakpoints[end]->right_jid == -1)){
          if(breakpoints[end]->left_jid == -1){
            breakpoints[end]->left_jid = start;
          }else{
            breakpoints[end]->right_jid = start;
          }
        }

        // add another adjacency from end to start
        int aid = adjacencies.rbegin()->first + 1;
        SV_type sv_type = set_var_type(breakpoints[end], breakpoints[start]);
        adjacency* adj = new adjacency(cell_ID, aid, p->id, end, start, VAR, NONTEL, sv_type);
        adj->is_inverted = true;

        p->edges.push_back(aid);
        adjacencies[aid] = adj;

        // p->nodes.push_back(start);
              
        breakpoints[end]->is_repaired = true;
        breakpoints[start]->is_repaired = true;
        
        p->is_circle = true;       
      }
    }
    if(verbose > 1) p->print();
  }


  // add one edge into path p, return success status
  bool update_path_by_adj(path* p, breakpoint* js, breakpoint* jn, map<int, adjacency*>& adjacencies, int verbose = 0){
    if(verbose > 1){
      cout << "adding edge from " << js->id << " to " << jn->id << " to path " << p->id + 1 << endl;
      p->print();
    }

    // find all the edges connected two breakpoints
    vector<int> aids;
    get_adjacency_ID(aids, js->id, jn->id, adjacencies);
    if(aids.size() == 0) return false;

    int aid = aids[0];
    adjacency* adj = adjacencies[aid];
    int aid2 = -1;
    adjacency* adj2;
    if(aids.size() > 1){ // a circle
      aid2 = aids[1];
      adj2 = adjacencies[aid2];
    }

    if(aids.size() == 1 && adj->path_ID == p->id){
      if(verbose > 1) cout << "adjacency " << aid << " has been visited!" << endl;
      return false;
    }else{ //aids.size() == 2
      if(adj->path_ID == p->id && adj2->path_ID == p->id){
        if(verbose > 1) cout << "adjacency " << aid << " and " << aid2 << " have been visited!" << endl;
        return false;
      }else{
        if(adj->path_ID == p->id && adj2->path_ID != p->id){
          aid = aid2;
          adj = adj2;
        }
      }
    }

    p->nodes.push_back(jn->id);
    p->edges.push_back(aid);
    // update adjacency direction based on current path
    if(js->id == adj->junc_id2){
      adj->is_inverted = true;
    }else{
      adj->is_inverted = false;
    }
    jn->path_ID = p->id;
    adj->path_ID = p->id;

    if((adj->type == INTERVAL) && (js->chr == jn->chr) &&
      ((js->pos <= CENT_STARTS[js->chr] && jn->pos >= CENT_ENDS[js->chr]) ||
      (jn->pos <= CENT_STARTS[js->chr] && js->pos >= CENT_ENDS[js->chr]))
      ){
    // if(adj->is_centromeric){
      // cout << "one more centromere" << endl;
      p->n_centromere += 1;
      adj->is_centromeric = true;
    }

    if(verbose > 1){
      cout << "\nupdate path " << p->id + 1 << " with adjacency " << aid << " with breakpoint " << js->id << ", " << jn->id << endl;
      adj->print();
      js->print();
      jn->print();
    }

    return true;
  }


  // a connected path of breakpoints starting from js
  // is_forward: true -- left to right (5' to 3', p-tel to q-tel), 1 -- right to left
  // check left breakpoint if the right breakpoint is the same as jp
  // check_interval: check whether the 1st adjacency is an interval
  // circular path will have the start node duplicated and removed during validation validate_path()
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
    while(nei != -1){
      if(verbose > 1) cout << "current breakpoint " << nei << endl;
      breakpoint* jn = breakpoints[nei];

      bool succ = update_path_by_adj(p, jp, jn, adjacencies, verbose);
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


  void validate_path(path* p){
    // check path is not in current paths
    vector<int> pids;
    for(auto p : paths){
      pids.push_back(p.first);
    }    
    if(find(pids.begin(), pids.end(), p->id) != pids.end()){
      cout << "path already included!" << endl;
      p->print();
      exit(1);
    }

    if(p->nodes[0] == p->nodes[p->nodes.size() - 1]){
      p->nodes.pop_back();
    }

    if(p->nodes.size() % 2 != 0){
      cout << "Incorrect number of nodes in path " << p->id + 1 << endl;
      p->print();
      write_path(p, cout);
      exit(1);
    }  

    if(!p->is_circle){
      if(p->nodes.size() != p->edges.size() + 1){
        cout << "Incorrect number of edges in linear path " << p->id + 1 << endl;
        p->print();
        write_path(p, cout);
        exit(1);        
      }  
    }else{
       if(p->nodes.size() != p->edges.size()){
        cout << "Incorrect number of edges in circular path " << p->id + 1 << endl;
        p->print();
        write_path(p, cout);
        exit(1);
      }
    }  
  }

  // # defines paths that start and end on a telomere prior to S phase
  // # i.e. these paths are fully connected upon completion of G1
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
        p->is_circle = true;
      }
      if(breakpoints[p->nodes[size-1]]->is_end){
        p->type = COMPLETE; // end nodes should also be telomere
      }else{
        p->type = PTEL; // pTel
      }

      validate_path(p);
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
        p->is_circle = true;
      }
      if(breakpoints[p->nodes[size-1]]->is_end){
        p->type = COMPLETE;
      }else{
        p->type = QTEL;
      }

      validate_path(p);
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
        p->is_circle = true;
      }
      assert(!breakpoints[p->nodes[size-1]]->is_end);
      
      validate_path(p);
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
        p->is_circle = true;
      }
      assert(!breakpoints[p->nodes[size-1]]->is_end);

      validate_path(p);
      if(verbose > 1){
        cout << "path after validation\n";
        p->print();
      }      
      paths[p->id] = p;
    }

    // all breakpoints should be in some path
    for(auto am : adjacencies){
      adjacency* a = am.second;
      if(a->path_ID < 0){
        cout << "\nwrong path connection!" << endl;
        a->print();
        exit(-1);
      }
    }
  }


// update breakpoint IDs at a duplicated adjacency with duplicated breakpoints
void update_adj_junc(adjacency*  adj_copy_prev, breakpoint* junc_copy_prev, breakpoint* junc_copy, int verbose = 0){
    if(adj_copy_prev->is_inverted){
      if(verbose > 1){
        cout << "Inverted adjacency" << endl;
      }
      adj_copy_prev->junc_id2 = junc_copy_prev->id;
      adj_copy_prev->junc_id1 = junc_copy->id;
    }else{
      adj_copy_prev->junc_id1 = junc_copy_prev->id;
      adj_copy_prev->junc_id2 = junc_copy->id;
    }
}


// update the other neighbor of last breakpoint in the duplicated path
void update_end_junc(int last_jid_orig, int last_jid) {
  if(!(breakpoints[last_jid_orig]->right_jid == breakpoints[last_jid]->right_jid ||
  breakpoints[last_jid_orig]->left_jid == breakpoints[last_jid]->left_jid)){
    cout << "Weird breakpoint neighbours!" << endl;
    breakpoints[last_jid_orig]->print();
    breakpoints[last_jid]->print();
    exit(1);
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
    adjacency* adj_var = new adjacency(cell_ID, aid, p.id, last_jid_orig, last_jid, VAR, NONTEL, H2HINV);
    adjacencies[adj_var->id] = adj_var;

    // one neighbor of the copied breakpoint must have pointed to its neighbour copy
    // only the neighbor not updated will be equal to original node's neighbor
    update_end_junc(last_jid_orig, last_jid);

    if((breakpoints[last_jid_orig]->left_jid < 0 || breakpoints[last_jid_orig]->right_jid < 0)){
      cout << "Weird breakpoint neighbours after updating!" << endl;
      breakpoints[last_jid_orig]->print();
      breakpoints[last_jid]->print();
      exit(1);
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
      adjacency* adj_var = new adjacency(cell_ID, aid, p.id, first_jid_orig, first_jid, VAR, NONTEL, T2TINV);
      adjacencies[adj_var->id] = adj_var;

      update_end_junc(first_jid_orig, first_jid);

      assert(breakpoints[first_jid_orig]->left_jid >= 0 && breakpoints[first_jid_orig]->right_jid >= 0);
      assert(breakpoints[first_jid]->left_jid >= 0 && breakpoints[first_jid]->right_jid >= 0);

      breakpoints[first_jid_orig]->is_repaired = true;
      breakpoints[first_jid]->is_repaired = true;

      p.edges.push_back(adj_var->id);
      p.is_circle = true;
    }    
  }else{
    assert(p.type == QTEL);  // qTel -- connect first breakpoint, path start from qTel
    int last_jid_orig = p.nodes[p.nodes.size() - 1];
    int last_jid = junc_copy->id;

    int aid = adjacencies.rbegin()->first + 1;
    adjacency* adj_var = new adjacency(cell_ID, aid, p.id, last_jid_orig, last_jid, VAR, NONTEL, T2TINV);
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


// "nDSB", "nDel", "nInv", "nIns", "nDup", 'cycleID', 'nBiasedChroms', 'nMitosisBreaks'
// "chr", "nDSB", "nOsc", "nDel", "nIns", "nInv", "nDup"
void get_summary_stats(/* arguments */) {
  // initialize for each chr
  vector<int> stat_total;
  vector<int> stat_chr;

  // writing SV deletions & duplications for shatterseek and circos
  for(int i = 0; i < NUM_CHR; i++){
    // write SV deletions - for Circos (by haplotype)
    // sv_type=DEL
    // write SV duplications - for shatterseek
    // write SV deletions - for shatterseek
  }

  // count oscillating CNs
  for(int i = 0; i < NUM_CHR; i++){

  }

}


void get_unique_interval(int verbose = 0){
  cn_by_chr_hap.clear();
  if(verbose > 1) cout << "Collecting all the unique genomic intervals\n";
  // get CN for each interval
  // chr, start, end, haplotype : cn
  map<haplotype_pos, int> cn_by_pos;  // default value is 0
  for(auto am : adjacencies){
    adjacency* a = am.second;
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

        if(intl_prev.end == intl_curr.start && intl_prev.cn == intl_curr.cn){
          interval intl_merged{intl_prev.start, intl_curr.end, intl_prev.cn};
          if(verbose > 1) cout << "interval after merging with previous one\t" << intl_merged.start << "\t" << intl_merged.end << "\t" << intl_merged.cn << endl;
          cn_by_chr_hap_merged[cnp.first].pop_back();
          cn_by_chr_hap_merged[cnp.first].push_back(intl_merged);
          intl_prev = intl_merged;
        }else{
          if(verbose > 1) cout << "interval without merging with previous one\t" << intl_curr.start << "\t" << intl_curr.end << "\t" << intl_curr.cn << endl;
          cn_by_chr_hap_merged[cnp.first].push_back(intl_curr);
          intl_prev = intl_curr;
        }
      }
    }

    if(verbose > 1){
      cout << "Merged intervals" << endl;
      cout << "cell\tchr\thaplotype\tstart\tend\tJID_start\tJID_end\tCN\n";
      for(auto cnp : cn_by_chr_hap_merged){
        int chr = cnp.first.first;
        int haplotype = cnp.first.second;
        // cout << cell_ID << "\t" << chr << "\t" << haplotype << endl;
        vector<interval> intls = cnp.second;
        for(int i = 0; i < intls.size(); i++){
          interval intl = intls[i];
          cout << cell_ID << "\t" << chr + 1 << "\t" << get_haplotype_string(haplotype) << "\t" << intl.start << "\t" << intl.end  << "\t" << intl.jid_start << "\t" << intl.end << "\t"  << intl.cn << "\n";
        }
      }
    }
  }


  // To avoid segments completely lost in one cell not shown, add chr end when needed (TODOs)
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

  // Compute allele-specific and total copy number
  void calculate_segment_cn(map<int, set<int>>& bps_by_chr, int verbose = 0){
    // verbose = 1;
    this->segments.clear();
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
      for (; it != bps.end(); it++){
        int pos2 = *it;
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
      }

      // find cns for each interval in each haplotype
      if(verbose > 1) cout << "CNs for haplotype A" << endl;
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
        this->segments.push_back(s);
        if(verbose > 1) s->print();
      }
    }
  }
};
