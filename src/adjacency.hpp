#ifndef ADJACENCY_HPP
#define ADJACENCY_HPP


#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

#include "util.hpp"


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
    cout << "Adjacency " << id << " in cell " << cell_ID << " at path " << path_ID + 1 << " with left breakpoint " << junc_id1 << " at " << breakpoints[junc_id1]->chr + 1 << "-" << get_haplotype_string(breakpoints[junc_id1]->haplotype) << "-" << breakpoints[junc_id1]->pos << " and right breakpoint " << junc_id2 << " at " << breakpoints[junc_id2]->chr + 1 << "-" << get_haplotype_string(breakpoints[junc_id2]->haplotype) << "-" << breakpoints[junc_id2]->pos << ", " << cent << ", " << invt << ", " << get_telomere_type_string(telomeric_type) << ", " << get_sv_type_string(sv_type) << endl;
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

// when j1 and j2 are on the same chromosome, need to ensure the coordinates of j1 is smaller than j2
SV_type set_var_type(breakpoint* j1, breakpoint* j2, int verbose = 0) {
    SV_type sv_type = NONE;

    if(j1->side == HEAD && j2->side == TAIL){
      assert(j1->left_jid >= 0);
      assert(j2->right_jid >= 0);
      j1->right_jid = j2->id;
      j2->left_jid = j1->id;
      // DEL (deletion-like; +/-), intervals also need to have CN < 2, updated later
      sv_type = DEL;
      // interval in the reference genome is lost
      // || j1->haplotype != j2->haplotype, ignore haplotype as it is hard to determine from real data
      if(j1->chr != j2->chr){
        sv_type = BND;
      }
    }else if(j1->side == HEAD && j2->side == HEAD){ //h2hINV
      assert(j1->left_jid >= 0);
      assert(j2->left_jid >= 0);
      j1->right_jid = j2->id;
      j2->right_jid = j1->id;
      // h2hINV (head-to-head inversion; +/+)
      sv_type = H2HINV;
      if(j1->chr != j2->chr){
        sv_type = BND;
      }
    }else if(j1->side == TAIL && j2->side == HEAD){
      assert(j1->right_jid >= 0);
      assert(j2->left_jid >= 0);
      j1->left_jid = j2->id;
      j2->right_jid = j1->id;
      // DUP (tandom duplication-like; -/+), intervals also need to have CN > 2, updated later
      // if(is_interval_overlap(j1, j2)){
      sv_type = DUP;
      // }
      if(j1->chr != j2->chr){
        sv_type = BND;
      }
    }else{  //if(j1->side == TAIL && j2->side == TAIL)
      assert(j1->right_jid >= 0);
      assert(j2->right_jid >= 0);
      j1->left_jid = j2->id;
      j2->left_jid = j1->id;
      // t2tINV (tail-to-tail inversion; -/-)
      sv_type = T2TINV;
      if(j1->chr != j2->chr){
        sv_type = BND;
      }
    }

    return sv_type;
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



#endif