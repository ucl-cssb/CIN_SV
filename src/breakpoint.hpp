#ifndef BREAKPOINT_HPP
#define BREAKPOINT_HPP


#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

#include "util.hpp"


// a breakpoint in the chromosome, introduced by a DSB, haplotype-specific,
// node in the genome graph
// one instance for each copy
class breakpoint{
public:
  int id;   // keep increasing when new breakpoints are added
  int cell_ID;
  int path_ID;  // a breakpoint can only belong to one path
  int chr;
  int pos;
  int side;   // 0 left or 1 right, strand in RCK, orientation
  int haplotype;
  // int cn;  // there may be multiple copies due to replication??
  bool is_repaired; // end breakpoint has only one connection (not necessarily telomere)
  // bool is_visited;  // whether in a path or not, used to determine circular paths
  int left_jid;
  int right_jid;
  bool is_end;  // whether it is at the end of genome or not, correspond to telomere when DSBs are not introduced to telomere regions
  // bool is_centromeric; // whether it overlaps with centromere, correspond to centromere when DSBs are not introduced to centromere regions

  breakpoint(){

  }


  breakpoint(int cell_ID, int id, int chr, int pos, int side, int haplotype, bool is_end, bool is_repaired){
    this->path_ID = -1;
    this->left_jid = -1;
    this->right_jid = -1;

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


// get position and IDs of breakpoints in the normal genome
// chr_bps:  breakpoints on a chromosome in a sorted order
// bp_jids: junction IDs for each bp, specified by chr, haplotype, loc, side
void intialize_chr_bp_jid(map<pair<int, int>, vector<int>>& chr_bps, map<pos_hap, int>& bp_jids){
  for(int chr = 0; chr < NUM_CHR; chr++){
    int haplotype = 0;
    pair<int, int> chp0(chr, haplotype);
    chr_bps[chp0].push_back(1);
    chr_bps[chp0].push_back(CHR_LENGTHS[chr]);

    int chr_jid = 4 * (chr);  // chr starts from 0
    pos_hap p0s = {chr, haplotype, 1, TAIL};
    bp_jids[p0s] = chr_jid;
    pos_hap p0e = {chr, haplotype, CHR_LENGTHS[chr], HEAD};
    bp_jids[p0e] = chr_jid + 1;

    haplotype = 1;
    pair<int, int> chp1(chr, haplotype);
    chr_bps[chp1].push_back(1);
    chr_bps[chp1].push_back(CHR_LENGTHS[chr]); 

    chr_jid = 4 * (chr) + 2;
    pos_hap p1s = {chr, haplotype, 1, TAIL};
    bp_jids[p1s] = chr_jid;
    pos_hap p1e = {chr, haplotype, CHR_LENGTHS[chr], HEAD};
    bp_jids[p1e] = chr_jid + 1;        
  }
}


void print_bp_id(const map<pair<int, int>, vector<int>>& chr_bps, const map<pos_hap, int>& bp_jids){
    cout << "sorted breakpoints on each chr each haplotype: \n";
    for(auto cbp : chr_bps){
      cout << cbp.first.first + 1<< " " << cbp.first.second;
      for(auto bp : cbp.second){
        cout << " " << bp;
      }
      cout << endl;
    }
    cout << "junction IDs for each breakpoint: " << endl;
    for(auto bid: bp_jids){
      pos_hap ph = bid.first;
      cout << ph.chr + 1 << " " << ph.haplotype << " " << ph.loc << " " << ph.side << " " << bid.second << endl;
    }   
}


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

    

#endif
