#ifndef SEGMENT_HPP
#define SEGMENT_HPP


#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

#include "util.hpp"



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


#endif