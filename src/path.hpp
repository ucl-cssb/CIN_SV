#ifndef PATH_HPP
#define PATH_HPP


#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

#include "util.hpp"


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


  // add one edge into path p, return success status
  bool update_path_by_adj(path* p, const breakpoint* js, breakpoint* jn, map<int, adjacency*>& adjacencies, int& prev_atype, int verbose = 0){
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

    if(prev_atype == adj->type){
      cout << "The adjacencies in a path must alternate with different types! " << endl;
      adj->print();
      exit(FAIL);
    }
    prev_atype = adj->type;

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
        if(verbose > 1){
          cout << "one more centromere" << endl;
        }
        p->n_centromere += 1;
        adj->is_centromeric = true;
    }

    if(verbose > 1){
      cout << "\nupdate path " << p->id + 1 << " with " << get_adj_type_string(adj->type) << " adjacency " << aid << " connected by breakpoint " << js->id << " and " << jn->id << endl;
      adj->print();
    }

    return true;
  }


#endif