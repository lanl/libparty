/*
 * Copyright 2014.  Los Alamos National Security, LLC. This material was produced 
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
 * for the U.S. Department of Energy. The U.S. Government has rights to use, 
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR 
 * LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, 
 * OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is 
 * modified to produce derivative works, such modified software should be 
 * clearly marked, so as not to confuse it with the version available from LANL.   
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not 
 * use this file except in compliance with the License. You may obtain a copy 
 * of the License at http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed 
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the 
 * specific language governing permissions and limitations under the License.
 * 
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *
 * LibParty, Version 1.x, Copyright number C114122, LA-CC-14-086
 *
 * Authors:
 *    Matthew Kinsey, kinsey@lanl.gov
 *    Verinder Rana, vrana@lanl.gov
 *    Bob Robey, XCP-2, LANL, brobey@lanl.gov
 *    Jeremy Sauer, EES-16, LANL, jsauer@lanl.gov
 *    Jon Reisner, XCP-4 LANL, reisner@lanl.gov
 */ 

#ifndef VECL
  #define VECL 4
#endif
#ifndef ALIGNL 
  #define ALIGNL 32
#endif


#ifndef PARTICLEAOSOA_H
#define PARTICLEAOSOA_H

#include <vector>
#include <unordered_map>
#include <cstdint>
//#include "util.hh"
#include <iostream>
#include "libparty.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>
#include <iostream>
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include <bitset>

struct part : CCpart
{
  public:
   // basic part info   
//  double pos[3];
//  double vel[3];
//  double mass;
//   double nn_h;
//   bool ghost;// = false;
//  int id;// = -1;
//   double lifetime = 0.0;
//  int proc;
//  bool del = false;

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
// NOTE!: YOU HAVE TO SERIALIZE _EVERY_ VARIABLE IN THE CLASS OR IT BREAKS
    ar & pos & vel & mass & id & del & proc & nn_h & rho;
  }
};

BOOST_IS_MPI_DATATYPE(part)

struct partvec {
  __declspec(align(ALIGNL)) double x[VECL], y[VECL], z[VECL];
  __declspec(align(ALIGNL)) double vx[VECL], vy[VECL], vz[VECL];
  __declspec(align(ALIGNL)) double mass[VECL], nn_h[VECL], rho[VECL];
//  int id[VECL], proc[VECL], del[VECL];
};
class part_sys_aosoa
{
  public:
  std::vector<partvec> partv;
  std::vector<std::bitset<VECL>> partvfill;
  void step( double dt );
  void build_hashtable(double bin_size);
  double current_bin_size;
  int table_size;
  tbb::concurrent_unordered_multimap<unsigned int,int> hash_table;
  tbb::concurrent_vector<tbb::concurrent_vector<int>> find_neighbors(int nmax);
  double nn_cutoff = pow(10.0,-6);
  double max_nn_h = pow((1/(4*3.1415*nn_cutoff)),(1.0/3.0));

  //IO
  int iter = 0;
  void output_h5part();
  
  void operator+=(const part &p) {
    if(partv.empty()) {
      partvec pv;
      pv.x[0] = p.pos[0];
      pv.y[0] = p.pos[1];
      pv.z[0] = p.pos[2];
      pv.vx[0] = p.vel[0];
      pv.vy[0] = p.vel[1];
      pv.vz[0] = p.vel[2];
      pv.nn_h[0] = p.nn_h;
      std::bitset<VECL> fill;
      fill[0]=1;
      partv.push_back(pv);
      partvfill.push_back(fill);
    } else {
      for(int i = 0; i < partv.size(); ++i) {
        for(int b = 0; b < VECL; ++b) {
          if(!partvfill[i][b]) {
            partv[i].x[b] = p.pos[0];
            partv[i].y[b] = p.pos[1];
            partv[i].z[b] = p.pos[2];
            partv[i].vx[b] = p.vel[0];
            partv[i].vy[b] = p.vel[1];
            partv[i].vz[b] = p.vel[2];
            partv[i].nn_h[b] = p.nn_h;
            partvfill[i][b] = 1;
            return;
          }
        }
      }
      //we didnt find an empty spot :(
      partvec pv;
      pv.x[0] = p.pos[0];
      pv.y[0] = p.pos[1];
      pv.z[0] = p.pos[2];
      pv.vx[0] = p.vel[0];
      pv.vy[0] = p.vel[1];
      pv.vz[0] = p.vel[2];
      pv.nn_h[0] = p.nn_h;
      std::bitset<VECL> fill;
      fill[0]=1;
      partv.push_back(pv);
      partvfill.push_back(fill);
    }
  }

  void operator-=(const int &ind) {
    int i, v;
    i = ind/VECL;
    v = ind%VECL;
    partvfill[i][v] = 0;
    if(partvfill[i].none()) {
      partvfill.erase(partvfill.begin()+i);
      partv.erase(partv.begin()+i);
    }
  }

};


typedef std::vector<double> bbox;
typedef std::vector<bbox> bboxv;

typedef std::vector<part>::iterator partIter;
typedef std::vector<part>::const_iterator c_partIter;
typedef std::vector<part*>::iterator p_partIter;

typedef std::pair<unsigned int,uintptr_t> keypart_p;
typedef std::pair<unsigned int,uintptr_t> ppart_p;
typedef std::pair<unsigned int,std::pair<unsigned int,uintptr_t>> keyprocpart_p;

#endif
