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

#include "libparty.hh"
#include "hashutil.hh"
#include <iostream>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

void part_sys::build_hashtable(double bin_size) {
  hash_table.clear();
  keypart_p tableval;
  //std::cout << "building table with h=" << bin_size << " " << std::endl;

  table_size = next_prime(2*(partv.size()+ghostv.size()));
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i = 0; i < partv.size(); ++i) {
    part* p = &partv[i];
    unsigned int key = hash_pos( p->pos, bin_size, table_size);
//    tableval = std::make_pair(key,(uintptr_t)&partv[i]);
    tableval = std::make_pair(key,i);
//    #pragma omp critical
    hash_table.insert(tableval);
  }
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i = 0; i < ghostv.size(); ++i) {
    part* p = &ghostv[i];
    unsigned int key = hash_pos( p->pos, bin_size, table_size);
    tableval = std::make_pair(key,i);
//    #pragma omp critical
    hash_table.insert(tableval);
  }
  current_bin_size = bin_size;
}

tbb::concurrent_vector<tbb::concurrent_vector<int>> part_sys::find_neighbors( 
//                                std::vector<part> parts,
//                                double search_dist)
                                std::vector<double> search_distv,
                                int nmax)
{

  tbb::concurrent_vector<tbb::concurrent_vector<int>> neighbors;
  neighbors.resize(partv.size());

#ifndef USE_MIC

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i = 0; i < partv.size(); ++i) {
    std::multimap<double, int> dist_map;
    part* p = &partv[i];
    int search_irange = (int)(ceil(p->nn_h/current_bin_size*sqrt(2)));
//    if (!search_irange) search_irange++;
    for( int zind = -search_irange; zind <= search_irange; zind++){
      for( int yind = -search_irange; yind <= search_irange; yind++){
        for( int xind = -search_irange; xind <= search_irange; xind++){
          double keypos[3];
          keypos[0] = p->pos[0] + xind*current_bin_size;
          keypos[1] = p->pos[1] + yind*current_bin_size;
          keypos[2] = p->pos[2] + zind*current_bin_size;
          unsigned int key = hash_pos( keypos, current_bin_size, table_size);
          auto range = hash_table.equal_range(key);
          for (auto it = range.first;  it != range.second;  ++it) {
//            #pragma omp critical
//            double dist2 = pow(p->pos[0] - ((part*)it->second)->pos[0],2.0)
//              + pow(p->pos[1] - ((part*)it->second)->pos[1],2.0)
//              + pow(p->pos[2] - ((part*)it->second)->pos[2],2.0);
            double dist2 = pow(p->pos[0] - partv[it->second].pos[0],2.0)
              + pow(p->pos[1] - partv[it->second].pos[1],2.0)
              + pow(p->pos[2] - partv[it->second].pos[2],2.0);
            dist_map.insert(std::make_pair(sqrt(dist2),it->second));
          }
        }
      }
    }
//  #pragma omp parallel for
    for (std::multimap<double, int>::iterator it=dist_map.begin(); it!=dist_map.end(); ++it) {
      if(it->first == 0.0) continue;
      if(neighbors[i].size() == nmax || it->first > max_nn_h) break;
      neighbors[i].push_back(it->second);
    }
  }
#endif
  return neighbors;
}

