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

void part_sys_soa::build_hashtable(double bin_size) {
  hash_table.clear();
  keypart_p tableval;
  //std::cout << "building table with h=" << bin_size << " " << std::endl;
  int np = x.size();
  keypart_p keyparts[np];
  table_size = next_prime(2*(x.size()));

#ifdef _OPENMP
  #pragma omp parallel for// schedule(dynamic,100)
#endif
//  #pragma simd  
  for(int i = 0; i < np; ++i) {
    unsigned int key = (unsigned int)((((int)(std::floor((x[i]+1000)/bin_size)) * 73856093)
                       ^ ((int)(std::floor((y[i]+1000)/bin_size)) * 19349663)
                       ^ ((int)(std::floor((z[i]+1000)/bin_size)) * 83492791))) % (2*table_size);
    keyparts[i] = std::make_pair(key,i);
    hash_table.insert(keyparts[i]);
  }
//  for(int i = 0; i < np; ++i) {
//    hash_table.insert(keyparts[i]);
//  }
  current_bin_size = bin_size;
}
tbb::concurrent_vector<tbb::concurrent_vector<int>> part_sys_soa::find_neighbors( 
//                                std::vector<part> parts,
//                                double search_dist)
//                                std::vector<double> search_distv,
                                int nmax)
{

  tbb::concurrent_vector<tbb::concurrent_vector<int>> neighbors;
  neighbors.resize(x.size());

#ifndef USE_MIC

//  if (!search_irange) search_irange++;
  #pragma ivdep
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int i = 0; i < x.size(); ++i) {
    if(ghost[i]) continue;
    std::multimap<double, int> dist_map;
    int search_irange = (int)(ceil(nn_h[i]/current_bin_size*sqrt(2)));
//    if (!search_irange) search_irange++;
//    querytable:
    for( int zind = -search_irange; zind <= search_irange; zind++){
      for( int yind = -search_irange; yind <= search_irange; yind++){
        for( int xind = -search_irange; xind <= search_irange; xind++){
          unsigned int key = (unsigned int)((((int)(std::floor((x[i]+1000)/current_bin_size)+xind) * 73856093)
                       ^ ((int)(std::floor((y[i]+1000)/current_bin_size)+yind) * 19349663)
                       ^ ((int)(std::floor((z[i]+1000)/current_bin_size)+zind) * 83492791))) % (2*table_size);

          auto range = hash_table.equal_range(key);
          for (auto it = range.first;  it != range.second;  ++it) {
//            #pragma omp critical
            double dist2 = pow(x[i] - x[it->second],2.0)
              + pow(y[i] - y[it->second],2.0)
              + pow(z[i] - z[it->second],2.0);
            dist_map.insert(std::make_pair(sqrt(dist2),it->second));
          }
        }
      }
    }
    for (std::multimap<double, int>::iterator it=dist_map.begin(); it!=dist_map.end(); ++it) {
      if(it->first == 0.0) continue;
      if(neighbors[i].size() == nmax || it->first > max_nn_h) break;
      neighbors[i].push_back(it->second);
    }
//    if( neighbors[i].size() != nmax && dist_map.end()->first < max_nn_h) {
//      dist_map.clear();
//      search_irange++;
//      goto querytable;
//    }
  }
#endif
  return neighbors;
}

