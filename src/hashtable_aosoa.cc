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

void part_sys_aosoa::build_hashtable(double bin_size) {
  hash_table.clear();
  keypart_p tableval;
  //std::cout << "building table with h=" << bin_size << " " << std::endl;
  int np;

  for(int i = 0; i < partvfill.size(); ++i) {
#pragma omp parallel for
  for(int v = 0; v < VECL; ++v) np += partvfill[i][v];
}
  keypart_p keyparts;
  __declspec(align(VECL)) unsigned int keys[VECL];
  table_size = next_prime(2*np);

  #pragma omp parallel for// schedule(dynamic,100)
//  #pragma simd  
  for(int i = 0; i < partv.size(); ++i) {
    keys[:] = (unsigned int)((((int)(std::floor((partv[i].x[:]+1000)/bin_size)) * 73856093)
                       ^ ((int)(std::floor((partv[i].y[:]+1000)/bin_size)) * 19349663)
                       ^ ((int)(std::floor((partv[i].z[:]+1000)/bin_size)) * 83492791))) % (2*table_size);
//#pragma omp parallel for                       
    for(int v = 0; v < VECL; ++v) {
    if(partvfill[i][v]){
      keypart_p keypart = std::make_pair(keys[v],i*VECL+v);
      hash_table.insert(keypart);
    }
  }
  }
//  for(int i = 0; i < np; ++i) {
//    keypart_p keypart = std::make_pair(keys[i],i);
//    hash_table.insert(keypart);
//  }
  current_bin_size = bin_size;
}

tbb::concurrent_vector<tbb::concurrent_vector<int>> part_sys_aosoa::find_neighbors( 
//                                std::vector<part> parts,
//                                double search_dist)
//                                std::vector<double> search_distv,
                                int nmax)
{

  tbb::concurrent_vector<tbb::concurrent_vector<int>> neighbors;
  neighbors.resize(partv.size()*VECL);

//  int search_irange = (int)(ceil(search_dist/current_bin_size*sqrt(2)));
//  if (!search_irange) search_irange++;
  #pragma ivdep
  #pragma omp parallel for
  for(int i = 0; i < partv.size(); ++i) {
    int search_irange[VECL];
    search_irange[:] = (int)(ceil(partv[i].nn_h[:]/current_bin_size*sqrt(2)));
//    if (!search_irange) search_irange++;
//    #pragma omp parallel for
//#pragma omp parallel for
    for(int v = 0; v < VECL; ++v) {
      std::multimap<double, int> dist_map;
      unsigned short int ncells = pow(2*search_irange[v]+1,3);
      __declspec(align(ALIGNL)) unsigned int keys[ncells];
      __declspec(align(ALIGNL)) int xind[ncells], yind[ncells], zind[ncells];
      unsigned short int ci = 0;

//#pragma omp parallel for collapse(3)
      for( int zi = -search_irange[v]; zi <= search_irange[v]; zi++){
        for( int yi = -search_irange[v]; yi <= search_irange[v]; yi++){
          for( int xi = -search_irange[v]; xi <= search_irange[v]; xi++){
            xind[ci] = xi;
            yind[ci] = yi;
            zind[ci] = zi;
            ci++;
          }
        }
      }

      keys[:] = (unsigned int)((((int)(std::floor((partv[i].x[v]+1000)/current_bin_size)+xind[:]) * 73856093)
                       ^ ((int)(std::floor((partv[i].y[v]+1000)/current_bin_size)+yind[:]) * 19349663)
                       ^ ((int)(std::floor((partv[i].z[v]+1000)/current_bin_size)+zind[:]) * 83492791))) % (2*table_size);
//#pragma omp parallel for
      for(int ki = 0; ki < ncells; ++ki) {
        auto range = hash_table.equal_range(keys[ki]);
//#pragma omp parallel for
        for (auto it = range.first;  it != range.second;  ++it) {
//            #pragma omp critical
          int ni, nv;
          ni = (int)(it->second/VECL);
          nv = it->second%VECL;
          double dist2 = pow(partv[i].x[v] - partv[ni].x[nv],2.0)
              + pow(partv[i].y[v] - partv[ni].y[nv],2.0)
              + pow(partv[i].z[v] - partv[ni].z[nv],2.0);
          dist_map.insert(std::make_pair(sqrt(dist2),it->second));
        }
      }
//#pragma omp parallel for
    for (std::multimap<double, int>::iterator it=dist_map.begin(); it!=dist_map.end(); ++it) {
//    for(int di = 0; di < dist_map.size(); ++di) {
      if(it->first == 0.0 || partvfill[i][v] ==0) continue;
      if(neighbors[i*VECL+v].size() == nmax || it->first > max_nn_h) break;
      neighbors[i*VECL+v].push_back(it->second);
    }
  }
  }
  return neighbors;
}

//tbb::concurrent_vector<int> part_sys_aosoa::find_neighbors(
      //short int range[2*search_irange[v]+1];
      ////      int range = 2*search_irange[v]+1;
      ////      int ri = 0;
      ////for( int ind = -search_irange[v]; ind <= search_irange[v]; ++ind){
      ////xind[ri:range*range*range:range] = ind;
      ////for( int r = 0; r < range; ++r) yind[r+range*ri:r + range*(ri+1)] = ind;
      ////zind[range*range*ri:range*range*(ri+1)] = ind;
      ////ri++;
      ////}
      ////      for(int ki = 0; ki < ncells; ++ki) {
      ////std::cout << xind[ki] << " " << yind[ki] << " " << zind[ki] << std::endl;
      ////}
      //
