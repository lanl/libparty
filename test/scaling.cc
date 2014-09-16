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

#include <mpi.h>
#include <iostream>
#include "libparty.hh"
#include <chrono>
#include <iomanip>
#include <fstream>

#define start_timer(id) auto id##_begin = std::chrono::high_resolution_clock::now();
#define stop_timer(id,name) auto id##_end = std::chrono::high_resolution_clock::now(); \
std::cout<< name << ": " << std::setprecision(9)<< std::chrono::duration_cast<std::chrono::nanoseconds>( id##_end - id##_begin).count()/1000. << std::endl;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  bbox lbox;
  lbox.resize(6);

  double xstride = 80./np;
  lbox[0] = -40+xstride*rank;
  lbox[1] = lbox[0]+xstride;
  lbox[2] = -40;
  lbox[3] = 40;
  lbox[4] = -40;
  lbox[5] = 40;

  std::ofstream faos;
  faos.open("nnaos.txt");
  std::ofstream fsoa;
  fsoa.open("nnsoa.txt");
  std::ofstream faosoa;
  faosoa.open("nnaosoa.txt");
int nstep = 1;

for(int n = 10; n < 45; ++n) {
  int nump = pow(10.0,(double)n/5.0);
  part_sys parts = part_sys(nump);
  parts.syncbboxes(lbox);
  std::vector<double> searchv;
double bin = 10*0.5/nump;
  for(int i = 0; i<nump; ++i) searchv.push_back(bin);
  part_sys_soa pdata;
  part_sys_aosoa pdata2;
  for(int i = 0; i < parts.partv.size(); ++i) {
    parts.partv[i].nn_h = bin;
    pdata += parts.partv[i];
    pdata2 += parts.partv[i];
  }

parts.build_hashtable(bin);
pdata.build_hashtable(bin);
pdata2.build_hashtable(bin);

auto begin = std::chrono::high_resolution_clock::now();
//parts.build_hashtable(bin);
parts.find_neighbors(searchv,1);
//for(int i = 0; i<nstep; ++i) parts.step(0.5);
auto end = std::chrono::high_resolution_clock::now();
faos << std::setprecision(9)<< nump << "," <<  std::chrono::duration_cast<std::chrono::nanoseconds>( end - begin).count()/1000000000./nstep << std::endl;

 begin = std::chrono::high_resolution_clock::now();
//pdata.build_hashtable(bin);
pdata.find_neighbors(1);
//for(int i = 0; i<nstep; ++i) pdata.step(0.5);
 end = std::chrono::high_resolution_clock::now();
fsoa << std::setprecision(9)<< nump << "," <<  std::chrono::duration_cast<std::chrono::nanoseconds>( end - begin).count()/1000000000./nstep << std::endl;

 begin = std::chrono::high_resolution_clock::now();
//pdata2.build_hashtable(bin);
pdata2.find_neighbors(1);
//for(int i = 0; i<nstep; ++i) pdata2.step(0.5);
 end = std::chrono::high_resolution_clock::now();
faosoa << std::setprecision(9)<< nump << "," <<  std::chrono::duration_cast<std::chrono::nanoseconds>( end - begin).count()/1000000000./nstep << std::endl;
}
  faos.close();
  fsoa.close();
  faosoa.close();

//  parts.build_hashtable(0.5);
//  parts.find_neighbors(searchv, 100);
//  parts.output_h5part();
//  start_timer(b)
//for(int i = 0; i<20; ++i) {
//  parts.calc_density();
//  parts.output_h5part();
//}
//stop_timer(b,"aos density calc")

//  start_timer(a)
//    for(int i = 0; i<100000; ++i) parts.step(0.5);
//  stop_timer(a,"aos step")

//  pdata.build_hashtable(0.5);
//  pdata.find_neighbors(100);
//  pdata.output_h5part();
//start_timer(d)
//  for(int i = 0; i<20; ++i) {
//  pdata.calc_density();
//  pdata.output_h5part();
//  }
//stop_timer(d,"soa density calc")
/*
start_timer(c)
  for(int i = 0; i<100; ++i) pdata.step(0.5);
stop_timer(c,"soa step")

//  pdata2.build_hashtable(0.5);
//  pdata2.find_neighbors(100);
  pdata2.output_h5part();
start_timer(f)
for(int i = 0; i<20; ++i) {
  pdata2.calc_density();
  pdata2.output_h5part();
}
stop_timer(f,"aosoa density calc")

start_timer(e)
  for(int i = 0; i<100; ++i) pdata2.step(0.5);
stop_timer(e,"aosoa step")
*/
//  std::cout << pdata;
  MPI_Finalize();
  return 0;
}
