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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef USE_MIC
#include <boost/mpi/operations.hpp>
#endif
#ifdef HAVE_H5PART
#define MPI_Comm h5part_mpi_comm
#include <H5Part.h>
#endif
#undef MPI_Comm
#include <sys/stat.h>

inline bool file_exists (const char* name) {
  struct stat buffer;   
  return (stat (name, &buffer) == 0); 
}

void part_sys_aosoa::output_h5part()
{
#ifdef HAVE_H5PART
  int totalnumpart;
  int localnumpart;// = partv.size();
  for(int i = 0; i < partvfill.size(); ++i) for(int v = 0; v < VECL; ++v) localnumpart += partvfill[i][v];
#if(0)
  localnumpart += ghostv.size();
#endif
#ifndef USE_MIC
  boost::mpi::all_reduce(world, localnumpart, totalnumpart, std::plus<int>());
#endif

  std::vector<int> numpart_per_proc;
  numpart_per_proc.resize(world.size());

#ifndef USE_MIC
  boost::mpi::all_gather(world, localnumpart, numpart_per_proc);
#endif
  int accum_part = 0;

//  H5PartFile* file;

  const char *filename = "particles_aosoa.h5part";//output_filename().c_str();

  for(int proc=0; proc < world.size(); ++proc) // Loop over procs in a known way to preserve order in th
  {
    if (0/*proc != mproc || localnumpart==0*/) {
      accum_part += numpart_per_proc[proc];
    } else {

      H5PartFile* file;
//      if(!file_exists(filename)) {
      if(iter==0&&proc==0) {
        file = H5PartOpenFile(filename,H5PART_WRITE);
//        H5PartSetStep(file, 0);
      } else {
        file = H5PartOpenFile(filename,H5PART_APPEND);
//        h5part_int64_t currentstep = H5PartGetNumSteps(file);
//        if(proc==0) currentstep++;
//        H5PartSetStep(file, currentstep-1);
      }
      H5PartSetStep(file,iter);
      //TODO something is wrong with the view
      H5PartSetNumParticles(file, totalnumpart);
      H5PartSetView(file,accum_part,accum_part+localnumpart-1);
      //if(debug)
      //std::cout << totalnumpart << " " << accum_part << " " << localnumpart << std::endl;
      int np = localnumpart;//partv.size();

      double p_x[np], p_y[np], p_z[np], p_vx[np], p_vy[np], p_vz[np], p_m[np], p_id[np], p_nn_h[np], p_rho[np];
//std::cout << "outputting " << np << std::endl;
      int arrayind = 0;
      for( int i = 0; i < partv.size(); ++i ) {
        for(int v = 0; v < VECL; ++v) {
          if(partvfill[i][v]) {
            p_x[arrayind] = partv[i].x[v];
            p_y[arrayind] = partv[i].y[v];
            p_z[arrayind] = partv[i].z[v];
            p_vx[arrayind] = partv[i].vx[v];
            p_vy[arrayind] = partv[i].vy[v];
            p_vz[arrayind] = partv[i].vz[v];
            p_m[arrayind] = partv[i].mass[v];
//          p_id[i*VECL:(i+1)*VECL] = partv[i].id[:];
            p_nn_h[arrayind] = partv[i].nn_h[v];
            p_rho[arrayind] = partv[i].rho[v];
            arrayind++;
          }
        }
//        }
      }
#if(0)
      for( partIter it = ghostv.begin(); it != ghostv.end(); ++it ) {
        int i = std::distance(ghostv.begin(),it) + partv.size();
        p_x[i] = it->pos[0];
        p_y[i] = it->pos[1];
        p_z[i] = it->pos[2];
        p_vx[i] = it->vel[0];
        p_vy[i] = it->vel[1];
        p_vz[i] = it->vel[2];
        p_m[i] = it->mass;
        p_id[i] = it->id;
        p_nn_h[i] = it->nn_h;
        p_rho[i] = it->rho;
      }
#endif
      H5PartWriteDataFloat64(file,"x", p_x);
      H5PartWriteDataFloat64(file,"y", p_y);
      H5PartWriteDataFloat64(file,"z", p_z);
      H5PartWriteDataFloat64(file,"v^x", p_vx);
      H5PartWriteDataFloat64(file,"v^y", p_vy);
      H5PartWriteDataFloat64(file,"v^z", p_vz);
      H5PartWriteDataFloat64(file,"m", p_m);
//      H5PartWriteDataFloat64(file,"id", p_id);
      H5PartWriteDataFloat64(file,"nn_h", p_nn_h);
      H5PartWriteDataFloat64(file,"rho", p_rho);
      H5PartSetView(file,-1,-1);
      H5PartCloseFile(file);
    }
    world.barrier();
  }
  iter++;
#endif
}
