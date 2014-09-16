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

#include "libparty_wrappers.h"
#include "libparty.hh"

#ifdef __cplusplus
extern "C" {
#endif
  Cpart_sys* part_sys__new (int n) {
    part_sys_soa *p = new part_sys_soa(n);
    return (Cpart_sys *)p;
  }
  void part_sys__initgauss (Cpart_sys *This) {
    ((part_sys_soa *)This)->initgauss();
  }
  int part_sys__npart (Cpart_sys *This) {
    return ((part_sys_soa *)This)->x.size();
  }
  void part_sys__delete (Cpart_sys *This) {
    delete ((part_sys_soa *)This);
  }
  void part_sys__step (Cpart_sys *This, double dt) {
    if(((part_sys_soa *)This)->have_bboxes) ((part_sys_soa *)This)->boostsend();
    if(((part_sys_soa *)This)->have_grid)((part_sys_soa *)This)->interp3dvelwithgrid(((part_sys_soa *)This)->vxptr,((part_sys_soa *)This)->vyptr, ((part_sys_soa *)This)->vzptr);
    ((part_sys_soa *)This)->step(dt);
  }
  void part_sys__print_sys (Cpart_sys *This) {
    std::cout << *(part_sys_soa *)This;
  }
  void part_sys__calc_density (Cpart_sys *This) {
    ((part_sys_soa *)This)->calc_density();
  }
  void part_sys__output_h5part (Cpart_sys *This) {
    ((part_sys_soa *)This)->output_h5part();
  }
  Cpart* part_sys__data_ptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->x[0]);
  }
  void * part_sys__xptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->x[0]);
  }
  void * part_sys__yptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->y[0]);
  }
  void * part_sys__zptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->z[0]);
  }
  void * part_sys__vxptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->vx[0]);
  }
  void * part_sys__vyptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->vy[0]);
  }
  void * part_sys__vzptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->vz[0]);
  }
  void * part_sys__massptr (Cpart_sys *This) {
    return (Cpart*)&(((part_sys_soa *)This)->mass[0]);
  }
  void part_sys__init_grid(Cpart_sys *This, double xmin, double ymin, double zmin, int nx, int ny, int nz, double xs, double ys, double zs) {
    ((part_sys_soa *)This)->grid.xmin = xmin;
    ((part_sys_soa *)This)->grid.ymin = ymin;
    ((part_sys_soa *)This)->grid.zmin = zmin;
    ((part_sys_soa *)This)->grid.nx = nx;
    ((part_sys_soa *)This)->grid.ny = ny;
    ((part_sys_soa *)This)->grid.nz = nz;
    ((part_sys_soa *)This)->grid.xstride = xs;
    ((part_sys_soa *)This)->grid.ystride = ys;
    ((part_sys_soa *)This)->grid.zstride = zs;
    bbox lbox;
    lbox.resize(6);
    lbox[0] = xmin;
    lbox[1] = xmin + nx*xs;
    lbox[2] = ymin;
    lbox[3] = ymin + ny*ys;
    lbox[4] = zmin;
    lbox[5] = zmin + nz*zs;
    ((part_sys_soa *)This)->have_grid = true;
    std::cout << lbox[0] << " " << lbox[1] << " " << lbox[2] << " " << lbox[3] << " " << lbox[4] << " " << lbox[5] << std::endl;
    ((part_sys_soa *)This)->syncbboxes(lbox);
  }
  void part_sys__syncbboxes(Cpart_sys *This, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax) {
    bbox lbox;
    lbox.resize(6);
    lbox[0] = xmin;
    lbox[1] = xmax;
    lbox[2] = ymin;
    lbox[3] = ymax;
    lbox[4] = zmin;
    lbox[5] = zmax;
    std::cout << xmin << " " << ymin << " " << zmin << " " << xmax << " " << ymax << " " << zmax << std::endl;
    ((part_sys_soa *)This)->syncbboxes(lbox);
  }
  void part_sys__give_velptr(Cpart_sys *This, double * vx, double * vy, double * vz) {
    ((part_sys_soa *)This)->vxptr = vx;
    ((part_sys_soa *)This)->vyptr = vy;
    ((part_sys_soa *)This)->vzptr = vz;
  }
#ifdef __cplusplus
}
#endif
