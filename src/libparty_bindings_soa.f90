!
! Copyright 2014.  Los Alamos National Security, LLC. This material was produced 
! under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
! Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
! for the U.S. Department of Energy. The U.S. Government has rights to use, 
! reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR 
! LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, 
! OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is 
! modified to produce derivative works, such modified software should be 
! clearly marked, so as not to confuse it with the version available from LANL.   
! 
! Licensed under the Apache License, Version 2.0 (the "License"); you may not 
! use this file except in compliance with the License. You may obtain a copy 
! of the License at http://www.apache.org/licenses/LICENSE-2.0
! 
! Unless required by applicable law or agreed to in writing, software distributed 
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
! CONDITIONS OF ANY KIND, either express or implied. See the License for the 
! specific language governing permissions and limitations under the License.
! 
!  Additionally, redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
!       National Laboratory, LANL, the U.S. Government, nor the names of its 
!       contributors may be used to endorse or promote products derived from 
!       this software without specific prior written permission.
!
! LibParty, Version 1.x, Copyright number C114122, LA-CC-14-086
!
! Authors:
!    Matthew Kinsey, kinsey@lanl.gov
!    Verinder Rana, vrana@lanl.gov
!    Bob Robey, XCP-2, LANL, brobey@lanl.gov
!    Jeremy Sauer, EES-16, LANL, jsauer@lanl.gov
!    Jon Reisner, XCP-4 LANL, reisner@lanl.gov
! 

! fort_test.f90
module libparty
  use, intrinsic :: ISO_C_Binding, only: C_int, C_double, C_ptr, C_NULL_ptr, &
                                         C_F_pointer, C_Loc
  implicit none
  private
  type part_sys_type
    private
    type(C_ptr) :: object = C_NULL_ptr
  end type part_sys_type
  type, BIND(C) :: part
    real(C_DOUBLE) :: pos(3)
    real(C_DOUBLE) :: vel(3)
    real(C_DOUBLE) :: m
    integer(C_INT) :: id
  end type part
  interface
    function C_part_sys__new (n) result(this) bind(C,name="part_sys__new")
      import
      type(C_ptr) :: this
      integer(C_int), value :: n
    end function C_part_sys__new
    subroutine C_part_sys__delete (this) bind(C,name="part_sys__delete")
      import
      type(C_ptr), value :: this
    end subroutine C_part_sys__delete
    function C_part_sys__npart (this) result(n_part) bind(C,name="part_sys__npart")
      import
      integer(C_int) :: n_part
      type(C_ptr), value :: this
    end function C_part_sys__npart
    subroutine C_part_sys__print_sys (this) bind(C,name="part_sys__print_sys")
      import
      type(C_ptr), value :: this
    end subroutine C_part_sys__print_sys
    subroutine C_part_sys__output_h5part (this) bind(C,name="part_sys__output_h5part")
      import
      type(C_ptr), value :: this
    end subroutine C_part_sys__output_h5part
    subroutine C_part_sys__step (this, dt) bind(C,name="part_sys__step")
      import
      REAL(C_double), value :: dt
      type(C_ptr), value :: this
    end subroutine C_part_sys__step
!    subroutine C_part_sys__returnpos(this, x, y, z) bind(C,name="part_sys__returnpos")
!      import
!      real(C_double), dimension(*), intent(out) :: x, y, z
!      type(C_ptr), value :: this
!    end subroutine C_part_sys__returnpos
!    subroutine C_part_sys__push_pos(this, x, y, z) bind(C,name="part_sys__push_pos")
!      import
!      real(C_double), dimension(*), intent(in):: x, y, z
!      type(C_ptr), value:: this
!    end subroutine C_part_sys__push_pos
!    subroutine C_part_sys__push_vel(this, vx, vy, vz) bind(C,name="part_sys__push_vel")
!      import
!      real(C_double), dimension(*), intent(in) :: vx, vy, vz
!      type(C_ptr), value:: this
!    end subroutine C_part_sys__push_vel
    function C_part_sys__data_ptr (this) result(dptr) bind(C,name="part_sys__data_ptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: dptr
    end function C_part_sys__data_ptr
    function C_part_sys__xptr (this) result(xptr) bind(C,name="part_sys__xptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: xptr
    end function C_part_sys__xptr
    function C_part_sys__yptr (this) result(yptr) bind(C,name="part_sys__yptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: yptr
    end function C_part_sys__yptr
    function C_part_sys__zptr (this) result(zptr) bind(C,name="part_sys__zptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: zptr
    end function C_part_sys__zptr
    function C_part_sys__vxptr (this) result(vxptr) bind(C,name="part_sys__vxptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: vxptr
    end function C_part_sys__vxptr
    function C_part_sys__vyptr (this) result(vyptr) bind(C,name="part_sys__vyptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: vyptr
    end function C_part_sys__vyptr
    function C_part_sys__vzptr (this) result(vzptr) bind(C,name="part_sys__vzptr")
      import
      type(C_ptr), value :: this
      type(C_ptr) :: vzptr
    end function C_part_sys__vzptr
  end interface
  interface psys_new
    module procedure part_sys__new
  end interface psys_new
  interface psys_delete
    module procedure part_sys__delete
  end interface psys_delete
  interface psys_npart
    module procedure part_sys__npart
  end interface psys_npart
  interface psys_step
    module procedure part_sys__step
  end interface psys_step
  interface psys_print_sys
    module procedure part_sys__print_sys
  end interface psys_print_sys
  interface psys_output_h5part
    module procedure part_sys__output_h5part
  end interface psys_output_h5part
  interface psys_xptr
    module procedure part_sys__xptr
  end interface psys_xptr
  interface psys_yptr
    module procedure part_sys__yptr
  end interface psys_yptr
  interface psys_zptr
    module procedure part_sys__zptr
  end interface psys_zptr
  interface psys_vxptr
    module procedure part_sys__vxptr
  end interface psys_vxptr
  interface psys_vyptr
    module procedure part_sys__vyptr
  end interface psys_vyptr
  interface psys_vzptr
    module procedure part_sys__vzptr
  end interface psys_vzptr
  interface psys_data_ptr
    module procedure part_sys__data_ptr
  end interface psys_data_ptr

  public :: psys_new, psys_delete, psys_npart, psys_print_sys, &
            psys_step, psys_output_h5part, psys_xptr, &
            psys_yptr, psys_zptr, psys_vxptr, psys_vyptr, psys_vzptr, &
            part_sys_type, part
contains
! Fortran wrapper routines to interface C wrappers
  subroutine part_sys__new(this,n)
    type(part_sys_type), intent(out) :: this
    integer :: n
    this%object = C_part_sys__new(int(n,C_int))
  end subroutine part_sys__new
  subroutine part_sys__delete(this)
    type(part_sys_type), intent(inout) :: this
    call C_part_sys__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine part_sys__delete
  function part_sys__npart(this) result(n_part)
    type(part_sys_type), intent(in) :: this
    integer :: n_part
    n_part = C_part_sys__npart(this%object)
  end function part_sys__npart
  subroutine part_sys__step(this,dt)
    type(part_sys_type), intent(in) :: this
    real*8 :: dt
    call C_part_sys__step(this%object, dt)
  end subroutine part_sys__step
  subroutine part_sys__print_sys(this)
    type(part_sys_type), intent(in) :: this
    call C_part_sys__print_sys(this%object)
  end subroutine part_sys__print_sys
  subroutine part_sys__output_h5part(this)
    type(part_sys_type), intent(in) :: this
    call C_part_sys__output_h5part(this%object)
  end subroutine part_sys__output_h5part

!  subroutine part_sys__returnpos(this, x, y, z)
!    type(part_sys_type), intent(in) :: this
!    real*8, intent(out), dimension(*):: x, y, z
!    call C_part_sys__returnpos(this%object, x, y, z)
!  end subroutine part_sys__returnpos
!  subroutine part_sys__push_pos(this, x, y, z)
!    type(part_sys_type), intent(in) :: this
!    real*8, dimension(*), intent(in) :: x, y, z
!    call C_part_sys__push_pos(this%object, x, y, z)
!  end subroutine part_sys__push_pos
!  subroutine part_sys__push_vel(this, vx, vy, vz)
!    type(part_sys_type), intent(in) :: this
!    real*8, dimension(*), intent(in) :: vx, vy, vz
!    call C_part_sys__push_vel(this%object, vx, vy, vz)
!  end subroutine part_sys__push_vel
  subroutine part_sys__data_ptr(this,dptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    type(part), pointer :: dptr(:)
!    np = C_part_sys__npart(this%object)
    firstp = C_part_sys__data_ptr(this%object)
    call c_f_pointer(firstp, dptr,[np])
  end subroutine part_sys__data_ptr
  subroutine part_sys__xptr(this,xptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    real(C_DOUBLE), pointer :: xptr(:)
    firstp = C_part_sys__xptr(this%object)
    call c_f_pointer(firstp, xptr,[np])
  end subroutine part_sys__xptr
  subroutine part_sys__yptr(this,yptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    real(C_DOUBLE), pointer :: yptr(:)
    firstp = C_part_sys__yptr(this%object)
    call c_f_pointer(firstp, yptr,[np])
  end subroutine part_sys__yptr
  subroutine part_sys__zptr(this,zptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    real(C_DOUBLE), pointer :: zptr(:)
    firstp = C_part_sys__zptr(this%object)
    call c_f_pointer(firstp, zptr,[np])
  end subroutine part_sys__zptr
  subroutine part_sys__vxptr(this,vxptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    real(C_DOUBLE), pointer :: vxptr(:)
    firstp = C_part_sys__vxptr(this%object)
    call c_f_pointer(firstp, vxptr,[np])
  end subroutine part_sys__vxptr
  subroutine part_sys__vyptr(this,vyptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    real(C_DOUBLE), pointer :: vyptr(:)
    firstp = C_part_sys__vyptr(this%object)
    call c_f_pointer(firstp, vyptr,[np])
  end subroutine part_sys__vyptr
  subroutine part_sys__vzptr(this,vzptr,np)
    type(part_sys_type), intent(in) :: this
    type(C_ptr) :: firstp
    integer ,intent(in) :: np
    real(C_DOUBLE), pointer :: vzptr(:)
    firstp = C_part_sys__vzptr(this%object)
    call c_f_pointer(firstp, vzptr,[np])
  end subroutine part_sys__vzptr

end module libparty
