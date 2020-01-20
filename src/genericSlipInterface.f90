!-----------------------------------------------------------------------
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
!
!   This file is part of NEXD 1D.
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but
!   WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with NEXD 1D. If not, see <http://www.gnu.org/licenses/>.
!---------------------------------
!-------------------------------------------------------------
!> \brief Generic interface 1d module
!
module genericSlipInterface
    use slipInterfaceMod
    implicit none
    interface dealloc; module procedure deallocGenericSlipInterface; end interface
    interface operator (.eta.); module procedure getComplianceGenericSlipInterface; end interface
    interface operator (.xi.); module procedure getFluidityGenericSlipInterface; end interface
   type generic_slip_interface
       type (transverse_slip_interface), pointer :: tsi => null()
    end type generic_slip_interface
!
contains
!-------------------------------------------------------------
!> \brief associate to transverse slip interface
!
    subroutine associateTransverseToGenericSlipInterface(this,tsi)
    type (generic_slip_interface) :: this
    type (transverse_slip_interface), target :: tsi
    call dealloc(this)
    this%tsi => tsi
    end subroutine associateTransverseToGenericSlipInterface
!-----------------------------------------------------------
!> brief deallocate generic interface object
!
    subroutine deallocGenericSlipInterface(this)
    type (generic_slip_interface) :: this
    if (associated(this%tsi)) nullify(this%tsi)
    end subroutine deallocGenericSlipInterface
!-----------------------------------------------------------------
!> \brief return selected compliance value
!
    function getComplianceGenericSlipInterface(this,j) result(res)
    type (generic_slip_interface), intent(in) :: this
    integer, intent(in) :: j
    real(kind=custom_real), dimension(size(this%tsi%eta(1,:))) :: res
    if (associated(this%tsi)) then
       res = getComplianceTransverseSlipInterface(this%tsi,j)
    endif
    end function getComplianceGenericSlipInterface

!-----------------------------------------------------------------
!> \brief return selected fluidity value
!
    function getFluidityGenericSlipInterface(this,j) result(res)
    type (generic_slip_interface), intent(in) :: this
    integer, intent(in) :: j
    real(kind=custom_real), dimension(size(this%tsi%xi(1,:))) :: res
    if (associated(this%tsi)) then
       res = getFluidityTransverseSlipInterface(this%tsi,j)
    endif
    end function getFluidityGenericSlipInterface
end module genericSlipInterface
