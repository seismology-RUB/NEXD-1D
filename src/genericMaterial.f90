!-----------------------------------------------------------------------
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
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
module genericMaterial
    !Module for handling different material descriptions
    use materials1DMod
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocGenericMaterial; end interface
    interface operator (.rho.); module procedure getDensityGenericMaterial; end interface
    interface operator (.vs.); module procedure getSVelocityGenericMaterial; end interface
    interface operator (.imp.); module procedure getImpedanceGenericMaterial; end interface
    interface operator (.A.); module procedure getMatrixAGenericMaterial; end interface
    interface operator (.AP.); module procedure getMatrixAPGenericMaterial; end interface
    interface operator (.AM.); module procedure getMatrixAMGenericMaterial; end interface
    interface operator (.E.); module procedure getMatrixEGenericMaterial; end interface
    type generic_material
        type (rovszimp_material), pointer :: rvzmat => null()
        type (porous_material), pointer :: poromat => null()
    end type generic_material

    contains

    subroutine associateRoVsZimpToGenericMaterial(this,rvzmat)
        !Function to associate to specific ro-vs-zimp-material
        type (generic_material) :: this
        type (rovszimp_material), target :: rvzmat

        call dealloc(this)
        this%rvzmat => rvzmat
    end subroutine associateRoVsZimpToGenericMaterial

    subroutine associatePorousToGenericMaterial(this,poromat)
        !Function to associate to specific ro-vs-zimp-material
        type (generic_material) :: this
        type (porous_material), target :: poromat

        call dealloc(this)
        this%poromat => poromat
    end subroutine associatePorousToGenericMaterial

    subroutine deallocGenericMaterial(this)
        !Function to deallocate object
        type (generic_material) :: this

        if (associated(this%rvzmat)) nullify(this%rvzmat)
        if (associated(this%poromat)) nullify(this%poromat)
    end subroutine deallocGenericMaterial

    function getDensityGenericMaterial(this,j) result(res)
        !Function to return selected density value
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real(kind=custom_real) :: res

        res = 0
        if (associated(this%rvzmat)) then
            res = getDensityRoVsZimpMaterial(this%rvzmat,j)
        endif
    end function getDensityGenericMaterial

    function getSVelocityGenericMaterial(this,j) result(res)
        !Function to return selected S-velocity value
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real(kind=custom_real) :: res

        res = 0
        if (associated(this%rvzmat)) then
            res = getSVelocityRoVsZimpMaterial(this%rvzmat,j)
        endif
    end function getSVelocityGenericMaterial

    function getImpedanceGenericMaterial(this,j) result(res)
        !Function to return selected impedance value
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real(kind=custom_real) :: res

        res = 0
        if (associated(this%rvzmat)) then
            res = getImpedanceRoVsZimpMaterial(this%rvzmat,j)
        endif
    end function getImpedanceGenericMaterial

    function getMatrixAGenericMaterial(this,j) result(res)
        !Function to return selected Matrix A
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real(kind=custom_real), dimension(:,:), allocatable :: res

        if (associated(this%poromat)) then
            if (associated(this%poromat%S2)) then
                allocate(res(6,6))
                res = getMatrixAPorousMaterial(this%poromat,j)
            else
                allocate(res(4,4))
                res = getMatrixAPorousMaterial(this%poromat,j)
            endif
        else
            allocate(res(2,2))
            res = getMatrixAElasticMaterial(this%rvzmat,j)
        endif
    end function getMatrixAGenericMaterial

    function getMatrixAPGenericMaterial(this,j) result(res)
        !Function to return selected Matrix AP
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res

        if (associated(this%poromat)) then
            if (associated(this%poromat%S2)) then
                allocate(res(6,6))
                res = getMatrixAPPorousMaterial(this%poromat,j)
            else
                allocate(res(4,4))
                res = getMatrixAPPorousMaterial(this%poromat,j)
            endif
        else
            allocate(res(2,2))
            res = getMatrixAPElasticMaterial(this%rvzmat,j)
        endif
    end function getMatrixAPGenericMaterial

    function getMatrixAMGenericMaterial(this,j) result(res)
        !Function to return selected Matrix AM
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res

        if (associated(this%poromat)) then
            if (associated(this%poromat%S2)) then
                allocate(res(6,6))
                res = getMatrixAMPorousMaterial(this%poromat,j)
            else
                allocate(res(4,4))
                res = getMatrixAMPorousMaterial(this%poromat,j)
            endif
        else
            allocate(res(2,2))
            res = getMatrixAMElasticMaterial(this%rvzmat,j)
        endif
    end function getMatrixAMGenericMaterial

    function getMatrixEGenericMaterial(this,j) result(res)
        !Function to return selected Matrix E
        type (generic_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res

        if (associated(this%poromat)) then
            if (associated(this%poromat%S2)) then
                allocate(res(6,6))
                res = getMatrixEPorousMaterial(this%poromat,j)
            else
                allocate(res(4,4))
                res = getMatrixEPorousMaterial(this%poromat,j)
            endif
        endif
    end function getMatrixEGenericMaterial

    function maximumVelocityGenericMaterial(this) result(res)
        !Function to return maximum velocity
        type (generic_material), intent(in) :: this
        real(kind=custom_real) :: res

        res = 0
        if (associated(this%rvzmat)) then
            res = maximumVelocityRoVsZimpMaterial(this%rvzmat)
        endif
    end function maximumVelocityGenericMaterial
end module genericMaterial
