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
module vandermonde1DMod
    ! Module to create the Vandermonde Matrix
    use constantsMod
    use parameterMod
    use jacobiMod
    use errorMessage
    implicit none

    contains

    subroutine vdm1d(mesh, r, Vdm)
        ! Purpose: Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:), intent(in) :: r
        !output
        real(kind=custom_real), dimension(size(r),mesh%np), intent(out) :: Vdm
        !local variables
        integer :: j
        real(kind=custom_real) :: alpha, beta

        alpha = 0.0
        beta = 0.0

        do j = 1, (mesh%N+1)
            call jacobiP(Vdm(:,j), r(:), alpha, beta, j-1)
        enddo
    end subroutine vdm1d

    subroutine gradvdm1D(mesh, DVr)
        ! Purpose: Initialize the gradient of the modal basis (i) at (r)
        !          at order N
        !input
        type(meshVar), intent(in) :: mesh
        !output
        real(kind=custom_real), dimension(size(mesh%r),mesh%N+1), intent(out) :: DVr
        !local variables
        integer :: j
        real(kind=custom_real) :: alpha, beta

        alpha = 0.0
        beta = 0.0

        ! Initialize matrix
        do j = 0, mesh%N
            call gradJacobiP(DVr(:,j+1), mesh%r(:), alpha, beta, j)
        enddo
    end subroutine gradvdm1D

    subroutine invVDM1D(mesh, Vdm, invV, errmsg)
        ! Purpose: Calulates the inverse of V, V^-1' if trans==0 and
        !          calculates the inverse transposed if trans == 1
        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: Vdm
        !output
        real(kind=custom_real), dimension(:,:), intent(out) :: invV
        !local variables
        real(kind=custom_real), dimension(2*mesh%Np) :: work !work array for LAPACK
        integer, dimension(mesh%Np) :: ipiv !pivot indices
        integer :: info, iwork
        character(len=8) :: myname = 'invVDM1D'
        type(error_message) :: errmsg

        invV = Vdm
        iwork = 2*mesh%Np

        ! Compute the LU factorization of W
        call addTrace(errmsg,myname)
        if (custom_real == size_real) then
            call sgetrf(mesh%Np, mesh%Np, invV, mesh%Np, ipiv, info)
        elseif (custom_real == size_double) then
            call dgetrf(mesh%Np, mesh%Np, invV, mesh%Np, ipiv, info)
        else
            call add(errmsg,2,'custom_real is neither size_real nor size_double',myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif
        if (info /= 0) then
            call add(errmsg,2,'Matrix is numerically singular!',myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        end if

        ! Compute the inverse of W
        call addTrace(errmsg,myname)
        if (custom_real == size_real) then
            call sgetri(mesh%Np, invV, mesh%Np, ipiv, work, iwork, info)
        elseif (custom_real == size_double) then
            call dgetri(mesh%Np, invV, mesh%Np, ipiv, work, iwork, info)
        else
            call add(errmsg,2,'custom_real is neither size_real nor size_double',myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif
        if (info /= 0) then
            call add(errmsg,2,'Matrix inversion failed!',myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        end if
    end subroutine invVDM1D
end module vandermonde1DMod
