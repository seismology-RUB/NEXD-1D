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
module localOperatorsMod
    use parameterMod
    use constantsMod
    use vandermonde1DMod
    use errorMessage
    implicit none

    type :: operator
        real(kind=custom_real), dimension(:,:), allocatable :: lift     !Array containing the first and last collum of the Mass Matrix
        real(kind=custom_real), dimension(:,:), allocatable :: Dr       !Differentiationsmatrix
        real(kind=custom_real), dimension(:,:), allocatable :: fscale   !Scaling factor ~2/elementwidth for the edges
        real(kind=custom_real), dimension(:,:), allocatable :: rx       !Scaling factor ~2/elementwidth for all grid points
        real(kind=custom_real), dimension(:,:), allocatable :: jacobian !Matrix containing the jacobi-polynomial
        real(kind=custom_real), dimension(:,:), allocatable :: invMass  !inverted Mass Matrix
    endtype operator

    contains

    subroutine lift1D(mesh, op, Vdm)
        !Purpose: Compute surface integral term in DG formulation
        !         i.e. pick first and last column of mass matrix

        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: Vdm
        !output
        type(operator) :: op
        !local variables
        real(kind=custom_real), dimension(mesh%Np, mesh%Nfaces*mesh%Nfp) :: emat, tmp

        emat = 0.0
        op%lift = 0.0

        !define Emat
        emat(1,1) = 1.0
        emat(mesh%Np, 2) = 1.0

        ! inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
        tmp = matmul(transpose(Vdm), emat)
        op%lift = matmul(Vdm, tmp)
    end subroutine lift1D

    subroutine dmatrix1d(mesh, op, Vdm, errmsg)
        !input
        type(error_message) :: errmsg
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: Vdm
        !output
        type(operator) :: op
        !local variables
        real(kind=custom_real), dimension(size(Vdm(:,1)),mesh%Np) :: Vr, Vinv
        character(len=9) :: myname = "dmatrix1d"

        call addTrace(errmsg, myname)

        call gradvdm1D(mesh, Vr)
        call invVDM1D(mesh, Vdm, Vinv, errmsg)

        op%Dr = matmul(Vr,Vinv)
    end subroutine dmatrix1d

    subroutine geometricFactors1D(mesh,op)
        ! Purpose: Compute the metrix elements for the local mappings
        !           of the 1D elements
        !input
        type(meshVar), intent(in) :: mesh
        !in/output
        type(operator) :: op
        !local
        real(kind=custom_real), dimension(mesh%Np,mesh%K) :: xr

        xr = matmul(op%Dr,mesh%x)
        op%jacobian = xr
        op%rx = 1/op%jacobian
    end subroutine geometricFactors1D

    subroutine fscale1D(mesh, op)
        !Build the surface normals and inverse metric at surface (Hesthaven & Warburton (2008), p.58)
        !input
        type(meshVar), intent(in) :: mesh
        !in/output
        type(operator) :: op
        !local
        integer :: i,j

        do i = 1, mesh%Nfaces
            do j = 1, mesh%K
                op%fscale(i,j) = 1.0/op%jacobian(mesh%ntof(i), j)
            enddo
        enddo
    end subroutine fscale1D

    subroutine invertedMassMatrix(mesh, op, vdm)
        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: Vdm
        !in/output
        type(operator) :: op
        allocate(op%invMass(mesh%np, mesh%np))

        op%invMass = matmul(vdm, transpose(vdm))
    end subroutine

    subroutine allocOperatorArrays(mesh, op)
        type(meshVar) :: mesh
        type(operator) :: op
        allocate(op%Dr(mesh%np, mesh%Np))
        allocate(op%lift(mesh%Np, mesh%Nfaces*mesh%Nfp))
        allocate(op%fscale(mesh%Nfaces, mesh%K))
        allocate(op%rx(mesh%Np,mesh%K))
    end subroutine

    subroutine deallocOperatorArrays(op)
        type(operator) :: op
        if (allocated(op%Dr)) deallocate(op%Dr)
        if (allocated(op%lift)) deallocate(op%lift)
        if (allocated(op%fscale)) deallocate(op%fscale)
        if (allocated(op%rx)) deallocate(op%rx)
    end subroutine


end module localOperatorsMod
