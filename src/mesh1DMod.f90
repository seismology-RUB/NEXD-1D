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
module mesh1DMod

    use constantsMod
    use gllMod
    use jacobiMod
    use matrixMod !obsolet?
    use parameterMod
    use vandermonde1DMod
    use localOperatorsMod
    use boundaryConditionsMod
    implicit none

    contains

    subroutine setup_mesh(mesh, op, vdm, invV, errmsg)
        !in-/output
        type(error_message) :: errmsg
        type(meshVar), intent(inout) :: mesh
        !output
        type(operator) :: op
        real(kind=custom_real), dimension(:,:), intent(out) :: vdm, invV
        character(len=10) :: myname = "setup_mesh"

        call addTrace(errmsg, myname)

        !Generate simple Mesh
        call MeshGen1D(mesh)

        !Determine gll-points
        mesh%r = -getGll(mesh)

        !Initialize grid and calculate the metric
        call gridAndMetric1D(mesh, op, vdm, invV, errmsg)
    end subroutine

    subroutine MeshGen1D(mesh)
        !input
        type(meshVar), intent(inout) :: mesh
        !local
        integer :: i, j
        real(kind=custom_real) :: dx

        !Calculate element width from input parameter (parfile)
        dx = (mesh%xmax - mesh%xmin)/mesh%ncell

        !Calculate vertex coordinates
        allocate(mesh%vx(mesh%nv))
        do i = 1, mesh%nv
            mesh%vx(i) = mesh%xmin + (i - 1)*dx - mesh%nghost*dx
        enddo

        !Define elements based on vertices
        allocate(mesh%etov(mesh%K,mesh%nfaces))
        do j = 1, mesh%K
            mesh%etov(j,1) = j
            mesh%etov(j,2) = j + 1
        enddo

        !Calculate elements' width (basically this is the same as dx, however, because of numerical precistion this
        ! is slightly different. In addition it is generally possible to have an irregular grid. Therefore, the
        ! elements' width is calculated again to ensure that it is easy to change the meshgeneration / import
        ! external meshes.
        allocate(mesh%h(mesh%k))
        do i = 1, mesh%k
            mesh%h(i) = mesh%vx(mesh%etov(i,2))-mesh%vx(mesh%etov(i,1))
        enddo

        !The following two loops replace the module Connect1D.m as presented in Hesthaven & Warburton (2008), p.59f
        !Create mapping matrix for neighbouring elements
        allocate(mesh%etoe(mesh%K,mesh%nfaces))
        do i = 1, mesh%K
            mesh%etoe(i,1) = max(i-1,1)             ! neighbouring element to the left, set to 1-neg for element 1-ng
            mesh%etoe(i,2) = min(i+1,mesh%K)        ! neighbouring element to the right, set to ner+neg for last element
        enddo

        if (trim(mesh%rbc) == "periodic" .and. trim(mesh%lbc) == "periodic") call periodicBC(mesh)

        !Calculate element-to-faces map, i.e (ne x 2) matrix
        allocate(mesh%etof(mesh%K,mesh%nfaces))
        do i = 1, mesh%k
            mesh%etof(i,1) = 2
            mesh%etof(i,2) = 1
        enddo

        !Construct normals (Normals1D.m from Hesthaven & Warburton (2008)
        !Normals are the local outward pointing normals
        allocate(mesh%nx(mesh%Nfp*mesh%Nfaces, mesh%K))
        mesh%nx(1,:) = -1
        mesh%nx(2,:) = 1
    end subroutine MeshGen1D

    subroutine gridAndMetric1D(mesh, op, vdm, invV, errmsg)
        !This is equivalent to StartUp1D.m described in Hesthaven & Warburton (2008)
        ! Purpose: Building operators, grid, metric and connectivity for 1D solver.
        !input
        type(meshVar), intent(in) :: mesh
        type(error_message) :: errmsg
        !output
        type(operator) :: op
        real(kind=custom_real), dimension(:,:), intent(out) :: Vdm, invV
        !local variables
        !real(kind=custom_real), dimension(:,:), allocatable :: jacobian
        real(kind=custom_real), dimension(:,:), allocatable :: fx !x-coordinates at faces of elements
        character(len=15) :: myname = "gridAndMetric1D"

        call addTrace(errmsg, myname)

        !Allocate arrays
        allocate(op%jacobian(mesh%np,mesh%K))
        allocate(fx(mesh%Nfaces*mesh%Nfp,mesh%K))

        !Build reference element matrices
        call vdm1d(mesh, mesh%r, Vdm)
        call invVDM1D(mesh, Vdm, invV, errmsg)
        call dmatrix1d(mesh, op, Vdm, errmsg)

        !Create Surface integral Terms
        call lift1D(mesh, op, Vdm)

        !Calculate inverted Mass Matrix
        call invertedMassMatrix(mesh, op, vdm)

        !Build coordinates for all nodes and calculate masks for edge nodes
        call nodes1D(mesh, fx)

        !Calculate geometric factors
        call geometricFactors1D(mesh, op)

        !Build the surface normals and inverse metric at surface (Hesthaven & Warburton (2008), p.58)
        call fscale1D(mesh, op)

        deallocate(fx)
    end subroutine gridAndMetric1D

    subroutine nodes1D(mesh, Fx)
        !input
        type(meshVar) :: mesh
        !output
        real(kind=custom_real), dimension(:,:), intent(out) :: Fx
        !local variables
        integer :: fmask1, fmask2, i, j

        fmask1 = 0
        fmask2 = 0

        allocate(mesh%x(mesh%Np,mesh%K))
        !Build coordinates of all the nodes using the affine mapping
        ! x E D^k: x(r) = x_l^k + (1+r)/2 * h^k, h^k = x_f^k - x_l^k (eq. 3.1 in Hesthaven & Warburton (2008))
        do i = 1, mesh%K
            do j = 1, mesh%Np
                mesh%x(j,i) = mesh%vx(i) + 0.5*(mesh%r(j) +1)*(mesh%vx(i+1) - mesh%vx(i))
            enddo
        enddo
        !Compute masks for edge nodes (Hesthaven & Warburton (2008), p.57)
        do i = 1, mesh%Np
            if (abs(mesh%r(i) + 1) < EPS) then
                fmask1 = i
            endif
        enddo
        do i = 1, mesh%Np
            if (abs(mesh%r(i) - 1) < EPS) then
                fmask2 = i
            endif
        enddo

        allocate(mesh%ntof(mesh%nfaces))
        mesh%ntof = (/fmask1, fmask2/)

        do i = 1, mesh%K
            do  j = 1, mesh%Nfaces
                Fx(j,i) = mesh%x(mesh%ntof(j),i)
            enddo
        enddo
    end subroutine nodes1D
end module mesh1DMod
