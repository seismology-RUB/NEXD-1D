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
module boundaryConditionsMod

    use parameterMod
    use constantsMod
    implicit none

    contains

    subroutine periodicBC(mesh)
        !in/out
        type(meshVar) :: mesh

        mesh%etoe(1+mesh%nghost,1) = mesh%K-mesh%nghost
        mesh%etoe(mesh%K-mesh%nghost,2) = 1+mesh%nghost
    end subroutine

    subroutine absorbingBC(mesh, solution)
        !input
        type(meshVar) :: mesh
        !in/out
        type(solutionvector), intent(inout) :: solution

        if (trim(mesh%lbc) == "absorbing") then
            solution%v(:,mesh%nghost-1)          = solution%v(2,mesh%nghost)
            solution%sigma(:,mesh%nghost-1)      = solution%sigma(2,mesh%nghost)
            if (allocated(solution%v1)) then
                solution%v1(:,mesh%nghost-1)     = solution%v1(2,mesh%nghost)
                solution%p1(:,mesh%nghost-1)     = solution%p1(2,mesh%nghost)
                if (allocated(solution%v2)) then
                    solution%v2(:,mesh%nghost-1) = solution%v2(2,mesh%nghost)
                    solution%p2(:,mesh%nghost-1) = solution%p2(2,mesh%nghost)
                endif
            endif
        endif
        if (trim(mesh%rbc) == "absorbing") then
            solution%v(:,mesh%K)          = solution%v(mesh%n,mesh%K-(mesh%nghost-1))
            solution%sigma(:,mesh%K)      = solution%sigma(mesh%n,mesh%K-(mesh%nghost-1))
            if (allocated(solution%v1)) then
                solution%v1(:,mesh%K)     = solution%v1(mesh%n,mesh%K-(mesh%nghost-1))
                solution%p1(:,mesh%K)     = solution%p1(mesh%n,mesh%K-(mesh%nghost-1))
                if (allocated(solution%v2)) then
                    solution%v2(:,mesh%K) = solution%v2(mesh%n,mesh%K-(mesh%nghost-1))
                    solution%p2(:,mesh%K) = solution%p2(mesh%n,mesh%K-(mesh%nghost-1))
                endif
            endif
        endif
    end subroutine

    subroutine reflectingBC(mesh, i, j, solution, delta)
        type(meshVar) :: mesh
        type(solutionvector), intent(in) :: solution
        integer, intent(in) :: i,j
        !output
        real(kind=custom_real), dimension(:), intent(out) :: delta

        if ((trim(mesh%lbc) == "reflecting" .and. (i == 1+mesh%nghost      .and. j == 1)) .or. &
            (trim(mesh%rbc) == "reflecting" .and. (i == mesh%k-mesh%nghost .and. j == 2))) then
            delta(1) = -2.*solution%sigma(mesh%ntof(j),i)
            delta(2) = 0.
            if (allocated(solution%v1)) then
                delta(3) = -2.*solution%p1(mesh%ntof(j),i)
                delta(4) = 0.
                if (allocated(solution%v2)) then
                    delta(5) = -2.*solution%p2(mesh%ntof(j),i)
                    delta(6) = 0.
                endif
            endif
        endif
    end subroutine
end module
