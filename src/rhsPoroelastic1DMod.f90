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
module rhsPoroelastic1DMod

    use parameterMod
    use constantsMod
    use matrixMod
    use riemannflux1Dmod
    use rhsSlipInterface1DMod
    use genericMaterial
    use slipInterfaceMod
    use boundaryConditionsMod

    implicit none

    contains

    subroutine rhsPoroelastic(par, mesh, genmat, op, solution, rhs)
        !Purpose: Evaluate RHS flux in 1D poroelastic wave propagation

        !input
        type (parameterVar), intent(in) :: par
        type (meshVar), intent(in) :: mesh
        type (generic_material), intent(in) :: genmat
        type (operator), intent(in) :: op
        type (solutionvector), intent(in) :: solution
        !output
        type (solutionvector), intent(out) :: rhs
        !local variables
        integer :: i, j, k, m, nif, nmf
        !integer :: l
        real(kind=custom_real), dimension(:,:), allocatable :: Ain,Ein!,Aout
        real(kind=custom_real), dimension(:,:), allocatable :: AP,AM
        type (solutionvector) :: flux
        type (solutionvector) :: tmp
        real(kind=custom_real), dimension(:), allocatable :: fluxR
        real(kind=custom_real), dimension(:), allocatable :: delta
        real(kind=custom_real), dimension(:), allocatable :: temp, temp2

        allocate(Ain(2+2*par%fluidn,2+2*par%fluidn))
        allocate(Ein(2+2*par%fluidn,2+2*par%fluidn))
        allocate(fluxR(2+2*par%fluidn))
        allocate(delta(2+2*par%fluidn))
        allocate(temp(2+2*par%fluidn))
        allocate(temp2(2+2*par%fluidn))

        allocate(flux%v(mesh%Nfp*mesh%Nfaces,mesh%K))
        allocate(flux%sigma(mesh%Nfp*mesh%Nfaces,mesh%K))
        allocate(tmp%v(mesh%Np,mesh%K))
        allocate(tmp%sigma(mesh%Np,mesh%K))
        allocate(flux%v1(mesh%Nfp*mesh%Nfaces,mesh%K))
        allocate(flux%p1(mesh%Nfp*mesh%Nfaces,mesh%K))
        allocate(tmp%v1(mesh%Np,mesh%K))
        allocate(tmp%p1(mesh%Np,mesh%K))
        if (par%fluidn == 2) then
            allocate(flux%v2(mesh%Nfp*mesh%Nfaces,mesh%K))
            allocate(flux%p2(mesh%Nfp*mesh%Nfaces,mesh%K))
            allocate(tmp%v2(mesh%Np,mesh%K))
            allocate(tmp%p2(mesh%Np,mesh%K))
        endif

        !calculate fluxes
        do i = 1, mesh%K                           !loop over elements
            do j = 1, mesh%Nfaces                  !loop over surfaces
                m = mesh%etoe(i,j)                 !neighboring element of element i at surface j
                nmf = mesh%ntof(mesh%etof(i,j))    !node at face of neighboring element
                nif = mesh%ntof(j)                 !node at own face

                !get difference of solutionvectors
                delta(1) = solution%sigma(nmf, m) - solution%sigma(nif, i)
                delta(2) = solution%v(nmf, m)     - solution%v(nif, i)
                delta(3) = solution%p1(nmf, m)    - solution%p1(nif, i)
                delta(4) = solution%v1(nmf, m)    - solution%v1(nif, i)
                if (par%fluidn == 2) then
                    delta(5) = solution%p2(nmf, m)    - solution%p2(nif, i)
                    delta(6) = solution%v2(nmf, m)    - solution%v2(nif, i)
                endif

                if (i == 1+mesh%nghost .or. i == mesh%K-mesh%nghost) then
                    if (trim(mesh%rbc) == "reflecting" .or. trim(mesh%lbc) == "reflecting") then
                        call reflectingBC(mesh, i, j, solution, delta)
                    endif
                endif

                select case (par%flux_type)
                    case (1:2) !Godunov or Rusanov
                        Ain = genmat.A.(mesh%ntom(nif,i))  !material properties at own face
                        allocate(AM(2+2*par%fluidn,2+2*par%fluidn))
                        select case (par%flux_type)
                            case (1) !Godunov
                                AM = genmat.AM.(mesh%ntom(nif,i))
                            case (2) !Rusanov
                                allocate(AP(2+2*par%fluidn,2+2*par%fluidn))
                                AP = genmat.AP.(mesh%ntom(nif,i))
                                AM = 0.
                                do k=1,2+2*par%fluidn
                                    AM(k,k) = -AP(1,1)
                                enddo
                                AM = (AM + Ain)/2.0
                                if (allocated(AP)) deallocate(AP)
                        end select

                        call GodunovFlux(delta, AM, mesh%nx(j,i), fluxR)
                        flux%sigma(j,i) = op%fscale(j,i) * fluxR(1)
                        flux%v(j,i)     = op%fscale(j,i) * fluxR(2)
                        flux%p1(j,i)    = op%fscale(j,i) * fluxR(3)
                        flux%v1(j,i)    = op%fscale(j,i) * fluxR(4)
                        if (par%fluidn == 2) then
                            flux%p2(j,i)    = op%fscale(j,i) * fluxR(5)
                            flux%v2(j,i)    = op%fscale(j,i) * fluxR(6)
                        endif
                        if (allocated(AM)) deallocate(AM)
                end select
            enddo
        enddo

        tmp%sigma = matmul(op%Dr, solution%sigma)
        tmp%v     = matmul(op%Dr, solution%v)
        tmp%p1    = matmul(op%Dr, solution%p1)
        tmp%v1    = matmul(op%Dr, solution%v1)
        if (par%fluidn == 2) then
            tmp%p2    = matmul(op%Dr, solution%p2)
            tmp%v2    = matmul(op%Dr, solution%v2)
        endif
        do i = 1, mesh%K                           !loop over elements
            do j = 1, mesh%Np                      !loop over nodes of element i
                Ain = genmat.A.(mesh%ntom(j,i))
                Ein = genmat.E.(mesh%ntom(j,i))
                select case (par%fluidn)
                    case (1)
                        temp  = (/ tmp%sigma(j,i), tmp%v(j,i), tmp%p1(j,i), tmp%v1(j,i)/)                            !Dr*Q
                        temp2 = (/ solution%sigma(j,i), solution%v(j,i), solution%p1(j,i), solution%v1(j,i)/)        !Q
                    case (2)
                        temp  = (/ tmp%sigma(j,i), tmp%v(j,i), tmp%p1(j,i), tmp%v1(j,i), tmp%p2(j,i), tmp%v2(j,i)/)  !Dr*Q
                        temp2 = (/ solution%sigma(j,i), solution%v(j,i), &
                         solution%p1(j,i), solution%v1(j,i), solution%p2(j,i), solution%v2(j,i) /)                   !Q
                end select
                ! r_x * A *  [Dr*Q](x_j) + E * Q(x_j)       (i.e. right-hand-side without flux.)
                temp = op%rx(j,i) * matmul(Ain,temp) - matmul(Ein,temp2)
                tmp%sigma(j,i) = temp(1)
                tmp%v(j,i)     = temp(2)
                tmp%p1(j,i)    = temp(3)
                tmp%v1(j,i)    = temp(4)
                if (par%fluidn == 2) then
                    tmp%p2(j,i)    = temp(5)
                    tmp%v2(j,i)    = temp(6)
                endif
            enddo
        enddo

        !compute right hand sides of the PDE
        rhs%v     = -tmp%v     - matmul(op%lift, flux%v)
        rhs%sigma = -tmp%sigma - matmul(op%lift, flux%sigma)
        rhs%v1    = -tmp%v1    - matmul(op%lift, flux%v1)
        rhs%p1    = -tmp%p1    - matmul(op%lift, flux%p1)
        if (par%fluidn == 2) then
            rhs%v2    = -tmp%v2    - matmul(op%lift, flux%v2)
            rhs%p2    = -tmp%p2    - matmul(op%lift, flux%p2)
        endif

        if (allocated(Ain)) deallocate(Ain)
        if (allocated(Ein)) deallocate(Ein)
        if (allocated(fluxR)) deallocate(fluxR)
        if (allocated(delta)) deallocate(delta)
        if (allocated(temp)) deallocate(temp)
        if (allocated(temp2)) deallocate(temp2)
        if (allocated(flux%v)) deallocate(flux%v)
        if (allocated(flux%sigma)) deallocate(flux%sigma)
        if (allocated(flux%v1)) deallocate(flux%v1)
        if (allocated(flux%p1)) deallocate(flux%p1)
        if (allocated(flux%v2)) deallocate(flux%v2)
        if (allocated(flux%p2)) deallocate(flux%p2)
        if (allocated(tmp%v)) deallocate(tmp%v)
        if (allocated(tmp%sigma)) deallocate(tmp%sigma)
        if (allocated(tmp%v1)) deallocate(tmp%v1)
        if (allocated(tmp%p1)) deallocate(tmp%p1)
        if (allocated(tmp%v2)) deallocate(tmp%v2)
        if (allocated(tmp%p2)) deallocate(tmp%p2)
    end subroutine rhsPoroelastic
end module
