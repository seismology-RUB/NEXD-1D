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
module rhsElastic1DMod

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

    subroutine rhsElastic(par, mesh, genmat, op, solution, lsi, rhs)
        !Purpose: Evaluate RHS flux in 1D elastic wave propagation

        !input
        type (parameterVar), intent(in) :: par
        type (meshVar), intent(in) :: mesh
        type (generic_material), intent(in) :: genmat
        type (operator), intent(in) :: op
        type (solutionvector), intent(in) :: solution
        !in-/output
        type (lsiVar), dimension(:) :: lsi
        !output
        type (solutionvector) :: rhs
        !local variables
        integer :: i, j, k, m, nif, nmf
        real(kind=custom_real) :: zin,zout,beta
        real(kind=custom_real), dimension(mesh%Nfp*mesh%Nfaces,mesh%K) :: fluxV
        real(kind=custom_real), dimension(mesh%Nfp*mesh%Nfaces,mesh%K) :: fluxSigma
        real(kind=custom_real), dimension(mesh%Np,mesh%K) :: tmpV
        real(kind=custom_real), dimension(mesh%Np,mesh%K) :: tmpSigma
        real(kind=custom_real), dimension(2) :: fluxR, fluxL
        real(kind=custom_real) :: dV
        real(kind=custom_real) :: dSigma
        real(kind=custom_real), dimension(2) :: delta

        real(kind=custom_real), dimension(2,2) :: Ain
        real(kind=custom_real), dimension(:,:), allocatable :: AP,AM

        !calculate fluxes
        do i = 1, mesh%K                           !loop over elements
            do j = 1, mesh%Nfaces                  !loop over surfaces
                m = mesh%etoe(i,j)                 !neighboring element of element i at surface j
                nmf = mesh%ntof(mesh%etof(i,j))    !node at face of neighboring element
                nif = mesh%ntof(j)                 !node at own face

                !get difference of solutionvectors
                delta(1) = solution%sigma(nmf, m) - solution%sigma(nif,i)
                delta(2) = solution%v(nmf, m)     - solution%v(nif,i)

                if (i == 1+mesh%nghost .or. i == mesh%K-mesh%nghost) then
                    if (trim(mesh%rbc) == "reflecting" .or. trim(mesh%lbc) == "reflecting") then
                        !call reflectingBC(mesh, i, j, V, Sigma, delta(2), delta(1))
                        call reflectingBC(mesh, i, j, solution, delta)
                    endif
                endif

                select case (par%flux_type)
                    case (0) !Riemann von Thomas
                        beta = genmat.vs.(mesh%ntom(nif,i))
                        zin = genmat.imp.(mesh%ntom(nif,i))  !material properties at own face
                        zout = genmat.imp.(mesh%ntom(nmf,m)) !material properties at neighboring face

                        dv     = delta(2)*mesh%nx(j,i)
                        dsigma = delta(1)*mesh%nx(j,i)

                        if (j == 1) then
                            !fluss nach rechts von der elementgrenze weg
                            call rightRiemannflux(dV, dSigma, zin, zout , beta, fluxR)
                            fluxV(j,i)     = mesh%nx(j,i) * op%fscale(j,i) * fluxR(1)
                            fluxSigma(j,i) = mesh%nx(j,i) * op%fscale(j,i) * fluxR(2)
                        else if (j == 2) then
                            !fluss nach links vom der elementgrenze weg
                            call leftRiemannflux(dV, dSigma, zin, zout, beta, fluxL)
                            fluxV(j,i)     = mesh%nx(j,i) * op%fscale(j,i) * fluxL(1)
                            fluxSigma(j,i) = mesh%nx(j,i) * op%fscale(j,i) * fluxL(2)
                        endif
                    case (1:2) !Godunov or Rusanov
                        Ain = genmat.A.(mesh%ntom(nif,i))    !material properties at own face
                        allocate(AM(2,2))
                        select case (par%flux_type)
                            case (1) !Godunov
                                AM = genmat.AM.(mesh%ntom(nif,i))
                            case (2) !Rusanov
                                allocate(AP(2,2))
                                AP = genmat.AP.(mesh%ntom(nif,i))
                                AM = 0.
                                do k=1,2
                                    AM(k,k) = -AP(1,1)
                                enddo
                                AM = (AM + Ain)/2.0
                                if (allocated(AP)) deallocate(AP)
                        end select

                        call GodunovFlux(delta, AM, mesh%nx(j,i), fluxR)
                        fluxSigma(j,i) = op%fscale(j,i) * fluxR(1)
                        fluxV(j,i)     = op%fscale(j,i) * fluxR(2)
                        if (allocated(AM)) deallocate(AM)
                end select
            enddo
        enddo

        !slip interface influence
        if (par%lsin /= 0) then
            do i = 1, par%lsin
                call rhsAddSlipInterfaceInfluence(mesh, genmat, op, lsi(i), fluxV, fluxSigma)
            enddo
        endif

        tmpSigma = matmul(op%Dr, solution%sigma)
        tmpV     = matmul(op%Dr, solution%v)
        do i = 1, mesh%K                           !loop over elements
            do j = 1, mesh%Np                      !loop over nodes of element i
                beta = genmat.vs.(mesh%ntom(j,i))
                zin = genmat.imp.(mesh%ntom(j,i))
                tmpV(j,i)     = op%rx(j,i) * beta*zin * tmpV(j,i)
                tmpSigma(j,i) = op%rx(j,i) * beta/zin * tmpSigma(j,i)
            enddo
        enddo

        !compute right hand sides of the PDE
        rhs%sigma = tmpV     - matmul(op%lift, fluxSigma)
        rhs%v     = tmpSigma - matmul(op%lift, fluxV)
    end subroutine rhsElastic
end module
