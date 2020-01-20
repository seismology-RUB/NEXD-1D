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
module rhsSlipInterface1DMod

    use parameterMod
    use constantsMod
    use matrixMod
    use riemannflux1Dmod
    use genericMaterial
    use slipInterfaceMod
    use localOperatorsMod
    use errorMessage
    implicit none

    contains

    subroutine rhsAddSlipInterfaceInfluence(mesh, genmat, op, lsi, fluxV, fluxSigma)
        !input
        type (meshVar), intent(in) :: mesh
        type (generic_material), intent(in) :: genmat
        type (operator), intent(in) :: op
        type(lsiVar) :: lsi
        !in-/output
        real(kind=custom_real), dimension(:,:), intent(inout) :: fluxV
        real(kind=custom_real), dimension(:,:), intent(inout) :: fluxSigma
        !local variables
        integer :: i, j, m, nif, nmf
        real(kind=custom_real) :: zin,zout,beta
        real(kind=custom_real), dimension(2) :: fluxR
        real(kind=custom_real), dimension(2) :: fluxL

        do j = 1, mesh%nfaces
            i = lsi%elements(j)
            m = mesh%etoe(i,j)
            nmf = mesh%ntof(mesh%etof(i,j))
            nif = mesh%ntof(j)
            beta = genmat.vs.(mesh%ntom(nif,i))
            zin = genmat.imp.(mesh%ntom(nif,i))
            zout = genmat.imp.(mesh%ntom(nmf,m))
            if (j == 1) then
                call rightRiemannFluxLSI(lsi%S(j), beta, zin, zout, fluxR)
                fluxV(j,i) = fluxV(j,i) + (mesh%nx(j,i) * op%fscale(j,i) * fluxR(1))
                fluxSigma(j,i) = fluxSigma(j,i) + (mesh%nx(j,i) * op%fscale(j,i) * fluxR(2))
            else if (j == 2) then
                call leftRiemannFluxLSI(lsi%S(j), beta, zin, zout, fluxL)
                fluxV(j,i) = fluxV(j,i) + (mesh%nx(j,i) * op%fscale(j,i) * fluxL(1))
                fluxSigma(j,i) = fluxSigma(j,i) + (mesh%nx(j,i) * op%fscale(j,i) * fluxL(2))
            endif
        enddo
    end subroutine

    subroutine activityFluxSlipInterface(rhsV, rhsSigma, genmat, mesh, lsi, lsi_spec)
        !input
        type (generic_material) :: genmat
        type (meshVar) :: mesh
        type (interface_spec), dimension(:) :: lsi_spec
        real(kind=custom_real), dimension(:,:), intent(in) :: rhsV      !time derivative of the velocity
        real(kind=custom_real), dimension(:,:), intent(in) :: rhsSigma  !time derivative of the stress.
        !input/output
        type (lsiVar) :: lsi
        !local
        integer :: i,j,m,nif, nmf, np
        real(kind=custom_real) :: zin
        real(kind=custom_real) :: zout
        real(kind=custom_real) :: nu
        real(kind=custom_real) :: dV

        do j = 1, mesh%nfaces
            i = lsi%elements(j)
            m = mesh%etoe(i,j)
            nif = mesh%ntof(j)
            nmf = mesh%ntof(mesh%etof(i,j))
            dV = (rhsV(nmf, m) - rhsV(nif,i))*mesh%nx(j,i)
            zin = genmat.imp.(mesh%etom(i))
            zout= genmat.imp.(mesh%etom(m))
            select case (trim(lsi_spec(lsi%prop)%type))
                case ("elastic")
                    nu = lsi_spec(lsi%prop)%nu
                    lsi%rhsS(j) = -nu*lsi%S(j) + ((zin+zout)/(2*zin*zout))*rhsSigma(nif,i) + 0.5*dV
            end select
        enddo
    end subroutine activityFluxSlipInterface
end module
