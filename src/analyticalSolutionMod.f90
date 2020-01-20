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
module analyticalSolutionMod
    use parameterMod
    use constantsMod
    use matrixMod
    use errorMessage
    use slipInterfaceMod
    use genericMaterial
    use sourcesMod
    implicit none

    type :: analyticalVar
        real(kind=custom_real) :: wavect
        real(kind=custom_real), dimension(:,:,:), allocatable :: misfit
        real(kind=custom_real), dimension(:,:,:), allocatable :: q
        real(kind=custom_real), dimension(:,:,:), allocatable :: h
    end type

    contains

    subroutine analyticalStep(tstep, t, mesh, genmat, lsi, V, ana, lsi_spec)
        !input
        type(meshVar) :: mesh
        type(lsiVar), dimension(:) :: lsi
        type (interface_spec), dimension(:) :: lsi_spec
        type(generic_material) :: genmat
        integer, intent(in) :: tstep
        real(kind=custom_real), intent(in) :: t
        !output
        type(analyticalVar) :: ana
        !local
        character(len=7) :: step
        integer :: i,j
        real(kind=custom_real) :: nu
        real(kind=custom_real) :: dx
        real(kind=custom_real) :: x
        real(kind=custom_real) :: xslip
        real(kind=custom_real) :: beta
        real(kind=custom_real) :: zimp
        real(kind=custom_real) :: expo
        real(kind=custom_real) :: ltor_amp
        real(kind=custom_real) :: rtol_amp
        real(kind=custom_real) :: hright
        real(kind=custom_real) :: hleft, etaz
        real(kind=custom_real), dimension(:,:) :: V

        dx = (mesh%xmax-mesh%xmin)/(mesh%ncell)
        etaz = genmat.imp.(mesh%etom(lsi(1)%prop))
        nu = lsi_spec(lsi(1)%prop)%nu
        xslip = lsi(1)%loc !mesh%xmin+lsiloc*dx
        ltor_amp = 1.
        rtol_amp = 0.

        write(step, "(i7.7)") tstep
        open(21, file = "output/misfit"//step)

        do i = 1, mesh%k
            do j = 1, mesh%np
                beta = genmat.vs.(mesh%etom(i))
                zimp = genmat.imp.(mesh%etom(i))
                x = mesh%x(j,i)
                if (t - (x-xslip)/beta > 0) then; hright = 1.0; else; hright = 0.; endif
                if (t + (x-xslip)/beta > 0) then; hleft = 1.0; else; hleft = 0.; endif
                if (x <= xslip) then
                    expo = exp(-nu*(t+(x-xslip)/beta))
                    ana%q(1,j,i) = ltor_amp*hright+rtol_amp*hleft-(rtol_amp-ltor_amp)*expo*hleft
                    ana%q(2,j,i) = -zimp*ltor_amp*hright+zimp*rtol_amp*hleft-zimp*(rtol_amp-ltor_amp)*expo*hleft
                else if (x > xslip) then
                    expo = exp(-nu*(t-(x-xslip)/beta))
                    ana%q(1,j,i) = rtol_amp*hleft+ltor_amp*hright+(rtol_amp-ltor_amp)*expo*hright
                    ana%q(2,j,i) = zimp*rtol_amp*hleft-zimp*ltor_amp*hright-zimp*(rtol_amp-ltor_amp)*expo*hright
                end if
                ana%misfit(tstep, j, i) = sqrt(((V(j,i) - ana%q(1,j,i))**2)/ana%q(1,j,i)**2)*100
                write(21,*) x, ana%misfit(tstep,j,i)
            enddo
        end do
    end subroutine

    subroutine analyticalSinSquared(tstep, time, mesh, src, genmat, lsi, V, ana, lsi_spec)
        !input
        type(meshVar) :: mesh
        type(lsiVar), dimension(:) :: lsi
        type (interface_spec), dimension(:) :: lsi_spec
        type(sourceVar) :: src
        type(generic_material) :: genmat
        integer, intent(in) :: tstep
        real(kind=custom_real), intent(in) :: time

        !output
        type(analyticalVar) :: ana
        !local
        character(len=7) :: step
        integer :: i,j
        real(kind=custom_real) :: nu
        real(kind=custom_real) :: dx
        real(kind=custom_real) :: x, x0, w
        real(kind=custom_real) :: xslip
        real(kind=custom_real) :: beta
        real(kind=custom_real) :: zimp
        real(kind=custom_real) :: etaz
        real(kind=custom_real) :: S
        real(kind=custom_real), dimension(:,:) :: V

        dx = (mesh%xmax-mesh%xmin)/(mesh%ncell)
        etaz = genmat.imp.(mesh%etom(lsi(1)%prop))
        nu = lsi_spec(lsi(1)%prop)%nu
        xslip = lsi(1)%loc
        w = src%width
        x0 = src%center

        write(step, "(i7.7)") tstep
        open(21, file = "output/misfit"//step)

        do i = 1, mesh%K
            do j = 1, mesh%np
                beta = genmat.vs.(mesh%etom(i))
                zimp = genmat.imp.(mesh%etom(i))
                x = mesh%x(j,i)
                if (-w/2 < x-x0-beta*time .and. x-x0-beta*time < w/2) then
                    ana%q(1,j,i) = src%amplitude * cos(pi/w*(x - x0 - beta*time))**2
                    ana%q(2,j,i) = -zimp*src%amplitude*cos(pi/w*(x - x0 - beta*time))**2
                else
                    ana%q(:,j,i) = 0.
                endif
                if (i <= lsi(1)%elements(2)) then
                     call analyticalS(src, beta, w, nu, time, x0, x, xslip, S)
                     ana%q(1,j,i) = ana%q(1,j,i) + S
                     ana%q(2,j,i) = ana%q(2,j,i) + zimp*S
                 else if (i > lsi(1)%elements(2)) then
                     call analyticalS(src, beta, w, nu, time, x0, x, xslip, S)
                     ana%q(1,j,i) = ana%q(1,j,i) - S
                     ana%q(2,j,i) = ana%q(2,j,i) + zimp*S
                 endif
                 ana%misfit(tstep, j, i) = V(j,i) - ana%q(1,j,i)
            enddo
        enddo
        ana%misfit(tstep,:,:) = 100*(abs(ana%misfit(tstep,:,:))/maxval(ana%q(1,:,:)))
        write(21,'(f8.3, e10.2)') ((mesh%x(j,i), ana%misfit(tstep,j,i), j=1,mesh%np), i = 1 + mesh%nghost,mesh%K-mesh%nghost)
        close(21)
    end subroutine

    subroutine analyticalS(src, beta, w, nu, time, x0, x, xslip, S)
        !input
        type(sourceVar) :: src
        real(kind=custom_real), intent(in) :: x
        real(kind=custom_real), intent(in) :: x0
        real(kind=custom_real), intent(in) :: w
        real(kind=custom_real), intent(in) :: nu
        real(kind=custom_real), intent(in) :: time
        real(kind=custom_real), intent(in) :: xslip
        real(kind=custom_real), intent(in) :: beta
        !output
        real(kind=custom_real), intent(out) :: S
        !local
        real(kind=custom_real) :: T1
        real(kind=custom_real) :: upper
        real(kind=custom_real) :: lower
        real(kind=custom_real) :: k
        real(kind=custom_real) :: b

        k = nu/beta
        b = 2*pi/w

        T1 = -b/(2*(k**2 + b**2)) * exp(-k*(beta*time + (x0-xslip) -abs(x-xslip)))
        if (abs(x-xslip)-(x0-xslip)-beta*time > w/2) then
            S = 0.
        else if (-w/2 < abs(x-xslip)-(x0-xslip)-beta*time .and. abs(x-xslip)-(x0-xslip)-beta*time < w/2) then
            upper = antiderivative(src, k, b, w/2)
            lower = antiderivative(src, k, b, abs(x-xslip)-(x0-xslip)-beta*time)
            S = T1*(upper-lower)
        else if (abs(x-xslip)-(x0-xslip)-beta*time < -w/2) then
            upper = antiderivative(src, k, b, w/2)
            lower = antiderivative(src, k, b, -w/2)
            S = T1*(upper-lower)
        end if
    end subroutine

    function antiderivative(src, k, b, boundary) result(anti)
        !input
        type(sourceVar) :: src
        real(kind=custom_real), intent(in) :: boundary
        real(kind=custom_real), intent(in) :: k
        real(kind=custom_real), intent(in) :: b
        !output
        real(kind=custom_real) :: anti

        anti = exp(-k*boundary) * src%amplitude * (k*sin(b*boundary) + b*cos(b*boundary))
    end function

    subroutine allocAnalyticalArrays(ana, par, mesh)
        type (parameterVar) :: par
        type (meshVar) :: mesh
        type (analyticalVar) :: ana
        allocate(ana%q(2,mesh%np,mesh%k))
        allocate(ana%h(2,mesh%np,mesh%k))
        allocate(ana%misfit(par%tsteps,mesh%np,mesh%k))
    end subroutine

    subroutine deallocAnalyticalArrays(ana)
        type (analyticalVar) :: ana
        if (allocated(ana%q)) deallocate(ana%q)
        if (allocated(ana%h)) deallocate(ana%h)
        if (allocated(ana%misfit)) deallocate(ana%misfit)
    end subroutine
end module
