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
program analyticalSolution

    use parameterMod
    use constantsMod
    use matrixMod
    use errorMessage
    use slipInterfaceMod
    use materials1DMod

    implicit none

    type(meshVar) :: mesh
    type(parameterVar) :: par
    type(materialVar) :: mat
    type(lsiVar) :: lsi
    type(error_message) :: errmsg
    type (interface_spec), dimension(:), allocatable :: lsi_spec
    type(rovszimp_material) :: rvzmat
    character(len=30) :: filename
    character(len=10) :: filenr
    integer :: i ,j
    integer :: lu
    integer :: lsiloc
    real(kind=custom_real) :: t
    real(kind=custom_real) :: t0
    real(kind=custom_real) :: dt
    real(kind=custom_real) :: nu
    real(kind=custom_real) :: dx
    real(kind=custom_real) :: xl
    real(kind=custom_real) :: xr
    real(kind=custom_real) :: xm
    real(kind=custom_real) :: xslip
    real(kind=custom_real) :: beta
    real(kind=custom_real) :: zimp
    real(kind=custom_real) :: expo
    real(kind=custom_real) :: ltor_amp
    real(kind=custom_real) :: rtol_amp
    real(kind=custom_real) :: hright
    real(kind=custom_real) :: hleft
    real(kind=custom_real), dimension(:,:,:), allocatable :: q

    call readParfile(par, mesh, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    call readSlipInterfaceProperties(par, lsi, mesh, 1, "input/lsi_eap" , errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    !call readTransverseSlipInterface(par, tsi, lsi, 1, "input/interfaces" ,errmsg)
    !Read the specifications of all possible interface types
    call readSlipInterfaceSpecs(lsi_spec, 1, "input/interfaces" ,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    call readMaterialProperties(mat, rvzmat, 1, "input/material", errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    allocate(q(2,mesh%K,par%tsteps))

    lu = 1
    dx = (mesh%xmax-mesh%xmin)/(mesh%ncell)
    beta = rvzmat%vs(1)
    zimp = rvzmat%zimp(1)
    if (size(lsi_spec) > 1) then
        call add(errmsg,2, "The analytical solution is only valid for a single interface and thus only a single interface type)",myname)
        call print(errmsg)
        stop
    else
        nu = lsi_spec(1)%nu
    end if
    lsiloc = 101-1
    xslip = mesh%xmin+lsiloc*dx
    ltor_amp = 1.
    rtol_amp = 0.
    t0 = 0.
    print *, dx
    dt = 0.5*par%cfl*dx/beta

    print *,'Time step: ',dt
    do j = 1, par%tsteps
        t = t0 + (j-1)*dt
        if (mod(j,100) == 0) write(6,'(i5,$)'), j
        write(filenr,"(I5.5)") j
        filename = "output/"//trim(filenr)
        open(lu,file = filename)
        open(2, file="output/S", position="append")
        do i = 1, mesh%K
            xl = mesh%xmin+(i-1)*dx
            xr = mesh%xmin+(i)*dx
            xm = 0.5*(xl+xr)
            if (t - (xm-xslip)/beta > 0) then; hright = 1.0; else; hright = 0.; endif
            if (t + (xm-xslip)/beta > 0) then; hleft = 1.0; else; hleft = 0.; endif
            if (i <= lsiloc) then
                expo = exp(-nu*(t+(xm-xslip)/beta))
                q(1,i,j) = ltor_amp*hright+rtol_amp*hleft-(rtol_amp-ltor_amp)*expo*hleft
                q(2,i,j) = -zimp*ltor_amp*hright+zimp*rtol_amp*hleft-zimp*(rtol_amp-ltor_amp)*expo*hleft
                if (i == lsiloc) write(2,*) t, (rtol_amp-ltor_amp)*expo
            else if (i > lsiloc) then
                expo = exp(-nu*(t-(xm-xslip)/beta))
                q(1,i,j) = rtol_amp*hleft+ltor_amp*hright+(rtol_amp-ltor_amp)*expo*hright
                q(2,i,j) = zimp*rtol_amp*hleft-zimp*ltor_amp*hright-zimp*(rtol_amp-ltor_amp)*expo*hright
            end if
            write(lu,*) xl, q(1,i,j), q(2,i,j)
            write(lu,*) xr, q(1,i,j), q(2,i,j)
        end do
        close(lu)
        close(2)
    end do

    deallocate(q)
    call deallocLsiArrays(lsi)
    call deallocMatArrays(mat)
    call deallocMeshArrays(mesh)
    call deallocTransverseSlipInterface(tsi)
    call deallocRoVsZimpMaterial(rvzmat)
    call deallocErrorMessage(errmsg)
    print *, ""
end program  analyticalSolution
