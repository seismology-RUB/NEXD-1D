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
program solverPro

    use constantsMod
    use mesh1DMod
    use parameterMod
    use matrixMod
    use outputMod
    use materials1DMod
    use rhsElastic1DMod
    use rhsPoroelastic1DMod
    use rhsSlipInterface1DMod
    use errorMessage
    use sourcesMod
    use slipInterfaceMod
    use genericMaterial
    use slopeLimiterMod
    use analyticalSolutionMod
    use receiverMod
    use matrixMod
    use boundaryConditionsMod

    implicit none

    type (parameterVar) :: par
    type (meshVar) :: mesh
    type (materialVar) :: mat
    type (error_message) :: errmsg
    type (rovszimp_material) :: rvzmat
    type (porous_material) :: poromat
    type (generic_material) :: genmat
    type (sourceVar), dimension(:), allocatable :: src
    type (operator) :: op
    type (limiterVar) :: lim
    type (analyticalVar) :: ana
    type (receiverVar), dimension(:), allocatable :: rec
    type (lsiVar), dimension(:), allocatable :: lsi
    type (interface_spec), dimension(:), allocatable :: lsi_spec
    type (timeStampVar) :: timestamp
    character(len=10) :: myname = "solver"
    character(len=10) :: filenr
    character(len=30) :: filename ,anasol
    character(len=80) :: outfile
    integer :: i, intrk
    integer :: tstep
    integer :: np
    integer :: out
    real(kind=custom_real) :: time
    real(kind=custom_real) :: vmax, vmin
    real(kind=custom_real) :: dt
    real(kind=custom_real) :: energy
    real(kind=custom_real) :: xmin
    real(kind=custom_real), dimension(:), allocatable :: timevec                                       !Vector that contains the timeindex for each timstep
    real(kind=custom_real), dimension(:), allocatable :: r
    real(kind=custom_real), dimension(:,:), allocatable :: Zimp      !Impedance array
    real(kind=custom_real), dimension(:,:), allocatable :: Vdm, invV    !Vandermonde- and inverted vandermonde Matrix

    type (solutionvector) :: solution   !solution vector
    type (solutionvector) :: src_vec    !source vector
    type (solutionvector) :: rhs        !right hand side of the DG equation
    type (solutionvector) :: res        !Runge Kutta residual
    type (solutionvector) :: rk1, rk2   !Runge-Kutta (SSP-RK) storage variables

    !Create new errormessagetrace
    call new(errmsg,myname)

    !! SETUP the environment (read inputfiles, create mesh and initialize wavelet)
    call readParfile(par, mesh, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    call allocSolverArrays()
    call allocOperatorArrays(mesh,op)
    !allocate necessary arrarys to comapre the numerical step function response to the analytical solution
    if (par%compare) call allocAnalyticalArrays(ana, par, mesh)
    call setup_mesh(mesh, op, vdm, invV, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    if (par%poroelastic) then
        call setup_materials(par, mat, mesh, poromat, errmsg)
        call associatePorousToGenericMaterial(genmat,poromat)
    else
        call setup_materials(par, mat, mesh, rvzmat, errmsg)
        call associateRoVsZimpToGenericMaterial(genmat,rvzmat)
    endif

    call compute_timestep()

    call init_sources(par, mesh, genmat, vdm, src, solution, dt, errmsg)

    if (par%recn > 0) call initReceiver(par, mesh, rec, vdm, 1, "input/receiver", errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    !LSI-only
    if (par%lsin /= 0) then
        allocate(lsi(par%lsin))
        call setup_lsi(par, mesh, genmat, lsi, lsi_spec, solution%v, errmsg)
    endif

    if (par%movie) call writeInitialSnapshot(par, mesh, ana, genmat, src, solution, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    if (par%log) write(*,"(a80)") "|                             starting timeloop...                             |"
    if (par%log) write(*,"(a80)") "|------------------------------------------------------------------------------|"
    if (par%logfile) then
        outfile = "output/logfile"//trim(par%title)
        outfile = trim(outfile)
        out = 16
        open(out, file = outfile, position="append")
        write(out,"(a80)") "|                             starting timeloop...                             |"
        write(out,"(a80)") "|------------------------------------------------------------------------------|"
    endif

    !create initial timestamp
    call initialTimestamp(timestamp)

    !! SOLVE
    call propagate()

    if (par%recn > 0) then
        call plotSeismo1D(par, src, rec, timevec)
    endif

    print *," "
    call deallocall()

    contains

    subroutine compute_timestep()
        !compute time step size
        ! find smallest distance between two nodes
        xmin = mesh%vx(2)-mesh%vx(1)
        do i = 3, size(mesh%vx)
            xmin = min(xmin, mesh%vx(i)-mesh%vx(i-1))
        enddo

        !!From Hesthaven & Warburton:
        !rmin = 1.
        !do i = 2, size(r)
        !   tmp = abs(r(i-1) - r(i))
        !   if (tmp <= rmin) rmin = tmp
        !enddo

        ! find maximum and minimum velocity
        vmax = 0.
        vmin = 300000. !;-)
        if (par%poroelastic) then
            do i = 1, size(mat%i)
                vmax = max(vmax,poromat%vmax(mat%i(i))) !this is the maximum eigenvalue of A (only valid for the inviscid case!)
                vmin = min(vmin,poromat%vmin(mat%i(i))) !this is the minimum eigenvalue of A (only valid for the inviscid case!)
            enddo
        else
            do i = 1, size(mat%i)
                if (rvzmat%vs(mat%i(i)) > vmax) vmax = rvzmat%vs(mat%i(i))
                if (rvzmat%vs(mat%i(i)) < vmin) vmin = rvzmat%vs(mat%i(i))
            enddo
        endif

        if (par%logfile) then
            outfile = "output/logfile"//trim(par%title)
            outfile = trim(outfile)
            out = 16
            open(out, file = outfile, position="append")
        endif
        ! calculate time step
        if (par%tint == 1) then
            if (par%log) write (*, "(a80)") "|                      Timeintegration: Euler                                  |"
            if (par%logfile) then
                write (out, "(a80)") "|                      Timeintegration: Euler                                  |"
            endif
            if (par%autodt) then
                dt = par%cfl * 0.5 * xmin/vmax/(mesh%np-1)
            else
                dt = par%dt
            endif
        else if (par%tint == 3) then
            if (par%log) write (*, "(a80)") "|                      Timeintegration: 3rd-Order SSP-Runge-Kutta              |"
            if (par%logfile) then
                write (out, "(a80)") "|                      Timeintegration: 3rd-Order SSP-Runge-Kutta              |"
            endif
            if (par%autodt) then
                dt = par%cfl * 0.5 * xmin/vmax/(mesh%np-1)
            else
                dt = par%dt
            endif
            allocate(rk1%v(mesh%Np,mesh%K))
            allocate(rk2%v(mesh%Np,mesh%K))
            allocate(rk1%sigma(mesh%Np,mesh%K))
            allocate(rk2%sigma(mesh%Np,mesh%K))
        else if (par%tint == 5)then
            if (par%log) write (*, "(a80)") "|                      Timeintegration: 4th-Order Runge-Kutta                  |"
            if (par%logfile) then
                write (out, "(a80)") "|                      Timeintegration: 4th-Order Runge-Kutta                  |"
            endif
            if (par%autodt) then
                dt = par%cfl * 0.5 * xmin/vmax/(mesh%np-1)
            else
                dt = par%dt
            endif
        else
            call add(errmsg, 2, "No accepted timeintegration. Abort!", myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif

        if (par%log) then
            write (*,"(a40, es12.5, a28)")   "|                            Time step: ", dt, "                           |"
            write (*,"(a40, es12.5, a28)")   "|                                DXMIN: ", xmin,&
                                             "                           |"
            write (*,"(a40, f10.5, a30)")   "|                     Maximum velocity: ", vmax, &
                                            "                             |"
            if (vmin >= 0000.00001) then
                write (*,"(a40, f10.5, a30)")   "|                     Minimum velocity: ", vmin, &
                                            "                             |"
            else
                write (*,"(a40, es10.3, a30)")   "|                     Minimum velocity: ", vmin, &
                                            "                             |"
            endif
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        endif

        if (par%logfile) then
            write (out,"(a40, es12.5, a28)")   "|                            Time step: ", dt, "                           |"
            write (out,"(a40, es12.5, a28)")   "|                                DXMIN: ", xmin,&
                                             "                           |"
            write (out,"(a40, f10.5, a30)")   "|                     Maximum velocity: ", vmax, &
                                            "                             |"
            if (vmin >= 0000.00001) then
                write (out,"(a40, f10.5, a30)")   "|                     Minimum velocity: ", vmin, &
                                            "                             |"
            else
                write (out,"(a40, es10.3, a30)")   "|                     Minimum velocity: ", vmin, &
                                            "                             |"
            endif
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            close(out)
        endif
    end subroutine

    subroutine propagate()
        !Propagation
        res%v = 0.
        res%sigma = 0.
        src_vec%v = 0.
        if (par%poroelastic) then
            res%v1 = 0.
            res%p1 = 0.
            src_vec%v1 = 0.
            if (par%fluidn == 2) then
                res%v2 = 0.
                res%p2 = 0.
                src_vec%v2 = 0.
            endif
        endif

        if (par%compare .and. src(1)%type == "sin2") ana%wavect = src(1)%center

        do tstep = 1, par%tsteps
            time = tstep*dt
            timevec(tstep) = (tstep-1)*dt
            call rhsSource(par, mesh, op, genmat, src, tstep, src_vec)
            select case(par%tint)
                case (1) !Euler
                    if (par%poroelastic) then
                        call rhsPoroelastic(par, mesh, genmat, op, solution, rhs)
                    else
                        call rhsElastic(par, mesh, genmat, op, solution, lsi, rhs)
                    endif
                    if (par%lsin /= 0) then
                        !LSI-only
                        do i = 1, par%lsin
                            call activityFluxSlipInterface(rhs%v, rhs%sigma, genmat, mesh, lsi(i), lsi_spec)
                            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                            lsi(i)%S = lsi(i)%S + dt*lsi(i)%rhsS
                        enddo
                    endif
                    solution%v              = solution%v     + dt*(rhs%v     + src_vec%v)
                    solution%sigma          = solution%sigma + dt*rhs%sigma
                    if (par%poroelastic) then
                            solution%v1     = solution%v1    + dt*(rhs%v1    + src_vec%v1)
                            solution%p1     = solution%p1    + dt*rhs%p1
                            if (par%fluidn == 2) then
                                solution%v2 = solution%v2    + dt*(rhs%v2    + src_vec%v2)
                                solution%p2 = solution%p2    + dt*rhs%p2
                            endif
                    endif
                case (3) !3rd-Order SSP-Runge-Kutta
                    !SSP-RK Stage 1
                    if (par%poroelastic) then
                        call rhsPoroelastic(par, mesh, genmat, op, solution, rhs)
                    else
                        call rhsElastic(par, mesh, genmat, op, solution, lsi, rhs)
                    endif
                    rk1%v     = solution%v     + dt*(rhs%v + src_vec%v)
                    rk1%sigma = solution%sigma + dt*rhs%sigma
                    if (par%poroelastic) then
                        rk1%v1     = solution%v1     + dt*(rhs%v1 + src_vec%v1)
                        rk1%p1     = solution%p1     + dt*rhs%p1
                        if (par%fluidn == 2) then
                            rk1%v2     = solution%v2     + dt*(rhs%v2 + src_vec%v2)
                            rk1%p2     = solution%p2     + dt*rhs%p2
                        endif
                    endif
                    call selectLimiter(par, mesh, vdm, invV, rk1%v, rk1%sigma, errmsg)
                    !SSP-RK Stage 2
                    if (par%poroelastic) then
                        call rhsPoroelastic(par, mesh, genmat, op, rk1, rhs)
                    else
                        call rhsElastic(par, mesh, genmat, op, rk1, lsi, rhs)
                    endif
                    rk2%v     = (3*solution%v     + rk1%v     + dt*(rhs%v + src_vec%v))/4.
                    rk2%sigma = (3*solution%sigma + rk1%sigma + dt*rhs%sigma)/4.
                    if (par%poroelastic) then
                        rk2%v1 = (3*solution%v1 + rk1%v1 + dt*(rhs%v1 + src_vec%v1))/4.
                        rk2%p1 = (3*solution%p1 + rk1%p1 + dt*rhs%p1)/4.
                        if (par%fluidn == 2) then
                            rk2%v2 = (3*solution%v2 + rk1%v2 + dt*(rhs%v2 + src_vec%v))/4.
                            rk2%p2 = (3*solution%p2 + rk1%p2 + dt*rhs%p2)/4.
                        endif
                    endif
                    call selectLimiter(par, mesh, vdm, invV, rk2%v, rk2%sigma, errmsg)
                    !SSP-RK Stage 3
                    if (par%poroelastic) then
                        call rhsPoroelastic(par, mesh, genmat, op, rk2, rhs)
                    else
                        call rhsElastic(par, mesh, genmat, op, rk2, lsi, rhs)
                    endif
                    solution%v     = (solution%v     + 2*rk2%v     + 2*dt*(rhs%v + src_vec%v))/3.
                    solution%sigma = (solution%sigma + 2*rk2%sigma + 2*dt*rhs%sigma)/3.
                    if (par%poroelastic) then
                        solution%v1 = (solution%v1 + 2*rk2%v1 + 2*dt*(rhs%v1 + src_vec%v1))/3.
                        solution%p1 = (solution%p1 + 2*rk2%p1 + 2*dt*rhs%p1)/3.
                        if (par%fluidn == 2) then
                            solution%v2 = (solution%v2 + 2*rk2%v2 + 2*dt*(rhs%v2 + src_vec%v2))/3.
                            solution%p2 = (solution%p2 + 2*rk2%p2 + 2*dt*rhs%p2)/3.
                        endif
                    endif
                    call selectLimiter(par, mesh, vdm, invV, solution%v, solution%sigma, errmsg)
                case (5) !4th-Order Runge-Kutta
                    do intrk = 1, par%tint
                        if (par%poroelastic) then
                            call rhsPoroelastic(par, mesh, genmat, op, solution, rhs)
                        else
                            call rhsElastic(par, mesh, genmat, op, solution, lsi, rhs)
                        endif
                        if (par%lsin /= 0) then
                            !LSI-only
                            do i = 1, par%lsin
                                call activityFluxSlipInterface(rhs%v, rhs%sigma, genmat, mesh, lsi(i), lsi_spec)
                                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                                lsi(i)%resS = rk4a(intrk)*lsi(i)%resS + dt*lsi(i)%rhsS
                                lsi(i)%S = lsi(i)%S + rk4b(intrk)*lsi(i)%resS
                            enddo
                        endif
                        res%v     = rk4a(intrk)*res%v     + dt*(rhs%v + src_vec%v)
                        solution%v     = solution%v     + rk4b(intrk)*res%v
                        res%sigma = rk4a(intrk)*res%sigma + dt*(rhs%sigma)
                        solution%sigma = solution%sigma + rk4b(intrk)*res%sigma
                        if (par%poroelastic) then
                            res%v1    = rk4a(intrk)*res%v1    + dt*(rhs%v1 + src_vec%v1)
                            solution%v1    = solution%v1    + rk4b(intrk)*res%v1
                            res%p1    = rk4a(intrk)*res%p1    + dt*rhs%p1
                            solution%p1    = solution%p1    + rk4b(intrk)*res%p1
                            if (par%fluidn == 2) then
                                res%v2    = rk4a(intrk)*res%v2    + dt*(rhs%v2 + src_vec%v2)
                                solution%v2    = solution%v2    + rk4b(intrk)*res%v2
                                res%p2    = rk4a(intrk)*res%p2    + dt*rhs%p2
                                solution%p2    = solution%p2    + rk4b(intrk)*res%p2
                            endif
                        endif
                    enddo
                case default
                    stop "No accepted timeintegration. Abort!"
            end select

            !check if "NaN" occured
            write (filenr,"(I7.7)") tstep
            if (any(isnan(solution%v))) then
                call add(errmsg,2,'NaN occurred in velocity solution at step '//trim(filenr)//'.',myname)
            elseif (any(isnan(solution%sigma))) then
                call add(errmsg,2,'NaN occurred in stress solution at step '//trim(filenr)//'.',myname)
            elseif (par%poroelastic .and. par%fluidn == 1) then
                if (any(isnan(solution%v1))) then
                    call add(errmsg,2,'NaN occurred in velocity solution of fluid 1 at step '//trim(filenr)//'.',myname)
                elseif (any(isnan(solution%p1))) then
                    call add(errmsg,2,'NaN occurred in pressure solution of fluid 1 at step '//trim(filenr)//'.',myname)
                endif
            elseif (par%poroelastic .and. par%fluidn == 2) then
                if (any(isnan(solution%v2))) then
                    call add(errmsg,2,'NaN occurred in velocity solution of fluid 2 at step '//trim(filenr)//'.',myname)
                elseif (any(isnan(solution%p2))) then
                    call add(errmsg,2,'NaN occurred in pressure solution of fluid 2 at step '//trim(filenr)//'.',myname)
                endif
            endif
            if (.level.errmsg == 2) then
                write(filenr,"(I7.7)") tstep
                filename = "output/"//trim(filenr)
                call writeSnapshot1d(solution, mesh, 1, filename)
                call print(errmsg)
                stop
            endif

            if (par%recn > 0) call recordSeismograms(par, mesh, rec, tstep, solution)

            if (par%compare) then
                select case(src(1)%type)
                    case ("box")
                        call analyticalStep(tstep, time, mesh, genmat, lsi, solution%v, ana, lsi_spec)
                    case ("sin2")
                        call analyticalSinSquared(tstep, time, mesh, src(1), genmat, lsi, solution%v, ana, lsi_spec)
                    case ("test")
                        call add(errmsg,2,'Analytical solution for the test wavelet does not exist!',myname)
                        return
                    case default
                        call add(errmsg,2,'No analytical solution is known for this problem. Abort...',myname)
                        return
                end select
            endif
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif

            !absorbing boundary conditions
            if (trim(mesh%rbc) == "absorbing" .or. trim(mesh%lbc) == "absorbing") then
                call absorbingBC(mesh, solution)
            endif

            !Write output
            if (par%movie) then
                if (mod(tstep,par%framestep) == 0 .or. tstep == par%tsteps) then
                    !write(6,'(i7,$)'), tstep
                    write(filenr,"(I7.7)") tstep
                    filename = "output/"//trim(filenr)
                    !print *, filename
                    call writeSnapshot1d(solution, mesh, 1, filename)
                    if (par%poroelastic) then
                        call writeEnergyPoroelastic(solution, mesh, poromat, time, energy, 1)
                    else
                        call writeEnergyElastic(solution, mesh, genmat, time, energy, 1)
                    endif
                    if (par%compare) then
                        anasol = "output/"//trim(filenr)//"ana"
                        call writeSnapshotAna1d(ana%q(1,:,:), mesh,2,anasol)
                        call plot1D(par, mesh, src, tstep,"velocity", "ana_velocity", ana%q(1,:,:))
                    endif
                    call plot1D(par, mesh, src, tstep, "velocity", "velocity", solution%v)
                    call plot1D(par, mesh, src, tstep, "stress", "stress", solution%sigma)
                    if (par%poroelastic) then
                        call plot1D(par, mesh, src, tstep, "v1", "v1", solution%v1)
                        call plot1D(par, mesh, src, tstep, "p1", "p1", solution%p1)
                        if (par%fluidn == 2) then
                            call plot1D(par, mesh, src, tstep, "v2", "v2", solution%v2)
                            call plot1D(par, mesh, src, tstep, "p2", "p2", solution%p2)
                        endif
                    endif
                endif
            endif
            if (mod(tstep,par%framestep) == 0 .or. tstep == par%tsteps) then
                call outputStamp(par, time, tstep, timestamp, solution, energy)
            endif
        enddo
        if (par%logfile) then
            outfile = "output/logfile"//trim(par%title)
            outfile = trim(outfile)
            out = 16
            open(out, file = outfile, position="append")
            write(out,"(a80)") "|                              timeloop finished                               |"
            write(out,"(a80)") "|------------------------------------------------------------------------------|"
            close(out)
        endif
        if (par%log) write(*,"(a80)") "|                              timeloop finished                               |"
        if (par%log) write(*,"(a80)") "|------------------------------------------------------------------------------|"

    end subroutine

    subroutine deallocall()
        call deallocSolverArrays()
        call deallocOperatorArrays(op)
        !call deallocLsiArrays(lsi, lsi_spec)
        call deallocMatArrays(mat)
        call deallocMeshArrays(mesh)
        call deallocGenericMaterial(genmat)
        if (par%lsin /= 0) then
            !LSI-only
            if (allocated(lsi)) deallocate(lsi)
            if (allocated(lsi_spec)) deallocate(lsi_spec)
        endif
        if (par%poroelastic) then
            call deallocPorousMaterial(poromat)
        else
            call deallocRoVsZimpMaterial(rvzmat)
        endif
        if (allocated(src)) deallocate(src)
        call deallocErrorMessage(errmsg)
        call deallocLimiterArrays(lim)
        if (allocated(rec)) deallocate(rec)
        if (par%compare) call deallocAnalyticalArrays(ana)
    end subroutine

    subroutine allocSolverArrays()
        allocate(r(mesh%np))
        allocate(solution%v(mesh%Np,mesh%K))
        allocate(solution%sigma(mesh%Np,mesh%K))
        allocate(rhs%v(mesh%Np,mesh%K))
        allocate(rhs%sigma(mesh%Np,mesh%K))
        allocate(res%v(mesh%Np,mesh%K))
        allocate(res%sigma(mesh%Np,mesh%K))
        allocate(src_vec%v(mesh%Np, mesh%K))
        allocate(timevec(par%tsteps))
        if (par%poroelastic) then
            allocate(solution%v1(mesh%Np,mesh%K))
            allocate(solution%p1(mesh%Np,mesh%K))
            allocate(rhs%v1(mesh%Np,mesh%K))
            allocate(rhs%p1(mesh%Np,mesh%K))
            allocate(res%v1(mesh%Np,mesh%K))
            allocate(res%p1(mesh%Np,mesh%K))
            allocate(src_vec%v1(mesh%Np, mesh%K))
            if (par%fluidn == 2) then
                allocate(solution%v2(mesh%Np,mesh%K))
                allocate(solution%p2(mesh%Np,mesh%K))
                allocate(rhs%v2(mesh%Np,mesh%K))
                allocate(rhs%p2(mesh%Np,mesh%K))
                allocate(res%v2(mesh%Np,mesh%K))
                allocate(res%p2(mesh%Np,mesh%K))
                allocate(src_vec%v2(mesh%Np, mesh%K))
            endif
        endif
        allocate(Zimp(mesh%Np,mesh%K))
        allocate(Vdm(mesh%np,mesh%n+1))
        allocate(invV(mesh%np,mesh%n+1))
    end subroutine

    subroutine deallocSolverArrays()
        if (allocated(r)) deallocate(r)
        if (allocated(solution%v)) deallocate(solution%v)
        if (allocated(solution%sigma)) deallocate(solution%sigma)
        if (allocated(rhs%v)) deallocate(rhs%v)
        if (allocated(rhs%sigma)) deallocate(rhs%sigma)
        if (allocated(res%v)) deallocate(res%v)
        if (allocated(res%sigma)) deallocate(res%sigma)
        if (allocated(solution%v1)) deallocate(solution%v1)
        if (allocated(solution%p1)) deallocate(solution%p1)
        if (allocated(solution%v2)) deallocate(solution%v2)
        if (allocated(solution%p2)) deallocate(solution%p2)
        if (allocated(rhs%v1)) deallocate(rhs%v1)
        if (allocated(rhs%p1)) deallocate(rhs%p1)
        if (allocated(rhs%v2)) deallocate(rhs%v2)
        if (allocated(rhs%p2)) deallocate(rhs%p2)
        if (allocated(res%v1)) deallocate(res%v1)
        if (allocated(res%p1)) deallocate(res%p1)
        if (allocated(res%v2)) deallocate(res%v2)
        if (allocated(res%p2)) deallocate(res%p2)
        if (allocated(Zimp)) deallocate(Zimp)
        if (allocated(Vdm)) deallocate(Vdm)
        if (allocated(invV)) deallocate(invV)
        if (allocated(rk1%v)) deallocate(rk1%v)
        if (allocated(rk2%v)) deallocate(rk2%v)
        if (allocated(rk1%sigma)) deallocate(rk1%sigma)
        if (allocated(rk2%sigma)) deallocate(rk2%sigma)
        if (allocated(timevec)) deallocate(timevec)
    end subroutine
end program solverPro
