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
module outputMod
    use parameterMod
    use constantsMod
    use analyticalSolutionMod
    use slipInterfaceMod
    use genericMaterial
    use sourcesMod
    use errorMessage
    use matrixMod
    use receiverMod
    use calendar
    use convert_time

    implicit none

    type :: timestampVar
        real(kind=SIZE_DOUBLE) :: start
        real(kind=SIZE_DOUBLE) :: current
        integer :: year
        integer :: month
        integer :: day
        integer :: hours
        integer :: minutes
        integer :: timestamp
    end type

    contains

    subroutine plot1D(par, mesh, src, step, folder, name, y)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in) :: mesh
        type(sourceVar), dimension(:), intent(in) :: src
        character(len=*), intent(in) :: name, folder
        integer, intent(in) :: step
        real(kind=custom_real), dimension(:,:), intent(in) :: y
        !local variables
        integer :: i,j,l
        real(kind=custom_real) :: maxamplitude
        real, dimension(:), allocatable :: vec_x
        real, dimension(:), allocatable :: vec_y
        character(len=7) :: filenr


        !variables for the plot
        integer :: idev, pgopen! pgid,
        character(len=30) :: xlab, ylab, title, plotfile
        real :: xmin, xmax, ymin, ymax

        allocate(vec_x(mesh%Np*(mesh%K-2*mesh%nghost)))
        allocate(vec_y(mesh%Np*(mesh%K-2*mesh%nghost)))

        !Transformation der Matrizen x und u in vectoren, damit die daten gepolottet werden können.
        l = 1
        do i = 1+mesh%nghost, mesh%K-mesh%nghost
            do j = 1, mesh%Np
                vec_x(l) = real(mesh%x(j,i))
                vec_y(l) = real(y(j,i))
                l = l+1
            enddo
        enddo

        maxamplitude = 0.
        do i = 1, par%srcn
            if (src(i)%amplitude > maxamplitude) then
                maxamplitude = src(i)%amplitude
            endif
        enddo

        xmin = real(mesh%xmin)!
        xmax = real(mesh%xmax)!
        ymin = real(-1.1*(max(abs(maxval(y(:,1+mesh%nghost:mesh%K-mesh%nghost))),&
                              abs(minval(y(:,1+mesh%nghost:mesh%K-mesh%nghost))))))
        ymax = real( 1.1*(max(abs(maxval(y(:,1+mesh%nghost:mesh%K-mesh%nghost))),&
                              abs(minval(y(:,1+mesh%nghost:mesh%K-mesh%nghost))))))
        if (abs(ymin) < epsilon(ymin) .and. (abs(ymax) < epsilon(ymax))) then
            ymin = -1.1*maxamplitude
            ymax =  1.1*maxamplitude
        elseif (abs(ymin) < epsilon(ymin)) then
            ymin = -ymax
        elseif (abs(ymax) < epsilon(ymax)) then
            ymax = -ymin
        endif

        write(filenr,"(I7.7)") step

        plotfile = trim(name)//"_"//trim(filenr)//".png"

        !How to add units to the png-files??? TM TM
        xlab = "x"
        ylab = name
        title = trim(name)//" (Timestep"//trim(filenr)//")"

        !pgid = pgopen('/XSERVE')
        idev = pgopen("output/"//trim(folder)//"/"//trim(plotfile)//"/png")
        if (idev <= 0) stop "something is wrong"

        call pgslct(idev)
        call PGSCR(0, 1., 1., 1.)    !black becomes white
        call PGSCR(1, 0., 0., 0.)    !white becomes black
        call pgenv(xmin, xmax, ymin, ymax, 0,0)
        call pglab(xlab, ylab, title)
        call pgsci(2)                       !set color index
        call pgslw(5)                       !line width
        call pgsch(1.2)                     !set character height
        call pgscf(2)                       !set chracter font
        call pgline(mesh%Np*(mesh%K-2*mesh%nghost), vec_x, vec_y)
        call pgclos

        deallocate(vec_x)
        deallocate(vec_y)
    end subroutine

    subroutine writeInitialSnapshot(par, mesh, ana, genmat, src, solution, errmsg)
        !input
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type(analyticalVar) :: ana
        type(generic_material) :: genmat
        type(sourceVar), dimension(:) :: src
        type(error_message) :: errmsg
        type (solutionvector) :: solution
        !real(kind=custom_real), dimension(:,:), intent(in) :: V, Sigma
        !local
        character(len=14) :: filename
        character(len=17) :: anasol
        character(len=21) :: myname = 'writeInitialSnapshot'
        type (solutionvector) :: anaq

        call addTrace(errmsg,myname)

        !Initial Plot
        filename = "output/0000001"
        call writeSnapshot1d(solution, mesh, 1, filename)
        if (par%compare) then
            allocate(anaq%v(mesh%Np,mesh%K))
            allocate(anaq%sigma(mesh%Np,mesh%K))

            anaq%v     = ana%q(1,:,:)
            anaq%sigma = ana%q(2,:,:)

            anasol = "output/0000001ana"
            select case (src(1)%type)
                case ("box")
                    call initRightpropBoxcar(par, mesh, src(1), genmat, anaq)
                case ("sin2")
                    call initRightpropSineSquared(par, mesh, src(1), genmat, anaq)
                case ("test")
                    call add(errmsg,2,'Analytical solution for the test wavelet does not exist!',myname)
                    return
                case default
                    call add(errmsg,2,'No vaild initial wavelet for selected',myname)
                    return
            end select
            call writeSnapshotAna1d(anaq%v, mesh,2,anasol)
            call plot1D(par, mesh, src, 1, "velocity", "ana_velocity", anaq%v)

            if (allocated(anaq%v)) deallocate(anaq%v)
            if (allocated(anaq%sigma)) deallocate(anaq%sigma)
        endif
        call plot1D(par, mesh, src, 1, "velocity", "velocity", solution%v)
        call plot1D(par, mesh, src, 1, "stress", "stress", solution%sigma)
        if (par%poroelastic) then
            call plot1D(par, mesh, src, 1, "v1", "v1", solution%v1)
            call plot1D(par, mesh, src, 1, "p1", "p1", solution%p1)
            if (par%fluidn == 2) then
                call plot1D(par, mesh, src, 1, "v2", "v2", solution%v2)
                call plot1D(par, mesh, src, 1, "p2", "p2", solution%p2)
            endif
        endif
    end subroutine

    subroutine writeSnapshot1d(solution, mesh, lu, filename)
        !input
        type (solutionvector) :: solution
        type (meshVar), intent(in) :: mesh
        integer :: lu
        character (len=*) :: filename
        !local
        integer :: i,j
        !real(kind=custom_real) :: energy_v, energy_s  # MB MB #

        open(lu,file = filename)

        if (allocated(solution%p2) .and. allocated(solution%v2)) then
            !poroelastic medium with two fluids
            do i = 1+mesh%nghost, mesh%k-mesh%nghost
                do j = 1,mesh%np
                    write(lu,*) mesh%x(j,i), solution%v(j,i), solution%sigma(j,i),&
                        solution%v1(j,i), solution%p1(j,i), solution%v2(j,i), solution%p2(j,i)
                enddo
            enddo
        elseif (allocated(solution%p1) .and. allocated(solution%v1)) then
            !poroelastic medium with one fluid
            do i = 1+mesh%nghost, mesh%k-mesh%nghost
                do j = 1,mesh%np
                    write(lu,*) mesh%x(j,i), solution%v(j,i), solution%sigma(j,i),&
                        solution%v1(j,i), solution%p1(j,i)
                enddo
            enddo
        else
            !elastic medium
            do i = 1+mesh%nghost, mesh%k-mesh%nghost
                do j = 1,mesh%np
                    write(lu,*) mesh%x(j,i), solution%v(j,i), solution%sigma(j,i)
                enddo
            enddo
        endif
        close(lu)
    end subroutine writeSnapshot1d

    subroutine writeSnapshotAna1d(q, mesh, lu, filename)
        !input
        type (meshVar), intent(in) :: mesh
        integer :: lu
        real(kind=custom_real), dimension(:,:), intent(in) :: q
        character (len=*) :: filename
        !local
        integer :: i,j

        open(lu,file = filename)

        do i = 1+mesh%nghost, mesh%k-mesh%nghost
            do j = 1,mesh%np
                write(lu,*) mesh%x(j,i),q(j,i)
            enddo
        enddo
        close(lu)
    end subroutine writeSnapshotAna1d

    subroutine outputStamp(par, localtime, tstep, timestamp, solution, energy)
        !input
        type(parameterVar) :: par
        type(timestampVar) :: timestamp
        type(solutionVector) :: solution
        real(kind=custom_real), intent(in) :: localtime, energy
        integer :: tstep
        !Local
        character(len=80) :: filename
        integer :: lu!, i
        !Variables to count elapsed time
        real(kind=SIZE_DOUBLE) :: cpu
        real(kind=SIZE_DOUBLE) :: remaining
        real(kind=SIZE_DOUBLE) :: total
        integer :: hours
        integer :: minutes
        integer :: seconds
        integer :: remaining_hours
        integer :: remaining_minutes
        integer :: remaining_seconds
        integer :: total_hours
        integer :: total_minutes
        integer :: total_seconds
        integer :: int_remaining
        integer :: int_total
        integer :: int_cpu
        !Variables to determine date and time at which the run will finish
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer, dimension(8) :: time_values
        character(len=3), dimension(12) :: month_name
        character(len=3), dimension(0:6) :: weekday_name
        data month_name /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
        data weekday_name /'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/
        integer :: julian_day_number
        integer :: day_of_week

        !strings for the flushed output:
        !character(len=80), dimension(13) :: line
        !character(len=58) :: line1
        !character(len=42) :: line2
        !character(len=53) :: line3
        !character(len=80) :: line4
        !character(len=80) :: line5
        !character(len=80) :: line6
        !character(len=80) :: line7
        !character(len=80) :: line8
        !character(len=80) :: line9
        !character(len=80) :: line10
        !character(len=80) :: line11
        !character(len=80) :: line12
        !character(len=80) :: line13
        !character(len=*) :: update

        !Create timestamp
        call date_and_time(date,time,zone,time_values)
        timestamp%year    = time_values(1)
        timestamp%month   = time_values(2)
        timestamp%day     = time_values(3)
        timestamp%hours   = time_values(5)
        timestamp%minutes = time_values(6)
        call convtime(timestamp%timestamp,timestamp%year,timestamp%month,timestamp%day,timestamp%hours,timestamp%minutes)

        !convert the current time to seconds
        timestamp%current = timestamp%timestamp*dble(60) + time_values(7) + time_values(8)/dble(1000)

        !Calculate elapsed time
        cpu     = timestamp%current - timestamp%start
        int_cpu = int(cpu)
        hours   = int_cpu/3600
        minutes = (int_cpu - 3600*hours)/60
        seconds = int_cpu - 3600*hours - 60*minutes

        !Calculate estimated remaining simulation time
        remaining         = (par%tsteps - tstep)*(cpu/dble(tstep))
        int_remaining     = int(remaining)
        remaining_hours   = int_remaining/3600
        remaining_minutes = (int_remaining - 3600*remaining_hours)/60
        remaining_seconds = int_remaining - 3600*remaining_hours-60*remaining_minutes

        !Calculate estimated total simulation time
        total         = remaining + cpu
        int_total     = int(total)
        total_hours   = int_total/3600
        total_minutes = (int_total - 3600*total_hours)/60
        total_seconds = int_total - 3600*total_hours - 60*total_minutes

        if (par%logfile) then
            filename = "output/logfile"//trim(par%title)
            filename = trim(filename)
            lu = 16
            open(lu, file = filename, position="append")
            write(lu,"('|           Time step number ',i7,' (t = ',es10.3,' s) out of ',i7, '         |')") tstep,localtime,par%tsteps
            write(lu,"('|           Elapsed time:                 ',i4,' h ',i2.2,' m ',i2.2,' s', '                     |')") hours,minutes,seconds
            write(lu,"('|           Estimated remaining time:     ',i4,' h ',i2.2,' m ',i2.2,' s', '                     |')")&
                remaining_hours, remaining_minutes, remaining_seconds
            write(lu,"('|           Estimated total time:         ',i4,' h ',i2.2,' m ',i2.2,' s', '                     |')") &
                 total_hours,total_minutes,total_seconds
            write(lu,"('|           Mean elapsed time per timestep:    ',es9.2,' s', '                     |')") cpu/dble(tstep)

            write (lu,"(a44, es14.7, a22)") "|           Current kinetic energy density: ", energy, "                    |"
            !calculate maximum Norms according to sqrt(value^2) -> absolute of that value. Output is the maxvalue of the array
            write (lu,"(a44, es14.7, a22)") "|           Current maximum norm of V     : ", maxval(abs(solution%v)), "                    |"
            write (lu,"(a44, es14.7, a22)") "|           Current maximum norm of Sigma : ", maxval(abs(solution%sigma)), "                    |"
            if (par%poroelastic) then
                write (lu,"(a44, es14.7, a22)") "|           Current maximum norm of v1    : ", maxval(abs(solution%v1)), "                    |"
                write (lu,"(a44, es14.7, a22)") "|           Current maximum norm of p1    : ", maxval(abs(solution%p1)), "                    |"
                if (par%fluidn == 2) then
                    write (lu,"(a44, es14.7, a22)") "|           Current maximum norm of v2    : ", maxval(abs(solution%v2)), "                    |"
                    write (lu,"(a44, es14.7, a22)") "|           Current maximum norm of p2    : ", maxval(abs(solution%p2)), "                    |"
                endif
            endif
        endif
        if (par%log) then
            write(*,"('|           Time step number ',i7,' (t = ',es10.3,' s) out of ',i7, '         |')") tstep,localtime,par%tsteps
            write(*,"('|           Elapsed time:                 ',i4,' h ',i2.2,' m ',i2.2,' s', '                     |')") hours,minutes,seconds
            write(*,"('|           Estimated remaining time:     ',i4,' h ',i2.2,' m ',i2.2,' s', '                     |')")&
                remaining_hours, remaining_minutes, remaining_seconds
            write(*,"('|           Estimated total time:         ',i4,' h ',i2.2,' m ',i2.2,' s', '                     |')") &
                 total_hours,total_minutes,total_seconds
            write(*,"('|           Mean elapsed time per timestep:    ',es9.2,' s', '                     |')") cpu/dble(tstep)

            write (*,"(a44, es14.7, a22)") "|           Current kinetic energy density: ", energy, "                     |"
            !calculate maximum Norms according to sqrt(value^2) -> absolute of that value. Output is the maxvalue of the array
            write (*,"(a44, es14.7, a22)") "|           Current maximum norm of V     : ", maxval(abs(solution%v)), "                     |"
            write (*,"(a44, es14.7, a22)") "|           Current maximum norm of Sigma : ", maxval(abs(solution%sigma)), "                     |"
            if (par%poroelastic) then
                write (*,"(a44, es14.7, a22)") "|           Current maximum norm of v1    : ", maxval(abs(solution%v1)), "                     |"
                write (*,"(a44, es14.7, a22)") "|           Current maximum norm of p1    : ", maxval(abs(solution%p1)), "                     |"
                if (par%fluidn == 2) then
                    write (*,"(a44, es14.7, a22)") "|           Current maximum norm of v2     : ", maxval(abs(solution%v2)), "                     |"
                    write (*,"(a44, es14.7, a22)") "|           Current maximum norm of p2     : ", maxval(abs(solution%p2)), "                     |"
                endif
            endif
        endif

        if(tstep <= par%tsteps) then
            ! compute date and time at which the run should finish (useful for long runs)
            ! add remaining minutes and get date and time of that future timestamp in minutes
            timestamp%timestamp = int((timestamp%current + remaining) / dble(60))
            call invtime(timestamp%timestamp,timestamp%year,timestamp%month,timestamp%day,timestamp%hours,timestamp%minutes)

            ! convert to Julian day to get day of the week
            call calndr(timestamp%day,timestamp%month,timestamp%year,julian_day_number)
            day_of_week = idaywk(julian_day_number)
            if (par%logfile) then
                write(lu,"('|           This run will finish on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' at ',i2.2,':',i2.2, '                 |')") &
                weekday_name(day_of_week),month_name(timestamp%month),&
                timestamp%day,timestamp%year,timestamp%hours,timestamp%minutes
            endif
            if (par%log) then
                write(*,"('|           This run will finish on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' at ',i2.2,':',i2.2, '                 |')") &
                weekday_name(day_of_week),month_name(timestamp%month),&
                timestamp%day,timestamp%year,timestamp%hours,timestamp%minutes
            endif
        endif
        if (par%logfile) then
            write(lu,'(a24, i3, a53)') "|                       ", 100*tstep/par%tsteps , "% of the simulation completed.                     |"
            write(lu,"(a80)") "|------------------------------------------------------------------------------|"
            close(lu)
        endif
        if (par%log) then
            write(*,'(a24, i3, a53)') "|                       ", 100*tstep/par%tsteps , "% of the simulation completed.                     |"
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
        endif
    end subroutine

    subroutine initialTimestamp(timestamp)
        !in/out
        type(timestampVar) :: timestamp
        !local
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer, dimension(8) :: time_values

        call date_and_time(date,time,zone,time_values)
        timestamp%year    = time_values(1)
        timestamp%month   = time_values(2)
        timestamp%day     = time_values(3)
        timestamp%hours   = time_values(5)
        timestamp%minutes = time_values(6)
        call convtime(timestamp%timestamp,timestamp%year,timestamp%month,timestamp%day,timestamp%hours,timestamp%minutes)

        !convert the current time to seconds
        timestamp%start = timestamp%timestamp*dble(60) + time_values(7) + time_values(8)/dble(1000)
    end subroutine

    subroutine plotSeismo1D(par, src, rec, timevec)
        !input
        type(parameterVar) :: par
        type(sourceVar), dimension(:) :: src
        type(receiverVar), dimension(:) :: rec
        real(kind=custom_real), dimension(:) :: timevec
        !local
        integer :: r, tstep
        real(kind=custom_real) :: shift                 !<insert comment>
        character(len=50) :: filename
        character(len=3) :: recno

        !This works only for 1 source!!! TM
        if (src(1)%type == "sin2") then
            shift = src(1)%width
        else
            shift = 0!src(1)%t0
        end if

        do r = 1, par%recn
            write(recno, "(i3.3)") r
            filename = "output/"//trim(par%title)//"_seismogram_v_rec"//trim(recno)
            open(unit=19, file = trim(filename))
            do tstep = 1, par%tsteps
                write(19,*) timevec(tstep)-shift, rec(r)%data(tstep)
            enddo
            close(19)
            filename = "output/"//trim(par%title)//"_seismogram_sigma_rec"//trim(recno)
            open(unit=19, file = trim(filename))
            do tstep = 1, par%tsteps
                write(19,*) timevec(tstep)-shift, rec(r)%dataSi(tstep)
            enddo
            close(19)
            if (par%poroelastic) then
                filename = "output/"//trim(par%title)//"_seismogram_v1_rec"//trim(recno)
                open(unit=19, file = trim(filename))
                do tstep = 1, par%tsteps
                    write(19,*) timevec(tstep)-shift, rec(r)%dataV1(tstep)
                enddo
                close(19)
                filename = "output/"//trim(par%title)//"_seismogram_p1_rec"//trim(recno)
                open(unit=19, file = trim(filename))
                do tstep = 1, par%tsteps
                    write(19,*) timevec(tstep)-shift, rec(r)%dataP1(tstep)
                enddo
                close(19)
                if (par%fluidn == 2) then
                    filename = "output/"//trim(par%title)//"_seismogram_v2_rec"//trim(recno)
                    open(unit=19, file = trim(filename))
                    do tstep = 1, par%tsteps
                        write(19,*) timevec(tstep)-shift, rec(r)%dataV2(tstep)
                    enddo
                    close(19)
                    filename = "output/"//trim(par%title)//"_seismogram_p2_rec"//trim(recno)
                    open(unit=19, file = trim(filename))
                    do tstep = 1, par%tsteps
                        write(19,*) timevec(tstep)-shift, rec(r)%dataP2(tstep)
                    enddo
                    close(19)
                endif
            endif
        enddo
    end subroutine

    subroutine writeEnergyPoroelastic(solution, mesh, poromat, time, energy, lu)
        !input
        type (solutionvector) :: solution
        type (porous_material) :: poromat
        type (meshVar), intent(in) :: mesh
        integer :: lu
        real(kind=custom_real) :: time
        !output
        real(kind=custom_real), intent(out) :: energy
        !local
        integer :: i,j
        logical :: isthere

        inquire(file = "output/energy", exist=isthere)
        if (isthere) then
            open(lu, file = "output/energy", status = "old", position = "append")
        else
            open(lu, file = "output/energy", status = "new")
        endif

        energy = 0.
        if (allocated(solution%p2) .and. allocated(solution%v2)) then
            !poroelastic medium with two fluids
            do i = 1+mesh%nghost, mesh%k-mesh%nghost
                do j = 1,mesh%np
                    energy = energy &
                     + .5 * (1-poromat%phi(mesh%ntom(j,i))) * &
                     poromat%rhos(mesh%ntom(j,i)) * solution%v(j,i)**2 &
                     + .5 * poromat%phi(mesh%ntom(j,i)) * poromat%S1(mesh%ntom(j,i)) * &
                     poromat%rho1(mesh%ntom(j,i)) * solution%v1(j,i)**2 &
                     + .5 * poromat%phi(mesh%ntom(j,i)) * poromat%S2(mesh%ntom(j,i)) * &
                     poromat%rho2(mesh%ntom(j,i)) * solution%v2(j,i)**2
                enddo
            enddo
        elseif (allocated(solution%p1) .and. allocated(solution%v1)) then
            !poroelastic medium with one fluid
            do i = 1+mesh%nghost, mesh%k-mesh%nghost
                do j = 1,mesh%np
                    energy = energy &
                     + .5 * (1-poromat%phi(mesh%ntom(j,i))) * &
                     poromat%rhos(mesh%ntom(j,i)) * solution%v(j,i)**2 &
                     + .5 * poromat%phi(mesh%ntom(j,i)) * &
                     poromat%rho1(mesh%ntom(j,i)) * solution%v1(j,i)**2
                enddo
            enddo
        else
            !elastic medium
            do i = 1+mesh%nghost, mesh%k-mesh%nghost
                do j = 1,mesh%np
                    energy = energy &
                     + .5 * poromat%rhos(mesh%ntom(j,i)) * solution%v(j,i)**2
                enddo
            enddo
        endif
        energy = energy / (mesh%k-2*mesh%nghost)*mesh%np
        write(lu,"(es15.8, a2, es15.8)") time, "  ", energy
        close(lu)
    end subroutine

    subroutine writeEnergyElastic(solution, mesh, genmat, time, energy, lu)
        !input
        type (solutionvector) :: solution
        type (generic_material) :: genmat
        type (meshVar), intent(in) :: mesh
        integer :: lu
        real(kind=custom_real) :: time
        !output
        real(kind=custom_real), intent(out) :: energy
        !local
        integer :: i,j
        logical :: isthere

        real(kind=custom_real) :: rho

        inquire(file = "output/energy", exist=isthere)
        if (isthere) then
            open(lu, file = "output/energy", status = "old", position = "append")
        else
            open(lu, file = "output/energy", status = "new")
        endif

        energy = 0.
        !elastic medium
        do i = 1+mesh%nghost, mesh%k-mesh%nghost
            do j = 1,mesh%np
                rho = genmat.rho.(mesh%ntom(j,i))
                energy = energy &
                 + .5 * rho * solution%v(j,i)**2
            enddo
        enddo
        energy = energy / (mesh%k-2*mesh%nghost)*mesh%np
        write(lu,"(es15.8, a2, es15.8)") time, "  ", energy
        close(lu)
    end subroutine
end module
