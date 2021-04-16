program main
    implicit none

    integer, parameter :: nMAX = 16000
    
    ! b = A * x
    double precision, dimension (nMAX,nMAX) :: A
    double precision, dimension (nMAX) :: x
    double precision, dimension (nMAX) :: b

    double precision :: s,t, start_ij,end_ij, start_ji,end_ji

    integer :: i,j,n,k

    call RANDOM_SEED()    
    call random_number(A)
    call random_number(x)

    open(unit=1, file="fortran_times.csv", status="unknown", action="write")
    write(1,"(a)") "n,Time_ij,Time_ji"
    
    do k=0, 16
        n = k * 1000

        ! mult IJ
        call cpu_time(start_ij)
        do i=1, n
            s = 0
            do j=1, n
                ! b(i) = b(i) + A(i,j) * x(j)
                s = s + A(i,j) * x(j)
            end do
            b(i) = s
        end do
        call cpu_time(end_ij)

        b(:n) = 0 ! reset array
        
        ! mult JI
        call cpu_time(start_ji)
        do j=1, n
            t = x(j)
            do i=1, n
                ! b(i) = b(i) + A(i,j) * x(j)
                b(i) = b(i) + A(i,j) * t
            end do
        end do
        call cpu_time(end_ji)
        
        write(*,*) n,(end_ij-start_ij),(end_ji-start_ji) ! write on console
        write(1,*) n,",",(end_ij-start_ij),",",(end_ji-start_ji) ! write on unit 1 (csv file)
    end do

    close(unit=1)
end program main
