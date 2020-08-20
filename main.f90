program main
implicit none

call readInput
call compute

end program


subroutine compute
!use para
implicit none
integer :: covMatResult

    call readCovMat(covMatResult)
    if ( covMatResult .lt. 0 ) call computeCovariance
    call computePC
    call projectPC

end subroutine


subroutine projectPC
use para
implicit none
integer :: pc, frame

    write(*,'(A,I3,A)') ' Project coordinates to first ', iPC, ' PCs. '
    allocate( PCA(Nframe, iPC) )
    
    ! project to PCs
    do pc = 1, iPC
        !! `matmul` is too slow
        PCA(1:Nframe, pc) = matmul(PCvector(1:Natom3,pc), coord1D(1:Natom3,1:Nframe))
        !call matrixMul(1, Natom3, Nframe, &
        !               PCvector(1:Natom3,pc), coord1D(1:Natom3,1:Nframe), &
        !               PCA(1:Nframe, pc))
    end do
    
    ! write to file
    write(*,*) 'Writing PCA to file...'
    open(727, file=trim(fileNamePure)//'_PCA', status='replace', action='write')
    do frame = 1, Nframe
        write(727, '(999999999E18.8)') PCA(frame,1:iPC) 
    end do
    close(727)

end subroutine


subroutine computePC
use para
implicit none
real(inum) :: PCvalueSum
integer :: i

    write(*,*) 'Diagonalizing corvariance matrix... '
    allocate( PCvalue(Natom3), PCvector(Natom3,Natom3) )
    call calcEigenvector(Natom3, cov, PCvalue, PCvector)

    !PCvalueSum = dsqrt(sum(PCvalue*PCvalue)) ! normalized
    PCvalueSum = sum(PCvalue) ! sum = 1
    !PCvalue = PCvalue / PCvalueSum
    write(*,*) 'PC value (normalized, first 10): '
    if ( Natom3 .ge. 10 ) then
        write(*,'(5F12.6)') PCvalue(1:10) / PCvalueSum
    else
        write(*,'(5F12.6)') PCvalue(1:Natom3) / PCvalueSum
    end if
    !# write(*,*) 'Sum of PC value: '
    !# write(*,*) PCvalueSum
    !write(*,*) 'PC vector 1: '
    !write(*,'(3F12.6)') PCvector(1:Natom3,1)
    
    ! write PC value
    !! PC value -> value -> val -> 825
    write(*,*) 'Writing PC value...'
    open(825, file=trim(fileNamePure)//'_PCvalue', status='replace', action='write')
    write(825,*) 'PC values with percentages'
    write(825, *) Natom3
    do i = 1, Natom3
        write(825, '(E18.8, 2X, F12.4)') PCvalue(i), PCvalue(i)*100d0/PCvalueSum
    end do
    close(825)
    
    ! write PC vector
    !! PC vector -> vector -> vec -> 832
    write(*,*) 'Writing first 5 PC vectors...'
    open(832, file=trim(fileNamePure)//'_PCvector', status='replace', action='write')
    write(832, *) Natom3
    do i = 1, Natom3
        write(832, '(99999999E18.8)') PCvector(1:5,i)
    end do
    close(832)

end subroutine


subroutine readCovMat(covMatResult)
use para
implicit none
integer, intent(out) :: covMatResult
logical covMatStatus
integer i

    inquire(file=trim(fileNamePure)//'_covMat', exist=covMatStatus)
    if (covMatStatus) then
        write(*,*) 'Covariance matrix exists. Reading... '
        open(268, file=trim(fileNamePure)//'_covMat', status='old', action='read') ! cov mat -> cmt -> 268
        read(268, *); read(268, *) Natom3
        
        if ( mod(Natom3, 3) .ne. 0 ) then
            call logError('Number of atoms*3 is wrong! ')
        else
            Natom = Natom3 / 3
        end if
        
        allocate( cov(Natom3, Natom3) )
        do i = 1, Natom3
            read(268, *) cov(1:Natom3,i)
        end do
        close(268)
        
        covMatResult = 100
        return
    else
        covMatResult = -100
        return
    end if

end subroutine


subroutine computeCovariance
use para
use omp_lib
implicit none
integer :: component, i, j
!integer :: constructCount, constructTotal

    write(*,*); write(*,*) 'Principal component analysis! '
    allocate( coordMean(Natom3), cov(Natom3, Natom3) )
    
    ! calculate the mean value of each components
    write(*,*) 'Calculating the mean value of each components...'
    !$OMP PARALLEL
    !$OMP DO
    do component = 1, Natom3
        coordMean(component) = sum(coord1D(component,1:Nframe)) / Nframe
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !write(*,*) 'Mean value of each components: '
    !call logMat1D( coordMean(1:Natom3) )
    call writeXYZ(trim(fileNamePure)//'_mean', 1, reshape(coordMean, (/3, Natom/)) )
    
    ! substract mean value of each components
    do component = 1, Natom3
        coord1D(component,1:Nframe) = coord1D(component,1:Nframe) - coordMean(component)
    end do
    
    ! construct covariance matrix
    write(*,*) 'Constructing a covariance matrix...'

    !constructTotal = Natom3 * (Natom3 + 1) / 2
    !constructCount = 0
    !$OMP PARALLEL
    !$OMP DO
    do i = 1, Natom3
        do j = i, Natom3
            !!! `sum` is too slow! 
            !cov(i,j) = sum( coord1D(i,1:Nframe) * coord1D(j,1:Nframe) ) / Nframe
            !write(*,*) i,j,OMP_get_thread_num()
            call matrixMul(1, Nframe, 1, &
                     coord1D(i,1:Nframe), coord1D(j,1:Nframe), cov(i,j) )
            !constructCount = constructCount + 1
            !if( mod(constructCount, int(constructTotal/100)) .eq. 0) write(*,*) '1%',constructCount
            !if (i .ne. j) cov(j,i) = cov(i,j)
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL
    !$OMP DO
    do i = 1, Natom3
        do j = i, Natom3
            if (i .ne. j) cov(j,i) = cov(i,j)
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    !write(*,*) 'Writing the covariance matrix...'
    !open(368, file=trim(fileNamePure)//'_covMat', status='replace', action='write')
    !write(368,*) 'Covariance matrix (F18.8)'
    !write(368,*) Natom3
    !do i = 1, Natom3
    !    write(368,'(999999999E18.8)') cov(1:Natom3,i)
    !end do

end subroutine


subroutine writeXYZ(fileName, frames, coordinates)
use para
implicit none
character(len=*), intent(in) :: fileName
integer, intent(in) :: frames
real(inum), intent(in) :: coordinates(3,Natom,frames)
integer :: frame, atom

    open(9999, file=trim(fileName)//'.xyz', status='replace', action='write')
    do frame = 1, frames
        write(9999,*) Natom; write(9999,*);
        do atom = 1, Natom
            write(9999,'(A4,2X,3F12.6)') eleName(atom), coordinates(1:3,atom,frame)
        end do
    end do
    close(9999)

end subroutine


subroutine logMat1D(mat)
use para
real(inum), intent(in) :: mat(Natom3)
integer :: atom

    write(*,*) ' Ele.           X           Y           Z'
    do atom = 1, Natom
        write(*,'(A4,2X,3F12.6)') trim(eleName(atom)), mat( (atom-1)*3+1:atom*3 )
    end do

end subroutine


subroutine readXYZ
use para
implicit none
integer :: fileStatus
character(icha) :: tempLine
integer :: atom, frame

    write(*,*) 'Reading .xyz file...'
    open(999, file=trim(inputFileName), status='old', action='read')
    
    ! count frames
    Nframe = 0
    do while (.true.)
        read(999, '(I)', iostat=fileStatus) Natom
        if (fileStatus /= 0) exit
        read(999, *, iostat=fileStatus) 
        do atom = 1, Natom
            read(999, *, iostat=fileStatus) tempLine
        end do
        Nframe = Nframe + 1
    end do
    write(*,*) 'Atoms  : ', Natom
    write(*,*) 'Frames : ', Nframe
    rewind(999)
    Natom3 = Natom * 3
    
    ! allocate array
    allocate( coord(3,Natom,Nframe), eleName(Natom) )
    ! memory needed: 
    write(*,*) 'Memory of coordinates needed (MB): '
    write(*,*) 3d0*Natom*Nframe * 8d0 / 1048576d0 ! 8 bytes, convert to megabytes ( mb, 1024 * 1024 )
    
    ! read coordinates
    do frame = 1, Nframe
        read(999, *); read(999, *)
        do atom = 1, Natom
            read(999, *) eleName(atom), coord(1:3,atom,frame)
            !# write(*,*) coord(1:3,atom,frame)
        end do
    end do 
    close(999)
    
    ! convert 2d coordinates (xyz, Natoms) to 1d (xyz * Natoms) 
    ! check the converted coordinate matrix, please. Because I'm not sure the `reshape` could deal with this problem correctly. Now the `ifort` could give the right result. ! 2020-08-19 12:31:18 Wenbin, FAN @ SHU
    allocate( coord1D(Natom3, Nframe) )
    coord1D = reshape( coord, (/Natom3, Nframe/) )
    !#write(*,'(3F12.6)') coord1D(:,1)

end subroutine


subroutine readInput
use para
implicit none
character(icha) :: fileExt
logical inputFileStatus

    call getarg(1, inputFileName)
    inquire(file=trim(inputFileName), exist=inputFileStatus)
    if (inputFileStatus) then
        call getExtension(inputFileName, fileNamePure, fileExt)
        write(*,*) 'Input file : ', trim(inputFileName)
    else
        call logError("Your input file `" // trim(inputFileName) // "` doesn't exist! ")
    end if
    
    if ( trim(fileExt) .eq. 'xyz') then
        call readXYZ
    else
        call logError('File extension '//trim(fileExt)//' is not supported! ')
    end if

end subroutine