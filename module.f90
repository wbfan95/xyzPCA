module constant
implicit none

! system precision
integer, parameter :: icha = 1024 ! character length
integer, parameter :: inum = 8 ! precision of real value

! constants
real(inum), parameter :: pi = dacos(-1d0)


end module


module para
use constant
implicit none

! file operation
character(icha) :: inputFileName
character(icha) :: fileNamePure

! trajectory
integer :: Natom, Nframe, Natom3 ! Natom3 = 3*Natom
real(inum), allocatable :: coord(:,:,:) ! 3 (xyz) * Natom * Nframe
real(inum), allocatable :: coord1D(:,:) ! Natom3 * Nframe
character(icha), allocatable :: eleName(:) ! Natom

! covariance matrix
real(inum), allocatable :: coordMean(:) ! Natom3
real(inum), allocatable :: cov(:,:) ! Natom3 * Natom3
real(inum), allocatable :: PCvalue(:), PCvector(:,:)

! project to PC
integer, parameter :: iPC = 3
real(inum), allocatable :: PCA(:,:)


end module


! get extension from a file name
! source : https://stackoverflow.com/a/36736660
subroutine getExtension(fileName, fileNameReal, extension)
use constant
implicit none
character(len=*), intent(in) :: fileName
character(icha), intent(out) :: fileNameReal, extension
integer :: ppos, fileNameLength

    fileNameLength = len(trim(fileName))
    ppos = scan(trim(fileName), '.', BACK=.true.)
    if ( ppos > 0 ) then
        fileNameReal = fileName(1:ppos-1)
        extension = fileName(ppos+1:fileNameLength)
    end if

end subroutine


! calculate eigenvalue and eigenvector
! N : size of matrix
! matIn : matrix(N,N)
! evReal : real eigenvalue(N)
! ev : eigenvector(N,N)
subroutine calcEigenvector(N, matIn, evReal, ev)
use lapack95, only: geev
integer, intent(in) :: N
real(8), intent(in) :: matIn(N,N)
real(8), intent(out) :: evReal(N)
real(8) :: evImag(N)
real(8) :: evl(N,N), evr(N,N)
real(8) :: work(N*N)
real(8), intent(out) :: ev(N,N)
integer :: info, i

    call dgeev('V','V',N,matIn,N,evReal,evImag,evl,N,evr,N,work,N*4,info)
    ev(1:N,1:N) = evl(1:N,1:N)

end subroutine


! matrix-matrix product
subroutine matrixMul(m, k, n, A, B, C)
!use lapack95, only: dgemm
!use blas95, only: gemm
use constant
implicit none
integer, intent(in) :: m, n, k
real(inum), intent(in) :: A(m,k), B(k,n)
real(inum), intent(inout) :: C(m,n)
integer :: maxSize

    maxSize = maxval( (/m,n,k/) )
    !write(*,*) maxSize
    call dgemm('N', 'N', m, n, k, 1d0, A, M, B, K, 0d0, C, M)
    !call gemm(A, B, C)

end subroutine

! log error and exit whole program
subroutine logError(content)
implicit none
character(len=*), intent(in) :: content

    print *, '[ERROR] ', content
    stop

end subroutine

subroutine logWarn(content)
implicit none
character(len=*), intent(in) :: content

    print *, '[WARNING] ', content

end subroutine