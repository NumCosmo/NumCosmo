!============================================================================
MODULE inversecov_selection 
! Parameters are defined in Highell_options module
! ===========================================================================

  USE highell_options

  REAL(8) :: inverse_s(1:datap_s,1:datap_s),inverse_e(1:datap_e,1:datap_e)
  REAL(8) :: covmat_s(1:datap_s,1:datap_s),covmat_e(1:datap_e,1:datap_e)
  
  PRIVATE
  public :: inverse_init
  public :: get_free_lun

contains
  
  ! ============================================================================
  SUBROUTINE inverse_init
  ! ============================================================================
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,k,l,lun,n,stat
    REAL(8)  :: dummy
    CHARACTER(LEN=240) :: binfilename,invcovfilenames,invcovfilenamee
    LOGICAL  :: good
    REAL(8)  :: dum_mat_s(datap_s,datap_s),dum_mat_e(datap_e,datap_e)
    INTEGER  :: lmin(nspec),lmax(nspec)
    REAL(8)  :: blmean(1:nspec,0:tbin-1),blmin(1:nspec,0:tbin-1),blmax(1:nspec,0:tbin-1),bmin(nspec),bmax(nspec)

    lmin(1)=lmin11
    lmax(1)=lmax11
    
    lmin(2)=lmin12
    lmax(2)=lmax12
    
    lmin(3)=lmin22
    lmax(3)=lmax22

    
   !------------------------------------------------- 
   !Read inverse covariance matrixes 
   !-------------------------------------------------

    invcovfilenames= trim(ACT_data_dir)//'south/Inverse_Realistic_Cov_Mat_SouthMulti.dat'
    invcovfilenamee= trim(ACT_data_dir)//'equa/Inverse_Realistic_Cov_Mat_Equa.dat'
    binfilename= trim(ACT_data_dir)//'south/binningFile.dat'

    call get_free_lun( lun )
    open(unit=lun,file=invcovfilenames,form='formatted',status='unknown',action='read')
    do i=1,datap_s
       read(lun,*) inverse_s(i,1:datap_s)
    enddo
    close(lun)

    call get_free_lun( lun )
    open(unit=lun,file=invcovfilenamee,form='formatted',status='unknown',action='read')
    do i=1,datap_e
       read(lun,*) inverse_e(i,1:datap_e)
    enddo
    close(lun)

   !------------------------------------------------- 
   !Invert to get covariance matrixes 
   !-------------------------------------------------
   
   dum_mat_s(1:datap_s,1:datap_s) = inverse_s(1:datap_s, 1:datap_s)
   n = size(inverse_s,1)
   call dpotrf( 'L', n, dum_mat_s, n, stat )
   call dpotri( 'L', n, dum_mat_s, n, stat )
   do i=1,datap_s
      dum_mat_s(i,1:datap_s)=dum_mat_s(1:datap_s,i)
   enddo
   covmat_s(1:datap_s,1:datap_s)=dum_mat_s(1:datap_s,1:datap_s)

   dum_mat_e(1:datap_e,1:datap_e) = inverse_e(1:datap_e, 1:datap_e)
   n = size(inverse_e,1)
   call dpotrf( 'L', n, dum_mat_e, n, stat )
   call dpotri( 'L', n, dum_mat_e, n, stat )
   do i=1,datap_e
      dum_mat_e(i,1:datap_e)=dum_mat_e(1:datap_e,i)
   enddo
   covmat_e(1:datap_e,1:datap_e)=dum_mat_e(1:datap_e,1:datap_e)


   !------------------------------------------------- 
   !Get ell ranges
   !-------------------------------------------------

   do j=1,nspec
      call get_free_lun( lun )
      open(unit=lun,file=binfilename,form='formatted',status='unknown',action='read')
      do i = 0,tbin-1
         read(lun,*) blmin(j,i), blmax(j,i), blmean(j,i)
      end do
      close(lun)

      do i=0,tbin-1
         if(blmean(j,i) <= lmin(j)) bmin(j)=i+1
         if(blmean(j,i) <= lmax(j)) bmax(j)=i
      enddo
   enddo

   bmin(1) = bmin(1)-4
   bmin(2) = bmin(2)-14
   bmin(3) = bmin(3)-14

   bmax(1) = bmax(1)-4
   bmax(2) = bmax(2)-14
   bmax(3) = bmax(3)-14

   !------------------------------------------------- 
   !Change covmat to select the ell range
   !-------------------------------------------------

   !-------------------------------------------------
   !south

   do i=1,nsp11_s
      do j=1,nbin11
         k= nbin11*(i-1)+j
         if( k .ge. nbin11*(i-1) .and. k .lt. nbin11*(i-1)+bmin(1)+1) then
             covmat_s(k,1:datap_s) = 0
             covmat_s(1:datap_s,k) = 0
             covmat_s(k,k) = 1E+10
         end if
         if( k .gt. nbin11*(i-1)+bmax(1)+1 .and. k .le. nbin11*(i+1)) then
             covmat_s(k,1:datap_s) = 0
             covmat_s(1:datap_s,k) = 0
             covmat_s(k,k) = 1E+10
         end if
      enddo 
   enddo

   do i=1,nsp12_s
      do j=1,nbin12
         k= nbin11*nsp11_s+nbin12*(i-1)+j
         if( k .ge. nbin11*nsp11_s+nbin12*(i-1) .and. k .lt. nbin11*nsp11_s+nbin12*(i-1)+bmin(2)+1) then
             covmat_s(k,1:datap_s) = 0
             covmat_s(1:datap_s,k) = 0
             covmat_s(k,k) = 1E+10
         end if
         if( k .gt. nbin11*nsp11_s+nbin12*(i-1)+bmax(2)+1 .and. k .le. nbin11*nsp11_s+nbin12*(i+1)) then
             covmat_s(k,1:datap_s) = 0
             covmat_s(1:datap_s,k) = 0
             covmat_s(k,k) = 1E+10
         end if
      enddo
   enddo

   do i=1,nsp22_s
      do j=1,nbin22
         k= nbin11*nsp11_s+nbin12*nsp12_s+nbin22*(i-1)+j
         if( k .ge. nbin11*nsp11_s+nbin12*nsp12_s+nbin22*(i-1) .and. k .lt. nbin11*nsp11_s+nbin12*nsp12_s+nbin22*(i-1)+bmin(3)+1) then
             covmat_s(k,1:datap_s) = 0
             covmat_s(1:datap_s,k) = 0
             covmat_s(k,k) = 1E+10
         end if
         if( k .gt. nbin11*nsp11_s+nbin12*nsp12_s+nbin22*(i-1)+bmax(3)+1 .and. k .le. nbin11*nsp11_s+nbin12*nsp12_s+nbin22*(i+1)) then
             covmat_s(k,1:datap_s) = 0
             covmat_s(1:datap_s,k) = 0
             covmat_s(k,k) = 1E+10 
         end if
      enddo
   enddo


   !-------------------------------------------------
   !same for equa

   do i=1,nsp11_e
      do j=1,nbin11
         k= nbin11*(i-1)+j
         if( k .ge. nbin11*(i-1) .and. k .lt. nbin11*(i-1)+bmin(1)+1) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
         if( k .gt. nbin11*(i-1)+bmax(1)+1 .and. k .le. nbin11*(i+1)) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
      enddo
   enddo

   do i=1,nsp12_e
      do j=1,nbin12
         k= nbin11*nsp11_e+nbin12*(i-1)+j
         if( k .ge. nbin11*nsp11_e+nbin12*(i-1) .and. k .lt. nbin11*nsp11_e+nbin12*(i-1)+bmin(2)+1) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
         if( k .gt. nbin11*nsp11_e+nbin12*(i-1)+bmax(2)+1 .and. k .le. nbin11*nsp11_e+nbin12*(i+1)) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
      enddo
   enddo

   do i=1,nsp22_e
      do j=1,nbin22
         k= nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1)+j
         if( k .ge. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1) .and. k .lt. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1)+bmin(3)+1) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
         if( k .gt. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1)+bmax(3)+1 .and. k .le. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i+1)) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
      enddo
   enddo





   !------------------------------------------------- 
   !Get new inverse covmats 
   !-------------------------------------------------

   dum_mat_s(1:datap_s,1:datap_s) = covmat_s(1:datap_s,1:datap_s)         
   n = size(covmat_s,1)
   call dpotrf( 'L', n, dum_mat_s, n, stat )
   call dpotri( 'L', n, dum_mat_s, n, stat )
   do i=1,datap_s
      dum_mat_s(i,1:datap_s)=dum_mat_s(1:datap_s,i)
   enddo
   inverse_s(1:datap_s,1:datap_s)=dum_mat_s(1:datap_s,1:datap_s)
   call get_free_lun( lun )
   open(unit=lun, file = trim(ACT_data_dir)//'south/Inverse_South.dat', action='write', status='unknown')
        write(lun,1001) inverse_s(1:datap_s,1:datap_s)
1001  format (e15.8)
   close(lun)

   dum_mat_e(1:datap_e,1:datap_e) = covmat_e(1:datap_e,1:datap_e)      
   n = size(covmat_e,1)
   call dpotrf( 'L', n, dum_mat_e, n, stat )
   call dpotri( 'L', n, dum_mat_e, n, stat )
   do i=1,datap_e
      dum_mat_e(i,1:datap_e)=dum_mat_e(1:datap_e,i)
   enddo
   inverse_e(1:datap_e,1:datap_e)=dum_mat_e(1:datap_e,1:datap_e)
   call get_free_lun( lun )
   open(unit=lun, file = trim(ACT_data_dir)//'equa/Inverse_Equa.dat', action='write', status='unknown')
        write(lun,1002) inverse_e(1:datap_e,1:datap_e)
1002  format (e15.8)
   close(lun)

   END SUBROUTINE inverse_init
  !================================================================================


  !================================================================================
  subroutine get_free_lun( lun )

  implicit none
  integer, intent(out) :: lun
  integer, save :: last_lun = 19
  logical :: used
  lun = last_lun

  do
    inquire( unit=lun, opened=used )
    if ( .not. used ) exit
        lun = lun + 1
  end do

  last_lun = lun
  end subroutine
  !================================================================================


END MODULE inversecov_selection
