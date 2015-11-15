module clik_plik
  use clik
contains  
  subroutine plik_get_vecsize(clikid,size)
    type(clik_object), intent(in) :: clikid
    integer, intent(out) :: size
    
    call fortran_plik_get_vecsize(clikid,size)  
  end subroutine

  subroutine plik_allocate_vec(clikid,vec)
    type(clik_object), intent(in) :: clikid
    real(8),allocatable,dimension(:),intent(out)::vec
    integer::size

    call plik_get_vecsize(clikid,size)
    allocate(vec(size))
  end subroutine

  subroutine plik_get_fg(clikid,cl_and_pars,vec)
    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:),intent(in) :: cl_and_pars
    real(8),dimension(:),intent(out)::vec
      
    call fortran_plik_get_fg(clikid,cl_and_pars,vec)  
  end subroutine

  subroutine plik_get_cal_beam(clikid,cl_and_pars,vec)
    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:),intent(in) :: cl_and_pars
    real(8),dimension(:),intent(out)::vec
      
    call fortran_plik_get_cal_beam(clikid,cl_and_pars,vec)  
  end subroutine

  subroutine plik_get_data(clikid,vec)
    type(clik_object), intent(in) :: clikid
    real(8),dimension(:),intent(out)::vec
      
    call fortran_plik_get_data(clikid,vec)  
  end subroutine

end module 