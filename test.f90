program test
   implicit none
   double complex :: comp
   double complex :: conjcomp

   comp = dcmplx(1.0d0,1.0d0)
   conjcomp = conjg(comp)

   write(0,*)'comp=',comp
   write(0,*)'conjcomp=',conjcomp

end program test
