      subroutine elemsubr ( npelm, x, y, nunk_pel, elem_mat,  &
                            elem_vec, elem_mass, prevsolution, itype )
!
!                       INPUT / OUTPUT PARAMETERS
!
      implicit none
      integer, intent(in) :: npelm, nunk_pel, itype
      double precision, intent(in) :: x(1:npelm), y(1:npelm),  &
                                      prevsolution(1:nunk_pel)
      double precision, intent(out) :: elem_mat(1:nunk_pel,1:nunk_pel),  &
                                       elem_vec(1:nunk_pel),  &
                                       elem_mass(1:nunk_pel)

! **********************************************************************
!
!                       LOCAL PARAMETERS
!
! ======================================================================
!

!
double precision :: beta(1:3), gamma(1:3), delta, h, s,e
integer :: i, j
! ======================================================================

s= 1d0/3*sqrt(2d0)
e=1d-6

if ( itype==1 ) then
!
!
!
!--- Type = 1: internal element
!Compute the factors beta, gamma and delta as defined in the NMSC
delta = (x(2)-x(1))*(y(3)-y(1))-(y(2)-y(1))*(x(3)-x(1))
beta(1) = (y(2)-y(3))/delta
beta(2) = (y(3)-y(1))/delta
beta(3) = (y(1)-y(2))/delta
gamma(1) = (x(3)-x(2))/delta
gamma(2) = (x(1)-x(3))/delta
gamma(3) = (x(2)-x(1))/delta
!--- Fill the element matrix as defined in the Lecture Notes
do j = 1, 3
do i = 1, 3
elem_mat(i,j) = 0.5d0 * abs(delta) * e*( beta(i)*beta(j) + gamma(i)*gamma(j) ) + &
0.5d0 * abs(delta) * s 	* (beta(j)+gamma(j))


end do
end do
elem_vec = 0d0

!
!



!--- The element vector is zero

else
!--- Type = 2: boundary element
!Compute Jacobian h

!h = sqrt ( (x(2)-x(1))**2 + (y(2)-y(1))**2 )
!
! --- The element matrix is zero
elem_mat = 0d0
do j= 1 , 3
elem_vec(j)= elem_vec(j)+  0.5d0 * abs(delta) *e*( beta(i)*beta(j) + gamma(i)*gamma(j) ) + &
0.5d0 * abs(delta) * s 	* (beta(j)+gamma(j))

!
 !--- Fill the element vector
!elem_vec(1:2) = h * 0.5d0 * x(1:2)
end do

end if
end subroutine elemsubr
