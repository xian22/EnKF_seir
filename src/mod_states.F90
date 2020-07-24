module mod_states
! Modelstate definition for ECLIPSE
   use mod_dimensions
 !  integer :: neq=3*na+9

   type states
! Solution data
      real S(na)        ! Susceptible age groups
      real E(na)        ! Exposed age groups
      real I(na)        ! Infectious age groups
      real Qm           ! Quarantened mild sickness
      real Qs           ! Quarantened severe sickness (will be hospitalized)
      real QfH          ! Quarantened fatal sickness (will die) H
      real QfR          ! Quarantened fatal sickness (will die) R
      real Hs           ! Hospitalized with severe sickness
      real HfH          ! Hospitalized with fatal sickness H
      real HfR          ! Hospitalized with fatal sickness R
      real CH           ! In Care home with fatal sickness H
      real CR           ! In Care home with fatal sickness R
      real Rm           ! Recovered from mild sickness
      real Rs           ! Recovered from severe sickness (from hospital)
      real DH           ! Dead H
      real DR           ! Dead R
   end type states


! Overloaded and generic operators
   interface operator(+)
      module procedure add_states
   end interface

   interface operator(-)
      module procedure subtract_states
   end interface

   interface operator(*)
      module procedure states_real_mult,&
                       real_states_mult,&
                       states_states_mult
   end interface

   interface assignment(=)
      module procedure assign_states
   end interface

   interface sqrt
      module procedure sqrt_states
   end interface

   interface sum
      module procedure sum_states
   end interface


contains

   function sum_states(A)
      real sum_states
      type(states), intent(in) :: A
      sum_states = sum(A%S) &
                 + sum(A%E) &
                 + sum(A%I) &
                 + A%Qm     &
                 + A%Qs     &
                 + A%QfH    &
                 + A%QfR    &
                 + A%Hs     &
                 + A%HfH    &
                 + A%HfR    &
                 + A%CH     &
                 + A%CR     &
                 + A%Rm     &
                 + A%Rs     &
                 + A%DH     &
                 + A%DR 
   end function sum_states

   function sqrt_states(A)
      type(states) sqrt_states
      type(states), intent(in) :: A
      real :: eps=0.1E-14
      sqrt_states%S       = sqrt(A%S+eps)
      sqrt_states%E       = sqrt(A%E+eps)
      sqrt_states%I       = sqrt(A%I+eps)
      sqrt_states%Qm      = sqrt(A%Qm+eps)
      sqrt_states%Qs      = sqrt(A%Qs+eps)
      sqrt_states%QfH     = sqrt(A%QfH+eps)
      sqrt_states%QfR     = sqrt(A%QfR+eps)
      sqrt_states%Hs      = sqrt(A%Hs+eps)
      sqrt_states%HfH     = sqrt(A%HfH+eps)
      sqrt_states%HfR     = sqrt(A%HfR+eps)
      sqrt_states%CH      = sqrt(A%CH+eps)
      sqrt_states%CR      = sqrt(A%CR+eps)
      sqrt_states%Rm      = sqrt(A%Rm+eps)
      sqrt_states%Rs      = sqrt(A%Rs+eps)
      sqrt_states%DH      = sqrt(A%DH+eps)
      sqrt_states%DR      = sqrt(A%DR+eps)
   end function sqrt_states

   function add_states(A,B)
      type(states) add_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
      add_states%S       = A%S  + B%S 
      add_states%E       = A%E  + B%E
      add_states%I       = A%I  + B%I
      add_states%Qm      = A%Qm + B%Qm
      add_states%Qs      = A%Qs + B%Qs
      add_states%QfH     = A%QfH + B%QfH
      add_states%QfR     = A%QfR + B%QfR
      add_states%Hs      = A%Hs + B%Hs
      add_states%HfH     = A%HfH + B%HfH
      add_states%HfR     = A%HfR + B%HfR
      add_states%CH      = A%CH  + B%CH 
      add_states%CR      = A%CR  + B%CR 
      add_states%Rm      = A%Rm + B%Rm
      add_states%Rs      = A%Rs + B%Rs
      add_states%DH      = A%DH  + B%DH
      add_states%DR      = A%DR  + B%DR
   end function add_states

   function subtract_states(A,B)
      type(states) subtract_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
      subtract_states%S       = A%S  - B%S 
      subtract_states%E       = A%E  - B%E
      subtract_states%I       = A%I  - B%I
      subtract_states%Qm      = A%Qm - B%Qm
      subtract_states%Qs      = A%Qs - B%Qs
      subtract_states%QfH     = A%QfH - B%QfH
      subtract_states%QfR     = A%QfR - B%QfR
      subtract_states%Hs      = A%Hs - B%Hs
      subtract_states%HfH     = A%HfH - B%HfH
      subtract_states%HfR     = A%HfR - B%HfR
      subtract_states%CH      = A%CH  - B%CH
      subtract_states%CR      = A%CR  - B%CR 
      subtract_states%Rm      = A%Rm - B%Rm
      subtract_states%Rs      = A%Rs - B%Rs
      subtract_states%DH      = A%DH  - B%DH
      subtract_states%DR      = A%DR  - B%DR
   end function subtract_states

   function states_real_mult(A,B)
      type(states) states_real_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
      states_real_mult%S       = B*A%S 
      states_real_mult%E       = B*A%E 
      states_real_mult%I       = B*A%I 
      states_real_mult%Qm      = B*A%Qm
      states_real_mult%Qs      = B*A%Qs
      states_real_mult%QfH     = B*A%QfH
      states_real_mult%QfR     = B*A%QfR
      states_real_mult%Hs      = B*A%Hs
      states_real_mult%HfH     = B*A%HfH
      states_real_mult%HfR     = B*A%HfR
      states_real_mult%CH      = B*A%CH
      states_real_mult%CR      = B*A%CR
      states_real_mult%Rm      = B*A%Rm
      states_real_mult%Rs      = B*A%Rs
      states_real_mult%DH      = B*A%DH
      states_real_mult%DR      = B*A%DR 
   end function states_real_mult

   function real_states_mult(B,A)
      type(states) real_states_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
      real_states_mult%S       = B*A%S 
      real_states_mult%E       = B*A%E 
      real_states_mult%I       = B*A%I 
      real_states_mult%Qm      = B*A%Qm
      real_states_mult%Qs      = B*A%Qs
      real_states_mult%QfH     = B*A%QfH
      real_states_mult%QfR     = B*A%QfR
      real_states_mult%Hs      = B*A%Hs
      real_states_mult%HfH     = B*A%HfH
      real_states_mult%HfR     = B*A%HfR
      real_states_mult%CH      = B*A%CH
      real_states_mult%CR      = B*A%CR
      real_states_mult%Rm      = B*A%Rm
      real_states_mult%Rs      = B*A%Rs
      real_states_mult%DH      = B*A%DH
      real_states_mult%DR      = B*A%DR 
   end function real_states_mult

   function states_states_mult(A,B)
      type(states) states_states_mult
      type(states), intent(in) :: A
      type(states), intent(in) :: B
      states_states_mult%S       = A%S  * B%S 
      states_states_mult%E       = A%E  * B%E 
      states_states_mult%I       = A%I  * B%I 
      states_states_mult%Qm      = A%Qm * B%Qm
      states_states_mult%Qs      = A%Qs * B%Qs
      states_states_mult%QfH     = A%QfH * B%QfH
      states_states_mult%QfR     = A%QfR * B%QfR
      states_states_mult%Hs      = A%Hs * B%Hs
      states_states_mult%HfH     = A%HfH * B%HfH
      states_states_mult%HfR     = A%HfR * B%HfR
      states_states_mult%CH      = A%CH  * B%CH
      states_states_mult%CR      = A%CR  * B%CR
      states_states_mult%Rm      = A%Rm * B%Rm
      states_states_mult%Rs      = A%Rs * B%Rs
      states_states_mult%DH      = A%DH  * B%DH
      states_states_mult%DR      = A%DR  * B%DR 
   end function states_states_mult


   subroutine assign_states(A,r)
      type(states), intent(out) :: A
      real, intent(in) :: r
      A%S       = r
      A%E       = r
      A%I       = r
      A%Qm      = r
      A%Qs      = r
      A%QfH     = r
      A%QfR     = r
      A%Hs      = r
      A%HfH     = r
      A%HfR     = r
      A%CH      = r
      A%CR      = r
      A%Rm      = r
      A%Rs      = r
      A%DH      = r
      A%DR      = r
   end subroutine assign_states


end module mod_states

