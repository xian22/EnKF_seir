module m_enkfprep
contains
subroutine enkfprep(ens,enspar)
   use mod_dimensions
   use mod_states
   use mod_params
   use mod_parameters
   use m_enkfini
   use m_agegroups
   use m_pseudo1D
   use m_fixsample1D
   use m_readinputs
   implicit none
   type(states), intent(in) :: ens(0:nt,nrens)
   type(params), intent(in) :: enspar(nrens)
   type(states) aveens
   type(params)  avepar
   real, allocatable :: obspertd(:,:)
   real, allocatable :: obsperth(:,:)
   real, allocatable :: obspertc(:,:)
   real, allocatable :: scaling(:)
   integer i,m,j
   logical, save :: lprt =.false.
 
   print *
   if (lprt) then
      print '(a)','Prior ensemble mean of parameters:'
      avepar=0.0
      do j=1,nrens
         avepar=avepar+enspar(j)*(1.0/real(nrens))
      enddo
      print '(100a10)',parnames
      print '(100f10.3)',avepar
      print *
      lprt=.false.
   endif

   print '(a)','Preparing for analysis computation'
   if (lmeascorr) then
      allocate(obspertd(0:nt,nrens))
      allocate(obsperth(0:nt,nrens))
      allocate(obspertc(0:nt,nrens))
      print '(2(a,i5))','nt+1=',nt+1, ' n=',nint(real(nt+1)*rh)
      call pseudo1D(obspertd,nt+1,nrens,rh,1.0,nint(real(nt+1)*rh))
      call pseudo1D(obsperth,nt+1,nrens,rh,1.0,nint(real(nt+1)*rh))
      call pseudo1D(obspertc,nt+1,nrens,rh,1.0,nint(real(nt+1)*rh))
      call fixsample1D(obspertd,nt+1,nrens)
      call fixsample1D(obsperth,nt+1,nrens)
      call fixsample1D(obspertc,nt+1,nrens)
   else
      call random(E,nrobs*nrens)
      call fixsample1D(E,nrobs,nrens)
   endif

   R=0.0
   do m=1,nrobs
      select case (cobs(m))
 !     case('dh')
 !       R(m,m)=real(nesmda)*min(maxerrdh,max(relerrdh*dobs(m),minerrdh))**2
 !       if (lmeascorr) E(m,:)=sqrt(R(m,m))*obspertd(iobs(m),:)
      case('d')
         R(m,m)=real(nesmda)*min(maxerrd,max(relerrd*dobs(m),minerrd))**2
         if (lmeascorr) E(m,:)=sqrt(R(m,m))*obspertd(iobs(m),:)
      case('h')
         R(m,m)=real(nesmda)*min(maxerrh,max(relerrh*dobs(m),minerrh))**2
         if (lmeascorr) E(m,:)=sqrt(R(m,m))*obsperth(iobs(m),:)
      case('c')
         R(m,m)=real(nesmda)*min(maxerrc,max(relerrc*dobs(m),minerrc))**2
         if (lmeascorr) E(m,:)=sqrt(R(m,m))*obspertc(iobs(m),:)
      case default
         stop 'Measurement type not found'
      end select
      if (.not.lmeascorr) E(m,:)=sqrt(R(m,m))*E(m,:)
      D(m,:)=dobs(m)+E(m,:)
   enddo

   if (lmeascorr) then
      R=matmul(E,transpose(E))/real(nrens)
      deallocate(obspertd,obsperth)
   endif


   do m=1,nrobs
! ensemble average of state for observation m
      aveens=0.0
      do j=1,nrens
         aveens= aveens + ens(iobs(m),j)
      enddo
      aveens=aveens*(1.0/real(nrens))

! ensemble average of state for observation m
      select case (cobs(m))
!      case('dh')
!        DH(m,:) = DH(m,:)-N*ens(iobs(m),:)%DH
!        S(m,:) = N*( ens(iobs(m),:)%DH - aveens%DH )
     case('d')
         D(m,:) = D(m,:)-N*ens(iobs(m),:)%DR - N*ens(iobs(m),:)%DH
         S(m,:) = N*( ens(iobs(m),:)%DH+ens(iobs(m),:)%DR - aveens%DR -aveens%DH )
      case('h')
!        DH(m,:) = DH(m,:)-N*(ens(iobs(m),:)%Hs + ens(iobs(m),:)%HfH)
         D(m,:) = D(m,:)-N*(ens(iobs(m),:)%Hs + ens(iobs(m),:)%HfH+ ens(iobs(m),:)%HfR)
         S(m,:) = N*( ens(iobs(m),:)%Hs - aveens%Hs &
                &    +ens(iobs(m),:)%HfH - aveens%HfH &
                &    +ens(iobs(m),:)%HfR - aveens%HfR )
      case('c')
         do j=1,nrens
            D(m,j) = D(m,j)- &
                        cfrac*N*(sum(ens(iobs(m),j)%I(1:na)) &
                                    +ens(iobs(m),j)%Qm       &
                                    +ens(iobs(m),j)%Qs       &
                                    +ens(iobs(m),j)%QfH      &
                                    +ens(iobs(m),j)%QfR      &
                                    +ens(iobs(m),j)%Hs       &
                                    +ens(iobs(m),j)%HfH      &
                                    +ens(iobs(m),j)%HfR      &
                                    +ens(iobs(m),j)%CH       &
                                    +ens(iobs(m),j)%CR       &
                                    +ens(iobs(m),j)%Rm       &
                                    +ens(iobs(m),j)%Rs       &
                                    +ens(iobs(m),j)%DR       &
                                    +ens(iobs(m),j)%DH        )
                                    

            S(m,j) =         N*( sum(ens(iobs(m),j)%I(:)) - sum(aveens%I(:))  &
                                    +ens(iobs(m),j)%Qm         -aveens%Qm     &
                                    +ens(iobs(m),j)%Qs         -aveens%Qs     &
                                    +ens(iobs(m),j)%QfH        -aveens%QfH    &
                                    +ens(iobs(m),j)%QfR        -aveens%QfR    &
                                    +ens(iobs(m),j)%Hs         -aveens%Hs     &
                                    +ens(iobs(m),j)%HfH        -aveens%HfH    &
                                    +ens(iobs(m),j)%HfR        -aveens%HfR    &
                                    +ens(iobs(m),j)%CH         -aveens%CH     &
                                    +ens(iobs(m),j)%CR         -aveens%CR     &
                                    +ens(iobs(m),j)%Rm         -aveens%Rm     &
                                    +ens(iobs(m),j)%Rs         -aveens%Rs     &
                                    +ens(iobs(m),j)%DH         -aveens%DH      &
                                    +ens(iobs(m),j)%DR         -aveens%DR      )
               enddo
      case default
         stop 'Measurement type not found'
      end select
!      print '(a,f12.2)','N=',N
!      print '(a,i3,10f10.2)','D',m,D(m,1:10)
!      print '(a,i3,10f10.2)','S',m,S(m,1:10)
      innovation(m)=sum(D(m,:))/real(nrens)
   enddo
! Scaling of matrices

   allocate(scaling(nrobs))
   do m=1,nrobs
      scaling(m)=1./sqrt(R(m,m))
      S(m,:)=scaling(m)*S(m,:)
      E(m,:)=scaling(m)*E(m,:)
      D(m,:)=scaling(m)*D(m,:)
      innovation(m)=scaling(m)*innovation(m)
   enddo

   do j=1,nrobs
   do i=1,nrobs
      R(i,j)=scaling(i)*R(i,j)*scaling(j)
   enddo
   enddo


end subroutine
end module
