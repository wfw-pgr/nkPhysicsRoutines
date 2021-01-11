! ====================================================== !
! === BField from Biot-Savart Law                    === !
! ====================================================== !
subroutine calc__biotsavartbfield( BField, coilS, I0, nBpt, nLpt )
  implicit none
  integer         , intent(in)    :: nBpt, nLpt
  double precision, intent(in)    :: I0
  double precision, intent(in)    :: coilS (3,nLpt)
  double precision, intent(inout) :: BField(6,nBpt)
  integer                         :: iB, iL
  double precision                :: dl(3), rvec(3), dlxrvec(3)
  double precision                :: r3Inv
  double precision, parameter     :: mu0_over_fpi = 1.d-7
  integer         , parameter     :: xp_=1, yp_=2, zp_=3, bx_=4, by_=5, bz_=6

  ! ====================================================== !
  ! === Biot-Savart Law Calculation                    === !
  ! ====================================================== !
  do iB=1, nBpt
     do iL=1, nLpt-1
        dl(xp_:zp_)        = coilS(xp_:zp_,iL+1) - coilS(xp_:zp_,iL)
        rvec(xp_:zp_)      = bfield(xp_:zp_,iB) &
             &               - 0.5d0 * ( coilS(xp_:zp_,iL+1) + coilS(xp_:zp_,iL) )
        r3Inv              = 1.d0 / ( ( rvec(xp_)**2 + rvec(yp_)**2 + rvec(zp_)**2 )**1.5d0 )
        dlxrvec(xp_)       = dl(yp_)*rvec(zp_) - dl(zp_)*rvec(yp_)
        dlxrvec(yp_)       = dl(zp_)*rvec(xp_) - dl(xp_)*rvec(zp_)
        dlxrvec(zp_)       = dl(xp_)*rvec(yp_) - dl(yp_)*rvec(xp_)
        BField(bx_:bz_,iB) = BField(bx_:bz_,iB) + mu0_over_fpi * I0 * dlxrvec(:) * r3Inv
     enddo
  enddo

  return
end subroutine calc__biotsavartbfield
