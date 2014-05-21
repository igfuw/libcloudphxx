module kessler_wrap

  use iso_c_binding, only: c_int, c_float
  use module_mp_kessler 

  implicit none

contains

  subroutine c_kessler(t, qv, qc, qr, rho, pii                  &
                      ,dt_in, z, xlv, cp                        &
                      ,EP2,SVP1,SVP2,SVP3,SVPT0,rhowater        &
                      ,dz8w                                     &
                      ,RAINNC, RAINNCV                          &
                      ,ids,ide, jds,jde, kds,kde                & 
                      ,ims,ime, jms,jme, kms,kme                & 
                      ,its,ite, jts,jte, kts,kte                & 
                      ) bind(c)
    integer(c_int), intent(in), value:: ids,ide, jds,jde, kds,kde, & 
                                        ims,ime, jms,jme, kms,kme, & 
                                        its,ite, jts,jte, kts,kte
    real(c_float),  intent(in), value:: xlv, cp,                   &
                                        EP2,SVP1,SVP2,SVP3,SVPT0,  &
                                        rhowater
    real(c_float),  intent(inout):: t(ims:ime,kms:kme,jms:jme),    &
                                   qv(ims:ime,kms:kme,jms:jme),    &
                                   qc(ims:ime,kms:kme,jms:jme),    &
                                   qr(ims:ime,kms:kme,jms:jme)
    real(c_float),  intent(in)::  rho(ims:ime,kms:kme,jms:jme),    &
                                  pii(ims:ime,kms:kme,jms:jme),    &
                                 dz8w(ims:ime,kms:kme,jms:jme),    &
                                    z(ims:ime,kms:kme,jms:jme)
    
    real(c_float),  intent(in), value:: dt_in
    real(c_float),  intent(inout):: RAINNC(ims:ime,jms:jme),       &
                                   RAINNCV(ims:ime,jms:jme)
    
    print*, "foo", ims, t(1,1,1), qv(1,1,1), qc(1,1,1)
    call kessler(t, qv, qc, qr, rho, pii                        &
                      ,dt_in, z, xlv, cp                        &
                      ,EP2,SVP1,SVP2,SVP3,SVPT0,rhowater        &
                      ,dz8w                                     &
                      ,RAINNC, RAINNCV                          &
                      ,ids,ide, jds,jde, kds,kde                & 
                      ,ims,ime, jms,jme, kms,kme                & 
                      ,its,ite, jts,jte, kts,kte                & 
                      )
    print*, "after calling", ims, t(1,1,1), qv(1,1,1), qc(1,1,1), xlv

    end subroutine c_kessler
end module kessler_wrap
