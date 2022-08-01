
!WRF:MEDIATION_LAYER:FIRE_MODEL
! Routines dealing with the atmosphere

module module_fr_fire_atm

use module_model_constants, only: cp,xlv
use module_fr_fire_util
use module_state_description, only: num_tracer 
use module_state_description, only: p_fire_smoke
use module_fr_fire_phys, only: mfuelcats, nfuelcats ! for emissions ME from sfire

! emission tables
REAL, dimension(mfuelcats), save:: oc1=0., bc1=0.
logical, save:: emis_read = .false.

contains

subroutine read_emissions_table
   implicit none
   logical, external:: wrf_dm_on_monitor
   external::wrf_dm_bcast_integer , wrf_dm_bcast_real
   integer, dimension(10)::compatible_chem_opt
   integer:: iounit,ierr,i
   character(len=128)::msg
   namelist/emissions/ oc1, bc1

IF ( wrf_dm_on_monitor() ) THEN
    ! we are the master task

    iounit=open_text_file('namelist.fire_emissions','read')

    read(iounit,emissions,iostat=ierr)
    if(ierr.ne.0)call crash('read_emissions_table: error reading namelist emissions in file namelist.fire_emissions')
    CLOSE(iounit)
    write(msg,'(a,i3,a)')'reading emissions table for',nfuelcats,' fuel categories'
    call message(msg,level=0)
ENDIF
!call wrf_dm_bcast_integer(printsums, 1)
call wrf_dm_bcast_real(oc1,  nfuelcats)
call wrf_dm_bcast_real(bc1,  nfuelcats)

emis_read=.true.

end subroutine read_emissions_table

! DME subroutine for passive tracers
subroutine add_fire_tracer_emissions(    &
       tracer_opt,dt,dx,dy,              &
       ifms,ifme,jfms,jfme,              &
       ifts,ifte,jtfs,jfte,              &
       ids,ide,kds,kde,jds,jde,          &
       ims,ime,kms,kme,jms,jme,          &
       its,ite,kts,kte,jts,jte,          &
       rho,dz8w,                         &
       burnt_area_dt,fgip,               &
       tracer,fire_tracer_smoke          &
)

implicit none
! arguments
integer,intent(in)::tracer_opt
real,intent(in)::fire_tracer_smoke
real,intent(in)::dt,dx,dy
integer,intent(in)::ifms,ifme,jfms,jfme,ifts,ifte,jtfs,jfte,ids,ide,kds,kde,jds,jde,ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte
real,intent(in)::rho(ims:ime,kms:kme,jms:jme),dz8w(ims:ime,kms:kme,jms:jme)
real,intent(in),dimension(ifms:ifme,jfms:jfme)::burnt_area_dt,fgip
real,intent(inout)::tracer(ims:ime,kms:kme,jms:jme,num_tracer)
! local
integer::isz1,jsz1,isz2,jsz2,ir,jr
integer::i,j,ibase,jbase,i_f,ioff,j_f,joff
real::avgw,emis,conv

isz1 = ite-its+1
jsz1 = jte-jts+1
isz2 = ifte-ifts+1
jsz2 = jfte-jtfs+1
ir=isz2/isz1
jr=jsz2/jsz1
avgw = 1.0/(ir*jr)

do j=max(jds+1,jts),min(jte,jde-2)
    jbase=jtfs+jr*(j-jts)
    do i=max(ids+1,its),min(ite,ide-2)
       ibase=ifts+ir*(i-its)
       do joff=0,jr-1
           j_f=joff+jbase
           do ioff=0,ir-1
               i_f=ioff+ibase
               if (num_tracer >0)then
                     emis=avgw*fire_tracer_smoke*burnt_area_dt(i_f,j_f)*fgip(i_f,j_f)*1000/(rho(i,kts,j)*dz8w(i,kts,j)) ! g_smoke/kg_air
                     tracer(i,kts,j,p_fire_smoke)=tracer(i,kts,j,p_fire_smoke)+emis
                endif
           enddo
       enddo
    enddo
enddo

end subroutine add_fire_tracer_emissions

!
!***
!

subroutine add_fire_aerosol_emissions(    &
       aer_fire_emit_opt,dt,dx,dy,     &
       ifms,ifme,jfms,jfme,              &
       ifts,ifte,jtfs,jfte,              &
       ids,ide,kds,kde,jds,jde,          &
       ims,ime,kms,kme,jms,jme,          &
       its,ite,kts,kte,jts,jte,          &
       rho,dz8w,                         &
       burnt_area_dt,fgip,nfuel_cat,     &
       qnfoc1_2d,qnfbc1_2d)
! ME subroutine for aerosols

implicit none
! arguments
integer,intent(in)::aer_fire_emit_opt
real,intent(inout),dimension(ims:ime, jms:jme)::qnfoc1_2d
real,intent(inout),dimension(ims:ime, jms:jme),optional::qnfbc1_2d
real,intent(in)::dt,dx,dy
integer,intent(in)::ifms,ifme,jfms,jfme,ifts,ifte,jtfs,jfte,ids,ide,kds,kde,jds,jde,ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte
real,intent(in)::rho(ims:ime,kms:kme,jms:jme),dz8w(ims:ime,kms:kme,jms:jme)
real,intent(in),dimension(ifms:ifme,jfms:jfme)::burnt_area_dt,fgip,nfuel_cat

! local
real,dimension(ims:ime, jms:jme)::mass_bc1, mass_oc1
integer::isz1,jsz1,isz2,jsz2,ir,jr
integer::i,j,ibase,jbase,i_f,ioff,j_f,joff,cat
real::avgw,emis,conv,fuel_burnt

if(.not.emis_read)call crash('add_fire_emissions: read_emissions_table must be called first') 

isz1 = ite-its+1
jsz1 = jte-jts+1
isz2 = ifte-ifts+1
jsz2 = jfte-jtfs+1
ir=isz2/isz1
jr=jsz2/jsz1
avgw = 1.0/(ir*jr)

do j=max(jds+1,jts),min(jte,jde-2)
    do i=max(ids+1,its),min(ite,ide-2)
       mass_bc1(i,j) = 0.
       mass_oc1(i,j) = 0.
    end do
end do

do j=max(jds+1,jts),min(jte,jde-2)
    jbase=jtfs+jr*(j-jts)
    do i=max(ids+1,its),min(ite,ide-2)
       ibase=ifts+ir*(i-its)
       do joff=0,jr-1
           j_f=joff+jbase
           do ioff=0,ir-1
               i_f=ioff+ibase

! avgw           :fire grid ratio                     (-)
! efactors(cat)  :emission factor for bc              (g/kg biomass burned)
! burnt_area_dt  :fraction of cell area burnt in dt   (-)
! fgip           :init mass of surface fuel           (kg/m^2)
! rho            :air density                         (kg/m^3)
! dz8w           :dz between full levels              (m)
! fire_qnxxxx    : # / kg

               fuel_burnt = fgip(i_f,j_f) * burnt_area_dt(i_f,j_f) ! kg/m^2 
               cat = nfuel_cat(i_f, j_f)                       ! usually 1 to 13

               ! Aerosols
               ! fuel_burnt kg/m^2 *  table g/kg 
               conv = avgw / (rho(i,kts,j) * dz8w(i,kts,j))  ! convert to g/kg dry air 
               ! BC mass
               emis = bc1(cat) * fuel_burnt
               mass_bc1(i,j) = mass_bc1(i,j) + emis * conv
               ! OC mass
               emis = oc1(cat) * fuel_burnt
               mass_oc1(i,j) = mass_oc1(i,j) + emis * conv
               
               if (aer_fire_emit_opt .eq. 11) then
                 call    calc_numc(16, mass_bc1, mass_oc1, rho, & 
                 ids,ide, jds,jde, kds,kde,       &
                 ims,ime, jms,jme, kms,kme,       &
                 its,ite, jts,jte, kts,kte,       &
                 qnfoc1_2d)
               elseif (aer_fire_emit_opt .eq. 12) then
                 call    calc_numc(16, mass_bc1, mass_oc1, rho, & 
                 ids,ide, jds,jde, kds,kde,       &
                 ims,ime, jms,jme, kms,kme,       &
                 its,ite, jts,jte, kts,kte,       &
                 qnfoc1_2d, qnfbc1_2d)
               endif
           enddo
       enddo
    enddo
enddo

end subroutine add_fire_aerosol_emissions

!
!***
!

SUBROUTINE fire_tendency( &
    ids,ide, kds,kde, jds,jde,   & ! dimensions
    ims,ime, kms,kme, jms,jme,   &
    its,ite, kts,kte, jts,jte,   &
    grnhfx,grnqfx,canhfx,canqfx, & ! heat fluxes summed up to  atm grid
    alfg,alfc,z1can,             & ! coeffients, properties, geometry
    zs,z_at_w,dz8w,mu,c1h,c2h,rho, &
    rthfrten,rqvfrten)             ! theta and Qv tendencies

! This routine is atmospheric physics
! it does NOT go into module_fr_fire_phys because it is not related to fire physical processes

! --- this routine takes fire generated heat and moisture fluxes and
!     calculates their influence on the theta and water vapor
! --- note that these tendencies are valid at the Arakawa-A location

   IMPLICIT NONE

! --- incoming variables

   INTEGER , INTENT(in) :: ids,ide, kds,kde, jds,jde, &
                           ims,ime, kms,kme, jms,jme, &
                           its,ite, kts,kte, jts,jte

   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: grnhfx,grnqfx  ! W/m^2
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: canhfx,canqfx  ! W/m^2
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: zs  ! topography (m abv sealvl)
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: mu  ! dry air mass (Pa)
   REAL, INTENT(in), DIMENSION( kms:kme         ) :: c1h, c2h ! Hybrid coordinate weights

   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: z_at_w ! m abv sealvl
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: dz8w   ! dz across w-lvl
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: rho    ! density

   REAL, INTENT(in) :: alfg ! extinction depth surface fire heat (m)
   REAL, INTENT(in) :: alfc ! extinction depth crown  fire heat (m)
   REAL, INTENT(in) :: z1can    ! height of crown fire heat release (m)

! --- outgoing variables

   REAL, INTENT(out), DIMENSION( ims:ime,kms:kme,jms:jme ) ::   &
       rthfrten, & ! theta tendency from fire (in mass units)
       rqvfrten    ! Qv tendency from fire (in mass units)
! --- local variables

   INTEGER :: i,j,k
   INTEGER :: i_st,i_en, j_st,j_en, k_st,k_en

   REAL :: cp_i
   REAL :: rho_i
   REAL :: xlv_i
   REAL :: z_w
   REAL :: fact_g, fact_c
   REAL :: alfg_i, alfc_i

   REAL, DIMENSION( its:ite,kts:kte,jts:jte ) :: hfx,qfx
   
!!   character(len=128)::msg

        do j=jts,jte
            do k=kts,min(kte+1,kde)
               do i=its,ite
                   rthfrten(i,k,j)=0.
                   rqvfrten(i,k,j)=0.
               enddo
            enddo
        enddo


! --- set some local constants
   

   cp_i = 1./cp     ! inverse of specific heat
   xlv_i = 1./xlv   ! inverse of latent heat
   alfg_i = 1./alfg
   alfc_i = 1./alfc

!!write(msg,'(8e11.3)')cp,cp_i,xlv,xlv_i,alfg,alfc,z1can
!!call message(msg)

   call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnhfx,'fire_tendency:grnhfx')
   call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnqfx,'fire_tendency:grnqfx')

! --- set loop indicies : note that

   i_st = MAX(its,ids+1)
   i_en = MIN(ite,ide-1)
   k_st = kts
   k_en = MIN(kte,kde-1)
   j_st = MAX(jts,jds+1)
   j_en = MIN(jte,jde-1)

! --- distribute fluxes

   DO j = j_st,j_en
      DO k = k_st,k_en
         DO i = i_st,i_en

            ! --- set z (in meters above ground)

            z_w = z_at_w(i,k,j) - zs(i,j) ! should be zero when k=k_st

            ! --- heat flux

            fact_g = cp_i * EXP( - alfg_i * z_w )
            IF ( z_w < z1can ) THEN
               fact_c = cp_i
            ELSE
               fact_c = cp_i * EXP( - alfc_i * (z_w - z1can) )
            END IF
            hfx(i,k,j) = fact_g * grnhfx(i,j) + fact_c * canhfx(i,j)

!!            write(msg,2)i,k,j,z_w,grnhfx(i,j),hfx(i,k,j)
!!2           format('hfx:',3i4,6e11.3)
!!            call message(msg)

            ! --- vapor flux

            fact_g = xlv_i * EXP( - alfg_i * z_w )
            IF (z_w < z1can) THEN
               fact_c = xlv_i
            ELSE
               fact_c = xlv_i * EXP( - alfc_i * (z_w - z1can) )
            END IF
            qfx(i,k,j) = fact_g * grnqfx(i,j) + fact_c * canqfx(i,j)
            
!!            if(hfx(i,k,j).ne.0. .or. qfx(i,k,j) .ne. 0.)then
!!                write(msg,1)i,k,j,hfx(i,k,j),qfx(i,k,j)
!!1               format('tend:',3i6,2e11.3)
!!                call message(msg)
!            endif

         END DO
      END DO
   END DO

! --- add flux divergence to tendencies
!
!   multiply by dry air mass (mu) to eliminate the need to
!   call sr. calculate_phy_tend (in dyn_em/module_em.F)

   DO j = j_st,j_en
      DO k = k_st,k_en-1
         DO i = i_st,i_en

            rho_i = 1./rho(i,k,j)

            rthfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (hfx(i,k+1,j)-hfx(i,k,j)) / dz8w(i,k,j)
            rqvfrten(i,k,j) = - (c1h(k)*mu(i,j)+c2h(k)) * rho_i * (qfx(i,k+1,j)-qfx(i,k,j)) / dz8w(i,k,j)

         END DO
      END DO
   END DO

   call print_3d_stats(its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme,rthfrten,'fire_tendency:rthfrten')
   call print_3d_stats(its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme,rqvfrten,'fire_tendency:rqvfrten')

   RETURN

END SUBROUTINE fire_tendency

!
!***
!

subroutine calc_numc(nbin_o, emism_bc1, emism_oc1, rho,   & 
          ids,ide, jds,jde, kds,kde,                    &
          ims,ime, jms,jme, kms,kme,                    &
          its,ite, jts,jte, kts,kte,                    &
          oc1_2d, bc1_2d)

    USE module_configure
    USE module_model_constants
    implicit none
    real, parameter :: hygro_msa_aer = 0.58
!
    INTEGER, INTENT(IN) :: its,ite, jts,jte, kts,kte, nbin_o
    INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme
    INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde
    REAL, DIMENSION(ims:ime, jms:jme),             &
          INTENT(IN ) ::  emism_bc1, emism_oc1
    REAL, DIMENSION(ims:ime, jms:jme), OPTIONAL,   &
          INTENT(OUT ) :: bc1_2d, oc1_2d
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  &
          INTENT(IN ) ::  rho!,relhum

!
! local variables
!

   REAL, DIMENSION(1:nbin_o):: number_bin

   integer i, j, k, l, m, n, isize

   real  dens_bc, dens_oc
   real  mass_bc, mass_oc
   real  mass_bci, mass_oci
   real  mass_bcj, mass_ocj
   real  conv1a, conv1b
   real  mass_fac

   integer ns
   real  dgnum_um,duma,dlo_um,dhi_um,dgmin,piosix,dtemp
   real, dimension(1:nbin_o) :: xmas_secti,xmas_sectj!,xmas_sectc
   real, dimension(1:nbin_o) :: xdia_cm
   REAL, PARAMETER :: FRAC2Aitken=0.25   ! Fraction of modal mass in Aitken mode - applied globally to each species

! 7/21/09  SAM variables needed to convert GOCART sectional dust and seasalt to MOZAIC sections
   real  dgnum, dhi,  dlo, xlo, xhi, dxbin, relh_frc
   real dlo_sect(nbin_o), dhi_sect(nbin_o)
   real*8 dlogoc, dhigoc
   integer istop
   integer, save :: kcall
   data  kcall / 0 /
   CHARACTER*256     :: outstring
   REAL, PARAMETER :: sginin=1.7     !  initial sigma-G for nucleimode
   REAL, PARAMETER :: sginia=2.0     !  initial sigma-G for accumulation mode
!   REAL, PARAMETER :: sginic=2.5     ! initial sigma-G for coarse mode 
   REAL, PARAMETER :: dginin=0.01E-6 !  initial mean diameter for nuclei mode [ m ]    
   REAL, PARAMETER :: dginia=0.07E-6 !  initial mean diameter for accumulation mode [ m ]
!   REAL, PARAMETER :: dginic=1.0E-6  ! initial mean diameter for coarse mode [ m ]  


      piosix=3.14159265359/6.0
!      dlo_um=0.0390625
      dlo_um=0.004
      dhi_um=10.0
!      drydens=1.8

      duma=1.0

      xlo = log(dlo_um * 1.0e-4)
      xhi = log(dhi_um * 1.0e-4)
      dxbin = (xhi - xlo)/nbin_o
      do isize=1,nbin_o
        dlo_sect(isize) = exp(xlo + dxbin * (isize-1))
        dhi_sect(isize) = exp(xlo + dxbin * isize)
        xdia_cm(isize)=(dlo_sect(isize)+dhi_sect(isize))/2.0
      enddo
!
! densities in g/cc !######!
!
      dens_oc    = 1.0        ! g/cc
! JCB, Feb. 20, 2008:  the density of BC is updated to reflect values
! published by Bond and Bergstrom, Light absorption by carboneceous
! particles: an investigative review 2006, Aerosol Sci. and Tech., 40:27-67.
      dens_bc    =  1.8        ! g/cc  midpoint of Bond and Bergstrom value

!
! Now divide mass into sections which is done by sect02:
! * xmas_secti is for aiken mode
! * xmas_sectj is for accumulation mode
! * sect02 expects input in um
! * pass in generic mass of 1.0 just to get a percentage distribution of mass among bins

      dgnum_um=dginin*1.E6
      call sect02(dgnum_um,sginin,duma,nbin_o,dlo_um,dhi_um, &
        xmas_secti)

      dgnum_um=dginia*1.E6
      call sect02(dgnum_um,sginia,duma,nbin_o,dlo_um,dhi_um, &
        xmas_sectj)

!
! units:
! * mass     - g/cc(air)
! * number   - #/cc(air)
! * volume   - cc(air)/cc(air)
! * diameter - cm
!
    
    do isize = 1, nbin_o
     ! Accumulation mode... (j)
     ! Aitken mode...       (i)
      mass_fac = ((1.-FRAC2Aitken)*xmas_sectj(isize) + FRAC2Aitken*xmas_secti(isize))
            ! Let's assume everything is dry in the beginning
      number_bin(isize) = mass_fac / (piosix*xdia_cm(isize)*xdia_cm(isize)*xdia_cm(isize))
    enddo  !isize

    do j=max(jds+1,jts),min(jte,jde-2)
      do i=max(ids+1,its),min(ite,ide-2)
        if (emism_bc1(i,j) .gt. 0.) then
  ! convert ug / kg dry air to g / kg
          conv1a =  1.0E-6
          if (present(bc1_2d)) then
            bc1_2d(i, j)= emism_bc1(i,j) * conv1a / dens_bc * sum(number_bin)
          endif
          oc1_2d(i, j)= emism_oc1(i,j) * conv1a / dens_oc * sum(number_bin)
        endif
      end do
    end do
    return
end subroutine calc_numc

!
!***
!

subroutine sect02(dgnum_um,sigmag,duma,nbin,dlo_um,dhi_um, &
                  xmas_sect)  ! Borrowed from CHEM (ME)
!
!   user specifies a single log-normal mode and a set of section boundaries
!   prog calculates mass for each section
!
    implicit none
    integer, INTENT(IN):: nbin
    REAL, DIMENSION(nbin), INTENT(OUT) :: xmas_sect
    integer n
    real &
      dgnum, dgnum_um, dhi, dhi_um, dlo, dlo_um,   &
       duma, dumfrac, dx,   &
      sigmag, sumnum, summas,   &
      sx, sxroot2, thi, tlo, vtot,   &
      x0, x3, xhi, xlo, xmtot, xvtot
    real dlo_sect(nbin), dhi_sect(nbin)
!       real erfc_num_recipes
    real pi
    parameter (pi = 3.1415926536)

    xmtot = duma

!   compute section boundaries
    dlo = dlo_um*1.0e-4
    dhi = dhi_um*1.0e-4
    xlo = log( dlo )
    xhi = log( dhi )
    dx = (xhi - xlo)/nbin
    do n = 1, nbin
        dlo_sect(n) = exp( xlo + dx*(n-1) )
        dhi_sect(n) = exp( xlo + dx*n )
    end do
!   compute modal "working" parameters including total num/vol/mass
    dgnum = dgnum_um*1.0e-4
    sx = log( sigmag )
    x0 = log( dgnum )
    x3 = x0 + 3.*sx*sx

!   compute mass for each section
    sxroot2 = sx * sqrt( 2.0 )
    summas = 0.

    do n = 1, nbin
        xlo = log( dlo_sect(n) )
        xhi = log( dhi_sect(n) )
        tlo = (xlo - x3)/sxroot2
        thi = (xhi - x3)/sxroot2
        if (tlo .le. 0.) then
            dumfrac = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
        else
            dumfrac = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
        end if
        xmas_sect(n) = xmtot*dumfrac
        summas = summas + xmas_sect(n)
    end do
end subroutine  sect02

!
!***
!
        real function erfc_num_recipes( x )
!
!   from press et al, numerical recipes, 1990, page 164
!
        implicit none
        real x
        double precision erfc_dbl, dum, t, z
        z = abs(x)
        t = 1.0/(1.0 + 0.5*z)
        dum =  ( -z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 +   &
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +   &
                                           t*(-1.13520398 +   &
          t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))
        erfc_dbl = t * exp(dum)
        if (x .lt. 0.0) erfc_dbl = 2.0d0 - erfc_dbl
        erfc_num_recipes = erfc_dbl
        return

        end function erfc_num_recipes
!-----------------------------------------------------------------------
!
!***
!

end module module_fr_fire_atm
