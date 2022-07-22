MODULE module_fire_emis

    CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Simple smoke aerosol emissions scheme for Thompson-Eidhammer Aerosol-Aware microphysics scheme  !
!                                                                                                  !
!  Implemented by T. W. Juliano (NCAR/RAL) on 4/13/2022                                            !
!                                                                                                  !
!  First implementation includes two options:                                                      !
!      1. wif_fire_inj = 0 --> Emit smoke aerosols in the lowest model grid cell                   !
!                          --> Recommended for LES configuration where dx is sufficiently small    !
!                              to resolve strong vertical velocities that can mix aerosols         !
!                              throughout the BL during convective conditions                      !
!      2. wif_fire_inj = 1 --> Default option                                                      !
!                          --> Emit smoke aerosols using simple plumerise scheme that evenly       !
!                              distributes throughout the BL column based on PBL height diagnosed  !
!                              by the PBL parameterization                                         !
!                          --> Recommended for mesoscale simulations where dx is not small enough  !
!                              to explicitly mix aerosols throughout the BL during convective      !
!                              conditions due to relatively weak updrafts                          !
!                                                                                                  !
!  Future implementation could include plumerise model by Freitas et al. (2007, ACP), which is     !
!  currently part of the WRF-Chem codebase                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE fire_emis_simple_plumerise(wif_fire_inj, aer_fire_emit_opt,                  &
              fire_aerosol_tau, z_at_mass, pblh,                                            &
              nwfa, nbca, noca1, nbca1, nocbb2d, nbcbb2d, nfbc1_2d, nfoc1_2d,               &
              dt_in, ids, ide, jds, jde, kds, kde,                                          &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      IMPLICIT NONE

      INTEGER                                   , INTENT(IN)    :: wif_fire_inj, aer_fire_emit_opt
      REAL, OPTIONAL                            , INTENT(IN)    :: fire_aerosol_tau
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN)    :: z_at_mass
      REAL, DIMENSION(ims:ime, jms:jme)         , INTENT(IN)    :: pblh
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: nwfa, nbca
      REAL, OPTIONAL, DIMENSION(ims:ime, jms:jme)         , INTENT(IN)    :: nocbb2d, nbcbb2d
      REAL, OPTIONAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: noca1, nbca1
      REAL, OPTIONAL, DIMENSION(ims:ime, jms:jme)  ,        INTENT(IN)    :: nfoc1_2d, nfbc1_2d
      REAL                                      , INTENT(IN)    :: dt_in
      INTEGER                                   , INTENT(IN)    :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local
      INTEGER  :: i, j, k, i_start, i_end, j_start, j_end, k_inj
      REAL     :: noc_emit, nbc_emit

      i_start = its
      j_start = jts
      i_end   = MIN(ite, ide-1)
      j_end   = MIN(jte, jde-1)

        ! Emit fire aerosols at surface (recommended for LES configuration)
      if (wif_fire_inj .eq. 0) then
        do j = j_start, j_end
          do i = i_start, i_end
            if (aer_fire_emit_opt .eq. 11) then
                ! Emit fresh aerosols 
              noca1(i,kts,j) = noca1(i,kts,j) + nfoc1_2d(i,j)
              nbca1(i,kts,j) = nbca1(i,kts,j) + nfbc1_2d(i,j)
            else
              nwfa(i,kts,j) = nwfa(i,kts,j) + nocbb2d(i,j)*dt_in
              if (aer_fire_emit_opt .eq. 2) then
                nbca(i,kts,j) = nbca(i,kts,j) + nbcbb2d(i,j)*dt_in
              end if
            end if 
          end do
        end do
        ! Distribute fire aerosols evenly through PBL (recommended for mesoscale simulations)
      else if (wif_fire_inj .eq. 1) then
        do j = j_start, j_end
          do i = i_start, i_end
              ! Find k index for PBL top
            k_inj = kts
            do while (z_at_mass(i,k_inj,j) .lt. pblh(i,j))
              k_inj = k_inj + 1
            end do
            if (aer_fire_emit_opt .eq. 11) then
                ! Compute num of OC aerosols to emit per grid cell
              noc_emit = nfoc1_2d(i,j)/k_inj
                ! Compute num of BC aerosols to emit per grid cell    
              nbc_emit = nfbc1_2d(i,j)/k_inj
                ! Emit aerosols
              do k = kts, k_inj
                noca1(i,k,j) = noca1(i,k,j) + noc_emit
                nbca1(i,k,j) = nbca1(i,k,j) + nbc_emit
              end do

            else
                ! Compute num of OC aerosols to emit per grid cell
              noc_emit = (nocbb2d(i,j)*dt_in)/k_inj
                ! Compute num of BC aerosols to emit per grid cell
              if (aer_fire_emit_opt .eq. 2) then
                nbc_emit = (nbcbb2d(i,j)*dt_in)/k_inj
              end if
                ! Emit aerosols
              do k = kts, k_inj
                nwfa(i,k,j) = nwfa(i,k,j) + noc_emit
                if (aer_fire_emit_opt .eq. 2) then
                  nbca(i,k,j) = nbca(i,k,j) + nbc_emit
                end if
              end do
            end if
          end do
        end do
        ! Option does not exist
      else
        call wrf_error_fatal('option wif_fire_inj = ', wif_fire_inj, ' does not exist. Please set =0 or =1')
      end if

        ! Convert the hydrophobic aerosols to hydrophilic
      if (aer_fire_emit_opt .eq. 11) then
        do j = j_start, j_end
          do i = i_start, i_end
            do k = kts, kte
              nwfa(i,k,j)  = nwfa(i,k,j) + MAX(0., dt_in / (fire_aerosol_tau * 3600.0) * (noca1(i,k,j)))
              nbca(i,k,j)  = nbca(i,k,j) + MAX(0., dt_in / (fire_aerosol_tau * 3600.0) * (nbca1(i,k,j)))
              ! Here we might add the nbca1 to nbca for the radiative effect? and not to the nwfa?
              nbca1(i,k,j) = MAX(0., nbca1(i,k,j) * (1.0 - dt_in / (fire_aerosol_tau * 3600.0)))
              noca1(i,k,j) = MAX(0., noca1(i,k,j) * (1.0 - dt_in / (fire_aerosol_tau * 3600.0)))
            end do 
          end do
        end do 
      end if
    END SUBROUTINE fire_emis_simple_plumerise

END MODULE module_fire_emis
