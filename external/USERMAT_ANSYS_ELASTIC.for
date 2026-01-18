c----------------------------------------------------------------------
c     USERMAT_ANSYS_ELASTIC.for
c     Sample Ansys USERMAT subroutine for isotropic linear elasticity
c     
c     This file is part of simcoon examples.
c     
c     Material properties:
c       prop(1) = E    (Young's modulus)
c       prop(2) = nu   (Poisson's ratio)
c       prop(3) = alpha (thermal expansion coefficient, optional)
c
c     Reference: Ansys Mechanical APDL Programmer's Reference
c----------------------------------------------------------------------

      subroutine usermat(
     &     matId, elemId, kDomIntPt, kLayer, kSectPt,
     &     ldstep, isubst, keycut,
     &     nDirect, nShear, ncomp, nStatev, nProp,
     &     Time, dTime, Temp, dTemp,
     &     stress, ustatev, dsdde,
     &     sedEl, sedPl, epseq,
     &     Strain, dStrain, epsPl,
     &     prop, coords, rotateM,
     &     defGrad_t, defGrad,
     &     tsstif, epsZZ,
     &     var1, var2, var3, var4, var5)
c
      implicit none
c
c     Arguments
      integer matId, elemId, kDomIntPt, kLayer, kSectPt
      integer ldstep, isubst, keycut
      integer nDirect, nShear, ncomp, nStatev, nProp
      integer var1, var2, var3, var4, var5
c
      double precision Time, dTime, Temp, dTemp
      double precision sedEl, sedPl, epseq, epsZZ
      double precision stress(ncomp), ustatev(nStatev)
      double precision dsdde(ncomp, ncomp)
      double precision Strain(ncomp), dStrain(ncomp)
      double precision epsPl(ncomp)
      double precision prop(nProp), coords(3)
      double precision rotateM(3,3)
      double precision defGrad_t(3,3), defGrad(3,3)
      double precision tsstif(2)
c
c     Local variables
      double precision E, nu, alpha
      double precision lambda, mu
      double precision Eel(6), T_init
      integer i, j
c
c     Material properties
      E = prop(1)
      nu = prop(2)
      if (nProp .ge. 3) then
         alpha = prop(3)
      else
         alpha = 0.0d0
      endif
c
c     Lame constants
      lambda = E * nu / ((1.0d0 + nu) * (1.0d0 - 2.0d0 * nu))
      mu = E / (2.0d0 * (1.0d0 + nu))
c
c     Initialize tangent modulus (isotropic elasticity)
c     Ansys Voigt notation: 11, 22, 33, 12, 23, 13
      do i = 1, ncomp
         do j = 1, ncomp
            dsdde(i,j) = 0.0d0
         enddo
      enddo
c
c     Diagonal terms (normal stresses)
      dsdde(1,1) = lambda + 2.0d0 * mu
      dsdde(2,2) = lambda + 2.0d0 * mu
      dsdde(3,3) = lambda + 2.0d0 * mu
c
c     Off-diagonal terms (coupling between normal stresses)
      dsdde(1,2) = lambda
      dsdde(1,3) = lambda
      dsdde(2,1) = lambda
      dsdde(2,3) = lambda
      dsdde(3,1) = lambda
      dsdde(3,2) = lambda
c
c     Shear terms
      dsdde(4,4) = mu
      dsdde(5,5) = mu
      dsdde(6,6) = mu
c
c     Initialize reference temperature from state variable
      if (nStatev .ge. 5) then
         T_init = ustatev(5)
      else
         T_init = Temp
      endif
c
c     First increment: initialize
      if (ldstep .eq. 1 .and. isubst .eq. 1) then
         T_init = Temp
         do i = 1, ncomp
            stress(i) = 0.0d0
         enddo
         sedEl = 0.0d0
         sedPl = 0.0d0
         if (nStatev .ge. 5) then
            ustatev(5) = T_init
         endif
      endif
c
c     Compute elastic strain (total strain minus thermal strain)
      do i = 1, 3
         Eel(i) = Strain(i) + dStrain(i) - alpha * (Temp + dTemp - T_init)
      enddo
      do i = 4, ncomp
         Eel(i) = Strain(i) + dStrain(i)
      enddo
c
c     Compute stress: sigma = C : epsilon_elastic
      do i = 1, ncomp
         stress(i) = 0.0d0
         do j = 1, ncomp
            stress(i) = stress(i) + dsdde(i,j) * Eel(j)
         enddo
      enddo
c
c     Update elastic strain energy density
      sedEl = 0.0d0
      do i = 1, ncomp
         sedEl = sedEl + 0.5d0 * stress(i) * Eel(i)
      enddo
c
c     No cutback requested
      keycut = 0
c
      return
      end
