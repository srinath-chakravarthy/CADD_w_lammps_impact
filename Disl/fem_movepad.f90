!!$!*==fem_move_pad.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!!$!
!!$! $Id: fem_movepad.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!!$!
SUBROUTINE FEM_MOVE_PAD(X,B,Ix,Bc,Prop,Cz,Id,Is_relaxed,E_out,&
     &                        F_out,Fullfield,Movedisl,Straine0,Numel,&
     &                        Avevirst,Mdtemp,Ifem,Moved)
!!$!
!!$!12345678901234567890123456789012345678901234567890123456789012345678901
!!$!         1         2         3         4         5         6         7
!!$!
!!$! if the node is on the interface u=b
!!$! if the node is on the boundary u=prop*bc
!!$! changes b
!!$! e_out for informational purposes only
!!$!
!!$!***********************************************************************
!!$!       Comment ---- here prop=time variable which is really
!!$!       the load increment and basically the total displacments/ forces
!!$!       are scaled with prop in all the fe_routines, but really time
!!$!       from the md/input files
!!$!       Correllation with md parameters
!!$!       cz --- z_length ... length of md box in the z direction
!!$!       prop = time
!!$!       f_out ---- dr (forces)
!!$!       e_out ---- eadd
!!$!       x, b, ix  are the same parameters
!!$!       f---- bc  boundary conditions (id=1, displacement, id=0, force)
!!$!       db---- displacement increment/step
!!$!***********************************************************************
      USE MOD_DD_SLIP
      USE MOD_FEM_PARAMETERS
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!!$!*--FEM_MOVE_PAD35

      DOUBLE PRECISION X(3,*) , B(3,*) , Bc(3,*) , F_out(3,*)
      DOUBLE PRECISION Avevirst(3,3,*) , tincr1
      INTEGER Numel
      INTEGER Id(3,*) , Ix(4,*) , Is_relaxed(*)
      DOUBLE PRECISION Prop , Cz
!!$!C--Jun Song: MDTemp=SysTemp defined in mod_global
      DOUBLE PRECISION E_out , Straine0 , Mdtemp
      LOGICAL Fullfield , Movedisl
      LOGICAL Moved

      DOUBLE PRECISION rhs(MAXEQS) , forces(MAXEQS) , e0 , eplus , eminus
      INTEGER i , Ifem
      INTEGER fe_locate , elem_old , elem_plus , elem_minus
      CHARACTER*80 error_message
      LOGICAL centraldiff , forwarddiff , backdiff
      LOGICAL changetime


!!$!      print *, 'Ifem in fem_move_pad is', iFem
      IF ( I_Flag/=1 ) THEN
         error_message = 'fem_solve: call fem_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
      IF ( I_Disl/=1 ) THEN
         error_message = 'fem_solve: call disl_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF

!!$!       solve FEM and DD problem
      TINcr = TINCR_SAV
!!$      IF ( Moved ) TINcr = 1.D-20
!!$      CALL ASSIGN_DISLOC_GLOBAL
      PRINT * , 'Total dislocations' , NDIsl
      PRINT * , 'Time entering fd_solve is' , Prop
      CALL FD_SOLVE(B,Bc,Prop,Cz,Id,Is_relaxed,rhs,forces,e0)

!!$!  compute the P.-K. force on dislocations
      IF ( Movedisl ) THEN
         CALL FD_PEACH_KOELLER(rhs)
!!$         IF ( NDIsl>0 ) THEN
!!$            IF ( Ifem==1 ) THEN
!!$               PRINT * , ' --- Entering Move Disl --'
!!$                 ! call fd_peach_koeller(rhs)
!!$               CALL MOVE_DISLOC(rhs)
!!$               CALL NUCLEATE_DISL(rhs)
!!$               DTIme_dd = DTIme_dd + TINcr
!!$               PRINT * , 'Dtime_dd =' , DTIme_dd
!!$            ELSE
!!$               TINcr = TINCR_SAV/100.0D0
!!$               PRINT * , ' --- Entering Move Disl --'
!!$                 ! call fd_peach_koeller(rhs)
!!$               CALL MOVE_DISLOC(rhs)
!!$               CALL NUCLEATE_DISL(rhs)
!!$               DTIme_dd = DTIme_dd + TINcr
!!$               PRINT * , 'Dtime_dd =' , DTIme_dd
!!$!
!!$            ENDIF
!!$         ELSEIF ( Ifem==1 ) THEN
!!$            CALL MOVE_DISLOC(rhs)
!!$            CALL NUCLEATE_DISL(rhs)
!!$            DTIme_dd = DTIme_dd + TINcr
!!$            PRINT * , 'Dtime_dd =' , DTIme_dd
!!$         ENDIF
      ENDIF

      DO i = 1 , NFIxed
         F_out(IFIx_hold(1,i),IMAp(IFIx_hold(2,i))) = F_out(IFIx_hold(1,i),IMAp(IFIx_hold(2,i))) - forces(IFIxed(i))*Cz
      ENDDO
      E_out = e0*Cz

!!$!       move dislocations based on P.-K. force
      if(MoveDisl) call move_dis(10.0,MDTemp)

!!$!       move pad atoms according to tilda and hat fields
      CALL FD_MOVEPAD(X,rhs,B)

!!$!       update b so b=u_hat+u_tilda
      IF ( Fullfield ) CALL FD_FULL_FIELD(rhs,B)

!!$!       compute total strain energy of fem region
      CALL FE_STRESS_CHECK(rhs,Straine0)

END SUBROUTINE FEM_MOVE_PAD
!!$!*==fd_movepad.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015



SUBROUTINE FD_MOVEPAD(X,Rhs,B)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!!$!*--FD_MOVEPAD131
      DOUBLE PRECISION X(3,*) , Rhs(*) , B(3,*)

      DOUBLE PRECISION u_out(3)
      INTEGER i , j , k , lmn , n

      DO i = 1 , NPAd
         n = PADmap(i)
         lmn = PADelmnt(i)
         DO j = 1 , NDOF
            B(j,n) = 0.D0
            DO k = 1 , KNODE
               B(j,n) = B(j,n) + Rhs((ICOnn(k,lmn)-1)*NDOF+j)&
                    &                  *PADtricoord(k,i)
            ENDDO
         ENDDO
         CALL DISL_DISPL(X(1,n),u_out)
         DO j = 1 , NDOF
            B(j,n) = B(j,n) + u_out(j)
         ENDDO
      ENDDO
END SUBROUTINE FD_MOVEPAD
!!$!*==fe_tricoord.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
                                !


SUBROUTINE FE_TRICOORD(R1,R2,R3,P,Coord)

!!$!  Yet another wrapper for intri from feaplib

      IMPLICIT NONE
!!$!*--FE_TRICOORD162
      DOUBLE PRECISION R1(3) , R2(3) , R3(3) , P(3) , Coord(3)

      DOUBLE PRECISION x1(2) , x2(2) , x3(2) , pt(2)
      LOGICAL INTRI , ontri
      INTEGER i
      CHARACTER*80 error_message

      DO i = 1 , 2
         x1(i) = R1(i)
         x2(i) = R2(i)
         x3(i) = R3(i)
         pt(i) = P(i)
      ENDDO
      IF ( .NOT.INTRI(x1,x2,x3,pt,Coord,ontri) ) THEN
         error_message = 'fe_tricoord: the atom is not in the element'
         CALL ERROR_HANDLER(error_message)
      ENDIF
END SUBROUTINE FE_TRICOORD
!!$!*==move_dis.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!!$ 
!!$!
!!$! $Log: fem_movepad.f,v $
!!$! Revision 1.1.1.1  2003/03/1220:09:00  shastry
!!$! vijay-   Initial import.
!!$!
!!$!




!!$!!!!!!!!dw/cs added subroutine!!!!!!!!!!!!!!!!!!
SUBROUTINE MOVE_DIS(Alpha,Temperature)
      USE MOD_DISL_PARAMETERS
      USE MOD_GLOBAL
      IMPLICIT NONE
!!$!*--MOVE_DIS197
      DOUBLE PRECISION Alpha , mobility , max_vel, max_ds
      DOUBLE PRECISION sf_f , aa , bb, deltas
      DOUBLE PRECISION min_pos , Temperature , time_step_con

      INTEGER FE_LOCATE , i , elem_old
      DOUBLE PRECISION b
      DOUBLE PRECISION ddis
      double precision :: ev_convert1
      CHARACTER*80 error_message
                                !
!!$!!!!    hacked parameters
      min_pos = -10.0
      time_step_con = 1000.0d0 * lammps_timestep

      !!--Jun Song: make sure temperature>0.0
      IF ( Temperature<0.0D0 ) THEN
         WRITE (*,*) "Temperature less than Zero!!!"
         STOP
      ENDIF
      ev_convert1 = 0.62415096471d-11 !> Convert from Pa to ev/A^3
      mobility = 5.0d-8 !> Actual value is Pa s
      mobility = mobility * 1.0d-15 * ev_convert1 !> Converted value is ev/A^3 fs
      max_vel = 2.0d3 !> Actual value in m/s
      max_vel = max_vel * 1.0d10/1.0d15 !> Converted to A/fs
      max_ds = max_vel * time_step_con
      sf_f = 0.0 ! zero stacking fault energy for hex al

      !!--JS: 6.242e-2 is unit conversion constant. Do not change
      !!--Change the stacking fault E for different materials
      !!mobility = time_step_con/(6.242E-2*5.0E-8*Temperature)
      !! 3rd parameter (5.0e-8) is damp coef from Olmstead paper
      !! "Atomistic simulations of dislocation mobility.."
!!$      max_vel = time_step_con*2000.0
!!$      sf_f = .089*6.242E-2 ! sf energy in J/m2,i.e., 0.089
!!!!    end of hacked parameters

      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) THEN
            PK_f(i) = (PK_force(1,i)*BURgers(1,i)+PK_force(2,i) *BURgers(2,i))/BURg_length(i)
            ddis = SQRT(R_Disl(1,i)**2+R_Disl(2,i)**2)
            PK_f(i) = PK_f(i) + sf_f
!!$        upper limit on velocities
            deltas = PK_f(i) * time_step_con/mobility
            if (deltas > max_ds) then
               deltas = max_ds
            end if

!!$            IF ( ABS(PK_f(i)*mobility)>max_vel ) then
!!$               PK_f(i) = max_vel/mobility*PK_f(i)/ABS(PK_f(i))
!!$            END IF
            WRITE (*,'(A,I7,2(1X,F15.6))') 'old disl pos = ' , i, R_Disl(1,i) , R_Disl(2,i)

!!$            IF ( R_Disl(2,i)<=min_pos .OR. PK_f(i)<0.0 ) THEN
            R_Disl(1,i) = R_Disl(1,i) + deltas*BURgers(1,i)/BURg_length(i)
            R_Disl(2,i) = R_Disl(2,i) + deltas*BURgers(2,i)/BURg_length(i)
!!$            ENDIF
            WRITE (*,'(A,I7,3(1X,F15.6))') 'new disl pos = ' , i, R_Disl(1,i) , R_Disl(2,i), deltas
            WRITE (*,*) ' '
            elem_old = ELEm_disl(i)
            ELEm_disl(i) = FE_LOCATE(R_Disl(1,i),ELEm_disl(i))
            IF ( ELEm_disl(i)==0 ) THEN
               IF ( elem_old>0 ) THEN
                  ELEm_disl(i) = -elem_old
               ELSE
                  ELEm_disl(i) = elem_old
               ENDIF
            ENDIF
         ENDIF
      ENDDO
END SUBROUTINE MOVE_DIS
