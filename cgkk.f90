!====================================================================
MODULE nrtype
	INTEGER, PARAMETER :: SP = KIND(1.0)
 	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
END MODULE
!====================================================================
MODULE parameters
USE nrtype

	! benchmark only or with experiment? General equilibrium?
	INTEGER(I4B),  PARAMETER :: RunExp=0, GE=1
	
	! general equilibrium indicator
	INTEGER(I4B),  PARAMETER :: max_Q_iter=50, LIP=0 ! 1=GE, 0=PE, LIP=0 means cubic-spline
	
	! parameters to play around
	REAL(DP), PARAMETER :: rho_z=0.80_DP, sigma_z_eps=0.95_DP	! parameters for z process
	REAL(DP), PARAMETER :: rho_m=0.80_DP, sigma_m_eps=0.15_DP	! parameters for z process
	REAL(DP), PARAMETER :: chi=0.0_DP ! bequest motives, chi*log(bequest)
	REAL(DP), PARAMETER :: beta=0.99_DP

	! production 
	REAL(DP), PARAMETER :: alpha=0.3_DP, Aprod=1.0_DP, delta=0.0_DP, hp=0.002_DP ! hp=home production (min cons)
	REAL(DP), PARAMETER :: rho=0.92_DP
	
	! grid size
	INTEGER(I4B),  PARAMETER  :: ne=5        ! labor idiosyncratic shocks
	INTEGER(I4B),  PARAMETER  :: ne_alpha=5  ! labor efficiency fixed effect
	INTEGER(I4B),  PARAMETER  :: nz=7        ! inv ability fixed component
	INTEGER(I4B),  PARAMETER  :: nm=1        ! inv ability idiosyncratic shocks over lifecycle
	INTEGER(I4B),  PARAMETER  :: ns=nz*ne_alpha*ne*nm
	REAL(DP),      PARAMETER  :: tauchen_m_z=3.0_DP ! how many std away from mean in tauchen
	
	! asset grid
	REAL(DP), PARAMETER      :: amax=10000.0_DP, amin=0.0_DP, a_theta=2.5_DP
	INTEGER(I4B),  PARAMETER :: na=200
	REAL(DP), DIMENSION(na)  :: agrid

	! use simulation panel
	INTEGER(I4B),  PARAMETER :: maxiter_panel=2100, numiter_fin_panel=100
	INTEGER(I4B), PARAMETER  :: panelN=10000 ! number of agents alive in last period
	INTEGER(I4B)             :: panel_pop    ! number of agents alive in the panel
	
	! tax rates
	REAL(DP), PARAMETER :: tauw0_initial=0.015_DP ! Initial wealth tax guess
	REAL(DP), PARAMETER :: tauw0_inc    =0.005_DP
	REAL(DP), PARAMETER :: taul0=0.30_DP ! labor income tax
	REAL(DP), PARAMETER :: tauk_bench=0.25_DP, tauk_exp=0.0_DP
	REAL(DP), PARAMETER :: tauw_bench=0.0_DP
	REAL(DP) :: tauw_exp
	
	! life-cycle elements
	INTEGER(I4B), PARAMETER       :: MaxAge=80, ReAge=45  ! if no retirement, put ReAge=maxAge+1
	REAL(DP), DIMENSION(MaxAge)   :: pr_cond_surv
	REAL(DP), PARAMETER           :: coeff_r=0.5 ! the scaling for retirement income
	REAL(DP), DIMENSION(ne_alpha) :: n_alpha, y_alpha
	REAL(DP)                      :: nr
	
	! labor efficiency shocks
	! log(y)=e_alpha+e_kappa+e where e=log(x')=rho*log(x)+e_epsilon'
	! e_alpha: inidividual fixed effect (fixed within generation)
	! e_kappa: life-cycle component
	REAL(DP), PARAMETER                      :: rho_e=0.9_DP, sigma_e_eps=0.20_DP
	REAL(DP), PARAMETER                      :: rho_alpha=0.5_DP, sigma_alpha_eps=0.1_DP		
	REAL(DP), DIMENSION(ne)                  :: egrid, Dist_e
	REAL(DP), DIMENSION(ne,ne)               :: pr_e
	REAL(DP), DIMENSION(ne_alpha)            :: e_alpha,Dist_alpha
	REAL(DP), DIMENSION(ne_alpha,ne_alpha)   :: pr_alpha
	REAL(DP), DIMENSION(MaxAge)              :: e_kappa
	
	! investment abilities
	REAL(DP), DIMENSION(nz)                  :: zgrid,Dist_z
	REAL(DP), DIMENSION(nz,nz)               :: pr_z
	REAL(DP), DIMENSION(nm)                  :: mgrid,Dist_m
	REAL(DP), DIMENSION(nm,nm)               :: pr_m
	
	! state variable mapping
	INTEGER(I4B), DIMENSION(ns) :: sz,sm,salpha,se
	INTEGER(I4B), DIMENSION(nz,nm,ne_alpha,ne) :: ssfun
	
	! tolerance
	REAL(DP), PARAMETER :: tolQ=1.0e-3_DP, tolbrent=1.0e-8_DP
	
END MODULE
!====================================================================
MODULE global
USE nrtype
USE parameters

	integer(i4b) :: ce_run
	
	REAL(DP), DIMENSION(ns,na,maxAge)    :: value_bench, value_exp, v2_bench, v2_exp
	REAL(DP), DIMENSION(ns,na,maxAge)    :: aprime_bench, aprime_exp
	REAL(DP) :: e_med, z_med
	
	REAL(DP), DIMENSION(nz,nm,na) :: x_product ! intermediate products households produce
	REAL(DP), DIMENSION(nz,nm,na,maxAge) :: x_product_rho
	REAL(DP), DIMENSION(ne_alpha,ne,maxAge) :: eff_labor
	
	REAL(DP) :: rr_bench, ww_bench, Gbar_bench, Q_bench, Nbar_bench, Nbar_exp
	REAL(DP) :: rr_exp, ww_exp, Gbar_exp, Q_exp, tauw0
	REAL(DP) :: tauw0_upper, tauw0_lower
	REAL(DP),DIMENSION(8) :: stats_bench, stats_exp

	INTEGER(I4B),  DIMENSION(maxAge)    :: size_by_age
	INTEGER(I4B),  DIMENSION(maxAge,nz) :: size_by_age_z
	
	INTEGER(I4B),  DIMENSION(panelN)          :: panel_s_init
	INTEGER(I4B),  DIMENSION(panelN,maxAge)   :: panel_bench_s, panel_bench_age
	INTEGER(I4B),  DIMENSION(panelN,maxAge)   :: panel_exp_s, panel_exp_age
	REAL(DP), DIMENSION(panelN,maxAge)   :: panel_bench_a, panel_exp_a
	REAL(DP), DIMENSION(maxAge) :: Coeff_age
	
	REAL(DP),DIMENSION(numiter_fin_panel) :: CE_fin, CE_pos_panel, CE_newborn


END MODULE
!====================================================================
PROGRAM main
USE GLOBAL
IMPLICIT NONE
! include "mpif.h"

	REAL(DP), DIMENSION(ns,na,maxAge) :: value, v2
	INTEGER(I4B)  ::  iterGbar
	REAL(DP), DIMENSION(2) :: Gbar_ab
	REAL(DP) :: time1, time2

	ce_run=0

	call cpu_time(time1)
	
	!----------------------------------------------
	! initialize asset grid and exogenous process
	!----------------------------------------------
	
	CALL initial

	!----------------------------------------------
	! benchmark model with capital income tax only
	!----------------------------------------------
	
	print*, '----------------------------------------------------------'
	print*, ' Benchmark                        unbundled               '
	print*, '----------------------------------------------------------'
	
	Q_bench =15.0_DP ! Nbar_bench from initial subroutine
	
	!----------------------------------------------
	! Calling solution for benchmark
	!----------------------------------------------
		
	CALL solution(rr_bench,ww_bench,tauk_bench,tauw_bench,taul0, &
		& Gbar_bench,value_bench,v2_bench,aprime_bench, &
		& Q_bench,Nbar_bench,panel_bench_s,panel_bench_a,panel_bench_age, &
		& 1,stats_bench)
		
	write(*,"('Gbar  = ',f15.8)") Gbar_bench
	write(*,"('rr    = ',f15.8)") rr_bench
	write(*,"('ww    = ',f15.8)") ww_bench
	write(*,"('Q     = ',f15.8)") Q_bench
	write(*,"('Nbar  = ',f15.8)") Nbar_bench
	
	write(*,"('beq/wealth     = ',f15.8)") stats_bench(1)
	write(*,"('wealth/output  = ',f15.8)") stats_bench(2)
	write(*,"('top 1 percent  = ',f15.8)") stats_bench(3)
	write(*,"('top 10 percent = ',f15.8)") stats_bench(4)
	write(*,"('top 20 percent = ',f15.8)") stats_bench(5)
	write(*,"('top 50 percent = ',f15.8)") stats_bench(6)
	write(*,"('top 90 percent = ',f15.8)") stats_bench(7)
	write(*,"('std earnings   = ',f15.8)") stats_bench(8)
	
	CALL saving
	
	!----------------------------------------------
	! Wealth tax experiment (run only if RunExp=1) 
	!----------------------------------------------
	
	IF (RunExp==1) THEN
	
	!----------------------------------------------
	! counterfactual with wealth tax
	!----------------------------------------------
	
	print*, '----------------------------------------------------------'
	print*, ' New Policy                                               '
	print*, '----------------------------------------------------------'
	Q_exp=Q_bench
	Nbar_exp=Nbar_bench
	
	tauw0=tauw0_initial
	Gbar_ab(1)=0.0_DP
	Gbar_ab(2)=0.0_DP
	Gbar_exp=-1.0_DP
	DO WHILE (Gbar_exp<Gbar_bench)
		CALL solution(rr_exp,ww_exp,tauk_exp,tauw0,taul0, &
			& Gbar_exp,value_exp,v2_exp,aprime_exp, &
			& Q_exp,Nbar_exp, panel_exp_s,panel_exp_a,panel_exp_age,0,stats_exp)
			! by assigning (Q,Nbar) to have intent (INOUT), this means in the while loop
			! they always take the values from previous tax rate loop. 
		print*, '-------------------------------------------------'
		print*, 'tauw0', tauw0,' Gbar_exp',Gbar_exp
		print*, '-------------------------------------------------'
		Gbar_ab(1)=Gbar_ab(2)
		Gbar_ab(2)=Gbar_exp
		tauw0=tauw0+tauw0_inc
	END DO
	
	tauw0_upper=tauw0-tauw0_inc
	tauw0_lower=tauw0-2*tauw0_inc   
	
	! gbar is between gbar_ab(1) and gbar_ab(2), tauw0 is between tauw0-2*tauw0_inc and tauw0-tauw0
	! bisection to get the precise tauw0
	tauw0=(Gbar_ab(2)-Gbar_bench)/(Gbar_ab(2)-Gbar_ab(1))*(tauw0-2.0_DP*tauw0_inc)+ &
		& (Gbar_bench-Gbar_ab(1))/(Gbar_ab(2)-Gbar_ab(1))*(tauw0-tauw0_inc)
	
	iterGbar=1
	
	DO WHILE ( abs(Gbar_exp/Gbar_bench-1.0_DP)>0.001_DP)
		IF (iterGbar>1) THEN
	        tauw0=(tauw0_upper+tauw0_lower)/2.0_DP
	END IF
	
	CALL solution(rr_exp,ww_exp,tauk_exp,tauw0,taul0, &
        & Gbar_exp,value_exp,v2_exp,aprime_exp, &
        & Q_exp,Nbar_exp,panel_exp_s,panel_exp_a,panel_exp_age,0,stats_exp)
        
	IF (Gbar_exp>Gbar_bench) THEN
		tauw0_upper = tauw0
		Gbar_ab(2)=Gbar_exp
	ELSE
		tauw0_lower = tauw0
		Gbar_ab(1)=Gbar_exp
	END IF
	print*,'Gbar_bench,Gbar_exp,tauw0,tauw0_lower,tauw0_upper'
	print*,Gbar_bench,Gbar_exp,tauw0,tauw0_lower,tauw0_upper
	
	iterGbar=iterGbar+1
	
	END DO
	tauw_exp=tauw0

	CALL solution(rr_exp,ww_exp,tauk_exp,tauw_exp,taul0, &
		& Gbar_exp,value_exp,v2_exp,aprime_exp, &
		& Q_exp,Nbar_exp,panel_exp_s,panel_exp_a,panel_exp_age,1,stats_exp)
	           
	write(*,"('Gbar  = ',f15.8)") Gbar_exp
	write(*,"('rr    = ',f15.8)") rr_exp
	write(*,"('ww    = ',f15.8)") ww_exp
	write(*,"('Q  = ',f15.8)") Q_exp
	write(*,"('Nbar  = ',f15.8)") Nbar_exp

	write(*,"('beq/wealth     = ',f15.8)") stats_exp(1)
	write(*,"('wealth/output  = ',f15.8)") stats_exp(2)
	write(*,"('top 1 percent  = ',f15.8)") stats_exp(3)
	write(*,"('top 10 percent = ',f15.8)") stats_exp(4)
	write(*,"('top 20 percent = ',f15.8)") stats_exp(5)
	write(*,"('top 50 percent = ',f15.8)") stats_exp(6)
	write(*,"('top 90 percent = ',f15.8)") stats_exp(7)
	
	write(*,"('**************************************')")
	print*,'Final tauw=   ', tauw_exp
	print*,'Final Gbar=   ', Gbar_exp
	print*,'Output change ', ((Aprod*Q_exp**alpha*Nbar_exp**(1.0_DP-alpha))- &
		& (Aprod*Q_bench**alpha*Nbar_bench**(1.0_DP-alpha))) &
		& /(Aprod*Q_bench**alpha*Nbar_bench**(1.0_DP-alpha))

	CALL saving
	
	!----------------------------------------------
	! welfare calculation (AGE ADJUSTED)
	!----------------------------------------------
	
	ce_run=1
	
	! need to set value and v2, otherwise program sometimes overwrite value_bench
	! and that will cause incorrect calculation of welfare
	value=0.0_DP
	v2=0.0_DP
	CALL solution(rr_bench,ww_bench,tauk_bench,tauw_bench,taul0, &
		& Gbar_bench,value,v2,aprime_bench, &
		& Q_bench,Nbar_bench,panel_bench_s,panel_bench_a,panel_bench_age, &
		& 0,stats_bench)

	END IF
	
	call cpu_time(time2)
	write(*,*) 'program completed in:  ',time2-time1,'  seconds.'
	
	OPEN (UNIT=1, FILE='computingtime', STATUS='replace')
	WRITE (UNIT=1, FMT=*) time2-time1
	CLOSE (UNIT=1)
	
END PROGRAM
!====================================================================
SUBROUTINE initial
USE global
IMPLICIT NONE

	INTEGER(I4B)  :: si, zi, zz, mi, ei, ee, ai, alphai, i, age, row, ISEED
	REAL(DP) :: m, RAN1
	REAL(DP), DIMENSION(81)         :: pop
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: panel_z, panel_m, panel_alpha, panel_e, rn_rank, rn_index
	REAL(DP),DIMENSION(:),ALLOCATABLE :: rn
	
	!----------------------------------------------
	! life-cycle component
	!----------------------------------------------
		
	! Population Numbers from Bell and Miller (2002)
	
	pop(1)=	197316.0_DP
	pop(2)=	197141.0_DP
	pop(3)=	196959.0_DP
	pop(4)=	196770.0_DP
	pop(5)=	196580.0_DP
	pop(6)=	196392.0_DP
	pop(7)=	196205.0_DP
	pop(8)=	196019.0_DP
	pop(9)=	195830.0_DP
	pop(10)=195634.0_DP
	pop(11)=195429.0_DP
	pop(12)=195211.0_DP
	pop(13)=194982.0_DP
	pop(14)=194739.0_DP
	pop(15)=194482.0_DP
	pop(16)=194211.0_DP
	pop(17)=193924.0_DP
	pop(18)=193619.0_DP
	pop(19)=193294.0_DP
	pop(20)=192945.0_DP
	pop(21)=192571.0_DP
	pop(22)=192169.0_DP
	pop(23)=191736.0_DP
	pop(24)=191271.0_DP
	pop(25)=190774.0_DP
	pop(26)=190243.0_DP
	pop(27)=189673.0_DP
	pop(28)=189060.0_DP
	pop(29)=188402.0_DP
	pop(30)=187699.0_DP
	pop(31)=186944.0_DP
	pop(32)=186133.0_DP
	pop(33)=185258.0_DP
	pop(34)=184313.0_DP
	pop(35)=183290.0_DP
	pop(36)=182181.0_DP
	pop(37)=180976.0_DP
	pop(38)=179665.0_DP
	pop(39)=178238.0_DP
	pop(40)=176689.0_DP
	pop(41)=175009.0_DP
	pop(42)=173187.0_DP
	pop(43)=171214.0_DP
	pop(44)=169064.0_DP
	pop(45)=166714.0_DP
	pop(46)=164147.0_DP
	pop(47)=161343.0_DP
	pop(48)=158304.0_DP
	pop(49)=155048.0_DP
	pop(50)=151604.0_DP
	pop(51)=147990.0_DP
	pop(52)=144189.0_DP
	pop(53)=140180.0_DP
	pop(54)=135960.0_DP
	pop(55)=131532.0_DP
	pop(56)=126888.0_DP
	pop(57)=122012.0_DP
	pop(58)=116888.0_DP
	pop(59)=111506.0_DP
	pop(60)=105861.0_DP
	pop(61)=99957.0_DP
	pop(62)=93806.0_DP
	pop(63)=87434.0_DP
	pop(64)=80882.0_DP
	pop(65)=74204.0_DP
	pop(66)=67462.0_DP
	pop(67)=60721.0_DP
	pop(68)=54053.0_DP
	pop(69)=47533.0_DP
	pop(70)=41241.0_DP
	pop(71)=35259.0_DP
	pop(72)=29663.0_DP
	pop(73)=24522.0_DP
	pop(74)=19890.0_DP
	pop(75)=15805.0_DP
	pop(76)=12284.0_DP
	pop(77)=9331.0_DP
	pop(78)=6924.0_DP
	pop(79)=5016.0_DP
	pop(80)=3550.0_DP
	pop(81)=2454.0_DP

	! Survival probabilities: surv(i)=prob(alive in i+1|alive in i)
	FORALL (age=1:maxAge-1) pr_cond_surv(age)= pop(age+1)/pop(age)
	pr_cond_surv(maxAge)=0.0_DP
! 	print*, 'pr_cond_surv',pr_cond_surv
! 	pause
	
	! the compounded discount coefficient for age (general with or without retirement)
	Coeff_age=1.0_DP
	DO age=1,maxAge
		DO i=age,maxAge-1
			Coeff_age(age)=Coeff_age(age)+beta**(i-age+1)*product(pr_cond_surv(age:i))
		END DO
	END DO
! 	print*, Coeff_age
	
	!----------------------------------------------
	! labor transitional matrices
	!----------------------------------------------
	IF (ne<2) THEN
		egrid=0.0_DP
		Dist_e=1.0_DP
		pr_e=1.0_DP
	ELSE
		CALL tauchen(rho_e,sigma_e_eps,ne,egrid,pr_e,Dist_e,3.0_DP)
	END IF 
	egrid=exp(egrid)
 	print*, 'egrid', egrid
! 	print*, 'egrid divided by 50'
	
	IF (ne_alpha<2) THEN
		e_alpha=0.0_DP
		Dist_alpha=1.0_DP
		pr_alpha=1.0_DP
	ELSE
		CALL tauchen(rho_alpha,sigma_alpha_eps,ne_alpha,e_alpha,pr_alpha,Dist_alpha,3.0_DP)
	END IF   
! 	print*, 'e_alpha', e_alpha
	
	! Age-Efficiency Units (peak at 30, 50% increase)
	FORALL (age=1:maxAge) e_kappa(age)=(60.0_DP*REAL(age-1,DP)-REAL(age-1,DP)**2)/1800.0_DP
! 	print*, 'e_kappa', e_kappa
	
	!----------------------------------------------
	! investment intergenerational transition
	!----------------------------------------------
	IF (nz<2) THEN
		zgrid=0.0_DP
		Dist_z=1.0_DP
		pr_z=1.0_DP
	ELSE
		CALL tauchen(rho_z,sigma_z_eps,nz,zgrid,pr_z,Dist_z,tauchen_m_z)
  	END IF
	zgrid=exp(zgrid)
	print*, 'zgrid', zgrid
 	print*, 'pr_z(nz,:)', pr_z(nz,:)
	print*, 'Gz', Dist_z
! 	pause

	IF (nm<2) THEN
		mgrid=0.0_DP
		Dist_m=1.0_DP
		pr_m=1.0_DP
	ELSE
		CALL tauchen(rho_m,sigma_m_eps,nm,mgrid,pr_m,Dist_m,3.0_DP)
  	END IF
	mgrid=exp(mgrid)
	print*, 'mgrid', mgrid
! 	print*, pr_z
	print*, 'Gm', Dist_m
! 	pause

	!----------------------------------------------
	! panel size conditional on (age,z)
	!----------------------------------------------
	FORALL (age=1:maxAge) size_by_age(age)=INT((pop(age)/pop(1)*REAL(panelN,DP)))
	FORALL (zi=1:nz,age=1:maxAge) size_by_age_z(age,zi)=INT(Dist_z(zi)*REAL(size_by_age(age),DP))
	
	! the following deals with population accounting issues
	! edited (3/26/2014), change nz to max(2,nz) so it works for nz=1 case (no z shock)
	size_by_age_z(1,max(2,nz)/2)=size_by_age_z(1,max(2,nz)/2)+size_by_age(1)-sum(size_by_age_z(1,:))
	DO age=2,maxAge	
		DO zi=max(2,nz)/2,1,-1 ! first downward allocate extra people
			size_by_age_z(age,zi)=min(size_by_age_z(age-1,zi),size_by_age_z(age,zi)+size_by_age(age)-sum(size_by_age_z(age,:)))
			IF (size_by_age(age) .eq. sum(size_by_age_z(age,:))) THEN
				EXIT
			END IF
		END DO
		DO zi=max(2,nz)/2+1,nz ! if still not enough, upward allocate
			size_by_age_z(age,zi)=min(size_by_age_z(age-1,zi),size_by_age_z(age,zi)+size_by_age(age)-sum(size_by_age_z(age,:)))
			IF (size_by_age(age) .eq. sum(size_by_age_z(age,:))) THEN
				EXIT
			END IF
		END DO
	END DO
	panel_pop=SUM(size_by_age)
	
	DO zi=1,nz
		DO age=maxAge,2,-1
			IF (size_by_age_z(age,zi)>size_by_age_z(age-1,zi)) THEN
				print*, age,zi,size_by_age_z(age,zi),size_by_age_z(age-1,zi)
				pause
			END IF
		END DO
	END DO
	
	DO age=1,maxAge
		IF (sum(size_by_age_z(age,:)) /= size_by_age(age)) Then
			pause
		END IF
	END DO
	
	i=size_by_age(maxAge)
	DO age=2,maxAge
		i=i+size_by_age(age-1)-size_by_age(age)
	END DO
	print*, 'death count', i, 'born count', panelN

	i=sum(size_by_age_z(maxAge,:))
	DO age=2,maxAge
		i=i+sum(size_by_age_z(age-1,:))-sum(size_by_age_z(age,:))
	END DO
	print*, 'death count', i, 'born count', panelN
	
!	DO age=1,maxAge
!        print*, age,size_by_age_z(age,:)
!    END DO
    
	IF (i /= panelN) THEN
		print*, 'error: death count not equal to born count', i, panelN
		pause
	END IF
! 	pause
	
	!----------------------------------------------
	! Calculating Nbar (exogenous, labor supply inelastic)
	! Agents supply labor up to (including) age ReAge-1 (no work >= ReAge)
	!----------------------------------------------
	Nbar_bench=0.0_DP
	DO ei=1,ne
		DO alphai=1,ne_alpha
			DO age=1,min(maxAge,ReAge-1)
				Nbar_bench=Nbar_bench+egrid(ei)*exp(e_alpha(alphai)+e_kappa(age))* &
					& Dist_e(ei)*Dist_alpha(alphai)*size_by_age(age)
			END DO
		END DO
	END DO
	Nbar_bench=Nbar_bench/REAL(panel_pop,DP)
	print*, 'Nbar', Nbar_bench
	
	!----------------------------------------------
	! Retirement incomes
	!----------------------------------------------
	! sum(log(egrid))/REAL(ne) is the median of egrid. 
	! the way we can calulate this is becuase we know log(egrid) is symmetric around 0. 
	IF (ReAge<=maxAge) THEN
		e_med=exp(sum(log(egrid))/REAL(ne,DP))
		z_med=exp(sum(log(zgrid))/REAL(nz,DP))
		DO alphai=1,ne_alpha
			n_alpha(alphai)=sum(exp(e_alpha(alphai)+e_kappa(1:min(maxAge,ReAge)-1))) &
				& /REAL(min(maxAge,ReAge)-1,DP)*e_med 
				! put min(maxAge,ReAge) to avoid out of bound if no retirement
		END DO
		nr=sum(n_alpha*Dist_alpha) ! average n_alpha
		y_alpha=n_alpha/nr
! 	print*, 'n_alpha', n_alpha
! 	print*, 'median e', exp(sum(log(egrid))/REAL(ne,DP))
! 	print*, 'nr', nr
! 	print*, 'y_alpha', y_alpha
! 	pause
	END IF
	
	!----------------------------------------------
	! Expanding asset grid 
	!----------------------------------------------
	m=(amax-amin)**(1.0_DP/a_theta)/REAL(na-1,DP)
	DO ai=1,na
		agrid(ai)=REAL(ai-1,DP)*m
	END DO
	agrid=amin+agrid**a_theta
	agrid(1)=1.0e-3_dp
! 	print*,'agrid',agrid

	!----------------------------------------------
	! Set state variable s
	!----------------------------------------------
	si=0
	DO zi=1,nz
		DO mi=1,nm
			DO alphai=1,ne_alpha
				DO ei=1,ne
					si=si+1
					sz(si)=zi
					sm(si)=mi
					salpha(si)=alphai
					se(si)=ei
					ssfun(zi,mi,alphai,ei)=si
				END DO
			END DO
		END DO
	END DO

	!----------------------------------------------
	! Initilize the simulation panel for newborns
	!----------------------------------------------

	allocate( panel_z(panelN), panel_m(panelN), panel_alpha(panelN), panel_e(panelN) ) 
	allocate( rn(panelN), rn_index(panelN), rn_rank(panelN) )

	panel_z=max(2,nz)/2 ! if nz=1,panel_z=1, if nz>1, take floor(nz/2)
	FORALL (i=1:size_by_age_z(1,1)) panel_z(i)=1
	DO zi=2,nz
		DO i=sum(size_by_age_z(1,1:zi-1))+1,sum(size_by_age_z(1,1:zi))
			panel_z(i)=zi
		END DO
	END DO
	
	ISEED=3626
	    		
	panel_m=1
	IF (nm>1) THEN
		DO row=1,panelN
			rn(row)=RAN1(ISEED)
		END DO
		CALL indexx(panelN,rn,rn_index)
		FORALL (i=1:panelN) rn_rank(rn_index(i))=i		    		
	    row=1
	    DO mi=1,nm
		    	i=1
		    	DO WHILE (i<INT(Dist_m(mi)*REAL(panelN,DP)) .and. row<=panelN)
		    		panel_m(rn_rank(row))=mi
		        	i=i+1
		        	row=row+1
		    	END DO
	    END DO
	END IF
		
	panel_alpha=1
	IF (ne_alpha>1) THEN
		DO row=1,panelN
			rn(row)=RAN1(ISEED)
		END DO
		CALL indexx(panelN,rn,rn_index)
		FORALL (i=1:panelN) rn_rank(rn_index(i))=i		    		
	    row=1
	    DO alphai=1,ne_alpha
		    	i=1
		    	DO WHILE (i<INT(Dist_alpha(alphai)*REAL(panelN,DP)) .and. row<=panelN)
		    		panel_alpha(rn_rank(row))=alphai
		        	i=i+1
		        	row=row+1
		    	END DO
	    END DO
	END IF
	
	panel_e=1
	IF (ne>1) THEN
		DO row=1,panelN
			rn(row)=RAN1(ISEED)
		END DO
		CALL indexx(panelN,rn,rn_index)
		FORALL (i=1:panelN) rn_rank(rn_index(i))=i		    		
	    row=1
	    DO ei=1,ne
		    	i=1
		    	DO WHILE (i<INT(Dist_e(ei)*REAL(panelN,DP)) .and. row<=panelN)
		    		panel_e(rn_rank(row))=ei
		        	i=i+1
		        	row=row+1
		    	END DO
	    END DO
	END IF
	
	DO row=1,panelN
		panel_s_init(row)=ssfun(panel_z(row),panel_m(row),panel_alpha(row),panel_e(row))
	END DO
		
	deallocate ( panel_z, panel_m, panel_alpha, panel_e ) 
	deallocate ( rn, rn_index, rn_rank )
	
	! intermediate goods that agents produce using capital, x=z*a
	FORALL (zi=1:nz,mi=1:nm,ai=1:na) x_product(zi,mi,ai)=zgrid(zi)*mgrid(mi)*agrid(ai)
	
	! effective labor
	DO age=1,maxAge
		x_product_rho(:,:,:,age)=x_product**rho
		
		IF (age>=ReAge) THEN
			DO alphai=1,ne_alpha
				eff_labor(alphai,:,age)=e_med*coeff_r*y_alpha(alphai)
			END DO
!			DO ai=1,na
!				x_product_rho(:,:,ai,age)=(z_med*agrid(ai))**rho
!			END DO
		ELSE
			FORALL (ei=1:ne,alphai=1:ne_alpha) eff_labor(alphai,ei,age)=egrid(ei)*exp(e_alpha(alphai)+e_kappa(age))
!			x_product_rho(:,:,:,age)=x_product**rho
		END IF
	END DO
	
	OPEN (UNIT=1, FILE='agrid', STATUS='replace')
	DO ai=1,na
		WRITE (UNIT=1, FMT=*) agrid(ai)
	END DO
	CLOSE (UNIT=1)

	OPEN (UNIT=1, FILE='zgrid', STATUS='replace')
	OPEN (UNIT=2, FILE='Dist_z', STATUS='replace')
	DO zi=1,nz
		WRITE (UNIT=1, FMT=*) zgrid(zi)
		WRITE (UNIT=2, FMT=*) Dist_z(zi)
	END DO
	CLOSE (UNIT=1)
	CLOSE (UNIT=2)
	
	OPEN (UNIT=1, FILE='pr_z', STATUS='replace')
	DO zi=1,nz
		DO zz=1,nz
			WRITE (UNIT=1, FMT=*) pr_z(zi,zz)
		END DO
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='mgrid', STATUS='replace')
	DO mi=1,nm
		WRITE (UNIT=1, FMT=*) mgrid(mi)
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='egrid', STATUS='replace')
	OPEN (UNIT=2, FILE='Dist_e', STATUS='replace')
	DO ei=1,ne
		WRITE (UNIT=1, FMT=*) egrid(ei)
		WRITE (UNIT=2, FMT=*) Dist_e(ei)
	END DO
	CLOSE (UNIT=1)
	CLOSE (UNIT=2)
	
	OPEN (UNIT=1, FILE='pr_e', STATUS='replace')
	DO ei=1,ne
		DO ee=1,ne
			WRITE (UNIT=1, FMT=*) pr_e(ei,ee)
		END DO
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='e_alpha', STATUS='replace')
	DO alphai=1,ne_alpha
		WRITE (UNIT=1, FMT=*) e_alpha(alphai)
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='e_kappa', STATUS='replace')
	DO age=1,maxAge
		WRITE (UNIT=1, FMT=*) e_kappa(age)
	END DO
	CLOSE (UNIT=1)
 
	OPEN (UNIT=1, FILE='size_by_age', STATUS='replace')
	DO age=1,maxAge
		WRITE (UNIT=1, FMT=*) size_by_age(age)
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='size_by_age_z', STATUS='replace')
	DO age=1,maxAge
		DO zi=1,nz
			WRITE (UNIT=1, FMT=*) size_by_age_z(age,zi)
		END DO
	END DO
	CLOSE (UNIT=1)
	  
	OPEN (UNIT=1, FILE='param_list', STATUS='replace')
	WRITE (UNIT=1,FMT=*) 'maxAge              ', maxAge
	WRITE (UNIT=1,FMT=*) 'ReAge               ', ReAge
	WRITE (UNIT=1,FMT=*) 'beta                ', beta
	WRITE (UNIT=1,FMT=*) 'rho_z               ', rho_z
	WRITE (UNIT=1,FMT=*) 'sigma_z_eps         ', sigma_z_eps
	WRITE (UNIT=1,FMT=*) 'rho_m               ', rho_m
	WRITE (UNIT=1,FMT=*) 'sigma_m_eps         ', sigma_m_eps
	WRITE (UNIT=1,FMT=*) 'a_theta             ', a_theta
	WRITE (UNIT=1,FMT=*) 'GE                  ', GE     
	WRITE (UNIT=1,FMT=*) 'LIP                 ', LIP
	WRITE (UNIT=1,FMT=*) 'max_Q_iter          ', max_Q_iter
	WRITE (UNIT=1,FMT=*) 'maxiter_panel       ', maxiter_panel
	WRITE (UNIT=1,FMT=*) 'panelN              ', panelN
	WRITE (UNIT=1,FMT=*) 'panel_pop           ', panel_pop
	WRITE (UNIT=1,FMT=*) 'tauw0_initial       ', tauw0_initial
	WRITE (UNIT=1,FMT=*) 'tauw0_inc           ', tauw0_inc
	WRITE (UNIT=1,FMT=*) 'taul0               ', taul0  
	WRITE (UNIT=1,FMT=*) 'alpha               ', alpha  
	WRITE (UNIT=1,FMT=*) 'Aprod               ', Aprod        
	WRITE (UNIT=1,FMT=*) 'delta               ', delta 
	WRITE (UNIT=1,FMT=*) 'hp                  ', hp
	WRITE (UNIT=1,FMT=*) 'rho                 ', rho
	WRITE (UNIT=1,FMT=*) 'rho_e               ', rho_e
	WRITE (UNIT=1,FMT=*) 'sigma_e_eps         ', sigma_e_eps
	WRITE (UNIT=1,FMT=*) 'rho_alpha           ', rho_alpha
	WRITE (UNIT=1,FMT=*) 'sigma_alpha_eps     ', sigma_alpha_eps        
	WRITE (UNIT=1,FMT=*) 'tolQ                ', tolQ         
	CLOSE (UNIT=1)

END SUBROUTINE
!====================================================================
SUBROUTINE solution(rr,ww,tauk,tauw,taul,Gbar,value,v2,aprime,Q,Nbar,panel_s,&
	& panel_a,panel_age,run_stats,stats_out)
USE global
IMPLICIT NONE

	REAL(DP), INTENT(INOUT)  :: Q, Nbar
	REAL(DP), INTENT(IN)     :: tauk, tauw, taul
	INTEGER(I4B), INTENT(IN) :: run_stats
	REAL(DP), DIMENSION(ns,na,maxAge), INTENT(OUT)  :: aprime, value, v2
	REAL(DP),DIMENSION(8),INTENT(OUT) :: stats_out
	REAL(DP), INTENT(OUT)   :: Gbar, rr, ww
	
	REAL(DP) :: rr_new, ww_new, tmp, tmp1, vprime
		
	REAL(DP) :: utility, Q_new, splint, distQ
	INTEGER(I4B)  :: zi, mi, mj, mm, alphai, alphaj, ei, ej, age, ai, i, j, ee, api, ss, si
	INTEGER(I4B)  :: iterQ, maxiterQ
	
	INTEGER(I4B),DIMENSION(2) :: xlim
	REAL(DP),    DIMENSION(2) :: ab

	REAL(DP), DIMENSION(maxiter_panel) :: panel_Q_seq, panel_wealth_seq, panel_xmax_seq, panel_N_seq
	INTEGER(I4B),  DIMENSION(panelN,maxAge) :: panel_s, panel_age
	REAL(DP),      DIMENSION(panelN,maxAge) :: panel_a, panel_x
	INTEGER(I4B),  DIMENSION(panelN) :: panel_alpha_end, panel_z0
	REAL(DP), DIMENSION(panelN)      :: panel_a_end, panel_z_end
	INTEGER(I4B) :: panel_z_ind, panel_m_ind, panel_alpha_ind, panel_e_ind
	
	REAL(DP) :: RAN1
	INTEGER(I4B)  :: row, ISEED, iter_panel, iter, final_run
	INTEGER(I4B), SAVE :: borncount=0
	REAL(DP) :: maxQ,  minQ, rtmp
	
	REAL(DP), DIMENSION(ns,na,maxAge)      :: consumption, asset, Bset, budget
	INTEGER(I4B),  DIMENSION(maxAge,nz)    :: livecount_age_z
	REAL(DP),     DIMENSION(panelN) :: rn_z, rn_age
	INTEGER(I4B), DIMENSION(panelN) :: rn_index, rn_rank
	REAL(DP) :: var1, var2, rsq, Q0, value0, value1
	
	REAL(DP),DIMENSION(numiter_fin_panel)   :: BeqWealthRatioRep, WealthOutputRatioRep
	REAL(DP),DIMENSION(numiter_fin_panel)   :: QRep, GbarRep, Std_logLaborY
	REAL(DP),DIMENSION(numiter_fin_panel,5) :: Top_wealth_Rep
	
	REAL(DP) :: cons_tilde, mybrent_EGa, coef1, coef2, coef3, power2
	REAL(DP) :: ax, bx, cx, tmp_k,tmp_w
	REAL(DP), EXTERNAL :: func_EGa
	
	REAL(DP) :: rr_coef, ww_coef 
	
	! indicator for whether this is the last run to compute stats
	final_run=0
	
	! initialize all the output variables
	value=0.0_DP
	aprime=0.0_DP
		
	distQ=1.0_DP
	iterQ=1
	maxiterQ=1+GE*max_Q_iter ! if GE=0, run once given (r,w). if GE=1, run up to max_Q_iter times
	
	! rearrange budget constraint into: coef1*a+coef2*a**power2+coef3=0
	! coef1 and power2 only depend on parameters. calculate outside of loop
	coef1=(1.0_DP-tauw)*(1.0_DP-(1.0_DP-tauk)*delta)
	power2=rho
			
	! Given input (Q,Nbar), calculate (r,w)
	! hold part of (r,w) that doesn't depend on Q in memory (rr_coef,ww_coef)
	rr_coef= alpha*Aprod*Nbar**(1.0_DP-alpha)
	ww_coef= (1.0_DP-alpha)*Aprod*Nbar**(-alpha)
	
	! (r, w)=(rr_coef*Q**(alpha-rho), ww_coef*Q**alpha)
	Q_new=Q
	rr_new = rr_coef*Q**(alpha-rho)
	ww_new = ww_coef*Q**alpha
	
	print*, rr_new, ww_new
	
 	maxQ  = agrid(na)
! 	maxQ  = 50.0_DP
	minQ  = 0.0_DP

	v2=0.0_DP ! initilize v2=0, if not updated, linear interpolation is implied.
		
	DO WHILE ( final_run==0 )
	
		IF (ce_run==1) THEN
			iterQ=maxiterQ+1
			aprime=aprime_bench
		END IF
		
		IF (distQ<tolQ .or. iterQ>maxiterQ) THEN
			final_run=1 ! for final run, (r,w) not updated
		ELSE
			! if not final run, given (r,w), calculate values and decisions
			Q0 = Q/5.0_DP ! assets given to newborns from moon	
			Q  = Q_new
			rr = rr_new
			ww = ww_new

			! total resources (=max possible consumption), c+a'<=budget
			DO si=1,ns
				zi=sz(si); mi=sm(si); alphai=salpha(si); ei=se(si)
				DO age=1,maxAge
					budget(si,:,age)=(1.0_DP-tauw)*(agrid+(1.0_DP-tauk)*(rr*x_product_rho(zi,mi,:,age))-delta*agrid)+ &
						& (1.0_DP-taul)*ww*eff_labor(alphai,ei,age)
				END DO
			END DO
			
			! value of last period if there is no retirement.
			! value is equal to per-period utility for last period of life
			aprime(:,:,maxAge)= (chi*beta*budget(:,:,maxAge))/(1.0_DP+chi*beta)
			consumption(:,:,maxAge)=budget(:,:,maxAge)-aprime(:,:,maxAge)			
			value(:,:,maxAge)=utility(consumption(:,:,maxAge))+beta*chi*log(aprime(:,:,maxAge)+hp)

			! starting from the next period to last, use backward induction to solve for values.
			DO age=maxAge-1,1,-1
!				print*, 'age', age
				IF (LIP==0) THEN	
					DO si=1,ns
						CALL spline(agrid,value(si,:,age+1),na,v2(si,:,age+1))
					END DO
				END IF
				
				DO api=1,na
					Bset(:,api,age)=beta*(1.0_DP-pr_cond_surv(age))*chi/(agrid(api)+hp) ! bequest term on RHS of Euler eq
					DO si=1,ns
						zi=sz(si); mi=sm(si); alphai=salpha(si); ei=se(si);
						DO mm=1,nm
!							IF (age<ReAge) THEN 
								tmp=(1.0_DP-tauw)*(1.0_DP+(1.0_DP-tauk)*(rho*rr*((zgrid(zi)*mgrid(mm))**rho)*(agrid(api)**(rho-1.0_DP))-delta))
!							ELSE
!								tmp=(1.0_DP-tauw)*(1.0_DP+(1.0_DP-tauk)*(rho*rr*(z_med**rho)*(agrid(api)**(rho-1.0_DP))-delta))				
!							END IF
							DO ee=1,ne
								ss=ssfun(zi,mm,alphai,ee)
								Bset(si,api,age)=Bset(si,api,age)+beta*pr_cond_surv(age)*tmp* &
									& (1.0_DP/(consumption(ss,api,age+1)+hp))*pr_e(ei,ee)*pr_m(mi,mm)
							END DO
						END DO
						cons_tilde=1.0_DP/Bset(si,api,age)-hp
						! rearrange budget constraint into: coef1*a+coef2*a**power2+coef3=0
						! coef1 and power2 already calculated outside loop. only need to compute coef2 and coef3
						coef2=(1.0_DP-tauw)*(1.0_DP-tauk)*rr*((zgrid(zi)*mgrid(mi))**rho)
						coef3=(1.0_DP-taul)*ww*eff_labor(alphai,ei,age)-cons_tilde-agrid(api)
!						IF (age>=ReAge) THEN
!							coef2=(1.0_DP-tauw)*(1.0_DP-tauk)*rr*(z_med**rho)
!						END IF
						
						IF (coef3>0.0_dp) THEN
							tmp=coef1*agrid(10)+coef2*agrid(10)**power2+coef3
							asset(si,api,age)=-agrid(10)*coef3/(tmp-coef3)
!							ax=-agrid(na)
!							cx=10.0_DP*agrid(na)
!							bx=agrid(api)
!							tmp=mybrent_EGa(ax,bx,cx,func_EGa,tolbrent,asset(si,api,age),coef1,coef2,coef3,power2)
!							print*, 'true sol', asset(si,api,age)
!							pause
						ELSE
						! solve for the asset position a that satisfies the budget constraint
						ax=0.0_dp
						cx=10.0_DP*agrid(na)
						bx=agrid(api)
						tmp=mybrent_EGa(ax,bx,cx,func_EGa,tolbrent,asset(si,api,age),coef1,coef2,coef3,power2)
						END IF
					END DO
				END DO
				
				DO ai=1,na
					DO si=1,ns
						zi=sz(si); mi=sm(si); alphai=salpha(si); ei=se(si)
						CALL bisec(asset(si,:,age),na,agrid(ai),xlim,ab)
						aprime(si,ai,age)=max(0.0_DP,ab(1)*agrid(xlim(1))+ab(2)*agrid(xlim(2)))
						consumption(si,ai,age)=budget(si,ai,age)-aprime(si,ai,age)
						vprime=0.0_DP
						DO ee=1,ne
							DO mm=1,nm
								ss=ssfun(zi,mm,alphai,ee)
								vprime=vprime+pr_e(ei,ee)*pr_m(mi,mm)* &
									& splint(agrid,value(ss,:,age+1),v2(ss,:,age+1),na,aprime(si,ai,age))
							END DO
						END DO
						value(si,ai,age)=utility(consumption(si,ai,age))+ &
							& beta*(pr_cond_surv(age)*vprime+(1.0_DP-pr_cond_surv(age))*chi*log(aprime(si,ai,age)+hp))
					END DO
				END DO
			END DO
		END IF
!		print*, 'min aprime', minval(aprime),minloc(aprime),minval(consumption),minloc(consumption)
!		print*, 'max aprime', maxval(aprime),maxloc(aprime),maxval(consumption),maxloc(consumption)

		!=============================================================
		! run simulation panel

		! initialize the seed for random numbers
		ISEED=314159265
			
		! initialize the first cohort with invariant distribution for (z,alpha,e)
		! start newborns with age 1 and assets Q
		
		panel_age=0 ! reset all age back to 0 before simulation
		panel_a=0.0_DP
		panel_s=0
		
		panel_age(:,1)=1
		panel_a(:,1)=Q0
		panel_s(:,1)=panel_s_init
		
		iter_panel=1
		DO WHILE (iter_panel<=maxiter_panel)
			! make random draw from a normal distribution N(0,sigma_z_eps)
			! modified from numerical recipes GASDEV subroutine
			DO row=1,panelN
				DO 
					var1=RAN1(ISEED)
					var2=RAN1(ISEED)
					var1=2.0_dp*var1-1.0_dp
					var2=2.0_dp*var2-1.0_dp
					rsq=var1**2+var2**2
					IF (rsq > 0.0 .and. rsq < 1.0) EXIT
				END DO
				rsq=sqrt(-2.0_dp*log(rsq)/rsq)
				rn_z(row)=sigma_z_eps*var1*rsq ! normal draw for z_eps~N(0,sigma_z_eps)
			END DO
			
			livecount_age_z=0
			borncount=0 ! counting number of newborns (need panelN of them)
			DO row=1,panelN
							
				IF (panel_age(row,maxAge)>0) THEN
				
					panel_z_ind=sz(panel_s(row,maxAge))
					panel_alpha_ind=salpha(panel_s(row,maxAge))
					
					borncount=borncount+1 ! all living agents die after maxAge
	                
					! make alpha draw from AR(1) for offspring
					tmp=RAN1(ISEED)
					alphaj=1
					rtmp=pr_alpha(panel_alpha_ind,1)
					DO WHILE (tmp>rtmp .and. alphaj<ne_alpha)
						alphaj=alphaj+1
						rtmp=SUM(pr_alpha(panel_alpha_ind,1:alphaj))
					END DO
					panel_alpha_end(borncount)=alphaj
					
					! make z draw from AR(1) for offspring
					panel_z_end(borncount)=exp(rho_z*log(zgrid(panel_z_ind))+rn_z(row))
	
					IF (chi==0.0_DP) THEN
						panel_a_end(borncount)=0.0_DP !without beq motive, leave zero assets.
					ELSE
						IF (panel_a(row,maxAge)<agrid(na)) THEN
							CALL bisec(agrid,na,panel_a(row,maxAge),xlim,ab)
						ELSE
							xlim(1)=na-1
							xlim(2)=na
							ab(1)=(agrid(na)-panel_a(row,maxAge))/(agrid(na)-agrid(na-1))
							ab(2)=(panel_a(row,maxAge)-agrid(na-1))/(agrid(na)-agrid(na-1))
						END IF
						panel_a_end(borncount)=max(0.0_dp,sum(ab*aprime(panel_s(row,maxAge),xlim,maxAge)))
					END IF
				END IF
				
			END DO
			
			DO age=maxAge-1,1,-1
			
	    		! survival counts by z, make draw for survival shocks
	    		! rank the shocks to determine whether this agent lives (conditional on age)
	    		DO row=1,panelN
	    			rn_age(row)=RAN1(ISEED)
	    		END DO
	    		CALL indexx(panelN,rn_age,rn_index)
	    		FORALL (i=1:panelN) rn_rank(rn_index(i))=i
    		
		    	DO i=1,panelN
		    	
		    		row=rn_rank(i) ! running through agents in the order of their survival draw
		    					   ! keep agents alive in draw order until counts equal to size_by_age_z(z,age+1)
		    		
		    		IF (panel_age(row,age)>0) THEN
		    			panel_z_ind=sz(panel_s(row,age))
		    			panel_m_ind=sm(panel_s(row,age))
						panel_e_ind=se(panel_s(row,age))
						panel_alpha_ind=salpha(panel_s(row,age))
		    			! make asset choice decisions
						IF (panel_a(row,age)<agrid(na)) THEN
							CALL bisec(agrid,na,panel_a(row,age),xlim,ab)
						ELSE
							xlim(1)=na-1
							xlim(2)=na
							ab(1)=(agrid(na)-panel_a(row,age))/(agrid(na)-agrid(na-1))
							ab(2)=(panel_a(row,age)-agrid(na-1))/(agrid(na)-agrid(na-1))
						END IF
						panel_a(row,age+1)=max(0.0_dp,sum(ab*aprime(panel_s(row,age),xlim,age)))
					
						IF (livecount_age_z(age+1,panel_z_ind)<size_by_age_z(age+1,panel_z_ind)) THEN
							livecount_age_z(age+1,panel_z_ind)=livecount_age_z(age+1,panel_z_ind)+1
	 						! make future draw for e
	 						tmp=RAN1(ISEED)
	 						ee=1
	 						rtmp=pr_e(panel_e_ind,1)
	 						DO WHILE (tmp>rtmp .and. ee<ne)
	 							ee=ee+1
	 							rtmp=SUM(pr_e(panel_e_ind,1:ee))
	 						END DO
	 						! make future draw for m
	 						tmp=RAN1(ISEED)
	 						mm=1
	 						rtmp=pr_m(panel_m_ind,1)
	 						DO WHILE (tmp>rtmp .and. mm<nm)
	 							mm=mm+1
	 							rtmp=SUM(pr_m(panel_m_ind,1:mm))
	 						END DO
							panel_s(row,age+1)=ssfun(panel_z_ind,mm,panel_alpha_ind,ee)
							panel_age(row,age+1)=age+1			
						ELSE
							borncount=borncount+1
							panel_age(row,age+1)=0 ! if not alive, mark 0 for panel_age
							! make z draw from AR(1) for offspring
		 					panel_z_end(borncount)=exp(rho_z*log(zgrid(panel_z_ind))+rn_z(row))
							! make alpha draw from AR(1) for offspring
							tmp=RAN1(ISEED)
							alphaj=1
							rtmp=pr_alpha(panel_alpha_ind,1)
							DO WHILE (tmp>rtmp .and. alphaj<ne_alpha)
								alphaj=alphaj+1
								rtmp=SUM(pr_alpha(panel_alpha_ind,1:alphaj))
							END DO					
							panel_alpha_end(borncount)=alphaj
							! pass wealth to offspring
							panel_a_end(borncount)=panel_a(row,age+1)
						END IF
					ELSE
						panel_age(row,age+1)=0 !! IMPORTANT, SO NOT TO TAKE VALUE FROM PREIVOUS COHORT! 
					END IF
				END DO
			END DO
					
			DO WHILE (borncount<panelN) ! creating offspring from moon
				borncount=borncount+1
				! make z draw from AR(1) with mean 0 last period		
				panel_z_end(borncount)=exp(rn_z(borncount))
				! make alpha draw from inv dist
				tmp1=RAN1(ISEED)
				alphaj=1
				rtmp=Dist_alpha(1)
				DO WHILE (tmp1>rtmp .and. alphaj<ne_alpha)
					alphaj=alphaj+1
					rtmp=SUM(Dist_alpha(1:alphaj))
				END DO
				panel_alpha_end(borncount)=alphaj
				! give Q0 as initial a
				panel_a_end(borncount)=Q0
			END DO
		
			! rank panel_z_end (real) into panel_z0 (INTEGER(I4B))
			! replace panel_z into z categories according to size_by_age_z
	 		CALL indexx(panelN,panel_z_end,panel_z0)
	 		FORALL (row=1:panelN,panel_z0(row)<=size_by_age_z(1,1)) panel_z0(row)=1
	 		DO row=1,panelN
	 		    DO zi=2,nz
	 		        IF (panel_z0(row)>sum(size_by_age_z(1,1:zi-1)) .and. &
						& panel_z0(row)<=sum(size_by_age_z(1,1:zi))) THEN
	            		panel_z0(row)=zi
	            		EXIT
            		END IF
        		END DO
			END DO
			
			DO row=1,panelN
				tmp=RAN1(ISEED)
				ej=1
				rtmp=Dist_e(1)
				DO WHILE (tmp>rtmp .and. ej<ne)
					ej=ej+1
					rtmp=SUM(Dist_e(1:ej))
				END DO
				
				tmp=RAN1(ISEED)
				mj=1
				rtmp=Dist_m(1)
				DO WHILE (tmp>rtmp .and. mj<nm)
					mj=mj+1
					rtmp=SUM(Dist_m(1:mj))
				END DO
				panel_s(row,1)=ssfun(panel_z0(row),mj,panel_alpha_end(row),ej)
			END DO
			
			panel_a(:,1)=panel_a_end
			panel_age(:,1)=1
	
			i=0
			panel_Q_seq(iter_panel)=0.0_DP
			panel_N_seq(iter_panel)=0.0_DP
			panel_x=0.0_dp
			DO row=1,panelN
				DO age=1,maxAge
					IF (panel_age(row,age)>0) THEN
						i=i+1
						panel_z_ind=sz(panel_s(row,age))
						panel_m_ind=sm(panel_s(row,age))
						panel_e_ind=se(panel_s(row,age))
						panel_alpha_ind=salpha(panel_s(row,age))
						IF (panel_a(row,age)>0.0_dp) THEN ! this if statement avoid program crashing when panel_a=0
!							IF (age<ReAge) THEN
								panel_Q_seq(iter_panel)=panel_Q_seq(iter_panel)+ &
									& (zgrid(panel_z_ind)*mgrid(panel_m_ind)*panel_a(row,age))**rho
								panel_x(row,age)=zgrid(panel_z_ind)*mgrid(panel_m_ind)*panel_a(row,age)
!							ELSE
!								panel_Q_seq(iter_panel)=panel_Q_seq(iter_panel)+(z_med*panel_a(row,age))**rho
!							END IF
						END IF
						
						IF (age<ReAge) THEN
							panel_N_seq(iter_panel)=panel_N_seq(iter_panel)+ &
								& egrid(panel_e_ind)*exp(e_alpha(panel_alpha_ind)+e_kappa(age))
						END IF
					END IF
				END DO
			END DO
			panel_Q_seq(iter_panel)=(panel_Q_seq(iter_panel)**(1.0_DP/rho))/REAL(i,DP)
			panel_wealth_seq(iter_panel)=sum(panel_a*real(min(1,panel_age),dp))/REAL(i,DP)
			panel_xmax_seq(iter_panel)=maxval(panel_x)
			panel_N_seq(iter_panel)=panel_N_seq(iter_panel)/REAL(i,DP)
! 			print*,iter_panel,panel_Q_seq(iter_panel),i
	
			IF (iter_panel>=maxiter_panel-numiter_fin_panel+1) THEN ! recording the last numiter_fin_panel rounds in simulation
				iter=iter_panel-(maxiter_panel-numiter_fin_panel)
				IF (ce_run==0) THEN
					QRep(iter)=panel_Q_seq(iter_panel)
					IF (final_run==1) THEN
						GbarRep(iter)=0; i=0
						DO row=1,panelN
							DO age=1,maxAge
								IF (panel_age(row,age)>0) THEN
									i=i+1
									panel_z_ind=sz(panel_s(row,age))
									panel_m_ind=sm(panel_s(row,age))
									panel_alpha_ind=salpha(panel_s(row,age))
									panel_e_ind=se(panel_s(row,age))
									! calculate tax base, tmp_k for tauk, and tmp_w for tauw
									IF (panel_a(row,age)>0.0_dp) THEN
!										IF (age<ReAge) THEN
											tmp_k=rr*((zgrid(panel_z_ind)*mgrid(panel_m_ind)*panel_a(row,age))**rho)-delta*panel_a(row,age)
!										ELSE
!											tmp_k=rr*((z_med*panel_a(row,age))**rho)-delta*panel_a(row,age)
!										END IF
										tmp_w=panel_a(row,age)+(1.0_DP-tauk)*tmp_k
									ELSE
										tmp_k=0.0_dp
										tmp_w=0.0_dp
									END IF
									GbarRep(iter)=GbarRep(iter)+tauk*tmp_k+tauw*tmp_w+ &
										& taul*ww*eff_labor(panel_alpha_ind,panel_e_ind,age)
								END IF
							END DO
						END DO
						GbarRep(iter)=GbarRep(iter)/REAL(i,DP)
						IF (run_stats==1) THEN
							CALL moment_panel(rr,ww,Q,Nbar,panel_s,panel_a,panel_age,stats_out)
							BeqWealthRatioRep(iter)=stats_out(1)
							WealthOutputRatioRep(iter)=stats_out(2)
							Top_wealth_Rep(iter,1:5)=stats_out(3:7)
							Std_logLaborY(iter)=stats_out(8)
						END IF						
					END IF
				ELSE
					CE_fin(iter)=0.0_DP
					CE_pos_panel(iter)=0.0_DP
					CE_newborn(iter)=0.0_DP
					i=0; j=0
					DO row=1,panelN
						DO age=1,maxAge
							IF (panel_age(row,age)>0) THEN
								i=i+1
								value0=splint(agrid,value_bench(panel_s(row,age),:,age), &
									& v2_bench(panel_s(row,age),:,age),na,panel_a(row,age))
								value1=splint(agrid,value_exp(panel_s(row,age),:,age), &
									& v2_bench(panel_s(row,age),:,age),na,panel_a(row,age))
								tmp=(exp((value1-value0)/Coeff_age(age))-1.0_DP)*100.0_DP
								CE_fin(iter)=CE_fin(iter)+tmp

								IF (tmp>0.0_DP) THEN
									CE_pos_panel(iter)=CE_pos_panel(iter)+1.0_DP
								END IF
								
								IF (panel_age(row,age)==1) THEN
									j=j+1
									CE_newborn(iter)=CE_newborn(iter)+tmp
								END IF
							END IF
						END DO
					END DO
					CE_fin(iter)=CE_fin(iter)/REAL(i,DP)
					CE_pos_panel(iter)=CE_pos_panel(iter)/REAL(i,DP)
					CE_newborn(iter)=CE_newborn(iter)/REAL(j,DP)	
				END IF
			END IF
			iter_panel=iter_panel+1
		END DO
		
		Q_new=sum(QRep)/REAL(numiter_fin_panel,DP)
		Gbar=sum(GbarRep)/REAL(numiter_fin_panel,DP)
			
		! end simulation panel
		!=============================================================
		
		maxQ=min(maxQ,max(Q_new,Q))
		minQ=max(minQ,min(Q,Q_new))
		
		distQ=ABS(1.0_DP-Q_new/Q)
		iterQ=iterQ+1
		
		!=======================================================
		! writtin by burhan, expand the bounds when not converging and hitting the bounds
		IF ((ABS(maxQ-minQ)<0.00001_DP) .and. (distQ>0.001_DP)) then
			IF (Q_new>Q) THEN 
				maxQ = maxQ+0.0001_DP
			ELSE
				minQ = minQ-0.0001_DP
			END IF
		END IF
		!=======================================================

		print*,'Q=' , Q,'Q_new=' , Q_new,'minQ=',minQ, 'maxQ=',maxQ, 'panel_N', panel_N_seq(maxiter_panel),Nbar
		
		Q_new  = (maxQ+minQ)/2.0_DP
		rr_new = rr_coef*(Q_new**(alpha-rho))
		ww_new = ww_coef*(Q_new**alpha)
	
	END DO
		
	IF (run_stats==0) THEN
		stats_out=0.0_DP
	ELSE
		stats_out(1)=sum(BeqWealthRatioRep)/REAL(numiter_fin_panel,DP)
		stats_out(2)=sum(WealthOutputRatioRep)/REAL(numiter_fin_panel,DP)
		stats_out(3)=sum(top_wealth_rep(:,1))/REAL(numiter_fin_panel,DP)
		stats_out(4)=sum(top_wealth_rep(:,2))/REAL(numiter_fin_panel,DP)
		stats_out(5)=sum(top_wealth_rep(:,3))/REAL(numiter_fin_panel,DP)
		stats_out(6)=sum(top_wealth_rep(:,4))/REAL(numiter_fin_panel,DP)
		stats_out(7)=sum(top_wealth_rep(:,5))/REAL(numiter_fin_panel,DP)
		stats_out(8)=sum(std_logLaborY)/REAL(numiter_fin_panel,DP)
	END IF

 	write(*,"('**************************************')")
 	write(*,"('Q  iter=',i4,' err =',f12.8,' conv?=',l1)") iterQ,distQ,(distQ<tolQ)
 	write(*,"('**************************************')")
	print*, 'max perc dev from Q', maxval(abs(QRep-Q)/Q)*100.0_DP
 	write(*,"('**************************************')")

	IF (ce_run==1) THEN
	
		OPEN (UNIT=1, FILE='CE_fin', STATUS='replace')
		DO i=1,numiter_fin_panel
		WRITE (UNIT=1, FMT=*) CE_fin(i)
		END DO
		CLOSE (UNIT=1)

		OPEN (UNIT=1, FILE='CE_pos_panel', STATUS='replace')
		DO i=1,numiter_fin_panel
		WRITE (UNIT=1, FMT=*) CE_pos_panel(i)
		END DO
		CLOSE (UNIT=1)
		
		OPEN (UNIT=1, FILE='CE_newborn', STATUS='replace')
		DO i=1,numiter_fin_panel
		WRITE (UNIT=1, FMT=*) CE_newborn(i)
		END DO
		CLOSE (UNIT=1)
		
		print*, 'avg CE       ', sum(CE_fin)/REAL(numiter_fin_panel,DP)
		print*, 'avg CE age 1 ', sum(CE_newborn)/REAL(numiter_fin_panel,DP)
		print*, 'frac pos CE  ', sum(CE_pos_panel)/REAL(numiter_fin_panel,DP)
		
	END IF
	
	IF (final_run==1 .and. run_stats==1) THEN
	IF (tauw==0.0_DP) THEN
		OPEN (UNIT=1, FILE='stats_bench_beqwealth', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) BeqWealthRatioRep(iter)
		END DO
		CLOSE (UNIT=1)
	
		OPEN (UNIT=1, FILE='stats_bench_Q', STATUS='replace')
		OPEN (UNIT=2, FILE='stats_bench_Gbar', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) QRep(iter)
			WRITE (UNIT=2, FMT=*) GbarRep(iter)
		END DO
		CLOSE (UNIT=1)
		CLOSE (UNIT=2)
		
		OPEN (UNIT=1, FILE='stats_bench_wealthoutput', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) WealthOutputRatioRep(iter)
		END DO
		CLOSE (UNIT=1)
	
		OPEN (UNIT=1, FILE='stats_bench_top_wealth', STATUS='replace')
		DO iter=1,numiter_fin_panel
			DO i=1,5
				WRITE (UNIT=1, FMT=*) Top_wealth_Rep(iter,i)
			END DO
		END DO
		CLOSE (UNIT=1)

		OPEN (UNIT=1, FILE='stats_bench_stdloglabory', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) std_logLaborY(iter)
		END DO
		CLOSE (UNIT=1)
		
		OPEN (UNIT=1, FILE='Q_seq_bench', STATUS='replace')
		DO iter=1,maxiter_panel
			WRITE (UNIT=1, FMT=*) panel_Q_seq(iter)
		END DO
		CLOSE (UNIT=1)

		OPEN (UNIT=1, FILE='N_seq_bench', STATUS='replace')
		DO iter=1,maxiter_panel
			WRITE (UNIT=1, FMT=*) panel_N_seq(iter)
		END DO
		CLOSE (UNIT=1)
		
		OPEN (UNIT=1, FILE='wealth_seq_bench', STATUS='replace')
		DO iter=1,maxiter_panel
			WRITE (UNIT=1, FMT=*) panel_wealth_seq(iter)
		END DO
		CLOSE (UNIT=1)	

		OPEN (UNIT=1, FILE='xmax_seq_bench', STATUS='replace')
		DO iter=1,maxiter_panel
			WRITE (UNIT=1, FMT=*) panel_xmax_seq(iter)
		END DO
		CLOSE (UNIT=1)	
	
	ELSE     
		OPEN (UNIT=1, FILE='stats_exp_beqwealth', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) BeqWealthRatioRep(iter)
		END DO
		CLOSE (UNIT=1)
	
		OPEN (UNIT=1, FILE='stats_exp_Q', STATUS='replace')
		OPEN (UNIT=2, FILE='stats_exp_Gbar', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) QRep(iter)
			WRITE (UNIT=2, FMT=*) GbarRep(iter)
		END DO
		CLOSE (UNIT=1)
		CLOSE (UNIT=2)
		
		OPEN (UNIT=1, FILE='stats_exp_wealthoutput', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) WealthOutputRatioRep(iter)
		END DO
		CLOSE (UNIT=1)
	
		OPEN (UNIT=1, FILE='stats_exp_top_wealth', STATUS='replace')
		DO iter=1,numiter_fin_panel
			DO i=1,5
				WRITE (UNIT=1, FMT=*) Top_wealth_Rep(iter,i)
			END DO
		END DO
		CLOSE (UNIT=1)

		OPEN (UNIT=1, FILE='stats_exp_stdloglabory', STATUS='replace')
		DO iter=1,numiter_fin_panel
			WRITE (UNIT=1, FMT=*) std_logLaborY(iter)
		END DO
		CLOSE (UNIT=1)
	
	END IF
	
	END IF
	
END SUBROUTINE
!====================================================================
REAL(DP) FUNCTION utility(c)
USE parameters
IMPLICIT NONE

	REAL(DP) :: c
	IF (c>=0.0_DP) THEN
		utility=LOG(c+hp)
	ELSE
		utility=-1.0e+16_DP
	END IF
	
	
END FUNCTION
!====================================================================
SUBROUTINE moment_panel(rr,ww,Q,Nbar,panel_s,panel_a,panel_age,stats_out)
USE parameters
IMPLICIT NONE

	REAL(DP), INTENT(IN) :: Q,Nbar, rr, ww
	REAL(DP),    DIMENSION(panelN,maxAge), INTENT(IN) :: panel_a
	INTEGER(I4B),DIMENSION(panelN,maxAge), INTENT(IN) :: panel_age, panel_s
	REAL(DP) :: BeqWealthRatio, Beq, WealthOutputRatio
	INTEGER(I4B)  :: row, i, age, tmp_ind, ei, alphai
	REAL(DP), DIMENSION(8) :: stats_out

	REAL(DP),    DIMENSION(panel_pop)  :: panelp_a, panelp_a_sort, logLaborY_panel
	INTEGER(I4B),DIMENSION(panel_pop)  :: panelp_a_index, panelp_a_rank
	REAL(DP),  DIMENSION(5),PARAMETER :: perc=(/0.01_DP,0.1_DP,0.2_DP,0.4_DP,0.6_DP/)
	REAL(DP),  DIMENSION(5),PARAMETER :: perc_wealth_data=(/0.34_DP,0.71_DP,0.83_DP,0.95_DP,0.99_DP/)
	REAL(DP),  DIMENSION(5) :: perc_wealth_model
!	REAL(DP),  DIMENSION(panel_pop) :: gini_s_sum, gini_s, gini_lag
!	REAL(DP) :: gini

	i=0
	DO row=1,panelN
		DO age=1,maxAge
			IF (panel_age(row,age)>0) THEN
				i=i+1
				panelp_a(i)=panel_a(row,age)
			END IF
		END DO
	END DO

	CALL indexx(panel_pop,panelp_a,panelp_a_index)
	FORALL (i=1:panel_pop) panelp_a_rank(panelp_a_index(i))=i
	FORALL (i=1:panel_pop) panelp_a_sort(panelp_a_rank(i))=panelp_a(i)

	DO i=1,5
		tmp_ind=NINT((1.0_DP-perc(i))*REAL(panel_pop,DP))
		perc_wealth_model(i)=SUM(panelp_a_sort(tmp_ind:panel_pop))/SUM(panelp_a)
!		print*, perc(i),perc_wealth_data(i),perc_wealth_model(i)
	END DO
	
	Beq=0.0_DP
	BeqWealthRatio=0.0_DP	
	DO row=1,panelN
		DO age=1,maxAge
			IF (panel_age(row,age)==1) THEN
	    			Beq=Beq+panel_a(row,age)
		    END IF
	    END DO
	END DO
	BeqWealthRatio=Beq/SUM(panelp_a)
	WealthOutputRatio=(SUM(panelp_a)/REAL(panel_pop,DP))/(Aprod*Q**alpha*Nbar**(1.0_DP-alpha))
!	write(*,"('**************************************')")
!	write(*,"(a15,1x,f12.8)" ) 'Bequest/wealth ',BeqWealthRatio
!	write(*,"(a15,1x,f12.8)" ) 'Wealth/Output  ',WealthOutputRatio
!	write(*,"('**************************************')")
	
!	FORALL (i=1:panel_pop) gini_s(i)=SUM(panelp_a_sort(1:i)/REAL(panel_pop,DP))
!
!	gini_lag(1)=0.0_DP
!	gini_lag(2:panel_pop)=gini_s(1:panel_pop-1)
!	gini_s_sum=(gini_lag+gini_s)/REAL(panel_pop,DP)
!
!	gini=1.0_DP-SUM(gini_s_sum)/gini_s(panel_pop)
! 	print*, 'gini', gini

	stats_out(1)=BeqWealthRatio
	stats_out(2)=WealthOutputRatio
	stats_out(3:7)=perc_wealth_model

	! std of cross sectional labor income
	logLaborY_panel=0.0_DP
	i=0
	DO row=1,panelN
		DO age=1,min(ReAge,maxAge)
			IF (panel_age(row,age)>0) THEN
				i=i+1
				ei=se(panel_s(row,age))
				alphai=salpha(panel_s(row,age))
				logLaborY_panel(i)=log(ww*egrid(ei)*exp(e_alpha(alphai)+e_kappa(age)))
			END IF
		END DO
	END DO
	stats_out(8)=sqrt(SUM((logLaborY_panel-SUM(logLaborY_panel)/ &
		& REAL(panel_pop,DP))**2)/REAL(panel_pop,DP))
	
END SUBROUTINE
!====================================================================
FUNCTION splint(xa,ya,y2a,n,x)
USE nrtype
USE global
IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), DIMENSION(n), INTENT(IN) :: xa,ya,y2a
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: splint
	INTEGER(I4B)  :: k,khi,klo
	REAL(DP) :: a,b,h, cspline_interp, linear_interp, weight_linear
	
	IF (x>xa(n)) THEN ! extrapolate above
		klo=n-1
		khi=n
		h=xa(khi)-xa(klo)
		a=(xa(khi)-x)/h
		b=(x-xa(klo))/h
		cspline_interp=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
		linear_interp=a*ya(klo)+b*ya(khi)
		IF (x>xa(n)+h) THEN
			weight_linear=1.0_DP ! if beyond one step past end point, ues linear interpolation
		ELSE
			weight_linear=(x-xa(n))/h ! if within one step past end point, use combo of linear and spline interpolation
		END IF
		splint=weight_linear*linear_interp+(1.0_DP-weight_linear)*cspline_interp
! 		print*, 'extrapolate here',x,splint,cspline_interp,linear_interp,ya(klo),ya(khi)
	ELSE
		klo=1
		khi=n
		DO while ((khi-klo)>1)
		     k=(khi+klo)/2
		     if (xa(k) .gt. x) then
		        khi=k
		     else
		        klo=k
		     end if
		end do
		h=xa(khi)-xa(klo)
		if (h == 0.0) pause 'bad xa input in splint'
		a=(xa(khi)-x)/h
		b=(x-xa(klo))/h
		splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
	END IF
	
END FUNCTION splint
!====================================================================
subroutine bisec(xa,n,x,xlim,ab)
USE nrtype
implicit none

INTEGER(I4B),  INTENT(IN) :: n
REAL(DP), DIMENSION(n), INTENT(IN) 	:: xa
REAL(DP), INTENT(IN) :: x

INTEGER(I4B),  DIMENSION(2), INTENT(OUT)	:: xlim
REAL(DP), DIMENSION(2), INTENT(OUT) :: ab

INTEGER(I4B)  :: klo, khi, k
REAL(DP) :: h


klo=1 						! We will find the right place in the table by means of bisection.
khi=n						! This is optimal if sequential calls to this routine are at random
							! values of x. If sequential calls are in order, and closely
							! spaced, one would do better to store previous values of
							! klo and khi and test if they remain appropriate on the next call.

DO while ( (khi-klo) > 1)
     k=(khi+klo)/2
     if (xa(k) .GT. x) then
        khi=k
     else
        klo=k
     endif
ENDdo    						! klo and khi now bracket the input value of x.
h=xa(khi)-xa(klo)

ab(1)=(xa(khi)-x)/h
ab(2)=(x-xa(klo))/h

xlim(1)=klo
xlim(2)=khi

END SUBROUTINE
!====================================================================
SUBROUTINE saving
USE global
IMPLICIT NONE

	INTEGER(I4B) :: si,ai,age,row

	OPEN (UNIT=1, FILE='rr_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) rr_bench
	CLOSE (UNIT=1)


	OPEN (UNIT=1, FILE='ww_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) ww_bench
	CLOSE (UNIT=1)


	OPEN (UNIT=1, FILE='Gbar_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Gbar_bench
	CLOSE (UNIT=1)

	OPEN (UNIT=1, FILE='Q_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Q_bench
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='Nbar_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Nbar_bench
	CLOSE (UNIT=1)
	
	
	OPEN (UNIT=1, FILE='tauk_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) tauk_bench
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='tauw_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) tauw_bench
	CLOSE (UNIT=1)

	OPEN (UNIT=1, FILE='taul', STATUS='replace')
	WRITE (UNIT=1, FMT=*) taul0
	CLOSE (UNIT=1)	

	OPEN (UNIT=1, FILE='beta', STATUS='replace')
	WRITE (UNIT=1, FMT=*) beta
	CLOSE (UNIT=1)	

! 	
	OPEN (UNIT=1, FILE='valuebench', STATUS='replace')
	DO si=1,ns
		DO ai=1,na
			DO age=1,maxAge
				WRITE (UNIT=1, FMT=*) value_bench(si,ai,age)
			END DO
		END DO
	END DO
	CLOSE (UNIT=1)

	OPEN (UNIT=2, FILE='aprimebench', STATUS='replace')
	DO si=1,ns
		DO ai=1,na
			DO age=1,maxAge
				WRITE (UNIT=2, FMT=*) aprime_bench(si,ai,age)
			END DO
		END DO
	END DO
	CLOSE (UNIT=2)
			
	OPEN (UNIT=1, FILE='Gbar_bench', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Gbar_bench
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='panel_bench_s', STATUS='replace')
	DO row=1,panelN
		DO age=1,maxAge
			WRITE (UNIT=1, FMT=*) panel_bench_s(row,age)
		END DO
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='panel_bench_age', STATUS='replace')
	DO row=1,panelN
		DO age=1,maxAge
			WRITE (UNIT=1, FMT=*) panel_bench_age(row,age)
		END DO
	END DO
	CLOSE (UNIT=1)

	OPEN (UNIT=1, FILE='panel_bench_a', STATUS='replace')
	DO row=1,panelN
		DO age=1,maxAge
			WRITE (UNIT=1, FMT=*) panel_bench_a(row,age)
		END DO
	END DO
	CLOSE (UNIT=1)
	
	
	! FOLLOWING POLICY VARIABLES SAVED ONLY IF EXPRUN=1
	IF (RunExp==1) THEN
	
	OPEN (UNIT=2, FILE='Gbar_exp', STATUS='replace')
	WRITE (UNIT=2, FMT=*) Gbar_exp
	CLOSE (UNIT=2)
	
	OPEN (UNIT=1, FILE='valueexp', STATUS='replace')
	DO si=1,ns
		DO ai=1,na
			DO age=1,maxAge
				WRITE (UNIT=1, FMT=*) value_exp(si,ai,age)
			END DO
		END DO
	END DO
	CLOSE (UNIT=1)

	OPEN (UNIT=2, FILE='aprimeexp', STATUS='replace')
	DO si=1,ns
		DO ai=1,na
			DO age=1,maxAge
				WRITE (UNIT=2, FMT=*) aprime_exp(si,ai,age)
			END DO
		END DO
	END DO
	CLOSE (UNIT=2)
	
	
	OPEN (UNIT=1, FILE='panel_exp_s', STATUS='replace')
	DO row=1,panelN
		DO age=1,maxAge
			WRITE (UNIT=1, FMT=*) panel_exp_s(row,age)
		END DO
	END DO
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='panel_exp_age', STATUS='replace')
	DO row=1,panelN
		DO age=1,maxAge
			WRITE (UNIT=1, FMT=*) panel_exp_age(row,age)
		END DO
	END DO
	CLOSE (UNIT=1)

	OPEN (UNIT=1, FILE='panel_exp_a', STATUS='replace')
	DO row=1,panelN
		DO age=1,maxAge
			WRITE (UNIT=1, FMT=*) panel_exp_a(row,age)
		END DO
	END DO
	CLOSE (UNIT=1)

	OPEN (UNIT=1, FILE='rr_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) rr_exp
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='ww_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) ww_exp
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='Gbar_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Gbar_exp
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='Q_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Q_exp
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='Nbar_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) Nbar_exp
	CLOSE (UNIT=1)		

	OPEN (UNIT=1, FILE='tauk_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) tauk_exp
	CLOSE (UNIT=1)
	
	OPEN (UNIT=1, FILE='tauw_exp', STATUS='replace')
	WRITE (UNIT=1, FMT=*) tauw_exp
	CLOSE (UNIT=1)

	END IF
										
END SUBROUTINE
!=======================================================================
REAL(DP) FUNCTION RAN1(ISEED)
  !--------------------------------------------------------------
  ! Returns a uniform random deviate between 0.0 and 1.0. 
  ! Based on: Park and Miller's "Minimal Standard" random number
  !           generator (Comm. ACM, 31, 1192, 1988)
  !--------------------------------------------------------------
    USE nrtype
    IMPLICIT NONE
    INTEGER(I4B), INTENT(INOUT) :: ISEED
    INTEGER(I4B), PARAMETER :: IM=2147483647, IA=16807, IQ=127773, IR= 2836
    REAL(DP),    PARAMETER :: AM=128.0/IM
    INTEGER(I4B) :: K

    K = ISEED/IQ
    ISEED = IA*(ISEED-K*IQ) - IR*K
    IF (ISEED < 0) ISEED = ISEED+IM
    RAN1 = AM*(ISEED/128)
    
END FUNCTION
!=======================================================================
REAL(DP) FUNCTION ran1_n(idum)
USE nrtype
IMPLICIT NONE

	INTEGER(I4B) :: idum
	INTEGER(I4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB
	REAL(DP), PARAMETER :: AM=1./IM,EPS=1.2e-7,RNMX=1.-EPS
	INTEGER(I4B) :: j,k,iv(NTAB),iy
	SAVE iv,iy
	DATA iv /NTAB*0/, iy /0/

	if (idum.le.0.or.iy.eq.0) then
		idum=max(-idum,1)
		do 11 j=NTAB+8,1,-1
			k=idum/IQ
			idum=IA*(idum-k*IQ)-IR*k
			if (idum.lt.0) idum=idum+IM
			if (j.le.NTAB) iv(j)=idum
11		continue
		iy=iv(1)
	endif
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	j=1+iy/NDIV
	iy=iv(j)
	iv(j)=idum
	ran1_n=min(AM*iy,RNMX)
    return
    END
!====================================================================
SUBROUTINE spline(x,y,n,y2)
USE nrtype
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: n
REAL(DP), DIMENSION(n), INTENT(IN) :: x,y
REAL(DP), DIMENSION(n), INTENT(OUT) :: y2
REAL(DP), DIMENSION(n) :: a,b,c,r
INTEGER(I4B)  :: j, NMAX
PARAMETER (NMAX=500)
REAL(DP) :: bet, gam(NMAX)
 
c(1:n-1)=x(2:n)-x(1:n-1)
r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
r(2:n-1)=r(2:n-1)-r(1:n-2)
a(2:n-1)=c(1:n-2)
b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
b(1)=1.0
b(n)=1.0

r(1)=0.0
c(1)=0.0

r(n)=0.0
a(n)=0.0

bet=b(1)
if (bet == 0.0)pause 'tridag: rewrite equations'
y2(1)=r(1)/bet
do j=2,n
	gam(j)=c(j-1)/bet
	bet=b(j)-a(j)*gam(j)
	if (bet == 0.0)pause 'tridag failed'
	y2(j)=(r(j)-a(j)*y2(j-1))/bet
end do
do j=n-1,1,-1
	y2(j)=y2(j)-gam(j+1)*y2(j+1)
end do
	
END SUBROUTINE spline
!====================================================================
SUBROUTINE tauchen(rho_z,sigma_e,nz,zgrid,pr_z,Gz,m)
USE nrtype
IMPLICIT NONE
	
	REAL(DP), INTENT(IN) :: rho_z, sigma_e
	INTEGER(I4B),  INTENT(IN) :: nz
	REAL(DP), INTENT(OUT), DIMENSION(nz)    :: zgrid, Gz
	REAL(DP), INTENT(OUT), DIMENSION(nz,nz) :: pr_z
	REAL(DP), DIMENSION(nz)    :: Gz_new
	REAL(DP) :: a, zstep, cdf_normal, mu_z
	REAL(DP), INTENT(IN) ::  m
	INTEGER(I4B)  :: i, j, k, zi, zz

	mu_z=0.0_DP
	zgrid=0.0_DP
	pr_z=0.0_DP
	a=(1.0_DP-rho_z)*mu_z;
	zgrid(nz)=m*sqrt(sigma_e**2.0_DP/(1.0_DP-rho_z**2))
	zgrid(1)=-zgrid(nz)
	zstep=(zgrid(nz)-zgrid(1))/REAL(nz-1,DP)
	DO i=2,nz-1
		zgrid(i)=zgrid(1)+zstep*REAL(i-1,DP)
	END DO
	zgrid=zgrid+a/(1.0_DP-rho_z)
	
	DO j=1,nz
		DO k=1,nz
			IF (k==1) THEN
				pr_z(j,k)=cdf_normal((zgrid(1)-a-rho_z*zgrid(j)+zstep/2.0_DP)/sigma_e)
			ELSE IF (k==nz) THEN
				pr_z(j,k)=1.0_DP-cdf_normal((zgrid(nz)-a-rho_z*zgrid(j)-zstep/2.0_DP)/sigma_e)
			ELSE
            	pr_z(j,k)=cdf_normal((zgrid(k)-a-rho_z*zgrid(j)+zstep/2.0_DP)/sigma_e)- &
					& cdf_normal((zgrid(k)-a-rho_z*zgrid(j)-zstep/2.0_DP)/sigma_e)
			END IF
		END DO
	END DO
	
	Gz(1)=cdf_normal((zgrid(1)+zstep/2.0_DP)/sigma_e)
	DO zi=2,nz-1
		Gz(zi)=cdf_normal((zgrid(zi)+zstep/2.0_DP)/sigma_e)- &
					& cdf_normal((zgrid(zi)-zstep/2.0_DP)/sigma_e)
	END DO
	Gz(nz)=1.0_DP-cdf_normal((zgrid(nz)-zstep/2.0_DP)/sigma_e)
! 	print*, 'Gz', Gz, 'sum', sum(Gz)

	DO i=1,1000
		Gz_new=0.0_DP
		DO zi=1,nz
			DO zz=1,nz
				Gz_new(zz)=Gz_new(zz)+Gz(zi)*pr_z(zi,zz)
			END DO
		END DO
		Gz=Gz_new
	END DO
	
! 	print*, 'Gz', Gz, 'sum', sum(Gz)
! 	pause
	
END SUBROUTINE
!====================================================================
SUBROUTINE rouwenhorst(rho_z,sigma_e,nz,zgrid,pr_z,Gz,m)
USE nrtype
! Rowenhurst method to approximate univariate AR(1) process by Markov chain
!      y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
! INPUTS: rho - serial correlation coefficient,
!         sigma - coefficient of variation
! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
!         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
!         s is an n-by-1 vector of stationary distribution (binomial)
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: rho_z, sigma_e
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: pr_z
	REAL(DP), DIMENSION(nz), INTENT(OUT) :: zgrid, Gz
	REAL(DP) :: ybar, q
	INTEGER :: zi, zz, i
	REAL(DP), INTENT(IN) ::  m
	INTEGER(I4B), INTENT(IN) :: nz
	REAL(DP) :: zstep, cdf_normal
	REAL(DP), DIMENSION(nz) :: Gz_new
	
	ybar=m*sqrt(sigma_e**2.0_DP/(1.0-rho_z**2))
	q=(1.0_DP+rho_z)/2.0_DP
	CALL rhmat(pr_z)
	
	zgrid(1)=-ybar
	FORALL (zi=2:nz) zgrid(zi)=zgrid(1)+REAL(zi-1,DP)*(2.0_DP*ybar)/REAL(nz-1,DP)

	Gz(1)=cdf_normal((zgrid(1)+zstep/2.0_DP)/sigma_e)
	DO zi=2,nz-1
		Gz(zi)=cdf_normal((zgrid(zi)+zstep/2.0_DP)/sigma_e)- &
					& cdf_normal((zgrid(zi)-zstep/2.0_DP)/sigma_e)
	END DO
	Gz(nz)=1.0_DP-cdf_normal((zgrid(nz)-zstep/2.0_DP)/sigma_e)
! 	print*, 'Gz', Gz, 'sum', sum(Gz)

	DO i=1,1000
		Gz_new=0.0_DP
		DO zi=1,nz
			DO zz=1,nz
				Gz_new(zz)=Gz_new(zz)+Gz(zi)*pr_z(zi,zz)
			END DO
		END DO
		Gz=Gz_new
	END DO
	
CONTAINS

	RECURSIVE SUBROUTINE rhmat(p)
		IMPLICIT NONE
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: p
		REAL(DP), DIMENSION(size(p,dim=1)-1,size(p,dim=2)-1) :: p1
		INTEGER :: h
		h=size(p,dim=1)
		IF (size(p,dim=2)/=h) STOP 'P must be a square matrix'
		IF (h<2) STOP 'P must be at least 2-by-2 matrix'
		IF (h==2) THEN
			p=reshape((/q,1-q,1-q,q/),(/2,2/))
		ELSE
			CALL rhmat(p1)
			p=0.0_DP
			p(1:h-1,1:h-1)=q*p1
			p(1:h-1,2:h)=(1.0_DP-q)*p1+p(1:h-1,2:h)
			p(2:h,1:h-1)=(1.0_DP-q)*p1+p(2:h,1:h-1)
			p(2:h,2:h)=q*p1+p(2:h,2:h)
			p(2:h-1,:)=p(2:h-1,:)/2.0_DP
		END IF
	END SUBROUTINE rhmat
	
END SUBROUTINE rouwenhorst
!====================================================================
REAL(DP) FUNCTION cdf_normal(x)
USE nrtype
IMPLICIT NONE

	real(DP), parameter :: a1 = 0.398942280444D+00
	real(DP), parameter :: a2 = 0.399903438504D+00
	real(DP), parameter :: a3 = 5.75885480458D+00
	real(DP), parameter :: a4 = 29.8213557808D+00
	real(DP), parameter :: a5 = 2.62433121679D+00
	real(DP), parameter :: a6 = 48.6959930692D+00
	real(DP), parameter :: a7 = 5.92885724438D+00
	real(DP), parameter :: b0 = 0.398942280385D+00
	real(DP), parameter :: b1 = 3.8052D-08
	real(DP), parameter :: b2 = 1.00000615302D+00
	real(DP), parameter :: b3 = 3.98064794D-04
	real(DP), parameter :: b4 = 1.98615381364D+00
	real(DP), parameter :: b5 = 0.151679116635D+00
	real(DP), parameter :: b6 = 5.29330324926D+00
	real(DP), parameter :: b7 = 4.8385912808D+00
	real(DP), parameter :: b8 = 15.1508972451D+00
	real(DP), parameter :: b9 = 0.742380924027D+00
	real(DP), parameter :: b10 = 30.789933034D+00
	real(DP), parameter :: b11 = 3.99019417011D+00
	real(DP) q
	real(DP) x
	real(DP) y
!
	!  |X| <= 1.28.
	!
	  if ( abs ( x ) <= 1.28D+00 ) then
	
	    y = 0.5D+00 * x * x
	
	    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
	      + a6 / ( y + a7 ) ) ) )
	!
	!  1.28 < |X| <= 12.7
	!
	  else if ( abs ( x ) <= 12.7D+00 ) then
	
	    y = 0.5D+00 * x * x
	
	    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
	      + b2 / ( abs ( x ) + b3 &
	      + b4 / ( abs ( x ) - b5 &
	      + b6 / ( abs ( x ) + b7 &
	      - b8 / ( abs ( x ) + b9 &
	      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
	!
	!  12.7 < |X|
	!
	  else
	
	    q = 0.0D+00
	
	  end if
	!
	!  Take account of negative X.
	!
	  if ( x < 0.0D+00 ) then
	    cdf_normal = q
	  else
	    cdf_normal = 1.0D+00 - q
	  end if
	
	  return
	end
!=======================================================================
SUBROUTINE indexx(n,arr,index)
USE nrtype
IMPLICIT NONE
INTEGER(I4B),  INTENT(IN) :: n
REAL(DP), DIMENSION(n), INTENT(IN)  :: arr
INTEGER(I4B),  DIMENSION(n), INTENT(OUT) :: index
INTEGER(I4B),  PARAMETER :: NN=15, NSTACK=50
REAL(DP) :: a
INTEGER(I4B) :: k,i,j,indext,jstack,l,r
INTEGER(I4B), DIMENSION(NSTACK) :: istack

FORALL (i=1:n) index(i)=i
jstack=0
l=1
r=n
do
	if (r-l < NN) then
		do j=l+1,r
			indext=index(j)
			a=arr(indext)
			do i=j-1,l,-1
				if (arr(index(i)) <= a) exit
				index(i+1)=index(i)
			end do
			index(i+1)=indext
		end do
		if (jstack == 0) RETURN
		r=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+r)/2
		call swap(index(k),index(l+1))
		call icomp_xchg(index(l),index(r))
		call icomp_xchg(index(l+1),index(r))
		call icomp_xchg(index(l),index(l+1))
		i=l+1
		j=r
		indext=index(l+1)
		a=arr(indext)
		do
			do
				i=i+1
				if (arr(index(i)) >= a) exit
			end do
			do
				j=j-1
				if (arr(index(j)) <= a) exit
			end do
			if (j < i) exit
			call swap(index(i),index(j))
		end do
		index(l+1)=index(j)
		index(j)=indext
		jstack=jstack+2
		if (jstack > NSTACK) then
			print*, 'indexx: NSTACK too small'
			pause
		end if
		
		if (r-i+1 >= j-l) then
			istack(jstack)=r
			istack(jstack-1)=i
			r=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
end do
CONTAINS
!BL
SUBROUTINE icomp_xchg(i,j)
INTEGER(I4B), INTENT(INOUT) :: i,j
INTEGER(I4B) :: swp
if (arr(j) < arr(i)) then
	swp=i
	i=j
	j=swp
end if
END SUBROUTINE icomp_xchg
!BL
SUBROUTINE swap(a,b)
INTEGER(I4B), INTENT(INOUT) :: a,b
INTEGER(I4B) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap
END SUBROUTINE indexx
!=======================================================================
REAL(DP) FUNCTION mybrent_EGa(ax,bx,cx,func_EGa,tol,xmin,coef1,coef2,coef3,power2)
USE nrtype
IMPLICIT NONE

	REAL(DP), INTENT(IN)  :: ax,bx,cx,tol
	REAL(DP), INTENT(OUT) :: xmin
	REAL(DP), EXTERNAL :: func_EGa
	INTEGER, PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
	INTEGER :: iter
	REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	
	REAL(DP), INTENT(IN)	 	:: coef1,coef2,coef3,power2

	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	
	!=========================================================
	fx=func_EGa(coef1,coef2,coef3,power2,x)
	!=========================================================	
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5_dp*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_sp*tol1
		if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
			xmin=x
			mybrent_EGa=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_sp*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		
		!=========================================================			
		fu=func_EGa(coef1,coef2,coef3,power2,u)
		!=========================================================	
		
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call shft(v,w,x,u)
			call shft(fv,fw,fx,fu)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
	PRINT*, 'brent: exceed maximum iterations'
	CONTAINS
	!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(DP), INTENT(OUT) :: a
	REAL(DP), INTENT(INOUT) :: b,c
	REAL(DP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
END FUNCTION mybrent_EGa
!=======================================================================
REAL(DP) FUNCTION func_EGa(coef1,coef2,coef3,power2,x)
USE nrtype
IMPLICIT NONE

	REAL(DP) :: coef1,coef2,coef3,power2,x
	func_EGa=(coef1*x+coef2*x**power2+coef3)**2
	
END FUNCTION
!=======================================================================
