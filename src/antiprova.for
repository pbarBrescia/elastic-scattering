      program antiproton
c In questo programma e` stata riscritta
c la sezione d'urto in termini di phase shift.

      implicit real*8 (a-h,o-z)
c parametri della somma finale
      parameter(lsum=120)
c par. relativi alla dstwav
      parameter(lmax=125,maxstep=1201,mxstep = maxstep+2)
c h e` lo step radiale elementare (h per maxstep-1 = raggio max).
c aq e` q in fermi inversi.
      parameter(h=2.5d-2)
      complex*16 esse,eye

      dimension onda(maxstep,lmax,2)
      dimension onda2(maxstep,lmax,2)
      dimension bessel0(maxstep,2),bessel1(maxstep,2)
      dimension efe(mxstep,3),ege(lmax,2),rr(mxstep)

c rr e` usata come raggio in sbr integrale. Successivamente la si 
c potrebbe usare in piu` routines

      common/parametri/r0r,r0i,a0r,a0i,w0d,r0c,r0mt,z0t,r0mp
      common/elle/lm
      common/egendre/ucost,pleg(130)
      common/scatt/amplr(125),ampli(125)

ccc Start of the code added in 2021
      CHARACTER(100) :: num1char
      CHARACTER(100) :: num2char
      CHARACTER(100) :: num3char
   !! IF(COMMAND_ARGUMENT_COUNT().NE.3)THEN
   !! WRITE(*,*)'ERROR, 3 COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
      ! example ./antip 50.0 40.078 20.0
   !! STOP
   !! ENDIF

   !! CALL GET_COMMAND_ARGUMENT(1,num1char)  
   !! CALL GET_COMMAND_ARGUMENT(2,num2char)
   !! CALL GET_COMMAND_ARGUMENT(3,num2char)
    
   !! READ(num1char,*)num1     ! lab momentum MeV/c              
   !! READ(num2char,*)num2     ! target mass
   !! READ(num2char,*)num3     ! target charge

!! Added in June 2021
!! Values for opt pot from terminal
      CHARACTER(100) :: numu0char
      CHARACTER(100) :: numw0char
      CHARACTER(100) :: numr0rchar
      CHARACTER(100) :: numr0ichar
      CHARACTER(100) :: numr0cchar
      CHARACTER(100) :: numa0rchar
      CHARACTER(100) :: numa0ichar
      REAL*8 :: num1
      REAL*8 :: num2
      REAL*8 :: num3
      REAL*8 :: numu0
      REAL*8 :: numw0
      REAL*8 :: numr0r
      REAL*8 :: numr0i
      REAL*8 :: numr0c
      REAL*8 :: numa0r
      REAL*8 :: numa0i
      IF(COMMAND_ARGUMENT_COUNT().NE.10)THEN
      WRITE(*,*)'ERROR, 10 COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
      !! example ./antip 50.0 40.078 20.0 30 150 1.25 1.25 1.2 0.6 0.5
      !! This number will be taken by a text file
      STOP
      ENDIF

      CALL GET_COMMAND_ARGUMENT(1,num1char)  
      CALL GET_COMMAND_ARGUMENT(2,num2char)
      CALL GET_COMMAND_ARGUMENT(3,num3char)
      CALL GET_COMMAND_ARGUMENT(4,numu0char)  
      CALL GET_COMMAND_ARGUMENT(5,numw0char)
      CALL GET_COMMAND_ARGUMENT(6,numr0rchar)
      CALL GET_COMMAND_ARGUMENT(7,numr0ichar)  
      CALL GET_COMMAND_ARGUMENT(8,numr0cchar)
      CALL GET_COMMAND_ARGUMENT(9,numa0rchar)
      CALL GET_COMMAND_ARGUMENT(10,numa0ichar)
      

      READ(num1char,*)num1     ! lab momentum MeV/c              
      READ(num2char,*)num2     ! target mass
      READ(num3char,*)num3     ! target charge
      READ(numu0char,*)numu0       !! real strength opt pot MeV              
      READ(numw0char,*)numw0       !! img strength opt pot MeV
      READ(numr0rchar,*)numr0r     !! real nuclear radius/A^(1/3) fm
      READ(numr0ichar,*)numr0i     !! img nuclear radius/A^(1/3) fm              
      READ(numr0cchar,*)numr0c     !! Coulomb radius/A^(1/3) fm              
      READ(numa0rchar,*)numa0r     !! real diffusness fm
      READ(numa0ichar,*)numa0i     !! img diffusness fm

ccc   End of the code added in 2021
ccc KEY PARAMETERS ///////////////////////////////////

c projectile charge
      z0p = -1.d0 

c target mass and radius 
      r0mt = num2 ! old value: 40.078d0 	      !! mass of target in amu
      z0t = num3  ! old value: 20.d0 		!! Z of target in electron charge

c lab momentum MeV/c (rescalable many times in do 7981, a few lines below)
      pmevc = num1 ! old value: 50.d0

c potential real and imaginary strength, radius/(A^1/3), diffuseness  
      u0 = numu0 				!! old value: 30 MeV
      w0 = numw0 				!! old value: 150 MeV
      w0d = 0.d0 			!! old value: 0
      r0r = numr0r			!! old value: 1.25 
      r0i = numr0i			!! old value: 1.25
      r0c = numr0c 			!! old value:1.2
      a0r = numa0r			!! old value: 0.6
      a0i = numa0i			!! old value: 0.5

ccc END PARAMETERS

      open(23,file='dwav.dat',status='old')

c inizializzazione (in parte interattiva in parte no)
      call iniziale
      int = 1200
      eye = cmplx(0.d0,1.d0)
      pi = 4.d0*atan(1.d0)

c momentum is converted into 1/fm
      pfm = pmevc/197.327

      do 7981, iij=1,1
      aq = pfm*iij
      aktr = 0.d0                                     !! k transverse (=0)
      aklon = aq                                      !! k longitudinal
      akp = sqrt(aklon*aklon+aktr*aktr)               !! |k|
      sigmaint = 0.d0
c      do 7345, icost = 100,-100,-1
c      ucost = 1.d-2*icost
      do 7345, ittt = 0,2000,1
         utheta = ittt*0.0005*pi              !! theta from 0 to pi rad ( step = (1/2000)*pi )  
         ucost = cos(utheta)
c risoluzione delle equazioni radiali (in "onda")
      call dstwav4(onda,maxstep,u0,w0,z0p,akp)    
c calcolo dei polinomi di legendre per un dato theta (ut in unita` pi)
      call legendre
      esser = 0.d0
      essei = 0.d0
      sig = 0.d0
      do 344, i=1,lsum
c      print *,i,amplr(i),ampli(i)
      esser = esser + (2*i-1)*(1.d0-amplr(i))*pleg(i)
      essei = essei + (2*i-1)*(-ampli(i))*pleg(i)
      sig = sig + (2*i-1)*((1.d0-amplr(i))**2+ampli(i)**2) 
344   continue
c fine do L 
      esse = cmplx(esser,essei)
c questa è l'ampiezza che al quadrato da la sigma diff senza coefficienti
      esse = esse*eye*0.5d0/akp
c sigma è la sez diff in fm2/sr
      sigma = abs(esse)**2
      sigtot = sig*pi/akp**2
      opt = 4.d0*pi*dimag(esse)/akp
c      print *,esse,sigma,sigtot,opt
ccc      print *,'       ',sigma,sigtot
      sigmaint = sigmaint+1.d-2*sigma*2.d0*pi
ccc      print *,ucost,sigmaint
	  !if(utheta.ne.0)then
      print *,utheta,dreal(esse),dimag(esse)
      !endif
7345  continue
c fine do angoli

7981   continue

      stop
      end

      subroutine iniziale
      implicit real*8 (a-h,o-z)
      common/parametri/r0r,r0i,a0r,a0i,w0d,r0c,r0mt,z0t,r0mp
      common/elle/lm
      common/wave/bind(3002)
c      open(25,file='snbound.dat',status='old')

      lm = 120
      
      do 525, i=0,599
      il = 5*i
c      read(25,*)bind(il+1),bind(il+2),bind(il+3),bind(il+4),bind(il+5) 
      bind(il+1) = 0;
      bind(il+2) = 0;
      bind(il+3) = 0;
      bind(il+4) = 0;
      bind(il+5) = 0;
 525  continue
      bind(3001) = 0.d0
      bind(3002) = 0.d0

c      do 625, i=1,3001
c      print *,bind(i),i
625   continue

      return
      end

      subroutine legendre

      implicit real*8 (a-h,o-z)
      common/egendre/ucost,pleg(130)
      
c calcolo polinomi di legendre in cos(u)
      pz = 1.d0
      p1 = ucost

      do 501, lp1=1,121
      l=lp1-1
      if (l.eq.0) then
      	pleg(1) = pz
      	pleg(2) = p1
      else
      	pleg(lp1+1) = ((2*l+1)*p1*pleg(lp1)
     *-l*pleg(lp1-1))/(l+1)
      end if      
501   continue

      return
      end


      subroutine dstwav4(dwav,max,u0,w0,z0p,r0k)

c questo e' un programma autosufficiente per far girare la subroutine
c dstwav in 1200 step. 

      implicit real*8 (a-h,o-z)


      parameter (maxstep=1201)
      dimension uffr(maxstep,125),uffi(maxstep,125)
      dimension dwav(max,125,2)

      COMMON /CONST/ PI,HXC,FM0,ESQ
      common/parametri/r0r,r0i,a0r,a0i,w0d,r0c,r0mt,z0t,r0mp
      common/elle/lm
C
      DATA PIMASS,PMASS,EMASS/.1498 ,1.00727647,0.00054858/ !! Masses of pion, proton and electron in amu -> old value: 1.007
C

      anfex = 1.d0
      ILC  =  -1
      PI   =  4.d0*ATAN(anfex)
      FM0  =  931.494                    !! old value: 931.478 -> amu
      HXC  =  197.327					 !! old value: 197.32
      ESQ  =  1.4399764					 !! old value: 1.43986

c nella chiamata della sbr tutte le "O" sono sost. da "0", 
c e altri "0" sono aggiunti per distinguere le stesse variabili 
c nel programma principale e nella sbr

      r0mp = pmass
      r0mu = num2*pmass/(num2+pmass)
      l0max = lm+1
      r0max = 30.d0
      irdim = maxstep
      h0 = r0max/(irdim-1)
      nonloc = 0
      bnldw = 0.d0
      iwrite = 1
      lctoff = l0max    

c iwrite = 1 annulla l'output
c irdim puo` essere minore di maxstep



      CALL DSTWAV(uffr,uffi,u0,w0,0.d0,0.d0,r0r,r0i,1.d0,1.d0,a0r,a0i,
     *1.d0,1.d0,w0d,r0c,r0mt,z0t,r0mp,z0p,0.d0,h0,r0max,l0max,ilc,
     *1,irdim,nonloc,bnldw,iwrite,lctoff,r0k,r0mu,0.d0,1.d0,1.d0)  !! r0mu not declared?

      kd = 1
      do 521, k = 2,1201,1
      kd = kd+1
      	x = (k-1)*h0*r0k
      do 521, lp1 = 1,121
      	dwav(kd,lp1,1) = uffr(k,lp1)
      	dwav(kd,lp1,2) = uffi(k,lp1)
521   continue      
c      print *,'xdw4=',x,k,kd

      return
      end


      SUBROUTINE DSTWAV (YDWR,YDWI,UO,WO,USO,WSO,ROR,ROI,RORS,ROIS,AR,
     *AI,ARS,AIS,WOD,ROC,RMT,ZT,RMP,ZP,SP,H,RMAX,LMAX,ILC,JDIM,
     *IRDIM,NONLOC,BNLDW,IWRITE,LCTOFF,RK,RMU,VE,ROE,AE)

C     THIS IS A real*8 VERSION WITH UP TO 125 PARTIAL WAVES
C     AND UP TO 1200 RADIAL STEPS
C     NO SPIN ORBIT OR EXCHANGE TERMS ARE PERMITTED
C     FOR ILC=0 USE LCTOFF=LMAX TO COMPUTE CORRECT NO OF COUL FUNCTIONS
C
      IMPLICIT real*8(A-H,O-Z)
C
      DIMENSION YDWR(IRDIM,125), YDWI(IRDIM,125), rextra(4,1201)
      COMMON /CONST/ PI,HXC,FM0,ESQ
      common/scatt/amplr(125),ampli(125)
      COMMON RDEPNL(4,1201), 
     1SLP1R(125),SLP1I(125),F(125),FD(125),G(125),
     1GD(125),SIGMA(125),UC(3),WC(3),RKR(3),RX(3),RERX(3),RERXSQ(3),
     2RKRT1(125),RKIT1(125),RKRT2(125),RKIT2(125),RKRT3(125),RKIT3(125)
      DATA XMIN,UPLIM,RESCAL/1.D-4,1.D15,1.D-15/
C
C     SPECIAL VERSION FOR 125 PARTIAL WAVES
C
      IF(JDIM.GT.1)CALL EXIT
      IF (IWRITE) 1,1,5
 1    WRITE (6,82)
      IF (NONLOC) 2,3,3
 2    WRITE (6,83)
      GO TO 4
 3    WRITE (6,84) BNLDW,NONLOC
 4    WRITE (6,85) RMT,ZT,RMP,ZP,SP
      WRITE (6,86)
      WRITE (6,87) UO,ROR,AR,WO,ROI,AI,WOD
      WRITE (6,90) ROC

c ------------------------------------
c             CINEMATICA
c ------------------------------------

 5    A11=RMT**(1./3.)
cc falied attempt change 2018, earlier used target frame mass. 
cc remark: rk is momentum in 1/fm
c      convert = rmt/(rmp+rmt);
c      RMU=RMT*RMP/(RMT+RMP)         
C      rk = rk*RMP/(RMT+RMP)
c      rmu = rmp 
c      rmu = rmp*convert
c      rk = rk*convert;
      rksq = rk*rk
c      EPCM=EP*RMT/(RMT+RMP)         
 401  NRMAX=RMAX/H+0.5                !! Max number of steps
      RMAX=NRMAX*H                    !! Max radius
      RMTCH=RMAX-3.0*H 			  !! Matching R -> to find delta?
      A1=2.*FM0*RMU/HXC**2            !! 2*931.478*1.007(?)/h_bar^2
      A12=ZP                          !! Projectile Charge (in electron unit)
      RC=A11*ROC                      !! A^(1/3)*Coulomb radius
      RR=A11*ROR                      !! A^(1/3)*Nuclear real-part radius 
      RI=A11*ROI                      !! A^(1/3)*Nuclear img-part radius              
cc remark: epcm is energy in MeV 
      epcm = rksq/a1
c      RKSQ=A1*EPCM
c      RK=sqrt(RKSQ)
      RKFO=RKSQ*RKSQ                  !! INDIFFERENT

      HSQ=H*H
      HFO=HSQ*HSQ                     !! INDIFFERENT
      ER=exp(DBLE(-H/AR))          !! e^(-step/a_r) -> Exponential for opt pot (real)
      EI=exp(DBLE(-H/AI))          !! e^(-step/a_i) -> Exponential for opt pot (img)
      FR=exp(  RR/AR)              !! e^(r_r/a_r) -> second exp for opt pot (real)
      FI=exp(  RI/AI)              !! e^(r_i/a_i) -> second exp for opt pot (img)
      GAMMA=ESQ*ZT*ZP/(2.d0*EPCM*RC) !! e^2*TargetCharge*ProjCharge/(2*Energy*Rcoul) -> CoulBarrier/Energy
      RROCSQ=1.0/(RKSQ*RC*RC) 	   !! 1/(k^2*rc^2) !! INDIFFERENT

c      print *,rk,epcm,rmp,rmt,convert

c ---------------------------------------
c        valori in 0,h,2h
c ---------------------------------------


      IF (ILC) 6,7,7
 6    LMINP1=1
      GO TO 8
 7    LMINP1=LMAX+1
 8    LMAXP1=LMAX+1
      IL=0
      DO 10 LP1=LMINP1,LMAXP1
      IL=IL+1
      YDWR(1,IL)=0.0
 10   YDWI(1,IL)=0.0
      DO 16 I=1,2 							!! I = # of h -> I = 1,2 => 1h,2h
      ISQ=I*I
      IFO=ISQ*ISQ
      FR=FR*ER
      FI=FI*EI
      UC(I)=UO*FR/(EPCM*(1.0+FR))					!! Real potential normalized for energy
      WC(I)=(WO*FI/(1.0+FI)+4.0*WOD*FI/(1.0+FI)**2)/EPCM 	!! N.B: W0D = 0 !!
      RKR(I)=1.d0-3.d0*GAMMA+UC(I) 					!! relative wave number (coulomb+real)
!!      write(23,*)gamma,uc(i),rkr(i)
      IL=0
      DO 16 LP1=LMINP1,LMAXP1 					!! LP1 = l (?)
      IL=IL+1
      L=LP1-1 								!! l-1 (?)
      A3=1.d0/(4.d0*L+6.d0) 						!! What are these factors? -> angular momentum?
      A5=1.0/(8.d0*L+20.d0)
      AR2=-RKR(I)*A3 							!! Add a factor to Uc -> Why?
      AI2=-WC(I)*A3 							!! Add a factor to Wc -> Why?
      AR4=-(RKR(I)*AR2-WC(I)*AI2+GAMMA*RROCSQ)*A5 		!! TO UNDERSTAND !!
      AI4=-(RKR(I)*AI2+WC(I)*AR2)*A5 				!! TO UNDERSTAND !!

      !! These are constants/coefficients to add
      !! in the definitions of YDWR and YDWI.

      if (lp1.le.101) then
      	if (i.eq.1) then
      	  alfaspeed = 1.d0
      	end if  
    	if (i.eq.2) then  
      	  alfaspeed = (2.d0)**lp1
      	end if
      end if 

      if (lp1.gt.101.and.lp1.le.151) then
      	if (i.eq.1) then
      	  alfaspeed = 1/(2.d0)**50
      	end if
      	if (i.eq.2) then
      	  alfaspeed = (2.d0)**(lp1-50)
      	end if
      end if

      if (lp1.gt.151) then
      	if (i.eq.1) then
      	  alfaspeed = 1/(2.d0)**100
      	end if
      	if (i.eq.2) then
      	  alfaspeed = (2.d0)**(lp1-100)
      	end if
      end if

      !! YDWR/YDWI
      

      YDWR(I+1,IL)= alfaspeed*(1.d0+AR2*RKSQ*ISQ*HSQ+AR4*RKFO*IFO       !! INDIFFERENT
     1*HFO)
 16   YDWI(I+1,IL)= alfaspeed*(AI2*RKSQ*ISQ*HSQ+AI4*RKFO*IFO*HFO)       !! INDIFFERENT
      C1=1.5d0/RC                   !! INDIFFERENT
      C2=-0.5d0/(RC**3)             !! INDIFFERENT
      COUL=ESQ*A1*ZT*ZP             !! COUL = 2m*e^2*Zp*Zt/h_bar^2
      A8=HSQ/12.d0
      DO 27 I=1,2
      RX(I)=I*H
      RERX(I)=1.0/RX(I)             !! INDIFFERENT
      RERXSQ(I)=RERX(I)**2          !! INDIFFERENT
      UC(I)=RKSQ*UC(I)              !! INDIFFERENT
      WC(I)=RKSQ*WC(I)              !! INDIFFERENT
      C3=COUL*F1BDS2(A12,RC,RX(I),RERX(I),C1,C2)      !! SEE SUBROUTINE F1BDS2      !! INDIFFERENT
      RKR(I)=UC(I)+RKSQ-C3          !! TOTAL K (REAL PART) -> TOTAL ENERGY/(2M/H_BAR^2) !! !! INDIFFERENT
      IF (NONLOC) 19,18,18          !! NONLOC <0 =0 >0
 18   rextra(1,I+1)=UC(I)
      rextra(2,I+1)=WC(I)
      rextra(3,I+1)=C3
      rextra(4,I+1)=RKR(I)
 19   IL=0                          !! NONLOC>0
      DO 27 LP1=LMINP1,LMAXP1
      IL=IL+1
      L=LP1-1
      A10=L*LP1*RERXSQ(I)           !! l(l+1)/Rx^2
      GO TO (24,25),I
 24   RKRT1(IL)=1.d0+A8*(RKR(1)-A10)
      RKIT1(IL)=A8*WC(1)
      GO TO 27
 25   RKRT2(IL)=1.d0+A8*(RKR(2)-A10)
      RKIT2(IL)=A8*WC(2)
 27   CONTINUE


c -------------------------------------
c	valori da 3h al matching
c -------------------------------------


      UOTA1=A1*UO
      WOTA1=A1*WO
      WODTA1=WOD*A1   !! SURFACE TERM = 0 !!
      RX(3)=RX(2)
      IUL=RMTCH/H+3.5d0

c do loop fondamentale

      DO 39 I=3,IUL                                         !! SIMILAR TO PREVIOUS PART...
      RX(3)=RX(3)+H
      rx45 = i*h
      RERX(3)=1.0/RX(3)
      RERXSQ(3)=RERX(3)**2
      IF (ILC) 29,29,32
 29   FR=FR*ER
      FI=FI*EI
      UC(3)=UOTA1*FR/(1.d0+FR)
      WC(3)=WOTA1*FI/(1.d0+FI)+4.d0*WODTA1*FI/(1.d0+FI)**2
      C3=COUL*F1BDS2(A12,RC,RX(3),RERX(3),C1,C2)            !! SEE SUBROUTINE F1BDS2
      RKR(3)=UC(3)+RKSQ-C3
      IF (ILC) 30,31,32
 30   IF (NONLOC) 33,31,31
 31   rextra(1,I+1)=UC(3)
      rextra(2,I+1)=WC(3)
      rextra(3,I+1)=C3
      rextra(4,I+1)=RKR(3)
      GO TO 33
 32   UC(3)=rextra(1,I+1)
      WC(3)=rextra(2,I+1)
      C3=rextra(3,I+1)
      RKR(3)=rextra(4,I+1)
 33   IL=0
      DO 39 LP1=LMINP1,LMAXP1
      IL=IL+1
      L=LP1-1
      A10=L*LP1*RERXSQ(3)
      RKRT3(IL)=1.d0+A8*(RKR(3)-A10)
      RKIT3(IL)=A8*WC(3)
      A81=12.d0-10.d0*RKRT2(IL)                                   !! NEW COEFFICIENT !!
      A82=10.d0*RKIT2(IL)                                         !! NEW COEFFICIENT !!
      RNR=A81*YDWR(I,IL)+A82*YDWI(I,IL)-RKRT1(IL)*YDWR(I-1,IL)+R  !! First time that it uses YDWR/I !!
     1KIT1(IL)*YDWI(I-1,IL)
      RNI=A81*YDWI(I,IL)-A82*YDWR(I,IL)-RKRT1(IL)*YDWI(I-1,IL)-R
     1KIT1(IL)*YDWR(I-1,IL)
      A9=1.d0/(RKRT3(IL)*RKRT3(IL)+RKIT3(IL)*RKIT3(IL))
      YDWR(I+1,IL)=(RNR*RKRT3(IL)+RNI*RKIT3(IL))*A9
      YDWI(I+1,IL)=(-RNR*RKIT3(IL)+RNI*RKRT3(IL))*A9
 	
c se le y sono > 1e15 vengono ridivise per 1e15
     
      IF(abs(YDWR(I+1,IL)).GE.UPLIM)GO TO 200
      IF(abs(YDWI(I+1,IL)).GE.UPLIM)GO TO 200
      GO TO 300
 200  IP1=I+1
      DO 301 K=1,IP1
      YDWR(K,IL)=YDWR(K,IL)*RESCAL
 301  YDWI(K,IL)=YDWI(K,IL)*RESCAL
 
 300  RKRT1(IL)=RKRT2(IL)
      RKIT1(IL)=RKIT2(IL)
      RKRT2(IL)=RKRT3(IL)
 39   RKIT2(IL)=RKIT3(IL)
c        print *,ydwr(199,30),ydwi(199,30) 
c      print *,rx(3),i,rx45,'*'

c -----------------------------------
c	sezione matching
c -----------------------------------

      IC=RMTCH/H+1.5
      IL=0
      DO 42 LP1=LMINP1,LMAXP1
      IL=IL+1
      RKRT1(IL)=((1.d0/6.d1)*(YDWR(IC+3,IL)-YDWR(IC-3,IL))
     1+0.15d0*(YDWR(IC-2,IL)-YDWR(IC+2,IL))+0.75d0*(YDWR(IC+1,IL)-
     2YDWR(IC-1,IL)))/H
 42   RKIT1(IL)=((1.d0/6.d1)*(YDWI(IC+3,IL)-YDWI(IC-3,IL))
     1+0.15d0*(YDWI(IC-2,IL)-YDWI(IC+2,IL))+0.75d0*(YDWI(IC+1,IL)-
     2YDWI(IC-1,IL)))/H
      ETA=RK*RC*GAMMA

c chiamata subroutine coulomb(r-match)

      IF (ILC) 45,43,46
 43   CALL COULFN (F,FD,G,GD,SIGMA,ETA,RK,RMTCH,LCTOFF,IWRITE)
      GO TO 46
 45   CALL COULFN (F,FD,G,GD,SIGMA,ETA,RK,RMTCH,LMAX,IWRITE)

 46   IL=0
      DO 49 LP1=LMINP1,LMAXP1
      L=LP1-1
      IL=IL+1
      A8=sin(SIGMA(LP1)-SIGMA(1))
      A9=cos(SIGMA(LP1)-SIGMA(1))
      ALP1=RKRT1(IL)*YDWR(IC,IL)+RKIT1(IL)*YDWI(IC,IL)            !! A_l
      BLP1=-RKRT1(IL)*YDWI(IC,IL)+RKIT1(IL)*YDWR(IC,IL)           !! B_l
      CLP1=YDWR(IC,IL)**2+YDWI(IC,IL)**2                          !! C_l
C------------------------------------------------------------------------------
C      TAKE RATIOS TO EASE OVERFLOW PROBLEMS ON THE VAX
C______________________________________________________________________________

      !! Here it uses continued fractions algorithm (I guess...)
      !! See article about RCWFN subroutine
      !! Barnett et al.
      ALP1=ALP1/CLP1
      BLP1=BLP1/CLP1
      A2=ALP1*F(LP1)
      A3=BLP1*G(LP1)
      A4=FD(LP1)
      A5=BLP1*F(LP1)
      A6=ALP1*G(LP1)
      A7=GD(LP1)
      VLP1=A2-A3-A4
      WLP1=A5+A6-A7
      XLP1=A2+A3-A4
      YLP1=A5-A6+A7
      XLP1=XLP1/YLP1
      VLP1=VLP1/YLP1
      WLP1=WLP1/YLP1
      YLP1=VLP1/WLP1
      DEN=1./(1.d0+XLP1**2)
      SLP1R(LP1)=-WLP1*(1.d0+XLP1*YLP1)*DEN
      SLP1I(LP1)=VLP1*(1.d0-XLP1/YLP1)*DEN
      amplr(lp1) = slp1r(lp1)                                     !! S_l (real)
      ampli(lp1) = slp1i(lp1)                                     !! S_l (img)
      !!write(23,*) lp1,slp1r(lp1),slp1i(lp1)
      slpabs = sqrt(slp1r(lp1)**2+slp1i(lp1)**2)                  !! |S_l|
c      if (lp1.lt.11)  print *,'**',lp1,amplr(lp1),ampli(lp1)
      A1=F(LP1)+F(LP1)*SLP1R(LP1)+G(LP1)*SLP1I(LP1)
      A2=G(LP1)-G(LP1)*SLP1R(LP1)+F(LP1)*SLP1I(LP1)
      UR=0.5d0*(A1*A9-A2*A8)
      UI=0.5d0*(A1*A8+A2*A9)
      A3=1.d0/CLP1
      RNR=(UR*YDWR(IC,IL)+UI*YDWI(IC,IL))*A3
      RNI=(-UR*YDWI(IC,IL)+UI*YDWR(IC,IL))*A3
      IULP1=IUL+1
      DO 49 K=1,IULP1
      A4=YDWR(K,IL)
      A5=YDWI(K,IL)
      YDWR(K,IL)=RNR*A4-RNI*A5
 49   YDWI(K,IL)=RNR*A5+RNI*A4
      iol = 50
c      print *,iol,ydwr(1200,iol),ydwi(1200,iol),'Ydwr,i'
c -----------------------------------
c		output
c -----------------------------------

      !! HERE ARE CALCULATED THE DSIGMAS !!
      IF (IWRITE) 60,60,110
 110  IF(IWRITE-2)63,111,63
 60   WRITE (6,94) EPCM,RMU,RK,ETA
      WRITE (6,95)
 111  REACT=0.
      DO 62 LP1=1,LMAXP1
      L=LP1-1
      REACT=REACT+(2*L+1)*(1.-SLP1R(LP1)**2-SLP1I(LP1)**2)        !! dSigma reaction !!
      elast=elast+(2*l+1)*((1-slp1r(lp1))**2+slp1i(lp1)**2)        !! dSigma elastic !!
      RJ=L
cc62   WRITE (6,96) L,RJ,SLP1R(LP1),SLP1I(LP1)
62    REACT=REACT*10.*PI/RK**2
      elast=elast*10.*PI/RK**2
      write(23,*)elast,react
      WRITE (6,97) REACT
 63   IF (NONLOC) 72,64,64
 64   A6=0.125*BNLDW**2
      IULP1=IUL+1
      DO 70 K=2,IULP1
      C3=-rextra(1,K)
      IF(NONLOC.EQ.0)GO TO 71
      RNI=-A6*rextra(2,K)
      GO TO (73,74), NONLOC
 73   C3=C3+rextra(3,K)
      GO TO 74
 71   RNI=0.
 74   RNR=A6*C3
      ERNR=exp(RNR)
      CRNI=cos(RNI)
      SRNI=sin(RNI)
      IL=0
      DO 70 LP1=LMINP1,LMAXP1
      IL=IL+1
      L=LP1-1
      A1=YDWR(K,IL)
      A2=YDWI(K,IL)
      YDWR(K,IL)=ERNR*(A1*CRNI-A2*SRNI)
 70   YDWI(K,IL)=ERNR*(A1*SRNI+A2*CRNI)
 72   RETURN

c -------------------------------------
c		format
c -------------------------------------

C
 82   FORMAT(//' DISTORTED WAVE SUBROUTINE')
 83   FORMAT (24H NO NON-LOCAL CORRECTION/ )
 84   FORMAT(' NON-LOCAL CORRECTION, BETA =',F9.4/' OPTION NUMBER =',
     1I4/)
 85   FORMAT (' TARGET- MASS',F7.3,' CHARGE',F7.3,5X,' PROJECTILE- ENERG
     1Y',1PE10.3,' MASS',E10.3,' CHARGE',0PF7.3,' SPIN',F7.3/ )
 86   FORMAT (27H NUCLEAR PARAMETERS (V,R,A))
 87   FORMAT(' CENTRAL POTL real',3F9.4,' IMAG VOL',3F9.4,' SURF',F9.4)
 90   FORMAT(' COULOMB   RADIUS ',F9.4/)
 94   FORMAT(' DERIVED DATA  PSQ/2*RMU, RMU, K, ETA =',4F10.4//
     1' S-MATRIX ELEMENTS')
 95   FORMAT (1X,4H  L ,5H  J  ,14H  real*8  PART  ,14H  IMAG  PART    )
 96   FORMAT (1X,I3,F5.1,1X,E12.4,2X,E12.4)
 97   FORMAT(//1X,'REACTION CROSS SECTION IN MB IS',1PE12.4)
 402  FORMAT(' *****  RELATIVISTIC KINEMATICS  *****')
      END
      
c ----------------------------------------------------------------
c 		END DSTWAV
c ----------------------------------------------------------------







      SUBROUTINE COULFN(F,FD,G,GD,SIGMA,ETA,AK,R,L,IWRITE)
      IMPLICIT real*8(A-H,O-Z)
      real*8 AK
      DIMENSION F(1),FD(1),G(1),GD(1),SIGMA(1)
C------------------------------------------------------------------------------
C      PRINT CONTROL: IWRITE = 1  DON'T PRINT
C                     IWRITE = 0  SOME  PRINTING
C                     IWRITE = -1 PRINT COULOMB FUNCTIONS IN FULL
C______________________________________________________________________________
      RHO=AK*R
      LP1=L+1
      IF (abs(ETA).GE.1.D-7) THEN
C------------------------------------------------------------------------------
C      COULOMB PHASE SHIFT FOR ETA NON-ZERO
C______________________________________________________________________________
      A = ETA**2
      B = A+16.
      C = ETA
      SIGMA(1)=-C+0.5*C*log(B)+3.5*atan(0.25*C)-atan(C)-atan(0.5*C)-
     1atan(C/3.)-C*(1.+(A-48.)/(30.*B**2)+(A**2-160.*A+1280.)/(105.*B**
     24))/(12.*B)
      DO 11 I=2,LP1
      II=I-1
 11   SIGMA(I)=SIGMA(I-1)+atan(ETA/II)
      GO TO 26
C------------------------------------------------------------------------------
C      COULOMB PHASE SHIFT FOR ETA = ZERO
C______________________________________________________________________________
      ELSE
      DO 25 LL=1,LP1
 25   SIGMA(LL)=0.0
      END IF
C------------------------------------------------------------------------------
C      COMPUTE COULOMB FUNCTIONS AND DERIVATIVES
C______________________________________________________________________________
 26   MINL  = 0
      ACCUR = 1.D-8
      STEP  = 999.
      CALL RCWFN(RHO,ETA,MINL,L,F,FD,G,GD,ACCUR,STEP)
C------------------------------------------------------------------------------
C      d/dr = k * d/d(rho)
C______________________________________________________________________________
      DO 32 LL=1,LP1
      FD(LL)=AK*FD(LL)
 32   GD(LL)=AK*GD(LL)
C
 56   continue

      RETURN
      END
 
      SUBROUTINE RCWFN(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP,ACCUR,STEP)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(1),FCP(1),GC(1),GCP(1)
C------------------------------------------------------------------------------
C      COULOMB WAVEFUNCTIONS CALCULATED AT R=RHO BY THE
C      CONTINUED FRACTION METHOD OF STEED
C      MINL,MAXL ARE ACTUAL L-VALUES
C      SEE BARNETT, FENG, STEED, AND GOLDFARB COMPUTER PHYSICS COMM 1974
       !! This paper is on Drive !!
C______________________________________________________________________________
      PACE = STEP
      ACC = ACCUR
      IF(PACE.LT.100.) PACE = 100.
      IF(ACC.LT.1.E-15.OR.ACC.GT.1.E-6) ACC = 1.E-6
      R    = RHO                                      !! RHO = k*r = (2*mu*E)^(1/2)*r/h_bar -> mu = reduced mass
      KTR  = 1
      LMAX = MAXL
      LMIN1= MINL + 1
      XLL1 = DFLOAT(MINL*LMIN1)
      ETA2 = ETA*ETA                                  !! ETA = mu*Z1*Z2*e^2/(hbar^2*k)
      TURN = ETA + sqrt(ETA2 + XLL1)
      IF(R.LT.TURN.AND.abs(ETA).GT.1.E-6) KTR = -1
      KTRP = KTR
      GO TO 2
 1    R    = TURN
      TF   = F
      TFP  = FP
      LMAX = MINL
      KTRP = 1
 2    ETAR = ETA*R
      RHO2 = R*R
      PL   = DFLOAT(LMAX + 1)
      PMX  = PL + 0.5
C------------------------------------------------------------------------------
C      CONTINUED FRACTION FOR FP(MAXL)/F(MAXL)
C______________________________________________________________________________
      FP   = ETA/PL + PL/R
      DK   = ETAR*2.
      DEL  = 0.
      D    = 0.
      F    = 1.
      K    = (PL*PL - PL + ETAR)*(2.*PL - 1.)
      IF(PL*PL+PL+ETAR.NE.0.) GO TO 3
      R    = R + 1.E-6
      GO TO 2
 3    H    = (PL*PL + ETA2)*(1. - PL*PL)*RHO2
      K    = K + DK + PL*PL*6.
      D    = 1./(D*H + K)
      DEL  = DEL*(D*K-1.)
      IF(PL.LT.PMX) DEL = -R*(PL*PL + ETA2)*(PL + 1.)*D/PL
      PL   = PL + 1.
      FP   = FP + DEL
      IF(D.LT.0.) F = -F
      IF(PL.GT.20000.) GO TO 11
      IF(abs(DEL/FP).GE.ACC) GO TO 3
      FP   = F*FP
      IF( LMAX.EQ.MINL) GO TO 5
      FC (LMAX+1) =  F*1.E-30
      FCP(LMAX+1) = FP*1.E-30
C------------------------------------------------------------------------------
C      DOWNWARD RECURSION TO MINL FOR F AND FP
C      ARRAYS GC, GCP ARE STORAGE
C______________________________________________________________________________
      L  = LMAX
      DO 4 LP = LMIN1,LMAX
      PL = DFLOAT(L)
      GC (L+1) = ETA/PL + PL/R
      GCP(L+1) = sqrt(ETA2 + PL*PL)/PL
      FC (L)   = (GC(L+1)*FC(L+1) + FCP(L+1))/GCP(L+1)
      FCP(L)   =  GC(L+1)*FC(L)   - GCP(L+1)*FC(L+1)
 4    L  = L - 1
C------------------------------------------------------------------------------
C      AFTER RECURSION      RENORMALISE FC AND FCP TO PREVENT OVERFLOW
C      IN CALCULATING WRONSKIAN ON VAX
C______________________________________________________________________________
      F  = abs(FC(LMIN1))
      DO 44 L=MINL,LMAX
      FC (L+1)  = FC (L+1)/F
 44   FCP(L+1)  = FCP(L+1)/F
      F  = FC (LMIN1)
      FP = FCP(LMIN1)
 5    IF(KTRP.EQ.-1) GO TO 1
C------------------------------------------------------------------------------
C      REPEAT FOR R = TURN IF RHO LT TURN
C      NOW OBTAIN P + I.Q FOR MINL FROM CONTINUED FRACTION (32)
C      real*8 ARITHMETIC TO FACILITATE CONVERSION TO real*8
C______________________________________________________________________________
      P  = 0.
      Q  = R - ETA
      PL = 0.
      AR = -(ETA2 + XLL1)
      AI =   ETA
      BR = 2.*Q
      BI = 2.
      WI = 2.*ETA
      DR =   BR/(BR*BR + BI*BI)
      DI =  -BI/(BR*BR + BI*BI)
      DP = -(AR*DI + AI*DR)
      DQ =  (AR*DR - AI*DI)
 6    P  =  P + DP
      Q  =  Q + DQ
      PL = PL + 2.
      AR = AR + PL
      AI = AI + WI
      BI = BI + 2.
      D  = AR*DR - AI*DI + BR
      DI = AI*DR + AR*DI + BI
      T  = 1./(D*D + DI*DI)
      DR =  T*D
      DI = -T*DI
      H  = BR*DR - BI*DI - 1.
      K  = BI*DR + BR*DI
      T  = DP*H  - DQ*K
      DQ = DP*K  + DQ*H
      DP = T
      IF(PL.GT.46000.) GO TO 11
      IF(abs(DP)+abs(DQ).GE.(abs(P)+abs(Q))*ACC) GO TO 6
      P  = P/R
      Q  = Q/R
C-----------------------------------------------------------------------------
C      SOLVE FOR FP, G, GP AND NORMALISE F AT L=MINL
C_____________________________________________________________________________
      G  = (FP - P*F)/Q
      GP = P*G - Q*F
      W  = 1./sqrt(FP*G - F*GP)
      G  = W*G
      GP = W*GP
      IF(KTR.EQ.1) GO TO 8
      F  = TF
      FP = TFP
      LMAX = MAXL
C-----------------------------------------------------------------------------
C      RUNGE-KUTTA INTEGRATION OF G(MINL) AND GP(MINL)
C            INWARDS FROM TURN
C            SEE FOX AND MAYERS 1968 PG 202
C_____________________________________________________________________________
      IF(RHO.LT.0.2*TURN) PACE = 999.
      R3 = 1./3.D0
      H  = (RHO - TURN)/(PACE + 1.)
      H2 = 0.5*H
      I2 =      PACE + 0.001
      ETAH = ETA*H
      H2LL = H2*XLL1
      S  = (ETAH + H2LL/R  )/R   - H2
 7    RH2= R + H2
      T  = (ETAH + H2LL/RH2)/RH2 - H2
      K1 = H2*GP
      M1 =  S*G
      K2 = H2*(GP + M1)
      M2 =  T*(G  + K1)
      K3 =  H*(GP + M2)
      M3 =  T*(G  + K2)
      M3 =     M3 + M3
      K4 = H2*(GP + M3)
      RH = R + H
      S  = (ETAH + H2LL/RH )/RH  - H2
      M4 =  S*(G + K3)
      G  = G  + (K1 + K2 + K2 + K3 + K4)*R3
      GP = GP + (M1 + M2 + M2 + M3 + M4)*R3
      R  = RH
      I2 = I2 - 1
      IF(abs(GP).GT.1.E36) GO TO 11
      IF(I2.GE.0) GO TO 7
      W  = 1./(FP*G - F*GP)
C------------------------------------------------------------------------------
C      UPWARD RECURSION FROM GC(MINL) AND GCP(MINL)
C      STORED VALUES ARE R,S
C      RENORMALISE FC,FCP FOR EACH L-VALUE
C______________________________________________________________________________
 8    GC (LMIN1) = G
      GCP(LMIN1) = GP
      IF(LMAX.EQ.MINL) GO TO 10
      DO  9  L = LMIN1,LMAX
      T        = GC(L+1)
      GC (L+1) = (GC(L)*GC (L+1) - GCP(L))/GCP(L+1)
      GCP(L+1) =  GC(L)*GCP(L+1) - GC (L+1)*T
      FC (L+1) = W*FC (L+1)
 9    FCP(L+1) = W*FCP(L+1)
      FC (LMIN1) = FC (LMIN1)*W
      FCP(LMIN1) = FCP(LMIN1)*W
C------------------------------------------------------------------------------
C      CHECK WRONSKIAN FOR L > MINL
C______________________________________________________________________________
      DO 45 L=LMIN1,LMAX
      W  = FCP(L+1)*GC(L+1)-FC(L+1)*GCP(L+1)
 45   IF(abs(W-1.).GT.1.D-10)WRITE (6,46) L
 46   FORMAT(' *** WRONSKIAN GT 1.E-10 FOR L = ',i4,'  ***')
      RETURN
 10   FC (LMIN1) = W*F
      FCP(LMIN1) = W*FP
      RETURN
 11   W  = 0.
      G  = 0.
      GP = 0.
      GO TO 8
      END

      !! FB1DS2 is the depencence on Rx 
      !!-> if Zp!=0 and Rx>Rc =>1/Rx (outside charged sphere)
      !!                Rx<Rc =>C1+C2*Rx^2 (inside charged sphere with charge volume distrib.)
      FUNCTION F1BDS2 (ZP ,RCOULE,RX,RECRX,A5,A6)
      IMPLICIT real*8 (A-H,O-Z)
      IF(abs(ZP).GT.1.D-6)GO TO 2
      F1BDS2=0.0D0
      RETURN
 2    IF (RCOULE-RX) 3,3,4
 3    F1BDS2=RECRX
      RETURN
 4    F1BDS2=A5+A6*RX**2
      RETURN
      END


