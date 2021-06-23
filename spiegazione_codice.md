# Spiegazione codice

Questo e' un tentativo di riassumere la logica utilizzata dal codice.
Mancano ancora alcune info da colmare.

Andro' per subroutine.

Nel codice i miei commenti sono segnati con `!!`.
Alcune subroutine potrebbero essere segnate con

```fortran
!! NOT CALLED IN THIS CASE! !!
```

e possono essere ignorate perche' non usate.

## MAIN

Righe: [1-194](./src/antiprova.for#L1) 

Qui vengono definite le variabili e costanti usate (fino a riga 75).
Nella subroutine `iniziale` (chiamata a riga 77) vengono settati dei valori di alcune variabili.
In particolare, vengono definite le costanti del **potenziale ottico**:

![equation](http://www.sciweavers.org/tex2img.php?eq=V%28r%29%20%3D%20%5Cfrac%7BU_%7B0%7D%20e%5E%7B%5Cfrac%7B%28r-RR%29%7D%7BAR%7D%7D%7D%7B1%2Be%5E%7B%5Cfrac%7B%28r-RR%29%7D%7BAR%7D%7D%7D%20%2B%20i%20%5C%5B%20%5Cfrac%7BW_%7B0%7D%20e%5E%7B%5Cfrac%7B%28r-RI%29%7D%7BAI%7D%7D%7D%7B1%2Be%5E%7B%5Cfrac%7B%28r-RI%29%7D%7BAI%7D%7D%7D%20%2B%204%20%5Cfrac%7BW_%7B0D%7D%20e%5E%7B%5Cfrac%7B%28r-RI%29%7D%7BAI%7D%7D%7D%7B%281%2Be%5E%7B%5Cfrac%7B%28r-RI%29%7D%7BAI%7D%7D%29%5E%7B2%7D%7D%20%5C%5D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

dove, in ordine, ci sono il termine reale (elastico), il termine immaginario di volume e il termine immaginario di superficie (non-elastici). Quest'ultimo, nel nostro caso, e' messo a zero.
Le costanti regolabili sono dunque `U0` (`u0`), `W0`(`w0`), `AR`(`ar`) e `AI`(`ai`), ovvero le intensita' del potenziale ottico (in MeV) e le *diffusness* (in fm).


Il ciclo `do` principale inizia a riga 93:

```fortran
      do 7345, ittt = 0,50,1
         utheta = ittt*2.d-2*3.14159265               !! theta from 0 to pi ( step = (1/50)*pi )  
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
      esse = -esse*eye*0.5d0/akp
c sigma è la sez diff in fm2/sr
      sigma = abs(esse)**2
      sigtot = sig*pi/akp**2
      opt = 4.d0*pi*dimag(esse)/akp
c      print *,esse,sigma,sigtot,opt
ccc      print *,'       ',sigma,sigtot
      sigmaint = sigmaint+1.d-2*sigma*2.d0*pi
ccc      print *,ucost,sigmaint

      print *,utheta,dreal(esse),dimag(esse)
7345  continue
```

La subroutine fondamentale e' `dstwav4`, che calcola tramite un'altra subroutine le matrici di scattering `S` (spesso diviso in parte reale e complessa, `esser` e `essei`).
`amplr` e `ampli` sono le ampiezze reale e immaginaria (ovvero le `f` usualmente usate sui libri) per ogni valore di `l` (il numero quantico angolare).

Questa e' sostanzialmente la parte "facile", quella che riportano i libri (seppur in forma di codice).
Questo ciclo stampa `utheta`, `dreal(esse)` (la parte reale di `f`) e `dimag(esse)` (la parte immaginaria di `f`).

## SUBROUTINE legendre 

Righe: [364-385](./src/antiprova.for#L364)

Calcola i polinomi di Legendre con l tra 1 e 120.

## SUBROUTINE dstwav4

Righe: [578-638](./src/antiprova.for#L578)

Subroutine ausiliaria: prende la subroutine [DSTWAV](#subroutine-dstwav) e la esegue in 1200 step.

Qui vengono definite alcune masse, ma l'unica usata e' `PMASS` (o `pmass`) che dovrebbe essere la massa del proiettile (antiprotone in questo caso).
Queste masse pero' sembrano essere distanti da questi valori (massa antiprotone: 938.27 MeV) o da quelli della massa ridotta (anche perche' questa dipende dalla massa del target e dalla massa del proiettile).

## SUBROUTINE COULFN e RCWFN 

Righe: [1037-1086](./src/antiprova.for#L1037) e [1088-1284](./src/antiprova.for#L1088)

In queste subroutine vengono calcolate le funzioni di Coulomb, ovvero le soluzioni esatte per lo scattering coulombiano.

Queste sono indispensabili poiche' il codice usa una Distorted-wave Born Approximation (DWBA) che usa le soluzioni Coulomb come onde entranti nel problema, risolvendo lo scattering con l'aggiunta del potenziale ottico definito prima con approssimazione di Born.

Queste subroutine vengono infatti chiamate subito dopo la sezione di matching della subroutine [DSTWAV](#subroutine-dstwav).

Per dettagli sul funzionamento della subroutine `RCWFN`, consultare il paper di [Barnett *et al*](https://www.sciencedirect.com/science/article/abs/pii/0010465574900137).

## SUBROUTINE DSTWAV ([641-1025](./src/antiprova.for#L641))

Questa e' la subroutine principale del codice.

Non tutte le variabili introdotte qui sono commentate o spiegate - anzi, quasi nessuna lo e' - quindi dobbiamo necessariamente fare speculazioni realistiche su cosa siano alcune di esse (ad esempio `YDWR/I`).

Alcuni limiti superiori vengono definiti insieme a costanti:

```fortran
NRMAX=RMAX/H+0.5            !! Max number of steps
RMAX=NRMAX*H                !! Max radius
RMTCH=RMAX-3.0*H            !! Matching R (?)
A1=2.*FM0*RMU/HXC**2        !! 2*931.478*1.007(?)/h_bar^2 
                            !! -> 1.007: conversion u to proton mass?
A12=ZP                      !! Projectile Charge (in electron unit)
RC=A11*ROC                  !! A^(1/3)*Coulomb radius
RR=A11*ROR                  !! A^(1/3)*Nuclear real-part radius 
RI=A11*ROI                  !! A^(1/3)*Nuclear img-part radius
```

Altre definizioni utili:

```fortran
ER=exp(DBLE(-H/AR))          !! e^(-step/a_r) -> Exponential for opt pot (real)
EI=exp(DBLE(-H/AI))          !! e^(-step/a_i) -> Exponential for opt pot (img)
FR=exp(  RR/AR)              !! e^(r_r/a_r) -> second exp for opt pot (real)
FI=exp(  RI/AI)              !! e^(r_i/a_i) -> second exp for opt pot (img)
GAMMA=ESQ*ZT*ZP/(2.*EPCM*RC) !! e^2*TargetCharge*ProjCharge/(2*Energy*Rcoul)
                             !!-> Coulomb barrier/kinetic energy
```

Le definizioni piu' problematiche sono quelle legate alle `YDW(R/I)`.
A prima vista sembrano delle sorta di armoniche sferiche, ma qualcosa non torna.
Altrimenti potrebbero essere funzioni radiali del tipo *u<sub>l</sub>(r)*, che dipendono dal numero quantico *l* su cui si sta facendo il loop. (Questa ipotesi e' rafforzata dalla chiamata in [dstwav4](#subroutine-dstwav4) di `uffr` e `uffi` nelle posizioni di `YDWR` e `YDWI` in `DSTWAV`).

In questo loop su `LP1` (*l*) viene definito il potenziale ottico (parte reale: `UC`, parte img:`WC`), diviso l'energia cinetica (`EPCM`). Questo potenziale e' quello definito all'inizio.
Viene inoltre definita una funzione `RKR(I)`:

```fortran
RKR(I)=1.0-3.0*GAMMA+UC(I)    !! relative wave number (coulomb+real)
```

che sostanzialmente pare essere il numero d'onda (`k`) del moto relativo della parte reale del potenziale (Coulomb + Nucleare reale).
La dipendenza da `I` e' dovuta al fatto che si sta calcolando queste funzioni per i primi step lungo il raggio (`H`). `I` indica dunque quanti step si stanno prendendo in considerazione. 
La prima parte del calcolo in `DSTWAV` e' infatti commentata con:

```fortran
c ---------------------------------------
c        valori in 0,h,2h
c ---------------------------------------
```

Una seconda definizione problematica e' quella dei valori `A3` e `A5`, e di conseguenza di `AR2`,`AI2`,`AR4` e `AI4`.
Quello che e' certo, e' che tutte queste funzioni vengono usate per calcolare le `YDWR/I`.

Queste `YDWR/I` sembrano avere la forma di un'espansione in `r` (in questo caso `I*H`) per r &rightarrow; 0 (ovvero piccoli `I`). In particolare si arriva a ordine r<sup>4</sup>:

```fortran
      YDWR(I+1,IL)= alfaspeed*(1.d0+AR2*RKSQ*ISQ*HSQ+AR4*RKFO*IFO*HFO)
 16   YDWI(I+1,IL)= alfaspeed*(AI2*RKSQ*ISQ*HSQ+AI4*RKFO*IFO*HFO)
```

dove `alfaspeed` sono dei coefficienti/funzioni definite/i precedentemente a seconda del valore di *l* raggiunto nel ciclo.
Queste funzioni dipendono dunque da k<sup>2</sup>, k<sub>r</sub> (`RKR`) e sviluppate ad ordine r<sup>2</sup> (`ISQ*HSQ`) e r<sup>4</sup> (`IFO*HFO`).
A zero, queste funzioni sono nulle.

A parte le precise definizioni delle singole variabili (fondamentali o ausiliarie che siano), il procedimento e' lo stesso usato nella subroutine `RCWFN`, ma usando come funzione d'onda base quella gia' passata attraverso il - aka distorta dal - campo coulombiano.

Dopo varie peripezie algoritmiche infatti si ottengono questi risultati:

```fortran
amplr(lp1) = slp1r(lp1)                     !! S_l (real)
ampli(lp1) = slp1i(lp1)                     !! S_l (img)
slpabs = sqrt(slp1r(lp1)**2+slp1i(lp1)**2)  !! |S_l|
```

dove `S_l` e' l'elemento di matrice di scattering dalla teoria.
Questo viene usato per il calcolo della **sezione d'urto differenziale di reazione** (a meno di fattori):

```fortran
REACT=REACT+(2*L+1)*(1.-SLP1R(LP1)**2-SLP1I(LP1)**2)    !! dSigma reaction !!
!! sum_l (2l+1)*(1-|S_l|^2) -> S_l = eta_l
```

In realta' il codice attuale stampa direttamente i valori di `amplr` e `ampli`, con le quali si possono calcolare sia le reazioni inelastiche che elastiche. 
Piu' precisamente, [nella prima parte del codice](#main) vengono stampate le f(&theta;) (parte reale: `dreal(esse)` e immaginaria: `dimag(esse)`), il cui modulo quadro da' la sezione d'urto differenziale totale.