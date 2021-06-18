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

## MAIN ([1-194](./src/antiprova.for#L1)) 

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

## SUBROUTINE legendre ([364-385](./src/antiprova.for#L364))

Calcola i polinomi di Legendre con l tra 1 e 120.

## SUBROUTINE dstwav4 ([578-638](./src/antiprova.for#L578))

## SUBROUTINE COULFN e RCWFN ([1037-1086](./src/antiprova.for#L1037) e [1088-1284](./src/antiprova.for#L1088)

Necessario prima queste di DSTWAV.

## SUBROUTINE DSTWAV ([641-1025](./src/antiprova.for#L641))

