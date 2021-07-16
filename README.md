# Elastic Scattering

Antiproton elastic scattering cross section using partial waves

... more explanation here (ask Andrea)...

## Download and compile

```bash
git clone https://github.com/pbarBrescia/elastic-scattering.git
cd elastic-scattering
./compile.sh
```

## Run

Input arguments are:

- Antiproton lab momentum (e.g. 50. MeV/c)
- Target mass in u.a.m. (e.g. 40.078 for Calcium)
- Target charge in e units (e.g. 20. for Calcium)

```bash
./run.sh 50.0 40.078 20.0
```

The printed result is organized in 5 columns: 

1. angle in degrees
2. angle in degrees again from previous code
3. dsigma/domega given by **strong** amplitudes
4. dsigma/domega given by **coulomb** amplitudes
5. dsigma/domega given by **strong+coulomb** amplitudes

Example:

```
[...]
7.2 3.6 4485.88 1.88047e+07 1.82919e+07
10.8 7.2 4352.16 3.72676e+06 3.77601e+06
14.4 10.8 4137.08 1.18462e+06 1.32648e+06
18 14.4 3851.73 488105 506526
[...]
```

## Plot example

Save output data for a given momentum values in a file... (the *sed* snippet is to skip the first line with theta=0)

```bash
./run.sh 50.0 40.078 20.0 | sed -n '2,$'p > out50.dat
```

... and plot it (ROOT example)

```bash
root -l plot.C
```

![](fig/example_fig_root.png)

## Scripts

A scan which saves a unique long output file

```bash
bash scripts/scan.sh
```

To create a scan + animation:

```bash
bash scripts/animated_gif.sh 
```

![](animationMom.gif)

## 2021 version

In the new version of the code the following features are introduced:
1. Control of the parameters from text file (e.g., `param.dat`)
1. 10 parameters in file requested: plab, A target, Z target and the parameters of optical potential (`u0`, `w0`, `r0r`, `r0i`, `r0c`, `a0r`, `a0i`)
1. New script `execute.sh`: it execute compile and run scripts and then plot the results.

For using `execute.sh`, do:

```bash
./execute.sh <input_file> <output_file>
```

- `<input_file>`: a text file containing the 10 parameters. If you insert a title line, please be sure to put a `#` before. Remember to insert a empty line below the parameters.
- `<output_file>`: a text file in which is saved the table to do the plot in `plot.C`

Here we provide an example of input file, also present in the repository - name: `param.dat`.

```plaintext
#PARAMETERS FOR EXECUTE.SH
#plab At Zt U0 W0 R0r R0i R0c A0r A0i
#MeV/c n n MeV MeV fm fm fm fm fm
#
#50.0 40.078 20.0 30.0 150.0 1.25 1.25 1.25 0.5 0.5
100.0 40.078 20.0 2.0 95.0 1.45 1.17 1.25 0.53 0.5
#
############################################################
#
#PARAMETERS FOR P_SCAN.SH
#plab At Zt U0 W0 R0r R0i R0c A0r A0i opt theta
#MeV/c n n MeV MeV fm fm fm fm fm [] deg
#
#50 40.078 20.0 30.0 150.0 1.25 1.25 1.25 0.5 0.5 mom 999
#50 40.078 20.0 30.0 150.0 1.25 1.25 1.25 0.5 0.5 ang 12.5
```

The lines that starts with a `#` are not read by the scripts - they are comments useful for the user or to store parameter already used.
In this case, only the sixth line is read by `execute.sh`. 

If this specific file is read by the other script (see next section), it gives an error due to the wrong number of parameters.

## p_scan.sh

`p_scan.sh` execute a scan on momentum or angle.

Two options are available:

1. `mom`: momentum scan

With this option, the script executes a scan on momentum from the momentum provided in the parameter file, adding steps of 10 MeV/c for 10 times. It uses the code `antip_scan.for` (in `src/`).

Then, it saves a figure in `fig/` with the **cross section** as a function of the momentum, using `gnuplot`. 

2. `ang`: angle scan

This this option, it executes the calculation of the **differential cross section** for a specific momentum and angle, printing the values on terminal.

The values printed are: 

angle, real(nuclear), imag(nuclear), real(nuc+coul), imag(nuc+coul) 

To execute the script, use:
```bash
./p_scan.sh param.dat
```

**WARNING**: look in the `param.dat` file - or your parameter file - and check if the number of parameters is correct (it must be 12), or check if you have commented the part you do not need for this script.
