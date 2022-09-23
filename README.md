# Python-FourierSpace-Modetomography-GraphCreation
Python script for the presentation of a fourier space image and a modetomography for the analysis of angular resolved spectroscopy measurement data

## Description
The script creates a plot of the fourier space (intensity distribution in dependence of the energy and the momentum value (k)) as well as
real space image (modetomography) (intensity distribution in dependence of two spatial coordinates (x,y)). 
The sample data shows the necessary structure of the input data. 

## Example graphs

Fourier space image <br/>
<img width="300" alt="fourier_space_image" src="https://user-images.githubusercontent.com/35634254/191931307-9dffe9dd-ff41-4d90-824a-04bd04197b49.png">

Real space image <br/>
<img width="300" alt="real_space_image" src="https://user-images.githubusercontent.com/35634254/191931320-ff75c8eb-6d8b-414c-b89d-5d2bf38d55dd.png">

## Excecution of the script

The following components must be installed locally:

- [python](https://www.python.org/downloads/) v3.9.7
- matplotlib, for example by typing:
```console
$ python -m pip install matplotlib
```

To run the project with the sample files locally, first change the path to where the python script is located in line 24:

<img width="620" alt="image" src="https://user-images.githubusercontent.com/35634254/191931042-e3009f42-0bf3-44c7-80ca-101f858b3821.png">

(in case of Windows also change / to \\ \\, for macOS just leave it like this)

and then enter the following in Commandline / Bash:

```console
$ git clone https://github.com/jessna96/Python-FourierSpace-Modetomography-GraphCreation.git
$ cd Python-FourierSpace-Modetomography-GraphCreation
$ python FourierAndRealSpaceGraphCreation.py
```
