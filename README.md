# Analysis for FSU's ResoNeut detector 

This new RN-analysis-et software that was based on code
that was written by Sean Kuvin a previous graduate student 
of Ingo Wiedenhover.

Best way to run if the number of data files is huge;
     time nohup ./Program &
     tail -f nohup.out 

cd build
cmake -G "Visual Studio 17 2022" -A Win32 -T host=x64 -DCMAKE_VERBOSE_MAKEFILE=ON ..\
cd ../
cmake --build .\build\ --target install --config Release

## Build
Make sure the ROOT-CERN framework is installed on your target computer and ensure the it is in the users path i.e. source /PATH-2-ROOT-SOURCE/bin/thisroot.sh

```
cmake -B ./build
cmake --build ./build --config Release
```


<!-- ## Details on the Simulation 
1) DWBA to Probability
2) Regular kinematics to Inverse kinematics 180-theta
3) from the n-theta create a neutron event using target and beam info
   randomize neutron phi
4) Calculate the break up rxn (p, ic-heavy)
   Use "experimental excitation value" and beam energy 
   Randomize proton phi, calculate what ic
5) Convert center of mass frame to lab frame
6) Check is the event is "inside" it's respective detector
7) Fill the _E(0) with the energy, _Pos[0] with the vector
8) Do physics, calculate mbars 
9) Plot, loop, do it again. -->


## Analysis Log
After speaking with ingo there are a few issues,
1) The cross section comparison between the 2+ state in the mirror nucleus, why is the transfer reaction's cross section significantly larger?
2) We are populating states from a 5/2+ ground state in 25Al, this causes a lot angular momentum coupling and makes it difficult to extract the l=0,l=2 contributions. 
3) It is possible that there is a unknown excited state in the continumm that decays down to the 1/2+ first exctited state at approximatedly 0.5 MeV. 


<!-- ## Energy Loss inputs -->
<!-- Energy (MeV/u)
0 - [He-base] F.Hubert et al, AD&ND Tables 46(1990)1	
Energy (MeV/u)	
1 - [H -base] J.F.Ziegler et al, Pergamon Press, NY (low energy)	
Energy (MeV/u)	
2 - ATIMA 1.2  LS-theory (recommended for high energy)	
Energy (MeV/u)	
3 - ATIMA 1.2  without LS-correction	
Energy (MeV/u)	
4 - electrical component of [1] - J.F.Ziegler et al	
Energy (MeV/u)	
5 - nuclear component of [1] - J.F.Ziegler et al	 -->

## Repo setup (just for me)
```
git init -b main
git add . && git commit -m "initial commit"
git submodule add https://github.com/eli-temanson/catima.git src/vendor/catima
git remote add origin https://github.com/eli-temanson/analysis_25Al-2H.git
git push origin main
```

## And Repo updates (just a reminder)
```
git status
git pull
git add . && git commit -m "comment"
git push origin main
```