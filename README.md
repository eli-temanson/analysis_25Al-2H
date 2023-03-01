# Analysis for FSU's RESONEUT Array 

This analysis software is an updated version based on the code that was written by Sean Kuvin, a previous graduate student of Ingo Wiedenhover.

Best way to run if the number of data files is huge;
```
time nohup ./Program &
tail -f nohup.out 
```

## Build
Make sure the ROOT-CERN framework is installed on your target computer and ensure the it is in the users path i.e. source /PATH-2-ROOT-SOURCE/bin/thisroot.sh

```
cmake -B ./build
cmake --build ./build --config Release
```

## Analysis Log
After speaking with ingo there are a few issues,
1) The cross section comparison between the 2+ state in the mirror nucleus, why is the transfer reaction's cross section significantly larger?
2) We are populating states from a 5/2+ ground state in 25Al, this causes a lot angular momentum coupling and makes it difficult to extract the l=0,l=2 contributions. 
3) It is possible that there is a unknown excited state in the continumm that decays down to the 1/2+ first exctited state at approximatedly 0.5 MeV. 

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