
# ZHadron code!

## Instructions 

First clone this repository somewhere

Make sure you have CERN `root` setup

Every time you start a new shell, this will need to be done: at the base directory, execute the `SetupAnalysis.sh` which will set the necessary environment variables
```
. SetupAnalysis.sh
```

Then go into the `CommonCode` directory and compile by typing `make` to compile the basic things

And have fun with things!


## Code structure

The folders for codes follow a two-layer structure.  The base folder is like the category, and the subfolder is the unit for specific piece of code.  The idea is that each subfolder is some code that specializes in one thing only, and we set things up so that all the parameters are passed in through command line for easier analysis steering.


| Base folder | Explanation |
|---|---|
| `CommonCode` | contains all the necessary shared header files |
| `SampleProcessing` | code related to sample processing |
| `Example` | contains some examples to run things |
| `BasicDistribution` | code to make basic distribution plots |
| `ZMass` | code related to Z mass and width studies |
| `MainAnalysis` | code for the main analysis |
| `Skims` | code to go from forest to internal skims |



