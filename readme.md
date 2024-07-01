
# ZHadron 2024 code

Simplified and modified from: https://github.com/pinchunchou/PhysicsHIZHadron2022 , 

which was forked from:
1. https://github.com/FHead/PhysicsHIZHadron2022
2. https://github.com/yenjie/Zhadron2022

Did not fork it because the commit history was already too messy.

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
| `BasicDistribution` | code to make basic distribution plots |
| `Skims` | code to go from forest to internal skims |



