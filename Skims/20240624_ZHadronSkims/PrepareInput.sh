#!/bin/bash
usage()
{
   echo "Usage: PrepareInput.sh InputPattern NumberPerExecute ScriptFolder OutputBase [the rest...]"
}
if [[ "$1" == "-h" || "$1" == "--help" || "$#" -lt 3 ]]; then
   usage
   exit
fi
InputPattern=$1
shift
Number=$1
shift
ScriptFolder=$1
shift
OutputBase=$1
shift
mkdir -p $ScriptFolder
mkdir -p $OutputBase
Condor=Submit.condor
echo "Universe              = vanilla"           > $Condor
echo "Executable            = $PWD/RunCondor.sh" >> $Condor
echo "should_transfer_files = NO"                >> $Condor
echo "+JobFlavour           = \"nextweek\""      >> $Condor
echo                                             >> $Condor
Count=0
batch=0
files=""

while IFS= read -r -d '' file
do
    if [ $((Count % Number)) -eq 0 ] && [ $Count -ne 0 ]; then
        echo "Arguments = $PWD $ProjectBase $CMSSW_BASE $OutputBase/Result${batch}.root $batch --Input ${files%,} $@" >> $Condor
        echo "Output    = $ScriptFolder/Part${batch}.out" >> $Condor
        echo "Error     = $ScriptFolder/Part${batch}.err" >> $Condor
        echo "Log       = $ScriptFolder/Part${batch}.log" >> $Condor
        echo "Queue" >> $Condor
        echo >> $Condor
        batch=$((batch + 1))
        files=""
    fi
    files+="$file,"
    Count=$((Count + 1))
done < <(find $InputPattern -name "*root" -print0)

# Handle the last batch if it exists
if [ -n "$files" ]; then
    echo "Arguments = $PWD $ProjectBase $CMSSW_BASE $OutputBase/Result${batch}.root $batch --Input ${files%,} $@" >> $Condor
    echo "Output    = $ScriptFolder/Part${batch}.out" >> $Condor
    echo "Error     = $ScriptFolder/Part${batch}.err" >> $Condor
    echo "Log       = $ScriptFolder/Part${batch}.log" >> $Condor
    echo "Queue" >> $Condor
fi