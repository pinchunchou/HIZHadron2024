#!/bin/bash

echo Yay

CurrentFolder=$1
shift

ProjectBase=$1
shift

CMSSWBASE=$1
shift

Output=$1
shift

ID=$1
shift

cd $CMSSWBASE
eval `scramv1 runtime -sh`
cd - > /dev/null
export PATH=$ProjectBase/CommonCode/binary/:$PATH

echo "Running job with argument = '--Output Temp${ID}.root $@'"

cp $CurrentFolder/Execute .
./Execute --RootOutput Temp${ID}.root --Output Temp${ID}.pdf $@
rm ./Execute

mv Temp${ID}.root $Output
mv Temp${ID}.pdf `dirname $Output`


