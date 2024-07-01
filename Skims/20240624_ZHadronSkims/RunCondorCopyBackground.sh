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

Background=$1
shift

cd $CMSSWBASE
eval `scramv1 runtime -sh`
cd - > /dev/null
export PATH=$ProjectBase/CommonCode/binary/:$PATH

echo "Running job with argument = '--Output Temp${ID}.root $@'"

BackgroundCount=0
LocalBackgroundList=
for i in `echo $Background | tr ',' ' '`
do
   cp $i Background$BackgroundCount.root
   LocalBackgroundList="Background$BackgroundCount.root,$LocalBackgroundList"
   BackgroundCount=$((BackgroundCount+1))
done

cp $CurrentFolder/Execute .
./Execute --Output Temp${ID}.root --Background $LocalBackgroundList $@
rm ./Execute

mv Temp${ID}.root $Output

rm Background*root


