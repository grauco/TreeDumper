# treeDumper

##to download and compile:
cmsrel CMSSW_8_0_20
cd src
cmsenv
git clone https://github.com/grauco/treeDumper.git

mkdir ttDM
cp -rf treeDumper/* ttDM/
rm -rf treeDumper
scram b -j 10

##to run:
cd test/
cmsRun topplusdmTrees_cfg.py maxEvts=N sample="mySample/sample.root" version="1"7 outputLabel="myoutput" 
