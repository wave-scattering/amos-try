#!/bin/bash
rm -rf FiNeR
git clone  git@github.com:kostyfisik/FiNeR.git
cd FiNeR
git checkout cmake
rm -rf .git

cd src/third_party
rm -rf *
git clone git@github.com:kostyfisik/BeFoR64.git
cd BeFoR64
git checkout cmake
cd ..

git clone git@github.com:szaghi/PENF.git
git clone git@github.com:szaghi/FACE.git
git clone git@github.com:pdebuyl/fortran_tester.git

git clone git@github.com:kostyfisik/FLAP.git
cd FLAP
git checkout cmake
cd ..

git clone git@github.com:kostyfisik/StringiFor.git
cd StringiFor
git checkout cmake


cd ../../../
find . -name ".git" -exec rm -rf "{}" \;
cd ..
./clean.sh
