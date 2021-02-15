#!/bin/bash
rm -rf FiNeR
git clone  git@github.com:szaghi/FiNeR.git
cd FiNeR
git checkout cmake
rm -rf .git

cd src/third_party
rm -rf *
git clone git@github.com:szaghi/BeFoR64.git
cd BeFoR64
git checkout cmake
cd ..

git clone git@github.com:szaghi/PENF.git
git clone git@github.com:szaghi/FACE.git
git clone git@github.com:pdebuyl/fortran_tester.git

git clone git@github.com:szaghi/FLAP.git
cd FLAP
git checkout cmake
cd ..

git clone git@github.com:szaghi/StringiFor.git
cd StringiFor
git checkout cmake

cd ../../../
find . -name ".git" -exec rm -rf "{}" \;
find . -name ".gitignore" -exec rm -rf "{}" \;
find . -name ".gitmodules" -exec rm -rf "{}" \;

cd ..
./clean.sh
