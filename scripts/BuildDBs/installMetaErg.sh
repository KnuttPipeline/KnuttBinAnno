#!/usr/bin/env bash
# Author: Lukas Jansen
# Licensed under the Academic Free License version 3.0
set -e




if (( $# != 1 )); then
    echo "Only one parameter allowed!"
    echo "installMetaErg.sh <installdir>"
    exit 1
fi

targetdir=$1
mkdir $targetdir
swissknife="https://sourceforge.net/projects/swissknife/files/swissknife/1.78/swissknife_1.78.tar.gz/download"
minpath="https://omics.informatics.indiana.edu/mg/get.php?justdoit=yes&software=minpath1.4.tar.gz"

mkdir -p ~/.cpanm

olddir="$(pwd)"
cd "$targetdir"
cpanm "$swissknife" -n

#
# MinPath
#
wget -qO- "$minpath" | tar -xzf -

mkdir bin

ln -sr $(echo "MinPath/MinPath?.?.py") "bin/MinPath.py"


#
# Metaerg
#
git clone https://github.com/xiaoli-dong/metaerg.git metaerg
# Current version is broken
git -C "metaerg" checkout 68907b737b12c7456275b028ea31314aebbe50e

ln -sr metaerg/bin/*.pl bin/

export MinPath=$PATH:$(realpath .)/MinPath
export PATH=$PATH:$(realpath .)/bin
check_tools.pl

cd $olddir

