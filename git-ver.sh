#!/bin/bash

GIT_HASH=$(git describe --match=nEvErMaTcH --always --abbrev=8 --dirty)
GIT_VERSION=($(git describe --dirty --always --tags | sed -nE 's/^v([0-9]+).([0-9]+).([0-9]+)(.*)/\1 \2 \3/p'))

VER_MAJOR=${GIT_VERSION[0]:-0}
VER_MINOR=${GIT_VERSION[1]:-0}
VER_PATCH=${GIT_VERSION[2]:-0}

case $1 in
    "hash" )
        echo $GIT_HASH ;;
    "major" )
        echo $VER_MAJOR ;;
    "minor" )
        echo $VER_MINOR ;;
    "patch" )
        echo $VER_PATCH ;;
    * )
        echo \"$GIT_HASH\" $VER_MAJOR $VER_MINOR $VER_PATCH ;;
esac
