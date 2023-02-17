#!/bin/bash

name=wfc
cc=gcc
src=main.c

flags=(
    -std=c99
    -O2
    -Wall
    -Wextra
    -pedantic
)

inc=(
    -I.
    -Iutopia
    -Iimgtool
    -Ispxe
)

lib=(
    -Llib
    -limgtool
    -lz
    -lpng
    -ljpeg
)

linux=(
    -lm
    -lpthread
    -D_POSIX_C_SOURCE=199309L
)

glmac=(
    -lglfw
    -framework OpenGL
)

gllin=(
    -lglfw
    -lGL
    -lGLEW
)

cc() {
    echo "$@" && $@
}

lib_build() {
    pushd $1/ && ./build.sh $2 && mv *.a ../lib/ && popd
}

build() {
    [ ! -d lib ] && mkdir lib && echo "mkdir lib"
    lib_build imgtool static
}

comp() {
    if echo "$OSTYPE" | grep -q "darwin"; then
        os=${glmac[*]}
    elif echo "$OSTYPE" | grep -q "linux"; then
        os=${linux[*]} ${gllin[*]}
    else
        echo "This OS not supported yet" && exit
    fi
    cc $cc $src -o $name ${flags[*]} ${inc[*]} ${lib[*]} ${os[*]}
}

cleanf() {
    [ -f $1 ] && rm $1 && echo "deleted $1"
}

cleand() {
    [ -d $1 ] && rm -r $1 && echo "deleted $1"
}

clean() {
    cleand lib
    cleanf $name
    return 0
}

case "$1" in
    "build")
        build;;
    "comp")
        comp;;
    "clean")
        clean;;
    "all")
        build && comp;;
    *)
        echo "Run with 'build' or 'comp' or to build."
        echo "Use 'clean' to remove local builds.";;
esac
