#!/bin/sh
make clean && make CFLAGS="-flto -fuse-linker-plugin -Ofast" CPPFLAGS="-flto -fuse-linker-plugin -Ofast" LDFLAGS="-flto -fuse-linker-plugin -Ofast" && cp src/yakmo ../encoder