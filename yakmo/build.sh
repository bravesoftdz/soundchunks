#!/bin/sh
flags="-flto -fuse-linker-plugin -O4"
make clean && make CFLAGS="$flags" CPPFLAGS="$flags" LDFLAGS="$flags" && make install-strip && cp /mingw64/bin/yakmo ../encoder
