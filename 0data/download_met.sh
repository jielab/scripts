#!/usr/bin/env bash



aria2c -c -x 4 -s 4 -j 8 \
  --auto-file-renaming=false \
  -d raw \
  -i ukb.EUR.lst \
  --save-session=saved.session \
  --save-session-interval=60