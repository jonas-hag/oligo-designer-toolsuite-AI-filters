#!/bin/sh
awk '{print NR "," length($0)}' "$1" > "$2"