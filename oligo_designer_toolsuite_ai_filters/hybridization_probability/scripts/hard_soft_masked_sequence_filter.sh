#!/bin/sh
awk -v pat='^(N|[[:lower:]]|[^>].*[N[:lower:]])' '
  NR==1 { prev=$0; next }
  {
    if ($0 ~ pat) {   # current line matches "N" -> drop previous and current
      prev = ""
      next
    } else {          # current does not match -> print previous (if any) and keep current as prev
      if (prev != "") print prev
      prev = $0
    }
  }
  END { if (prev != "") print prev }
' $1 > $2