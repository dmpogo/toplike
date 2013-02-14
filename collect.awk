#!/bin/awk
BEGIN {
     t = 0
}

FNR == 1           { t = 0 }
FNR == 1 && NR > 1 { print s }

$2 ~ /D/ {
     gsub(/D/,"E")
     if ( $2 < t ) {
        t = $2
        s = $0
     } 
}


END {
    print s
}
