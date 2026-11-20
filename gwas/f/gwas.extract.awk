# Header-based PLINK2 dosage counts -> biallelic variant IDs, including indels.
function numeric(x) { return x ~ /^[0-9]+([.][0-9]*)?([eE][+-]?[0-9]+)?$/ }
NR == 1 {
    for (i=1; i<=NF; i++) col[$i]=i
    if (!("ID" in col) || !("ALT" in col) || !("ALT_CTS" in col) || !("OBS_CT" in col)) {
        print "Unsupported .acount header" > "/dev/stderr"; bad=1; exit 2
    }
    next
}
{
    total++
    id=$(col["ID"]); alt=$(col["ALT"]); ac=$(col["ALT_CTS"]); obs=$(col["OBS_CT"])
    if (alt ~ /,/) next
    if (!numeric(ac) || !numeric(obs) || obs+0 <= 0 || ac+0 > obs+0) next
    mac=(ac+0 < obs-ac ? ac+0 : obs-ac)
    if (id != "." && mac/(obs+0) >= maf) { print id; kept++ }
}
END {
    if (!bad) {
        print "TOTAL\tRETAINED\tMIN_MAF" > summary
        print total+0 "\t" kept+0 "\t" maf > summary
        print "Shared MAF >= " maf ": " kept+0 " / " total+0 " variants retained" > "/dev/stderr"
    }
}
