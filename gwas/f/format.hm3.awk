# Filter an already standardized GWAS without changing retained rows.
# Match std_hm3() in gwas_post_perf.f.sh: HM3 rsID OR GRCh38 position OR P < 1e-3.
# Optional first input is an existing uncompressed thin file: validate, do not emit.
BEGIN {
    FS=OFS="\t"
    while ((getline line < hm3)>0) {
        gsub(/\r/, "", line); split(line,a,/[ \t]+/)
        if (a[1]!="" && a[1]!="SNP") ids[a[1]]=1
    }
    close(hm3)
    while ((getline line < hm3pos)>0) {
        gsub(/\r/, "", line); split(line,a,/[ \t]+/)
        ch=normchr(a[1]); bp=a[4]
        if (ch~/^[0-9]+$/ && bp~/^[0-9]+$/ && bp+0>0)
            positions[ch SUBSEP (bp+0)]=1
    }
    close(hm3pos)
    if (length(ids)==0 || length(positions)==0) fail("empty HM3 reference")
}
function fail(message) {
    print "ERROR: " message > "/dev/stderr"; failed=1; exit 2
}
function normchr(x) {
    sub(/^chr/, "", x); x=toupper(x)
    if (x=="X") x=23; else if (x=="Y") x=24; else if (x=="MT" || x=="M") x=25
    sub(/^0+/, "", x); return x
}
function numeric(x) {return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
FNR==1 {
    if ($0!="SNP\tCHR\tPOS\tEA\tNEA\tEAF\tN\tBETA\tSE\tP\tLOG10P")
        fail("unexpected standardized header: " FILENAME)
    if (FILENAME=="-") {print; have_header=1}
    next
}
{
    if (NF!=11) fail("expected 11 columns: " FILENAME ":" FNR)
    ch=normchr($2); bp=$3
    valid=(ch~/^[0-9]+$/ && ch+0>=1 && ch+0<=25 && bp~/^[0-9]+$/ && bp+0>0)
    keep=valid && (($1 in ids) || ((ch SUBSEP (bp+0)) in positions) || (numeric($10) && $10+0<0.001))
    if (FILENAME!="-") {
        if (!keep) fail("existing thin row would be lost: " $1)
        thin_rows++; next
    }
    input_rows++
    if (keep) {print; output_rows++}
}
END {
    if (!failed) {
        if (!have_header || !input_rows || !output_rows) fail("empty GWAS input/output")
        print "input_rows\toutput_rows\tthin_rows" > audit
        print input_rows,output_rows,thin_rows+0 >> audit
    }
}
