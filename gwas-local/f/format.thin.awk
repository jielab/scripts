# Display candidates: HM3 (rsID or build-specific position) or P < 1e-3.
BEGIN {
  FS=OFS="\t"
  while ((getline line < hm3)>0) {split(line,a,/[ \t]+/); snps[a[1]]=1}
  close(hm3)
  if (hm3pos!="") while ((getline line < hm3pos)>0) {
    split(line,a,/[ \t]+/); if (a[4]~/^[0-9]+$/) positions[normchr(a[1]) SUBSEP (a[4]+0)]=1
  }
  if (hm3pos!="") close(hm3pos)
}
function normchr(x) {
  x=toupper(x);sub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25
  sub(/^0+/,"",x);return x
}
function numeric(x) {return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
NR==1 {
  for(i=1;i<=NF;i++)c[$i]=i
  if (!("SNP" in c && "CHR" in c && "POS" in c && "P" in c)) {
    print "ERROR: thin requires standardized SNP/CHR/POS/P columns" > "/dev/stderr";exit 2
  }
  print;next
}
{
  p=$(c["P"]);ch=normchr($(c["CHR"]));pos=$(c["POS"])
  if (numeric(p) && p+0>=0 && p+0<=1 && ch~/^[0-9]+$/ && ch+0>=1 && ch+0<=25 && pos~/^[0-9]+$/ && pos+0>0 &&
      (p+0<0.001 || $(c["SNP"]) in snps || (ch SUBSEP (pos+0)) in positions)) print
}
