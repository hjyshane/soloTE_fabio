# build_te.awk : convert UCSC RepeatMasker .txt to SoloTE BEDs
BEGIN {
  FS = "\t"; OFS = "\t";
  four_out = "te4.txt";   # minimal 4-col BED
  seven_out = "te7.txt";  # BED6+1 (includes score,strand,.)
}

# skip header lines if any
NR==1 && ($1=="bin" || $6=="genoName") { next }

function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }

{
  chrom = trim($6)
  start = trim($7)
  end   = trim($8)
  strand = trim($10)
  subf  = trim($11)
  classf = trim($12)
  fam    = trim($13)

  # handle Class/Family combined in repClass (e.g. LINE/L1)
  cls = classf; fam2 = ""
  n = split(classf, toks, "/")
  if (n >= 2) { cls = toks[1]; fam2 = toks[2] }

  # fill missing family
  if (fam=="" || fam==".") fam = fam2

  # sanity checks
  if (chrom=="" || start=="" || end=="" || subf=="" || cls=="" || strand=="") next
  if (start >= end) next

  id = chrom "|" start "|" end "|" subf ":" fam ":" cls "|" strand

  # output
  print chrom, start, end, id >> four_out
  print chrom, start, end, id, 0, strand, "." >> seven_out
}
END { }