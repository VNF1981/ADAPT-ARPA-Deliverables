import os
import shutil
import sys
import tempfile
from pathlib import Path

path = Path(sys.argv[1])
backup = Path(str(path) + ".before_AD.vcf")

if not backup.exists():
    shutil.copy2(path, backup)

fd, temporary = tempfile.mkstemp(dir=path.parent, suffix=".vcf")
records = 0
modified = 0

try:
    contigs = set()
    declared_contigs = set()
    with backup.open() as scan:
        for line in scan:
            if line.startswith("##contig=<ID="):
                declared_contigs.add(line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
            elif not line.startswith("#"):
                contigs.add(line.split("\t", 1)[0])

    with backup.open() as src, os.fdopen(fd, "w") as dst:
        has_ad_header = False

        for line in src:
            if line.startswith("##FORMAT=<ID=AD,"):
                has_ad_header = True

            if line.startswith("#CHROM"):
                for contig in sorted(contigs - declared_contigs):
                    dst.write(f"##contig=<ID={contig}>\n")
                if not has_ad_header:
                    dst.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths derived from FreeBayes RO and AO">\n')
                dst.write(line)
                continue

            if line.startswith("#"):
                dst.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            fmt = fields[8].split(":")
            records += 1

            if "AD" not in fmt:
                if "RO" not in fmt or "AO" not in fmt:
                    raise ValueError(f"Missing RO/AO at {fields[0]}:{fields[1]}")

                ro_index = fmt.index("RO")
                ao_index = fmt.index("AO")
                alt_count = len(fields[4].split(","))

                for i in range(9, len(fields)):
                    sample = fields[i].split(":")
                    ro = sample[ro_index]
                    ao = sample[ao_index]

                    if ro == "." or ao == ".":
                        raise ValueError(f"Missing allele counts at {fields[0]}:{fields[1]}")

                    values = [ro] + ao.split(",")

                    if len(values) != alt_count + 1:
                        raise ValueError(f"Invalid allele count at {fields[0]}:{fields[1]}")

                    sample.append(",".join(values))
                    fields[i] = ":".join(sample)

                fields[8] += ":AD"
                modified += 1

            dst.write("\t".join(fields) + "\n")

    os.replace(temporary, path)

except Exception:
    if os.path.exists(temporary):
        os.unlink(temporary)
    raise

print(f"VCF records: {records}; AD added: {modified}", flush=True)
