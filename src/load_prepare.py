import csv
import gzip
import os
import shutil
import sys


raw_folder = "raw_data"
summary_file = "results/sample_summary.tsv"

adapter = "AGATCGGAAGAGC"   
scan_reads = 100000         
pair_checks = 10000         
max_errors = 5              

samples = {
    "SRR12192530": {"condition": "WT", "replicate": 1, "instrument": "HiSeq", "layout": "interleaved", "gsm": "GSM4665362"},
    "SRR12192531": {"condition": "WT", "replicate": 1, "instrument": "MiSeq", "layout": "interleaved", "gsm": "GSM4665363"},
    "SRR12192532": {"condition": "WT", "replicate": 2, "instrument": "HiSeq", "layout": "interleaved", "gsm": "GSM4665364"},
    "SRR12192533": {"condition": "WT", "replicate": 2, "instrument": "MiSeq", "layout": "paired", "gsm": "GSM4665365"},
    "SRR12192534": {"condition": "WT", "replicate": 3, "instrument": "HiSeq", "layout": "paired", "gsm": "GSM4665366"},
    "SRR12192535": {"condition": "WT", "replicate": 4, "instrument": "HiSeq", "layout": "paired", "gsm": "GSM4665367"},
    "SRR12192536": {"condition": "R168C", "replicate": 1, "instrument": "HiSeq", "layout": "paired", "gsm": "GSM4665368"},
    "SRR12192537": {"condition": "R168C", "replicate": 1, "instrument": "MiSeq", "layout": "paired", "gsm": "GSM4665369"},
    "SRR12192538": {"condition": "R168C", "replicate": 2, "instrument": "HiSeq", "layout": "paired", "gsm": "GSM4665370"},
    "SRR12192539": {"condition": "R168C", "replicate": 2, "instrument": "MiSeq", "layout": "paired", "gsm": "GSM4665371"},
    "SRR12192540": {"condition": "R168C", "replicate": 3, "instrument": "HiSeq", "layout": "paired", "gsm": "GSM4665372"},
    "SRR12192541": {"condition": "R168C", "replicate": 4, "instrument": "HiSeq", "layout": "paired", "gsm": "GSM4665373"},
}


def find_fastq(sample_id, layout):
    
    endings = [""] if layout == "interleaved" else ["_1", "_2"]
    files = []
    for ending in endings:
        found = None
        for ext in (".fastq", ".fastq.gz"):
            path = os.path.join(raw_folder, sample_id + ending + ext)
            if os.path.exists(path):
                found = path
                break
        if found is None:
            return []
        files.append(found)
    return files


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_id(header):
    return header.split()[0].lstrip("@")


def count_and_validate(path):
    
    
    reads = 0
    errors = []
    with open_fastq(path) as f:
        lines = iter(f)
        for header in lines:
            try:
                seq = next(lines)
                plus = next(lines)
                qual = next(lines)
            except StopIteration:
                errors.append("file ends in the middle of read %d" % (reads + 1))
                break
            reads += 1
            if not header.startswith("@") or not plus.startswith("+"):
                errors.append("bad FASTQ format at read %d" % reads)
            elif len(seq.strip()) != len(qual.strip()):
                errors.append("sequence/quality length mismatch at read %d" % reads)
            if len(errors) >= max_errors:
                errors.append("(stopped checking after %d errors)" % max_errors)
                break
    return reads, errors


def check_pair_ids(files, layout):
    
    
    checked = 0
    mismatched = 0
    if layout == "paired":
        f1 = open_fastq(files[0])
        f2 = open_fastq(files[1])
        for i in range(pair_checks):
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            for j in range(3):
                f1.readline()
                f2.readline()
            checked += 1
            if read_id(h1) != read_id(h2):
                mismatched += 1
        f1.close()
        f2.close()
    else:
        f = open_fastq(files[0])
        for i in range(pair_checks):
            h1 = f.readline()
            if not h1:
                break
            for j in range(3):
                f.readline()
            h2 = f.readline()
            if not h2:
                mismatched += 1
                checked += 1
                break
            for j in range(3):
                f.readline()
            checked += 1
            if read_id(h1) != read_id(h2):
                mismatched += 1
        f.close()
    return checked, mismatched


def header_condition(path):
    
    
    condition = "UNKNOWN"
    checked = 0
    with open_fastq(path) as f:
        n = 0
        for line in f:
            if n % 4 == 0:
                checked += 1
                fields = line.split()
                if len(fields) >= 2:
                    parts = fields[1].split("_")
                    if len(parts) >= 2 and parts[1].endswith("-R168C"):
                        condition = "R168C"
                    elif len(parts) >= 2 and parts[1] == "KHDRBS2":
                        condition = "WT"
                if checked >= pair_checks:
                    break
            n += 1
    return condition


def count_adapter(path, interleaved):
    
    
    hits = 0
    scanned = 0
    record = 0
    with open_fastq(path) as f:
        n = 0
        for line in f:
            if n % 4 == 1:
                if not interleaved or record % 2 == 0:
                    scanned += 1
                    if adapter in line:
                        hits += 1
                    if scanned >= scan_reads:
                        break
                record += 1
            n += 1
    return hits, scanned


def main():
    print("=" * 70)
    print("Module 1: Load and Prepare Data")
    print("KHDRBS2 WT vs R168C easyCLIP-seq (GSE154168)")
    print("=" * 70)

    if not os.path.isdir(raw_folder):
        print("ERROR: raw data folder not found:", raw_folder)
        sys.exit(1)

    results = []
    problems = []

    for sample_id in sorted(samples):
        info = samples[sample_id]
        print("\nSample %s (%s, replicate %d, %s)"
              % (sample_id, info["condition"], info["replicate"],
                 info["instrument"]))

        files = find_fastq(sample_id, info["layout"])
        if len(files) == 0:
            print("  ERROR: expected FASTQ file(s) not found in", raw_folder)
            problems.append(sample_id + ": missing FASTQ files")
            results.append([sample_id, info["condition"], info["replicate"],
                            info["instrument"], info["layout"], info["gsm"],
                            "", "", 0, 0, "MISSING FILES"])
            continue

        
        print("  counting and validating reads ...")
        counts = []
        errors = []
        for path in files:
            try:
                n, errs = count_and_validate(path)
            except OSError as e:
                n, errs = 0, ["cannot open file: %s" % e]
            counts.append(n)
            errors.extend(errs)

        if info["layout"] == "paired":
            reads_r1, reads_r2 = counts
            print("  reads: R1=%d  R2=%d" % (reads_r1, reads_r2))
        else:
            total = counts[0]
            if total % 2 != 0:
                errors.append("interleaved file has %d reads (odd number)" % total)
            reads_r1 = reads_r2 = total // 2
            print("  reads: %d records = %d pairs (interleaved R1/R2)"
                  % (total, reads_r1))

        for err in errors:
            print("  WARNING:", err)
            problems.append(sample_id + ": " + err)

        
        checked, mismatched = check_pair_ids(files, info["layout"])
        if mismatched:
            msg = "%d/%d checked read pairs have mismatched IDs" % (mismatched, checked)
            print("  ERROR:", msg)
            problems.append(sample_id + ": " + msg)
        else:
            print("  mate read IDs match (first %d pairs checked)" % checked)

        
        header_cond = header_condition(files[0])
        print("  headers say condition=%s (table says %s)"
              % (header_cond, info["condition"]))
        if header_cond != info["condition"]:
            msg = ("header condition %s does not match assignment %s"
                   % (header_cond, info["condition"]))
            print("  ERROR:", msg)
            problems.append(sample_id + ": " + msg)

        
        hits, scanned = count_adapter(files[0], info["layout"] == "interleaved")
        pct = 100.0 * hits / scanned if scanned else 0.0
        print("  adapter %s in %.2f%% of the first %d R1 reads"
              % (adapter, pct, scanned))

        if info["layout"] == "paired":
            r1_file, r2_file = files
        else:
            r1_file = r2_file = files[0]

        status = "OK"
        if errors or mismatched or reads_r1 == 0 or reads_r2 == 0:
            status = "ERROR"

        results.append([sample_id, info["condition"], info["replicate"],
                        info["instrument"], info["layout"], info["gsm"],
                        r1_file, r2_file, reads_r1, reads_r2, status])

    
    os.makedirs("results", exist_ok=True)
    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "condition", "replicate", "instrument",
                         "layout", "gsm", "r1_file", "r2_file",
                         "reads_r1", "reads_r2", "status"])
        writer.writerows(results)
    print("\nSummary table written to", summary_file)

    
    print("\nSummary by condition")
    for condition in ("WT", "R168C"):
        rows = [r for r in results if r[1] == condition]
        total = sum(r[8] for r in rows)
        print("  %s: %d sample(s), %d read pairs total"
              % (condition, len(rows), total))
        for r in rows:
            print("    %s  rep%d  %-6s  %-11s  %9d pairs  [%s]"
                  % (r[0], r[2], r[3], r[4], r[8], r[10]))

    
    print("\nPreprocessing needs:")
    print("  UMI extraction:   not needed - the 6-nt sample barcode and the")
    print("                    7-nt random barcode (UMI) are already in the headers")
    print("  deduplication:    later, after alignment, using the header barcode")
    print("  deinterleaving:   needed for SRR12192530/31/32 before alignment")
    print("  adapter trimming: recommended (adapter found in R1 reads), cutadapt",
          "installed" if shutil.which("cutadapt") else "NOT installed")

    if problems:
        print("\nModule 1 finished with %d problem(s)" % len(problems))
        sys.exit(1)
    print("\nModule 1 finished successfully. Ready for Module 2 (alignment).")


if __name__ == "__main__":
    main()
