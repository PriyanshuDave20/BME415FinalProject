import csv
import os
import math
import subprocess
import sys


bin_size = 4000
min_mapq = 10
min_reads = 5
alpha = 0.001

align_file = "results/alignment_summary.tsv"
sample_file = "results/sample_summary.tsv"
out_folder = "results/binding_regions"

use_wsl = True
wsl_distro = "Ubuntu"
samtools = "/root/miniconda3/envs/bme415/bin/samtools"


def to_wsl_path(path):
    path = os.path.abspath(path)
    drive, rest = os.path.splitdrive(path)
    return "/mnt/" + drive[0].lower() + rest.replace("\\", "/")


def get_cmd(args):
    if use_wsl:
        return ["wsl.exe", "-d", wsl_distro, "-u", "root", "--",
                samtools] + args
    return ["samtools"] + args


def run_cmd(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    text = result.stdout.decode(errors="replace")
    text += result.stderr.decode(errors="replace")
    return result.returncode, text


def load_data():
    if not os.path.exists(align_file):
        print("Alignment file not found. Run Module 2 first.")
        sys.exit(1)

    with open(align_file, newline="") as f:
        samples = list(csv.DictReader(f, delimiter="\t"))

    if os.path.exists(sample_file):
        with open(sample_file, newline="") as f:
            old_samples = {
                x["sample_id"]: x["condition"]
                for x in csv.DictReader(f, delimiter="\t")
            }

        for s in samples:
            if s["sample_id"] not in old_samples:
                print("Sample missing from Module 1:", s["sample_id"])
                sys.exit(1)

            if old_samples[s["sample_id"]] != s["condition"]:
                print("Condition mismatch:", s["sample_id"])
                sys.exit(1)

    for s in samples:
        if s["status"] != "OK":
            print("Sample failed:", s["sample_id"])
            sys.exit(1)

        if not os.path.exists(s["bam_file"]):
            print("BAM file missing:", s["bam_file"])
            sys.exit(1)

        if not os.path.exists(s["bai_file"]):
            print("BAM index missing:", s["bai_file"])
            sys.exit(1)

    return samples


def count_reads(bam_file):
    bam = to_wsl_path(bam_file) if use_wsl else os.path.abspath(bam_file)

    cmd = get_cmd([
        "view", "-F", "0x904", "-f", "0x40",
        "-q", str(min_mapq), bam
    ])

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    counts = {}
    strands = {}
    total = 0

    for line in process.stdout:
        parts = line.decode(errors="replace").split("\t", 4)

        flag = int(parts[1])
        chrom = parts[2]
        pos = int(parts[3]) - 1

        bin_num = pos // bin_size

        if chrom not in counts:
            counts[chrom] = {}

        counts[chrom][bin_num] = counts[chrom].get(bin_num, 0) + 1

        if chrom not in strands:
            strands[chrom] = {}

        if bin_num not in strands[chrom]:
            strands[chrom][bin_num] = [0, 0]

        if flag & 16:
            strands[chrom][bin_num][1] += 1
        else:
            strands[chrom][bin_num][0] += 1

        total += 1

    process.stdout.close()
    error = process.stderr.read().decode(errors="replace")
    code = process.wait()

    if code != 0:
        raise RuntimeError(error)

    return counts, strands, total


def get_threshold(counts):
    values = []

    for chrom in counts:
        for b in counts[chrom]:
            values.append(counts[chrom][b])

    if len(values) == 0:
        return min_reads, 0

    values.sort()

    
    cutoff = values[int(0.99 * (len(values) - 1))]
    values = [x for x in values if x <= cutoff]

    avg = sum(values) / len(values)

    if avg <= 0:
        return min_reads, 0

    prob = math.exp(-avg)
    total_prob = prob
    k = 0

    while total_prob < 1 - alpha:
        k += 1
        prob = prob * avg / k
        total_prob += prob

    return max(min_reads, k), avg


def find_regions(counts, strands, threshold):
    regions = []

    for chrom in counts:
        good_bins = []

        for b in counts[chrom]:
            if counts[chrom][b] >= threshold:
                good_bins.append(b)

        good_bins.sort()

        groups = []

        for b in good_bins:
            if groups and b == groups[-1][-1] + 1:
                groups[-1].append(b)
            else:
                groups.append([b])

        for group in groups:
            start = group[0] * bin_size
            end = (group[-1] + 1) * bin_size

            support = sum(counts[chrom][b] for b in group)

            plus = sum(strands[chrom][b][0] for b in group)
            minus = sum(strands[chrom][b][1] for b in group)

            total = plus + minus

            strand = "."

            if total > 0:
                if plus / total >= 0.75:
                    strand = "+"
                elif minus / total >= 0.75:
                    strand = "-"

            regions.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "reads": support,
                "bins": len(group),
                "strand": strand
            })

    return regions


def save_regions(sample, regions):
    filename = os.path.join(
        out_folder, sample["sample_id"] + "_regions.tsv"
    )

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        writer.writerow([
            "chrom", "start", "end", "read_support",
            "n_enriched_bins", "strand"
        ])

        for r in regions:
            writer.writerow([
                r["chrom"], r["start"], r["end"],
                r["reads"], r["bins"], r["strand"]
            ])


def main():
    os.makedirs(out_folder, exist_ok=True)

    samples = load_data()

    
    if "--samples" in sys.argv:
        i = sys.argv.index("--samples")
        if i + 1 >= len(sys.argv):
            print("Give a sample ID after --samples")
            sys.exit(1)

        wanted = sys.argv[i + 1].split(",")
        samples = [s for s in samples if s["sample_id"] in wanted]

        if len(samples) == 0:
            print("No samples matched:", sys.argv[i + 1])
            sys.exit(1)

    all_regions = []
    summary = []
    failed = []

    for sample in samples:
        name = sample["sample_id"]
        condition = sample["condition"]
        replicate = sample["replicate"]

        print("Processing", name)

        
        code, output = run_cmd(
            get_cmd(["quickcheck", "-v", to_wsl_path(sample["bam_file"])])
        )

        if code != 0:
            print("BAM problem:", name, output.strip())
            failed.append(name)
            summary.append([name, condition, replicate, 0, 0, 0, 0,
                            0, 0, 0, "FAILED"])
            continue

        try:
            counts, strands, total = count_reads(sample["bam_file"])
        except RuntimeError as e:
            print("Could not count reads:", name, e)
            failed.append(name)
            summary.append([name, condition, replicate, 0, 0, 0, 0,
                            0, 0, 0, "FAILED"])
            continue

        if total == 0:
            print("No reads found:", name)
            failed.append(name)
            summary.append([name, condition, replicate, 0, 0, 0, 0,
                            0, 0, 0, "FAILED"])
            continue

        threshold, avg = get_threshold(counts)

        regions = find_regions(counts, strands, threshold)

        save_regions(sample, regions)

        for r in regions:
            all_regions.append([
                name, condition, replicate,
                r["chrom"], r["start"], r["end"],
                r["reads"], r["bins"], r["strand"]
            ])

        enriched = sum(r["bins"] for r in regions)
        support = sum(r["reads"] for r in regions)
        occupied = len([b for c in counts for b in counts[c]])

        summary.append([
            name, condition, replicate, total,
            occupied,
            round(avg, 2), threshold,
            enriched, len(regions), support, "OK"
        ])

        print("  Reads:", total)
        print("  Background:", round(avg, 2))
        print("  Threshold:", threshold)
        print("  Candidate regions:", len(regions))

    
    with open(os.path.join(out_folder, "all_samples_regions.tsv"),
              "w", newline="") as f:

        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "condition", "replicate", "chrom",
            "start", "end", "read_support", "n_enriched_bins",
            "strand"
        ])
        writer.writerows(all_regions)

    
    with open(os.path.join(out_folder, "binding_regions_summary.tsv"),
              "w", newline="") as f:

        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "condition", "replicate",
            "fragments_counted", "occupied_bins",
            "background_lambda", "threshold_used",
            "n_enriched_bins", "n_regions", "total_support",
            "status"
        ])
        writer.writerows(summary)

    print("\nFinished Module 3")

    if failed:
        print("Samples that failed:", ", ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
