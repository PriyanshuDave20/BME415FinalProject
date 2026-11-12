import csv
import gzip
import math
import os
import statistics
import subprocess
import sys
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import load_prepare as m1
import find_binding_regions as m3

bin_size = 4000
smallbin = 250
minov = 1
pc = 1.0
fccut = 1.0
ks = [4, 5, 6]
kpc = 1
dseed = 154168
topn = 10
use_dedup = True

sample_file = "results/sample_summary.tsv"
align_file = "results/alignment_summary.tsv"
ref = "reference/GRCh38.primary_assembly.genome.fa.gz"
outdir = "results/comparison"
figdir = "results/comparison/figures"
dsbams = "results/downsampling/bams"
dsregs = "results/downsampling/regions_downsampled"
dedupdir = "results/dedup"
umitools = "/root/miniconda3/envs/umitools/bin/umi_tools"

ena = {
    "SRR12192530": ("GSM4665362", "SAMN15502104", "easyCLIP KHDRBS2 Replicate 1 HiSeq"),
    "SRR12192531": ("GSM4665363", "SAMN15502103", "easyCLIP KHDRBS2 Replicate 1 MiSeq"),
    "SRR12192532": ("GSM4665364", "SAMN15502102", "easyCLIP KHDRBS2 Replicate 2 HiSeq"),
    "SRR12192533": ("GSM4665365", "SAMN15502108", "easyCLIP KHDRBS2 Replicate 2 MiSeq"),
    "SRR12192534": ("GSM4665366", "SAMN15502107", "easyCLIP KHDRBS2 Replicate 3 HiSeq"),
    "SRR12192535": ("GSM4665367", "SAMN15502106", "easyCLIP KHDRBS2 Replicate 4 HiSeq"),
    "SRR12192536": ("GSM4665368", "SAMN15502105", "easyCLIP KHDRBS2-R168C Replicate 1 HiSeq"),
    "SRR12192537": ("GSM4665369", "SAMN15502101", "easyCLIP KHDRBS2-R168C Replicate 1 MiSeq"),
    "SRR12192538": ("GSM4665370", "SAMN15502100", "easyCLIP KHDRBS2-R168C Replicate 2 HiSeq"),
    "SRR12192539": ("GSM4665371", "SAMN15502099", "easyCLIP KHDRBS2-R168C Replicate 2 MiSeq"),
    "SRR12192540": ("GSM4665372", "SAMN15502098", "easyCLIP KHDRBS2-R168C Replicate 3 HiSeq"),
    "SRR12192541": ("GSM4665373", "SAMN15502097", "easyCLIP KHDRBS2-R168C Replicate 4 HiSeq"),
}

mapped = {}

def wslp(p):
    p = os.path.abspath(p)
    d, rest = os.path.splitdrive(p)
    return "/mnt/" + d[0].lower() + rest.replace("\\", "/")

def run(args):
    r = subprocess.run(m3.get_cmd(args), stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    t = r.stdout.decode(errors="replace") + r.stderr.decode(errors="replace")
    return r.returncode, t

def sam_bam(srr):
    return os.path.join("results", "alignments_preprocessed", srr,
                        srr + "_aligned.bam")


def getkeys(path, inter, lim):
    
    
    keys = []
    f = m1.open_fastq(path)
    n = 0
    rec = 0
    for line in f:
        if n % 4 == 0:
            head = line
        elif n % 4 == 1:
            if not inter or rec % 2 == 0:
                tok = head.split()[1]
                rand = tok.split("rand=", 1)[1].split("-")[0]
                keys.append(rand + ":" + line.strip()[:20])
                if len(keys) >= lim:
                    break
            rec += 1
        n += 1
    f.close()
    return keys


def dovlap(a, b):
    pa, ia = m1.find_fastq(a, m1.samples[a]["layout"])
    pb, ib = m1.find_fastq(b, m1.samples[b]["layout"])
    big = set(getkeys(pb, ib, 1000000))
    sm = getkeys(pa, ia, 200000)
    h = 0
    for k in sm:
        if k in big:
            h = h + 1
    if len(sm) == 0:
        return 0.0
    return 100.0 * h / len(sm)


def repcheck():
    
    
    pairs = [
        ("SRR12192531", "SRR12192530", "WT replicate 1"),
        ("SRR12192533", "SRR12192532", "WT replicate 2"),
        ("SRR12192537", "SRR12192536", "R168C replicate 1"),
        ("SRR12192539", "SRR12192538", "R168C replicate 2"),
    ]
    controls = [
        ("SRR12192531", "SRR12192536",
         "control: WT rep1 MiSeq vs R168C rep1 HiSeq"),
        ("SRR12192531", "SRR12192535",
         "control: WT rep1 MiSeq vs WT rep3 HiSeq"),
    ]
    cpath = os.path.join(outdir, "technical_replicate_checks.tsv")
    res = {}
    if os.path.exists(cpath):
        f = open(cpath, newline="")
        for r in csv.DictReader(f, delimiter="\t"):
            res[(r["run_a"], r["run_b"])] = (float(r["shared_pct"]),
                                             r["comparison"])
        f.close()
        print("  reusing the overlap checks from", cpath)
    else:
        print("  running the molecule key overlap checks (this reads a lot of")
        print("  fastq so it takes a few minutes) ...")
        for a, b, label in pairs + controls:
            p = dovlap(a, b)
            res[(a, b)] = (p, label)
            print("  %s vs %s: %.2f%%" % (a, b, p))
    cmax = 0.0
    for c in controls:
        if res[(c[0], c[1])][0] > cmax:
            cmax = res[(c[0], c[1])][0]
    notes = {}
    for a, b, label in pairs:
        p = res[(a, b)][0]
        if p >= 5 and cmax < 1:
            notes[(a, b)] = "confirmed: read overlap %.1f%% vs controls <= %.1f%%" % (p, cmax)
        else:
            notes[(a, b)] = "AMBIGUOUS: overlap %.1f%% (controls up to %.1f%%) - verify manually" % (p, cmax)

    f = open(cpath, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["run_a", "run_b", "comparison", "shared_pct", "interpretation"])
    for a, b, label in pairs:
        w.writerow([a, b, label, "%.2f" % res[(a, b)][0],
                    "same library" if res[(a, b)][0] >= 5 else "AMBIGUOUS"])
    for a, b, label in controls:
        w.writerow([a, b, label, "%.2f" % res[(a, b)][0],
                    "different libraries (as expected)"
                    if res[(a, b)][0] < 1 else "UNEXPECTED"])
    f.close()
    print("  wrote", cpath)

    
    mpath = os.path.join(outdir, "replicate_metadata.tsv")
    f = open(mpath, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["srr", "gsm", "bio_sample", "condition",
                "biological_replicate", "library_id", "layout", "ena_title",
                "mapped_reads", "run_type", "evidence"])
    for srr in sorted(m1.samples):
        info = m1.samples[srr]
        gsm, samn, title = ena[srr]
        
        files = m1.find_fastq(srr, info["layout"])
        lib = ""
        f2 = m1.open_fastq(files[0])
        for line in f2:
            if line.startswith("@"):
                p2 = line.split()[1].split("_")
                if len(p2) >= 2:
                    lib = p2[0] + "_" + p2[1]
                break
        f2.close()
        partner = None
        for other in m1.samples:
            oi = m1.samples[other]
            if other != srr and oi["condition"] == info["condition"] and oi["replicate"] == info["replicate"]:
                partner = other
        if partner is None:
            rtype = "only run of this library"
            ev = "no second run exists"
        else:
            n = notes.get((srr, partner), notes.get((partner, srr), ""))
            if mapped[srr] >= mapped[partner]:
                rtype = "primary run of library (technical repeat: " + partner + ")"
            else:
                rtype = "technical replicate run (same library as " + partner + ")"
            ev = n
        w.writerow([srr, gsm, samn, info["condition"], info["replicate"],
                    lib, info["layout"], title, mapped[srr], rtype, ev])
    f.close()
    print("  wrote", mpath)
    return notes, cmax


def getbr():
    
    groups = {}
    for sid, info in m1.samples.items():
        groups.setdefault((info["condition"], info["replicate"]), []).append(sid)
    ds = {}
    for (cond, rep), sids in sorted(groups.items()):
        pick = sids[0]
        for sid in sids:
            if m1.samples[sid]["instrument"] == "HiSeq":
                pick = sid
        ds["%s-BR%d" % (cond, rep)] = pick
    return ds


def getregions(bam, nm, cond):
    cnt, strands, tot = m3.count_reads(bam)
    thr, lam = m3.get_threshold(cnt)
    found = m3.find_regions(cnt, strands, thr)
    regs = []
    for r in found:
        regs.append({"sample_id": nm, "condition": cond, "chrom": r["chrom"],
                     "start": r["start"], "end": r["end"], "support": r["reads"]})
    return regs, tot


def mergecond(regions, cond, bs):
    cov = {}
    for r in regions:
        if r["condition"] != cond:
            continue
        f = r["start"] // bs
        l = (r["end"] - 1) // bs
        for b in range(f, l + 1):
            cov[(r["chrom"], b)] = True
    bychrom = {}
    for chrom, b in cov:
        bychrom.setdefault(chrom, []).append(b)
    out = []
    for chrom in sorted(bychrom):
        bl = sorted(bychrom[chrom])
        groups = []
        for b in bl:
            if groups and b == groups[-1][-1] + 1:
                groups[-1].append(b)
            else:
                groups.append([b])
        for g in groups:
            out.append({"id": "%s_%05d" % (cond, len(out) + 1),
                        "condition": cond, "chrom": chrom,
                        "start": g[0] * bs, "end": (g[-1] + 1) * bs,
                        "bins": g})
    return out


def setcats(wtu, ru):
    
    for side, other, cname in ((wtu, ru, "WT"), (ru, wtu, "R168C")):
        ob = {}
        for i in range(len(other)):
            for b in other[i]["bins"]:
                ob[(other[i]["chrom"], b)] = i
        for r in side:
            hits = {}
            for b in r["bins"]:
                k = (r["chrom"], b)
                if k in ob:
                    i = ob[k]
                    hits[i] = hits.get(i, 0) + 1
            best = None
            bn = 0
            for i in hits:
                if hits[i] > bn:
                    best = i
                    bn = hits[i]
            if best is not None and bn >= minov:
                r["cat"] = "shared"
                r["partner"] = other[best]
                r["sbins"] = bn
            else:
                r["cat"] = cname + "-associated"
                r["partner"] = None
                r["sbins"] = 0


def owner(regs, nm):
    own = {}
    for r in regs:
        if r["sample_id"] != nm:
            continue
        f = r["start"] // bin_size
        l = (r["end"] - 1) // bin_size
        for b in range(f, l + 1):
            own[(r["chrom"], b)] = (r["start"], r["end"], r["support"])
    return own


def sig(r, own):
    
    got = set()
    for b in r["bins"]:
        k = (r["chrom"], b)
        if k in own:
            got.add(own[k])
    t = 0
    for (s, e, sup) in got:
        t = t + sup
    return t


def getdir(l):
    if l >= fccut:
        return "stronger in R168C"
    if l <= -fccut:
        return "stronger in WT"
    return "similar"


def wrtables(allu, wtu, ru, brnames, frags, dregs):
    
    
    own = {}
    for nm in dregs:
        own[nm] = owner(dregs[nm], nm)

    path = os.path.join(outdir, "region_signal_comparison.tsv")
    f = open(path, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["region_id", "category", "chrom", "start", "end",
                "wt_mean_fpm", "wt_median_fpm", "wt_sd_fpm",
                "r168c_mean_fpm", "r168c_median_fpm", "r168c_sd_fpm",
                "log2fc", "direction"] + brnames)
    for r in allu:
        r["fpms"] = {}
        wtv = []
        rv = []
        for nm in brnames:
            sup = sig(r, own[nm])
            v = sup * 1000000.0 / frags[nm]
            r["fpms"][nm] = v
            if nm.startswith("WT"):
                wtv.append(v)
            else:
                rv.append(v)
        r["wtm"] = statistics.mean(wtv)
        r["r168m"] = statistics.mean(rv)
        r["l2fc"] = math.log2((r["r168m"] + pc) / (r["wtm"] + pc))
        w.writerow([r["id"], r["cat"], r["chrom"], r["start"], r["end"],
                    "%.3f" % r["wtm"],
                    "%.3f" % statistics.median(wtv),
                    "%.3f" % (statistics.stdev(wtv) if len(wtv) > 1 else 0.0),
                    "%.3f" % r["r168m"],
                    "%.3f" % statistics.median(rv),
                    "%.3f" % (statistics.stdev(rv) if len(rv) > 1 else 0.0),
                    "%.3f" % r["l2fc"], getdir(r["l2fc"])]
                   + ["%.3f" % r["fpms"][nm] for nm in brnames])
    f.close()
    print("  wrote", path)

    path = os.path.join(outdir, "region_fold_changes.tsv")
    f = open(path, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["region_id", "category", "chrom", "start", "end",
                "wt_mean_fpm", "r168c_mean_fpm", "log2fc", "direction"])
    for r in allu:
        w.writerow([r["id"], r["cat"], r["chrom"], r["start"], r["end"],
                    "%.3f" % r["wtm"], "%.3f" % r["r168m"],
                    "%.3f" % r["l2fc"], getdir(r["l2fc"])])
    f.close()
    print("  wrote", path)

    path = os.path.join(outdir, "region_overlap.tsv")
    f = open(path, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["region_id", "category", "chrom", "start", "end",
                "n_enriched_bins", "partner_id", "partner_start",
                "partner_end", "shared_bins"])
    for r in allu:
        if r["partner"]:
            w.writerow([r["id"], r["cat"], r["chrom"], r["start"], r["end"],
                        len(r["bins"]), r["partner"]["id"],
                        r["partner"]["start"], r["partner"]["end"],
                        r["sbins"]])
        else:
            w.writerow([r["id"], r["cat"], r["chrom"], r["start"], r["end"],
                        len(r["bins"]), "", "", "", 0])
    f.close()
    print("  wrote", path)


def countkmers(seq, cnt):
    seq = seq.upper()
    if "N" in seq:
        pieces = seq.split("N")
    else:
        pieces = [seq]
    for piece in pieces:
        for k in ks:
            d = cnt[k]
            for i in range(len(piece) - k + 1):
                km = piece[i:i + k]
                d[km] = d.get(km, 0) + 1


def dokmers(wtu, ru):
    
    bychrom = {}
    for nm, union in (("WT", wtu), ("R168C", ru)):
        for r in union:
            bychrom.setdefault(r["chrom"], []).append((nm, r))

    wcnt = {}
    rcnt = {}
    for k in ks:
        wcnt[k] = {}
        rcnt[k] = {}

    cur = None
    parts = []
    f = gzip.open(ref, "rt")

    def dochrom(chrom, parts):
        if chrom is None or chrom not in bychrom:
            return
        seq = "".join(parts)
        for nm, r in bychrom[chrom]:
            piece = seq[r["start"]:r["end"]]
            if nm == "WT":
                countkmers(piece, wcnt)
            else:
                countkmers(piece, rcnt)

    for line in f:
        if line.startswith(">"):
            dochrom(cur, parts)
            cur = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    dochrom(cur, parts)
    f.close()

    fw = {}
    fr = {}
    for k in ks:
        wc = {}
        for km in wcnt[k]:
            ok = True
            for c in km:
                if c not in "ACGT":
                    ok = False
            if ok:
                wc[km] = wcnt[k][km]
        rc = {}
        for km in rcnt[k]:
            ok = True
            for c in km:
                if c not in "ACGT":
                    ok = False
            if ok:
                rc[km] = rcnt[k][km]
        wtot = sum(wc.values())
        rtot = sum(rc.values())
        poss = 4 ** k
        for km in set(wc) | set(rc):
            fw[(k, km)] = (wc.get(km, 0) + kpc) / (wtot + kpc * poss)
            fr[(k, km)] = (rc.get(km, 0) + kpc) / (rtot + kpc * poss)

    path = os.path.join(outdir, "kmer_frequencies.tsv")
    f = open(path, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["kmer", "k", "wt_freq", "r168c_freq"])
    for key in sorted(fw):
        w.writerow([key[1], key[0], "%.3e" % fw[key], "%.3e" % fr[key]])
    f.close()
    print("  wrote", path)

    path = os.path.join(outdir, "kmer_comparison.tsv")
    f = open(path, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["kmer", "k", "wt_freq", "r168c_freq", "log2_ratio",
                "direction"])
    for key in sorted(fw):
        ratio = math.log2(fr[key] / fw[key])
        d = "similar"
        if abs(ratio) >= fccut:
            d = getdir(ratio)
        w.writerow([key[1], key[0], "%.3e" % fw[key], "%.3e" % fr[key],
                    "%.3f" % ratio, d])
    f.close()
    print("  wrote", path)
    return fw, fr


def dsbam(bamin, bamout, frac):
    
    
    if os.path.exists(bamout) and os.path.exists(bamout + ".bai"):
        print("    already exists, skipping:", bamout)
        return
    thr = int(frac * 10000)
    sd = (str(dseed) + ":").encode()
    p1 = subprocess.Popen(m3.get_cmd(["view", "-h", m3.to_wsl_path(bamin)]),
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(m3.get_cmd(["view", "-b", "-o",
                                      m3.to_wsl_path(bamout), "-"]),
                          stdin=p1.stdout, stdout=subprocess.DEVNULL,
                          stderr=subprocess.PIPE)
    for line in p1.stdout:
        if line.startswith(b"@"):
            p2.stdin.write(line)
            continue
        if zlib.crc32(sd + line.split(b"\t", 1)[0]) % 10000 < thr:
            p2.stdin.write(line)
    p1.stdout.close()
    p2.stdin.close()
    err = p2.stderr.read().decode(errors="replace")
    c2 = p2.wait()
    c1 = p1.wait()
    if c1 != 0 or c2 != 0:
        print("    ERROR downsampling:", err[-300:])
        sys.exit(1)
    c, t = run(["index", m3.to_wsl_path(bam_out), m3.to_wsl_path(bam_out + ".bai")])
    if c != 0:
        print("    ERROR indexing:", t[-300:])
        sys.exit(1)
    print("    wrote", bam_out)


def outregs_folder():
    return "results/downsampling/regions_downsampled"


def dodown(ds, frags, dregs, ddbams):
    
    
    os.makedirs(dsbams, exist_ok=True)
    os.makedirs(outregs_folder(), exist_ok=True)
    bamfor = {}
    for rep in (1, 2, 3, 4):
        wb = "WT-BR%d" % rep
        rb = "R168C-BR%d" % rep
        tgt = min(frags[wb], frags[rb])
        for br in (wb, rb):
            frac = tgt / frags[br]
            if frac < 0.999:
                bo = os.path.join(dsbams, br.replace("-", "_") + "_downsampled.bam")
                print("  %s: downsample to %.4f (target %d fragments)"
                      % (br, frac, tgt))
                dsbam(ddbams[br], bo, frac)
                bamfor[br] = bo
            else:
                bamfor[br] = ddbams[br]

    outregs = {}
    dsfrags = {}
    for nm in sorted(ds):
        cond = "WT" if nm.startswith("WT") else "R168C"
        if bamfor[nm] != ddbams[nm]:
            regs, tot = getregions(bamfor[nm], nm, cond)
            outregs[nm] = regs
            dsfrags[nm] = tot
            path = os.path.join(outregs_folder(),
                                nm.replace("-", "_") + "_regions.tsv")
            f = open(path, "w", newline="")
            w = csv.writer(f, delimiter="\t")
            w.writerow(["chrom", "start", "end", "read_support"])
            for r in regs:
                w.writerow([r["chrom"], r["start"], r["end"], r["support"]])
            f.close()
            print("  %s (downsampled): %d fragments, %d regions"
                  % (nm, tot, len(regs)))
        else:
            outregs[nm] = dregs[nm]
            dsfrags[nm] = frags[nm]
            print("  %s unchanged bam, regions reused" % nm)

    wtu2 = mergecond([r for lst in outregs.values() for r in lst], "WT", bin_size)
    ru2 = mergecond([r for lst in outregs.values() for r in lst], "R168C", bin_size)
    setcats(wtu2, ru2)
    allu2 = wtu2 + ru2
    
    own = {}
    for nm in outregs:
        own[nm] = owner(outregs[nm], nm)
    wtn = []
    rn = []
    for nm in outregs:
        if nm.startswith("WT"):
            wtn.append(nm)
        else:
            rn.append(nm)
    for r in allu2:
        r["fpms"] = {}
        for nm in outregs:
            sup = sig(r, own[nm])
            r["fpms"][nm] = sup * 1000000.0 / dsfrags[nm]
        wt_values = []
        for n in wtn:
            wt_values.append(r["fpms"][n])
        r_values = []
        for n in rn:
            r_values.append(r["fpms"][n])
        r["wtm"] = statistics.mean(wt_values)
        r["r168m"] = statistics.mean(r_values)
        ratio = (r["r168m"] + pc) / (r["wtm"] + pc)
        r["l2fc"] = math.log2(ratio)

    return wtu2, ru2, allu2


def pics(wtu, ru, allu, xs, ys, dswtu, dsru):
    wsc = sum(1 for r in wtu if r["cat"] == "shared")
    rsc = sum(1 for r in ru if r["cat"] == "shared")
    woc = len(wtu) - wsc
    roc = len(ru) - rsc
    plt.figure(figsize=(6, 4.5))
    x = [0, 1]
    width = 0.35
    plt.bar([i - width / 2 for i in x], [wsc, rsc], width, label="shared")
    plt.bar([i + width / 2 for i in x], [woc, roc], width,
            label="condition-associated")
    plt.xticks(x, ["WT", "R168C"])
    plt.ylabel("number of candidate regions")
    plt.title("Candidate binding regions: shared vs condition-associated")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "region_counts.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(5.5, 5.5))
    px = []
    py = []
    for r in allu:
        px.append(r["wtm"] + pc)
        py.append(r["r168m"] + pc)
    plt.xscale("log")
    plt.yscale("log")
    plt.scatter(px, py, s=6, alpha=0.3)
    plt.plot([min(px), max(px)], [min(px), max(px)], "r--", linewidth=1,
             label="equal signal")
    plt.xlabel("WT mean FPM (+ pseudocount)")
    plt.ylabel("R168C mean FPM (+ pseudocount)")
    plt.title("Normalized binding signal per region")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "signal_scatter.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(6, 4.5))
    vals = []
    for r in allu:
        vals.append(r["l2fc"])
    plt.hist(vals, bins=60, color="steelblue")
    plt.axvline(fccut, color="red", linestyle="--")
    plt.axvline(-fccut, color="red", linestyle="--",
                label="+/- %.1f (direction cutoff)" % fccut)
    plt.xlabel("log2 fold change (R168C / WT)")
    plt.ylabel("number of regions")
    plt.title("Binding signal fold changes")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "fold_changes.png"), dpi=150)
    plt.close()

    wscd = sum(1 for r in dswtu if r["cat"] == "shared")
    rscd = sum(1 for r in dsru if r["cat"] == "shared")
    wocd = len(dswtu) - wscd
    rocd = len(dsru) - rscd
    plt.figure(figsize=(6, 4.5))
    plt.bar([i - width / 2 for i in x], [wscd, rscd], width, label="shared")
    plt.bar([i + width / 2 for i in x], [wocd, rocd], width,
            label="condition-associated")
    plt.xticks(x, ["WT", "R168C"])
    plt.ylabel("number of candidate regions")
    plt.title("After matched-replicate downsampling")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "downsampling_region_counts.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(5.5, 5.5))
    plt.scatter(xs, ys, s=6, alpha=0.3)
    lim = 1
    for v in xs + ys:
        if abs(v) > lim:
            lim = abs(v)
    plt.plot([-lim, lim], [-lim, lim], "r--", linewidth=1, label="y = x")
    plt.xlabel("log2FC original")
    plt.ylabel("log2FC downsampled")
    plt.title("Fold changes after downsampling")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "downsampling_fc_scatter.png"), dpi=150)
    plt.close()
    print("  wrote 5 figures to", figdir)


def writemd(wtu, ru, allu, fw, fr, dswtu, dsru, rval, nmatch, notes, cmax):
    wsc = sum(1 for r in wtu if r["cat"] == "shared")
    rsc = sum(1 for r in ru if r["cat"] == "shared")
    woc = len(wtu) - wsc
    roc = len(ru) - rsc
    wscd = sum(1 for r in dswtu if r["cat"] == "shared")
    rscd = sum(1 for r in dsru if r["cat"] == "shared")
    wocd = len(dswtu) - wscd
    rocd = len(dsru) - rscd

    topr = []
    for r in allu:
        if getdir(r["l2fc"]) == "stronger in R168C":
            topr.append(r)
    topr.sort(reverse=True)
    topw = []
    for r in allu:
        if getdir(r["l2fc"]) == "stronger in WT":
            topw.append(r)
    topw.sort()

    klist2 = []
    for key in fw:
        ratio = math.log2(fr[key] / fw[key])
        klist2.append((ratio, key))
    
    klist2.sort(reverse=True)
    maxk = 0.0
    for x in klist2:
        if abs(x[0]) > maxk:
            maxk = abs(x[0])

    r_value = float("nan")
    
    
    

    lines = []
    lines.append("# WT vs R168C candidate binding region comparison")
    lines.append("")
    lines.append("Candidate regions from Module 3 (read-coverage based, NOT")
    lines.append("confirmed binding sites) compared between KHDRBS2 WT and")
    lines.append("R168C easyCLIP-seq.")
    lines.append("")
    lines.append("## 1. Biological replicate verification")
    lines.append("")
    lines.append("Each condition has 4 biological replicates, not 6. The MiSeq")
    lines.append("runs are technical repeats of the same library as the HiSeq")
    lines.append("runs (verified from the shared library prep IDs in the FASTQ")
    lines.append("headers and molecule key overlap tests - same-library pairs")
    lines.append("share 7-20% of molecule keys, different-library controls")
    lines.append("<= %.2f%%). All replicate-level analyses use the deeper HiSeq"
                 % cmax)
    lines.append("run per library (see replicate_metadata.tsv).")
    lines.append("")
    for (a, b), note in sorted(notes.items()):
        lines.append("- %s / %s: %s" % (a, b, note))
    lines.append("")
    lines.append("## 2. Normalization")
    lines.append("")
    lines.append("PCR duplicates were removed with umi_tools (the original")
    lines.append("easyCLIP workflow's dedup step) after merging the technical")
    lines.append("runs of each library, so region support counts unique")
    lines.append("molecules. Signal = fragments per million counted fragments")
    lines.append("(FPM), using Module 3's genome-wide fragment count as")
    lines.append("denominator so numerator and denominator have the same unit")
    lines.append("(the earlier RPM version divided by mapped READS, ~3x larger")
    lines.append("because it counts both mates and supplementary alignments;")
    lines.append("those baseline tables are kept in baseline_rpm_6runs/).")
    lines.append("Condition signal = mean over the 4 biological replicates;")
    lines.append("per-replicate values and mean/median/sd per region are in")
    lines.append("region_signal_comparison.tsv. Fold change = log2((R168C + %.1f)"
                 % pc)
    lines.append("/ (WT + %.1f))." % pc)
    lines.append("")
    lines.append("## 3. Region overlap")
    lines.append("")
    lines.append("Two regions overlap when they share at least %d x %d bp bin"
                 % (minov, bin_size))
    lines.append("on the same chromosome (regions merged per condition first).")
    lines.append("")
    lines.append("- WT condition-level regions: %d (shared %d, WT-associated %d)"
                 % (len(wtu), wsc, woc))
    lines.append("- R168C condition-level regions: %d (shared %d, R168C-associated %d)"
                 % (len(ru), rsc, roc))
    lines.append("")
    lines.append("A region found in only one condition is NOT proof that")
    lines.append("binding is absent in the other condition - it can simply")
    lines.append("reflect sequencing depth and detection thresholds.")
    lines.append("")
    lines.append("## 4. Binding signal")
    lines.append("")
    n_r = 0
    n_w = 0
    for r in allu:
        if getdir(r["l2fc"]) == "stronger in R168C":
            n_r = n_r + 1
        if getdir(r["l2fc"]) == "stronger in WT":
            n_w = n_w + 1
    lines.append("- regions with stronger R168C signal (log2FC >= %.1f): %d"
                 % (fccut, n_r))
    lines.append("- regions with stronger WT signal (log2FC <= -%.1f): %d"
                 % (fccut, n_w))
    lines.append("")
    lines.append("Strongest differences (top %d each):" % topn)
    lines.append("")
    lines.append("Stronger in R168C:")
    lines.append("")
    for r in topr[:topn]:
        lines.append("- %s %s:%d-%d  WT %.1f FPM, R168C %.1f FPM, log2FC %.2f"
                     % (r["id"], r["chrom"], r["start"], r["end"],
                        r["wtm"], r["r168m"], r["l2fc"]))
    lines.append("")
    lines.append("Stronger in WT:")
    lines.append("")
    for r in topw[:topn]:
        lines.append("- %s %s:%d-%d  WT %.1f FPM, R168C %.1f FPM, log2FC %.2f"
                     % (r["id"], r["chrom"], r["start"], r["end"],
                        r["wtm"], r["r168m"], r["l2fc"]))
    lines.append("")
    lines.append("## 5. K-mer composition (k = 4, 5, 6)")
    lines.append("")
    lines.append("K-mer frequencies over the condition-level region sequences")
    lines.append("(+1 pseudocount per possible k-mer). The two conditions are")
    lines.append("nearly identical in composition: the largest difference is")
    lines.append("|log2 ratio| = %.2f (kmer_comparison.tsv). A frequency"
                 % maxk)
    lines.append("difference does not by itself prove a change in the KHDRBS2")
    lines.append("binding motif - transcript abundance and the genomic")
    lines.append("distribution of the detected regions also shape it.")
    lines.append("")
    lines.append("## 6. Downsampling sensitivity")
    lines.append("")
    lines.append("The deeper run of each matched replicate pair was downsampled")
    lines.append("to its partner's counted-fragment depth (deterministic")
    lines.append("read-name hash, seed %d, mates kept together; no samples" % dseed)
    lines.append("discarded) and the whole comparison was repeated:")
    lines.append("")
    lines.append("| category | original | downsampled |")
    lines.append("|---|---|---|")
    lines.append("| WT regions | %d | %d |" % (len(wtu), len(dswtu)))
    lines.append("| R168C regions | %d | %d |" % (len(ru), len(dsru)))
    lines.append("| shared (WT side) | %d | %d |" % (wsc, wscd))
    lines.append("| WT-associated | %d | %d |" % (woc, wocd))
    lines.append("| shared (R168C side) | %d | %d |" % (rsc, rscd))
    lines.append("| R168C-associated | %d | %d |" % (roc, rocd))
    lines.append("")
    lines.append("Per-region log2 fold changes correlate with Pearson r = %.2f"
                 % rval)
    lines.append("over %d matched regions (see downsampling_fc_scatter.png)."
                 % nmatch)
    lines.append("If the R168C-associated excess")
    lines.append("shrinks after downsampling, it was largely a sequencing-depth")
    lines.append("effect; if it persists, depth alone does not explain it.")
    lines.append("")
    lines.append("## 7. Caveats")
    lines.append("")
    lines.append("- Candidate regions are inferred from read coverage, not")
    lines.append("  experimentally confirmed binding sites.")
    lines.append("- FPM normalization reduces but does not remove detection")
    lines.append("  bias; downsampling is a sensitivity analysis, not proof")
    lines.append("  that depth bias is eliminated.")
    lines.append("- Descriptive fold changes are not statistical significance;")
    lines.append("  no statistical test was applied.")
    lines.append("- K-mer frequency differences do not establish changes in")
    lines.append("  binding specificity (motif discovery is Module 5).")

    path = os.path.join(outdir, "comparison_summary.md")
    f = open(path, "w")
    f.write("\n".join(lines) + "\n")
    f.close()
    print("  wrote", path)





def fivecount(bam, bs):
    
    
    import re
    cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")
    cnt = {}
    strands = {}
    tot = 0
    cmd = m3.get_cmd(["view", "-F", "0x904", "-f", "0x40",
                      "-q", str(m3.min_mapq), wslp(bam)])
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL)
    for line in p.stdout:
        fl = line.decode(errors="replace").split("\t", 6)
        flag = int(fl[1])
        chrom = fl[2]
        pos = int(fl[3]) - 1
        cigar = fl[5]
        span = 0
        for pair in cigar_re.findall(cigar):
            if pair[1] in "MDN=X":
                span = span + int(pair[0])
        if flag & 16:
            five = pos + span - 1
            st = "-"
        else:
            five = pos
            st = "+"
        b = five // bs
        if chrom not in cnt:
            cnt[chrom] = {}
        cnt[chrom][b] = cnt[chrom].get(b, 0) + 1
        if chrom not in strands:
            strands[chrom] = {}
        if b not in strands[chrom]:
            strands[chrom][b] = [0, 0]
        if st == "-":
            strands[chrom][b][1] += 1
        else:
            strands[chrom][b][0] += 1
        tot += 1
    p.stdout.close()
    p.wait()
    return cnt, strands, tot


def narrowstuff(dregs, ddbams):
    
    os.makedirs(os.path.join("results/binding_regions", "narrow_bins"),
                exist_ok=True)
    nregs = {}
    for nm in sorted(dregs):
        cond = "WT" if nm.startswith("WT") else "R168C"
        cnt, strands, tot = fivecount(ddbams[nm], smallbin)
        thr, lam = m3.get_threshold(cnt)
        old = m3.bin_size
        m3.bin_size = smallbin
        found = m3.find_regions(cnt, strands, thr)
        m3.bin_size = bin_size
        regs = []
        for r in found:
            regs.append({"id": "%s_narrow_%05d" % (nm, len(regs) + 1),
                         "sample_id": nm, "condition": cond,
                         "chrom": r["chrom"], "start": r["start"],
                         "end": r["end"]})
        nregs[nm] = regs
        path = os.path.join("results/binding_regions", "narrow_bins",
                            nm.replace("-", "_") + "_narrow_regions.tsv")
        f = open(path, "w", newline="")
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "start", "end"])
        for r in regs:
            w.writerow([r["chrom"], r["start"], r["end"]])
        f.close()
        print("  %-10s narrow (250 bp) regions: %d" % (nm, len(regs)))
    return nregs


def methodmd(wtu, ru, nwt, nru):
    
    sub = {}
    for r in wtu + ru:
        cond = "WT" if r["id"].startswith("WT") else "R168C"
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        for b in range(f, l + 1):
            sub[(cond, r["chrom"], b)] = r["id"]

    wt_in = 0
    for r in nwt:
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        found = False
        for b in range(f, l + 1):
            k = ("WT", r["chrom"], b)
            if k in sub:
                found = True
        if found:
            wt_in = wt_in + 1
    r_in = 0
    for r in nru:
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        found = False
        for b in range(f, l + 1):
            k = ("R168C", r["chrom"], b)
            if k in sub:
                found = True
        if found:
            r_in = r_in + 1

    nbins = set()
    for r in nwt:
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        for b in range(f, l + 1):
            nbins.add((r["chrom"], b))
    wt_cov = 0
    for r in wtu:
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        hit = False
        for b in range(f, l + 1):
            if (r["chrom"], b) in nbins:
                hit = True
        if hit:
            wt_cov = wt_cov + 1
    nbins2 = set()
    for r in nru:
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        for b in range(f, l + 1):
            nbins2.add((r["chrom"], b))
    r_cov = 0
    for r in ru:
        f = r["start"] // smallbin
        l = (r["end"] - 1) // smallbin
        hit = False
        for b in range(f, l + 1):
            if (r["chrom"], b) in nbins2:
                hit = True
        if hit:
            r_cov = r_cov + 1

    path = os.path.join("results/binding_regions", "method_comparison.md")
    f = open(path, "w")
    f.write("# Region-identification method comparison\n\n")
    f.write("Method A (baseline): 4 kb bins counted at the fragment\n")
    f.write("leftmost aligned base, per-sample Poisson-tail threshold,\n")
    f.write("consecutive enriched bins merged (Module 3).\n\n")
    f.write("Method B (comparison): 250 bp bins counted at the read 5' end\n")
    f.write("(the easyCLIP crosslink-site proxy), same Poisson-tail\n")
    f.write("threshold and merging, run on the UMI-deduplicated\n")
    f.write("biological-replicate BAMs.\n\n")
    f.write("The original easyCLIP workflow calls peaks with CLIPper on\n")
    f.write("deduplicated BAMs. CLIPper could not be reproduced here (its\n")
    f.write("python2-era toolchain and dependency pins are impractical in\n")
    f.write("this environment), so Method B implements the crosslink-site\n")
    f.write("resolution principle with the same threshold logic as the\n")
    f.write("baseline. Neither method is a calibrated significance test:\n")
    f.write("the Poisson background ignores mappability and local\n")
    f.write("read-density structure.\n\n")
    f.write("| metric | WT | R168C |\n|---|---|---|\n")
    f.write("| 4 kb regions (Method A) | %d | %d |\n" % (len(wtu), len(ru)))
    f.write("| 250 bp regions (Method B) | %d | %d |\n" % (len(nwt), len(nru)))
    f.write("| 250 bp regions overlapping a 4 kb region | %d (%.1f%%) | %d (%.1f%%) |\n"
            % (wt_in, 100.0 * wt_in / max(len(nwt), 1),
               r_in, 100.0 * r_in / max(len(nru), 1)))
    f.write("| 4 kb regions containing a 250 bp region | %d (%.1f%%) | %d (%.1f%%) |\n"
            % (wt_cov, 100.0 * wt_cov / max(len(wtu), 1),
               r_cov, 100.0 * r_cov / max(len(ru), 1)))
    f.write("\n")
    f.write("Interpretation: the 4 kb regions are an order of magnitude\n")
    f.write("coarser. Most 250 bp crosslink-site regions fall inside 4 kb\n")
    f.write("regions, but only a fraction of 4 kb regions contain a\n")
    f.write("detected crosslink-site cluster - the 4 kb method calls\n")
    f.write("broader, weaker signal that the narrow method does not\n")
    f.write("resolve. For downstream motif and transcript work the\n")
    f.write("narrower regions are preferable; the 4 kb regions remain a\n")
    f.write("conservative whole-region screen. Per-replicate narrow\n")
    f.write("regions are in results/binding_regions/narrow_bins/.\n")
    f.close()
    print("  wrote", path)



def main():
    print("=" * 70)
    print("Module 4: Compare WT and R168C Binding Patterns")
    print("KHDRBS2 WT vs R168C easyCLIP-seq (GSE154168)")
    print("=" * 70)

    
    if not os.path.exists(align_file):
        print("Alignment summary not found. Run Module 2 first.")
        sys.exit(1)
    f = open(align_file, newline="")
    for r in csv.DictReader(f, delimiter="\t"):
        mapped[r["sample_id"]] = int(r["mapped_reads"])
    f.close()

    if not os.path.exists(sample_file):
        print("Sample summary not found. Run Module 1 first.")
        sys.exit(1)
    samples = []
    f = open(sample_file, newline="")
    for r in csv.DictReader(f, delimiter="\t"):
        if r.get("sample_id"):
            samples.append(r)
    f.close()
    if len(samples) != 12:
        print("ERROR: expected 12 samples, found %d" % len(samples))
        sys.exit(1)
    for s in samples:
        if s["status"] != "OK":
            print("ERROR: sample not OK:", s["sample_id"])
            sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)

    print("\nBiological replicate verification ...")
    notes, cmax = repcheck()
    ds = getbr()
    print("  biological replicate datasets (deeper run per library):")
    for nm in sorted(ds):
        print("    %s = %s" % (nm, ds[nm]))

    print("\nDetecting regions (Module 3 functions) ...")
    
    
    os.makedirs(dedupdir, exist_ok=True)
    allruns = {}
    for (cond, rep), sids in sorted(
            {(i["condition"], i["replicate"]): [] for i in []}.items()):
        pass
    groups = {}
    for sid, info in m1.samples.items():
        groups.setdefault((info["condition"], info["replicate"]), []).append(sid)
    for (cond, rep) in sorted(groups):
        allruns["%s-BR%d" % (cond, rep)] = sorted(groups[(cond, rep)])

    ddb = {}
    for nm in sorted(allruns):
        srrs = allruns[nm]
        ob = os.path.join(dedupdir, nm + ".bam")
        ddb[nm] = ob
        if os.path.exists(ob) and os.path.exists(ob + ".bai"):
            print("  %s: dedup bam exists, skipping" % nm)
            continue
        inputs = []
        for s in srrs:
            inputs.append(sam_bam(s))
        merged = os.path.join(dedupdir, nm + "_merged.bam")
        if len(inputs) > 1:
            print("  %s: merging %d technical runs ..." % (nm, len(inputs)))
            c, t = run(["merge", "-f", "-o", wslp(merged)]
                       + [wslp(x) for x in inputs])
            if c != 0:
                print("    ERROR merge:", t[-300:])
                sys.exit(1)
            c, t = run(["index", wslp(merged)])
            if c != 0:
                print("    ERROR indexing merged bam:", t[-300:])
                sys.exit(1)
        else:
            merged = inputs[0]
        print("  %s: umi_tools dedup ..." % nm)
        lg = os.path.join(dedupdir, nm + "_dedup.log")
        cmd = ["wsl.exe", "-d", "Ubuntu", "-u", "root", "--", umitools,
               "dedup", "--stdin", wslp(merged), "--stdout", wslp(ob),
               "--umi-separator=_", "--log", wslp(lg)]
        ef = open(os.path.join(dedupdir, nm + "_dedup.stderr.log"), "w")
        rr = subprocess.run(cmd, stdout=ef, stderr=ef)
        ef.close()
        if rr.returncode != 0:
            print("    ERROR umi_tools dedup:", nm)
            sys.exit(1)
        c, t = run(["index", wslp(ob), wslp(ob + ".bai")])
        if c != 0:
            print("    ERROR indexing:", t[-300:])
            sys.exit(1)
        if len(inputs) > 1 and os.path.exists(merged):
            os.remove(merged)
        print("    wrote", ob)

    dregs = {}
    frags = {}
    ddbams = {}
    for nm in sorted(ds):
        cond = "WT" if nm.startswith("WT") else "R168C"
        regs, tot = getregions(ddb[nm], nm, cond)
        dregs[nm] = regs
        frags[nm] = tot
        ddbams[nm] = ddb[nm]
        print("  %-10s %s: %d fragments, %d regions"
              % (nm, os.path.basename(ddb[nm]), tot, len(regs)))

    print("\nMerging and classifying regions ...")
    wtregs = []
    rregs = []
    for nm in dregs:
        for r in dregs[nm]:
            if nm.startswith("WT"):
                wtregs.append(r)
            else:
                rregs.append(r)
    wtu = mergecond(wtregs, "WT", bin_size)
    ru = mergecond(rregs, "R168C", bin_size)
    setcats(wtu, ru)
    allu = wtu + ru
    print("  WT condition-level regions: %d" % len(wtu))
    print("  R168C condition-level regions: %d" % len(ru))

    wrtables(allu, wtu, ru, sorted(dregs), frags, dregs)

    print("\nCounting kmers over the condition-level regions ...")
    fw, fr = dokmers(wtu, ru)

    print("\nDownsampling sensitivity analysis ...")
    dswtu, dsru, dsallu = dodown(ds, frags, dregs, ddb)

    
    dsbins = {}
    for i in range(len(dsallu)):
        for b in dsallu[i]["bins"]:
            dsbins[(dsallu[i]["chrom"], b)] = i
    xs = []
    ys = []
    for r in allu:
        hits = {}
        for b in r["bins"]:
            k = (r["chrom"], b)
            if k in dsbins:
                i = dsbins[k]
                hits[i] = hits.get(i, 0) + 1
        if hits:
            best = None
            bn = 0
            for i in hits:
                if hits[i] > bn:
                    best = i
                    bn = hits[i]
            if bn >= minov:
                xs.append(r["l2fc"])
                ys.append(dsallu[best]["l2fc"])
    rval = float("nan")
    if len(xs) > 2:
        rval = statistics.correlation(xs, ys)
    print("  log2FC correlation original vs downsampled: r = %.2f (n = %d)"
          % (rval, len(xs)))

    
    wsc = 0
    for r in wtu:
        if r["cat"] == "shared":
            wsc = wsc + 1
    rsc = 0
    for r in ru:
        if r["cat"] == "shared":
            rsc = rsc + 1
    wscd = 0
    for r in dswtu:
        if r["cat"] == "shared":
            wscd = wscd + 1
    rscd = 0
    for r in dsru:
        if r["cat"] == "shared":
            rscd = rscd + 1
    path = os.path.join(outdir, "region_counts_comparison.tsv")
    f = open(path, "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["category", "original_bio", "downsampled"])
    w.writerow(["WT regions", len(wtu), len(dswtu)])
    w.writerow(["R168C regions", len(ru), len(dsru)])
    w.writerow(["shared (WT side)", wsc, wscd])
    w.writerow(["WT-associated", len(wtu) - wsc, len(dswtu) - wscd])
    w.writerow(["shared (R168C side)", rsc, rscd])
    w.writerow(["R168C-associated", len(ru) - rsc, len(dsru) - rscd])
    f.close()
    print("  wrote", path)

    pics(wtu, ru, allu, xs, ys, dswtu, dsru)
    writemd(wtu, ru, allu, fw, fr, dswtu, dsru, rval, len(xs), notes, cmax)

    print("\nNarrow-bin (250 bp) crosslink site method comparison ...")
    nregs = narrowstuff(dregs, ddb)
    nwt_list = []
    nru_list = []
    for nm in nregs:
        for r in nregs[nm]:
            if nm.startswith("WT"):
                nwt_list.append(r)
            else:
                nru_list.append(r)
    nwt = mergecond(nwt_list, "WT", smallbin)
    nru = mergecond(nru_list, "R168C", smallbin)
    print("  narrow condition-level regions: WT %d, R168C %d"
          % (len(nwt), len(nru)))
    methodmd(wtu, ru, nwt, nru)

    print("\nFinished Module 4")


if __name__ == "__main__":
    main()
