// SPDX-License-Identifier: GPL-3.0-only
// GU multi-reference genotype preparation. The merge semantics follow
// PrincetonUniversity/IBDmix generate_gt.cc (6a9d4bc67dcaf602e65cac0fa6191e1e7a1291ef).
// See ibdmix.LICENSE. Modern GT is decoded once, with no text-VCF round trip.
#include <htslib/bgzf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Intervals = std::vector<std::pair<int64_t, int64_t>>;
static void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
static std::string chrom_name(std::string s) {
    if (s.rfind("chr", 0) == 0) s.erase(0, 3);
    return s == "X" ? "23" : s;
}

struct Reader {
    std::unique_ptr<bcf_srs_t, decltype(&bcf_sr_destroy)> sr{bcf_sr_init(), bcf_sr_destroy};
    std::string path, chrom, genotypes, blank;
    const Intervals& intervals;
    int32_t* gt = nullptr;
    int capacity = 0, samples = 0;
    int64_t pos = 0, previous = -1;
    uint64_t records = 0;
    bool valid = false;
    char ref = '.', alt = '.';

    Reader(const std::string& file, const std::string& region, const std::string& chr,
           const Intervals& allowed) : path(file), chrom(chr), intervals(allowed) {
        require(bool(sr), "Cannot allocate VCF reader");
        if (!region.empty()) require(bcf_sr_set_regions(sr.get(), region.c_str(), 0) == 0,
                                     "Invalid region: " + region);
        require(bcf_sr_add_reader(sr.get(), file.c_str()) == 1,
                "Cannot open " + file + ": " + bcf_sr_strerror(sr->errnum));
        samples = bcf_hdr_nsamples(header());
        require(samples > 0, "No samples in " + file);
        genotypes.resize(2 * samples);
        blank.resize(2 * samples);
        for (int i = 0; i < samples; ++i) {
            genotypes[2*i+1] = blank[2*i+1] = '\t';
            blank[2*i] = '0';
        }
    }
    ~Reader() { free(gt); }
    bcf_hdr_t* header() { return bcf_sr_get_header(sr.get(), 0); }

    bool next(bool archaic = false) {
        while (bcf_sr_next_line(sr.get())) {
            auto* r = bcf_sr_get_line(sr.get(), 0);
            require(r && !r->errcode, "Invalid VCF record in " + path);
            require(chrom_name(bcf_hdr_id2name(header(), r->rid)) == chrom,
                    "Unexpected chromosome in " + path);
            pos = r->pos + 1;
            require(pos >= previous, "Unsorted VCF in " + path);
            previous = pos;
            if (!intervals.empty() && !std::any_of(intervals.begin(), intervals.end(),
                [&](auto p) { return p.first <= r->pos && r->pos < p.second; })) continue;
            ++records;
            require(bcf_unpack(r, BCF_UN_STR) == 0, "Cannot decode alleles in " + path);
            // Retain modern non-SNPs as merge barriers. Dropping them invents
            // homozygous-reference modern calls at archaic-only positions.
            valid = r->n_allele >= 1 && std::string(r->d.allele[0]).size() == 1 &&
                    (r->n_allele == 1 || (r->n_allele == 2 && std::string(r->d.allele[1]).size() == 1));
            if (!valid) { if (archaic) continue; return true; }
            ref = r->d.allele[0][0];
            alt = r->n_allele == 1 ? '.' : r->d.allele[1][0];
            int n = bcf_get_genotypes(header(), r, &gt, &capacity);
            require(n > 0 && n % samples == 0, "Missing/invalid GT in " + path + ":" + std::to_string(pos));
            int width = n / samples;
            bool informative = false;
            for (int i = 0; i < samples; ++i) {
                const auto* g = gt + i * width;
                int ploidy = 0;
                bool missing = false;
                while (ploidy < width && g[ploidy] != bcf_int32_vector_end) {
                    missing |= bcf_gt_is_missing(g[ploidy]);
                    ++ploidy;
                }
                char code = '9';
                if (!missing && ploidy > 0) {
                    // Match +setGT -t . -n . followed by +fixploidy -f 2:
                    // any missing allele -> missing; haploid -> duplicate;
                    // higher ploidy -> first two alleles.
                    int a = bcf_gt_allele(g[0]), b = bcf_gt_allele(g[ploidy > 1 ? 1 : 0]);
                    // Do not format diagnostics in the per-sample hot loop.
                    if (a < 0 || a > 1 || b < 0 || b > 1)
                        throw std::runtime_error("Unsupported GT allele in " + path + ":" + std::to_string(pos));
                    code = char('0' + a + b);
                    informative = true;
                }
                genotypes[2*i] = code;
            }
            if (!archaic || informative) return true;
        }
        require(!sr->errnum, "VCF read failed: " + path + ": " + bcf_sr_strerror(sr->errnum));
        return false;
    }
};

struct Output {
    Reader archaic;
    std::string name, base, buffer;
    BGZF* file = nullptr;
    bool recheck = false;
    uint64_t rows = 0, hom_ref = 0;
    int64_t last = -1;

    Output(const std::string& label, const std::string& path, const std::string& region,
           const std::string& chr, const Intervals& intervals, const std::string& stem,
           Reader& modern, hts_tpool* pool) : archaic(path, region, chr, intervals), name(label), base(stem) {
        require(archaic.samples == 1, "Expected exactly one archaic sample: " + path);
        std::set<std::string> ids;
        std::string header = "chrom\tpos\tref\talt";
        for (Reader* reader : {&archaic, &modern}) {
            for (int i = 0; i < reader->samples; ++i) {
                std::string id = reader->header()->samples[i];
                require(ids.insert(id).second, "Duplicate genotype sample ID: " + id);
                header += '\t' + id;
            }
        }
        file = bgzf_open((base + ".gt.txt.gz.part").c_str(), "w1");
        require(file, "Cannot create genotype output: " + base);
        // All five compressors share a bounded pool; no per-reference worker
        // explosion and no complete chromosome resident in memory.
        if (bgzf_thread_pool(file, pool, 8) != 0) {
            bgzf_close(file); file = nullptr;
            throw std::runtime_error("Cannot initialize genotype compressor");
        }
        buffer = header + '\n';
        buffer.reserve(1 << 20);
    }
    ~Output() { if (file) bgzf_close(file); }
    void flush() {
        require(bgzf_write(file, buffer.data(), buffer.size()) == static_cast<ssize_t>(buffer.size()),
                "Genotype write failed: " + base);
        buffer.clear();
    }
    void emit(char alt, const std::string& modern) {
        require(archaic.pos > last, "Genotype positions must be strictly increasing: " + name);
        last = archaic.pos;
        ++rows;
        hom_ref += archaic.genotypes[0] == '0';
        buffer += archaic.chrom + '\t' + std::to_string(archaic.pos) + '\t' + archaic.ref + '\t' + alt + '\t';
        buffer += archaic.genotypes;
        buffer += modern;
        buffer += '\n';
        if (buffer.size() >= (1 << 20)) flush();
    }
    void consume(const Reader& modern) {
        // Preserve upstream EOF, allele compatibility and one-to-one equal-POS
        // handling, including archaic-only non-reference sites before modern EOF.
        while (recheck || archaic.next(true)) {
            recheck = false;
            if (archaic.pos < modern.pos) {
                if (archaic.genotypes[0] != '0') emit(archaic.alt, modern.blank);
            } else if (archaic.pos == modern.pos) {
                if (modern.valid && archaic.ref == modern.ref &&
                    (archaic.alt == '.' || archaic.alt == modern.alt)) emit(modern.alt, modern.genotypes);
                break;
            } else { recheck = true; break; }
        }
    }
    void finish(int modern_samples, const std::string& build, const Intervals& intervals) {
        require(rows > 0, "No informative genotype records: " + name);
        flush();
        auto* closing = file;
        file = nullptr;
        require(bgzf_close(closing) == 0, "Genotype close failed: " + base);
        // QC is computed while emitting, replacing a full gzip decode/scan.
        std::ofstream qc(base + ".qc.json.part");
        qc << "{\n  \"rows\": " << rows << ",\n  \"archaic_hom_ref_rows\": " << hom_ref
           << ",\n  \"modern_samples\": " << modern_samples << ",\n  \"scoring_build\": "
           << (build.empty() ? "null" : "\"" + build + "\"") << ",\n  \"allowed_intervals_bed0\": ";
        if (intervals.empty()) qc << "null";
        else {
            qc << '[';
            for (size_t i = 0; i < intervals.size(); ++i)
                qc << (i ? ", " : "") << '[' << intervals[i].first << ", " << intervals[i].second << ']';
            qc << ']';
        }
        qc << "\n}\n";
        qc.close();
        require(bool(qc), "Genotype QC write failed: " + base);
    }
    void commit() {
        fs::rename(base + ".gt.txt.gz.part", base + ".gt.txt.gz");
        fs::rename(base + ".qc.json.part", base + ".qc.json");
    }
};

int main(int argc, char** argv) try {
    std::string modern_file, modern_region, outdir, unit, chrom, bed, build;
    struct Archaic { std::string name, file, region; };
    std::vector<Archaic> refs;
    int threads = 4;
    auto value = [&](int& i) { require(i + 1 < argc, "Missing argument value"); return std::string(argv[++i]); };
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--version") { std::cout << "ibdmix-prepare shared-hts-v1\n"; return 0; }
        else if (a == "--modern") modern_file = value(i);
        else if (a == "--modern-region") modern_region = value(i);
        else if (a == "--out") outdir = value(i);
        else if (a == "--unit") unit = value(i);
        else if (a == "--chrom") chrom = chrom_name(value(i));
        else if (a == "--threads") threads = std::stoi(value(i));
        else if (a == "--nonpar-bed") bed = value(i);
        else if (a == "--build") build = value(i);
        else if (a == "--archaic") { auto name = value(i), file = value(i), region = value(i); refs.push_back({name, file, region}); }
        else throw std::runtime_error("Unknown argument: " + a);
    }
    require(!modern_file.empty() && !outdir.empty() && !unit.empty() && !chrom.empty() && !refs.empty(),
            "Required: --modern VCF --out DIR --unit UNIT --chrom CHR --archaic NAME VCF REGION [...]");
    require(threads >= 1 && threads <= 64, "--threads must be 1..64");
    auto safe = [](const std::string& s) { return !s.empty() && s != "." && s != ".." && s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-") == std::string::npos; };
    require(safe(unit), "Invalid unit name");
    std::set<std::string> names;
    for (const auto& r : refs) require(safe(r.name) && names.insert(r.name).second, "Invalid/duplicate reference name");
    Intervals intervals;
    if (!bed.empty()) {
        require(chrom == "23" && (build == "37" || build == "38"), "Non-PAR filtering requires X and --build 37|38");
        std::ifstream input(bed);
        require(bool(input), "Cannot open non-PAR BED");
        std::string c;
        int64_t lo, hi;
        while (input >> c >> lo >> hi) {
            require(chrom_name(c) == chrom && lo >= 0 && hi > lo && (intervals.empty() || lo >= intervals.back().second), "Invalid non-PAR BED");
            intervals.emplace_back(lo, hi);
        }
        require(input.eof() && !intervals.empty(), "Malformed/empty non-PAR BED");
    } else require(build.empty(), "--build requires --nonpar-bed");
    fs::create_directories(outdir);
    Reader modern(modern_file, modern_region, chrom, intervals);
    std::unique_ptr<hts_tpool, decltype(&hts_tpool_destroy)> pool(hts_tpool_init(threads), hts_tpool_destroy);
    require(bool(pool), "Cannot create compression pool");
    std::vector<std::unique_ptr<Output>> outputs;
    for (const auto& r : refs) outputs.push_back(std::make_unique<Output>(r.name, r.file, r.region, chrom,
        intervals, outdir + '/' + r.name + '.' + unit, modern, pool.get()));
    auto start = std::chrono::steady_clock::now();
    while (modern.next()) for (auto& out : outputs) out->consume(modern);
    for (auto& out : outputs) out->finish(modern.samples, build, intervals);
    // Nothing becomes reusable until every input and every compressor succeeded.
    for (auto& out : outputs) out->commit();
    auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cerr << "[IBDMIX PREPARE] modern_passes=1 modern_records=" << modern.records
              << " samples=" << modern.samples << " references=" << refs.size() << " elapsed_s=" << seconds << '\n';
    for (auto& out : outputs) std::cerr << "[IBDMIX PREPARE] ref=" << out->name << " rows=" << out->rows << '\n';
    return 0;
} catch (const std::exception& e) {
    std::cerr << "ERROR: IBDmix genotype preparation: " << e.what() << '\n';
    return 1;
}
