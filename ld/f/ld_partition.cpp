// McManus et al. (2023), Methods S1, linear LD partition objective.
// Exact DP *subject to* every block's last POS - first POS <= D.
// Q is reported for unordered pairs: sum_blocks sum_i<j (r2_ij - penalty).
// This differs from a symmetric/ordered-pair objective only by a factor of two.
// D is a search restriction, not an r2 threshold or proof of an unrestricted optimum.
// O(n * variants-in-D) time and O(n) RAM for one D. A grid of k distances
// maintains k independent DP states, sharing a single pass through the LD table.
// The PLINK table is read backwards in bounded buffers; it must be in ascending
// PVAR-index (ID_A, ID_B) order.
#include <algorithm>
#include <charconv>
#include <cerrno>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

static vector<string> fields(const string& s) {
    istringstream in(s);
    vector<string> v;
    string x;
    while (in >> x) v.push_back(x);
    return v;
}

static string normalized_chrom(string c) {
    if (c.rfind("chr", 0) == 0) c.erase(0, 3);
    return c;
}

static long long integer(const string& value, const string& label) {
    size_t used = 0;
    long long x = stoll(value, &used);
    if (used != value.size()) throw runtime_error("Invalid integer in " + label);
    return x;
}

static long double real_number(const string& value, const string& label) {
    size_t used = 0;
    long double x = stold(value, &used);
    if (used != value.size() || !isfinite(x))
        throw runtime_error("Nonfinite or invalid value in " + label);
    return x;
}


// Fast path only for the exact numeric-ID, three-column table emitted by
// ld_pipeline.py. Floating-point conversion uses the same strtold function
// underlying std::stold; DP arithmetic, order and precision are unchanged.
static bool ascii_space(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v';
}
static int numeric_id(const char* begin, const char* end, int count) {
    if (begin == end || *begin++ != 'v' || begin == end ||
        *begin < '0' || *begin > '9' || (*begin == '0' && end - begin != 1))
        throw runtime_error("Unknown VCOR variant ID");
    int index = -1;
    auto parsed = from_chars(begin, end, index);
    if (parsed.ec != errc() || parsed.ptr != end || index < 0 || index >= count)
        throw runtime_error("Unknown VCOR variant ID");
    return index;
}
static void numeric_row(const string& line, int count, int& a, int& b, long double& r) {
    const char* cursor = line.c_str();
    const char* limit = cursor + line.size();
    const char* starts[3];
    const char* ends[3];
    for (int i = 0; i != 3; ++i) {
        while (cursor != limit && ascii_space(*cursor)) ++cursor;
        starts[i] = cursor;
        while (cursor != limit && !ascii_space(*cursor)) ++cursor;
        ends[i] = cursor;
        if (starts[i] == ends[i]) throw runtime_error("Malformed VCOR row");
    }
    while (cursor != limit && ascii_space(*cursor)) ++cursor;
    if (cursor != limit) throw runtime_error("Malformed VCOR row");
    a = numeric_id(starts[0], ends[0], count);
    b = numeric_id(starts[1], ends[1], count);
    char* end = nullptr;
    errno = 0;
    r = strtold(starts[2], &end);
    if (end == starts[2] || end != ends[2] || errno == ERANGE || !isfinite(r))
        throw runtime_error("Nonfinite or invalid value in PHASED_R2");
}

// LDZST001: independently checksummed Zstandard frames, a protected index,
// and a completion footer. All integers are little endian. Only storage changes;
// decompressed bytes, line order and the DP arithmetic below remain unchanged.
#include <zstd.h>
#include <memory>
#include <array>
namespace ldz {
constexpr uint64_t block_size = 4 * 1024 * 1024;
struct Entry { uint64_t offset, compressed, raw; };
static void require(bool ok, const char* message) { if (!ok) throw runtime_error(message); }
static uint64_t get64(const char* s) {
    uint64_t v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | static_cast<unsigned char>(s[i]);
    return v;
}
static void put64(string& s, uint64_t v) {
    for (int i = 0; i < 8; ++i) { s.push_back(static_cast<char>(v & 255)); v >>= 8; }
}
static uint64_t hash_bytes(const string& s, uint64_t h = 14695981039346656037ULL) {
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}
static string read_at(ifstream& in, uint64_t offset, size_t size) {
    in.clear(); in.seekg(static_cast<streamoff>(offset));
    string s(size, '\0'); in.read(s.data(), static_cast<streamsize>(size));
    require(in && static_cast<size_t>(in.gcount()) == size, "Truncated LDZ or read failure");
    return s;
}
static void zcheck(size_t code) {
    if (ZSTD_isError(code)) throw runtime_error(string("Zstandard: ") + ZSTD_getErrorName(code));
}
struct Input {
    ifstream in;
    bool compressed;
    uint64_t bytes;
    vector<Entry> entries;
    // Batch adjacent compressed frames into bounded reads. On WSL/DrvFS each
    // seek/read has significant latency; one tiny read per frame loses most of
    // the storage benefit. This cache holds bytes verbatim, in original order.
    string packed_cache;
    uint64_t cache_at = 0;
    explicit Input(const string& path) : in(path, ios::binary), compressed(false), bytes(0) {
        require(bool(in), "Cannot open LD input");
        in.seekg(0, ios::end); auto end = in.tellg();
        require(end >= 0, "LD input must be seekable"); bytes = static_cast<uint64_t>(end);
        if (bytes < 8 || read_at(in, 0, 8) != "LDZST001") return;
        compressed = true;
        require(bytes >= 56, "Truncated LDZ header/footer");
        string header = read_at(in, 0, 16), footer = read_at(in, bytes - 40, 40);
        require(get64(header.data() + 8) == block_size, "Unsupported LDZ block size");
        require(footer.substr(0, 8) == "LDZEND01", "Missing LDZ completion footer");
        uint64_t index_at = get64(footer.data() + 8), count = get64(footer.data() + 16);
        uint64_t total = get64(footer.data() + 24);
        require(index_at >= 16 && index_at <= bytes - 40 && count <= (bytes - 56) / 24,
                "Invalid LDZ index bounds");
        require(count * 24 == bytes - 40 - index_at && count <= 10000000,
                "Invalid or excessive LDZ index size");
        string index = read_at(in, index_at, static_cast<size_t>(count * 24));
        auto hash = hash_bytes(footer.substr(0, 32), hash_bytes(index, hash_bytes(header)));
        require(hash == get64(footer.data() + 32), "LDZ index/footer checksum mismatch");
        uint64_t offset = 16, raw = 0;
        for (uint64_t i = 0; i < count; ++i) {
            const char* p = index.data() + i * 24;
            Entry e{get64(p), get64(p + 8), get64(p + 16)};
            require(e.offset == offset && e.compressed > 0 && e.compressed <= ZSTD_compressBound(block_size) &&
                    e.compressed <= index_at - offset && e.raw > 0 && e.raw <= block_size &&
                    (i + 1 == count || e.raw == block_size), "Invalid LDZ block index");
            offset += e.compressed; raw += e.raw; entries.push_back(e);
        }
        require(offset == index_at && raw == total, "LDZ length/count mismatch");
    }
    string block(size_t i) {
        const auto& e = entries.at(i);
        if (packed_cache.empty() || e.offset < cache_at ||
            e.offset + e.compressed > cache_at + packed_cache.size()) {
            uint64_t end = e.offset + e.compressed;
            if (packed_cache.empty()) {
                cache_at = e.offset; // Header inspection needs only its own frame.
            } else if (e.offset < cache_at) {
                cache_at = end > 16 + 32 * 1024 * 1024 ? end - 32 * 1024 * 1024 : 16;
            } else {
                cache_at = e.offset; // Also keep forward --unpack reads linear.
                const auto& last = entries.back();
                end = min<uint64_t>(last.offset + last.compressed, cache_at + 32 * 1024 * 1024);
            }
            packed_cache = read_at(in, cache_at, static_cast<size_t>(end - cache_at));
        }
        string packed = packed_cache.substr(static_cast<size_t>(e.offset - cache_at),
                                            static_cast<size_t>(e.compressed));
        // The frame descriptor's checksum flag is mandatory for this format.
        require(packed.size() >= 5 && (static_cast<unsigned char>(packed[4]) & 4),
                "LDZ frame checksum missing");
        require(ZSTD_getFrameContentSize(packed.data(), packed.size()) == e.raw,
                "LDZ frame raw length mismatch");
        size_t frame = ZSTD_findFrameCompressedSize(packed.data(), packed.size()); zcheck(frame);
        require(frame == packed.size(), "LDZ block contains trailing/multiple frames");
        string raw(static_cast<size_t>(e.raw), '\0');
        size_t n = ZSTD_decompress(raw.data(), raw.size(), packed.data(), packed.size()); zcheck(n);
        require(n == raw.size(), "LDZ decompressed length mismatch");
        return raw;
    }
    string header() {
        string s;
        if (compressed) {
            if (entries.empty()) throw runtime_error("Empty VCOR");
            s = block(0);
        } else s = read_at(in, 0, static_cast<size_t>(min<uint64_t>(bytes, block_size)));
        auto end = s.find('\n');
        require(end != string::npos || s.size() < block_size, "VCOR header exceeds 4 MiB");
        return s.substr(0, end);
    }
};
static void pack(const string& source, const string& dest) {
    require(!filesystem::exists(dest), "Refusing to overwrite packed output");
    ifstream in(source, ios::binary); require(bool(in), "Cannot open compression input");
    ofstream out(dest, ios::binary); require(bool(out), "Cannot create compression output");
    string header = "LDZST001"; put64(header, block_size); out.write(header.data(), header.size());
    unique_ptr<ZSTD_CCtx, decltype(&ZSTD_freeCCtx)> ctx(ZSTD_createCCtx(), ZSTD_freeCCtx);
    require(bool(ctx), "Cannot allocate Zstandard context");
    zcheck(ZSTD_CCtx_setParameter(ctx.get(), ZSTD_c_compressionLevel, 1));
    zcheck(ZSTD_CCtx_setParameter(ctx.get(), ZSTD_c_checksumFlag, 1));
    string raw(block_size, '\0'), packed(ZSTD_compressBound(block_size), '\0'), index;
    uint64_t offset = 16, total = 0, count = 0;
    while (true) {
        in.read(raw.data(), raw.size()); size_t n = static_cast<size_t>(in.gcount());
        require(!in.bad() && (in.eof() || !in.fail()), "LD compression input read failure");
        if (!n) break;
        size_t c = ZSTD_compress2(ctx.get(), packed.data(), packed.size(), raw.data(), n); zcheck(c);
        out.write(packed.data(), c); require(bool(out), "LD compressed output write failure");
        put64(index, offset); put64(index, c); put64(index, n);
        offset += c; total += n; ++count;
    }
    string footer = "LDZEND01"; put64(footer, offset); put64(footer, count); put64(footer, total);
    put64(footer, hash_bytes(footer, hash_bytes(index, hash_bytes(header))));
    out.write(index.data(), index.size()); out.write(footer.data(), footer.size()); out.close();
    require(bool(out), "LD compressed output close failure");
}
static void unpack(const string& source, const string& dest) {
    Input input(source); require(input.compressed, "Expected LDZ input");
    require(!filesystem::exists(dest), "Refusing to overwrite unpacked output");
    ofstream out(dest, ios::binary); require(bool(out), "Cannot create unpacked output");
    for (size_t i = 0; i < input.entries.size(); ++i) {
        string raw = input.block(i); out.write(raw.data(), raw.size());
        require(bool(out), "Unpacked output write failure");
    }
    out.close(); require(bool(out), "Unpacked output close failure");
}
} // namespace ldz

struct ReverseLines {
    ldz::Input input;
    streamoff pos;
    size_t next_block;
    string buf;
    explicit ReverseLines(const string& path) : input(path), pos(input.bytes), next_block(input.entries.size()) {
        if (input.compressed) pos = next_block;
    }
    bool next(string& s) {
        for (;;) {
            auto k = buf.rfind('\n');
            if (k != string::npos) {
                s = buf.substr(k + 1);
                buf.resize(k);
                if (!s.empty() && s.back() == '\r') s.pop_back();
                if (!s.empty()) return true;
                continue;
            }
            if (!pos) {
                if (buf.empty()) return false;
                s = move(buf);
                buf.clear();
                if (!s.empty() && s.back() == '\r') s.pop_back();
                return !s.empty();
            }
            string chunk;
            if (input.compressed) {
                chunk = input.block(--next_block);
                pos = next_block;
            } else {
                auto n = min<streamoff>(4 * 1024 * 1024, pos);
                pos -= n;
                chunk = ldz::read_at(input.in, pos, static_cast<size_t>(n));
            }
            buf = chunk + buf;
            if (buf.size() > 16 * 1024 * 1024)
                throw runtime_error("VCOR line exceeds 16 MiB; expected tabular text");
        }
    }
};

struct PartitionState {
    long long distance;
    string out_path, qc_path;
    vector<long double> score, best;
    vector<int> block_end, window_end;
    int hi;
    uint64_t pairs_used = 0, pairs_expected = 0;
    PartitionState(long long d, const string& out, const string& qc, int n)
        : distance(d), out_path(out), qc_path(qc), score(n, 0), best(n + 1, 0),
          block_end(n), window_end(n), hi(n - 1) {}
};

int main(int argc, char** argv) {
    try {
        if (argc == 4 && string(argv[1]) == "--pack") { ldz::pack(argv[2], argv[3]); return 0; }
        if (argc == 4 && string(argv[1]) == "--unpack") { ldz::unpack(argv[2], argv[3]); return 0; }
        const unordered_set<string> allowed = {
            "--chrom", "--max-distance-bp", "--penalty", "--pvar", "--vcor",
            "--out", "--qc", "--ld-mode", "--distance-grid-kb", "--out-dir",
            "--coincident-policy"
        };
        unordered_map<string, string> opt;
        for (int i = 1; i < argc; i += 2) {
            if (i + 1 == argc) throw runtime_error("Missing option value");
            if (!allowed.count(argv[i])) throw runtime_error("Unknown option " + string(argv[i]));
            if (!opt.emplace(argv[i], argv[i + 1]).second)
                throw runtime_error("Repeated option " + string(argv[i]));
        }
        auto req = [&](const string& k) {
            if (!opt.count(k)) throw runtime_error("Missing " + k);
            return opt[k];
        };
        const string chrom = normalized_chrom(req("--chrom"));
        vector<long long> distances;
        const bool grid_mode = opt.count("--distance-grid-kb");
        if (grid_mode) {
            if (opt.count("--max-distance-bp") || opt.count("--out") || opt.count("--qc"))
                throw runtime_error("Grid mode cannot be combined with single-distance/output options");
            req("--out-dir");
            const string grid = req("--distance-grid-kb");
            if (grid.empty() || grid.back() == ',') throw runtime_error("Invalid distance grid");
            istringstream input(grid);
            string token;
            while (getline(input, token, ',')) {
                const long long kb = integer(token, "distance grid");
                if (kb <= 0 || kb > numeric_limits<long long>::max() / 1000 ||
                    (!distances.empty() && kb * 1000 <= distances.back()))
                    throw runtime_error("Distance grid must contain increasing positive integer kb values");
                distances.push_back(kb * 1000);
            }
            if (distances.size() > 64) throw runtime_error("Distance grid has too many states (maximum 64)");
        } else {
            if (opt.count("--out-dir")) throw runtime_error("--out-dir requires --distance-grid-kb");
            distances.push_back(integer(req("--max-distance-bp"), "D"));
            req("--out");
            req("--qc");
        }
        const long double pen = real_number(req("--penalty"), "penalty");
        const string mode = opt.count("--ld-mode") ? opt.at("--ld-mode") : "strict";
        const string coincident = opt.count("--coincident-policy") ? opt.at("--coincident-policy") : "group";
        if (distances.empty() || distances[0] < 0 || pen < 0 || pen >= 1)
            throw runtime_error("Invalid D or penalty");
        if (mode != "strict" && mode != "sparse")
            throw runtime_error("--ld-mode must be strict or sparse");
        if (coincident != "group" && coincident != "error")
            throw runtime_error("--coincident-policy must be group or error");

        ifstream pv(req("--pvar"));
        if (!pv) throw runtime_error("Cannot open PVAR");
        vector<long long> pos;
        vector<string> names;
        unordered_map<string, int> ids;
        string line;
        bool pvar_header = false;
        uint64_t duplicate_position_groups = 0, extra_colocated_records = 0;
        while (getline(pv, line)) {
            if (line.empty()) continue;
            if (line.rfind("##", 0) == 0) continue;
            if (line[0] == '#') {
                auto v = fields(line);
                if (v.size() < 5 || v[0] != "#CHROM" || v[1] != "POS" ||
                    v[2] != "ID" || v[3] != "REF" || v[4] != "ALT" || pvar_header)
                    throw runtime_error("Invalid or repeated PVAR header");
                pvar_header = true;
                continue;
            }
            if (!pvar_header) throw runtime_error("PVAR #CHROM/POS/ID/REF/ALT header is required");
            auto v = fields(line);
            if (v.size() < 5) throw runtime_error("Invalid PVAR row");
            if (normalized_chrom(v[0]) != chrom) throw runtime_error("Unexpected PVAR chromosome");
            const long long p = integer(v[1], "PVAR POS");
            if (p < 1 || (!pos.empty() && p < pos.back())) throw runtime_error("Unsorted or invalid PVAR POS");
            if (!pos.empty() && p == pos.back()) {
                ++extra_colocated_records;
                if (pos.size() < 2 || pos[pos.size() - 2] != p) ++duplicate_position_groups;
                if (coincident == "error")
                    throw runtime_error("Co-located variants violate the paper's injective-coordinate assumption at POS " +
                        to_string(p) + ". --coincident-policy group retains all records and disallows cuts within a coordinate.");
            }
            if (pos.size() >= static_cast<size_t>(numeric_limits<int>::max()))
                throw runtime_error("Too many variants for this build");
            if (v[2] == "." || !ids.emplace(v[2], static_cast<int>(pos.size())).second)
                throw runtime_error("Missing or duplicate PVAR variant ID");
            if (v[4] == "." || v[4].find(',') != string::npos)
                throw runtime_error("Expected retained biallelic variants in PVAR");
            pos.push_back(p);
            names.push_back(v[2]);
        }
        if (pv.bad()) throw runtime_error("PVAR read failed");
        const int n = static_cast<int>(pos.size());
        if (!n) throw runtime_error("No retained variants");

        ReverseLines rev(req("--vcor"));
        line = rev.input.header();
        if (line.empty()) throw runtime_error("Empty VCOR");
        auto header = fields(line);
        if (header.empty() || header[0].empty() || header[0][0] != '#')
            throw runtime_error("VCOR must begin with a PLINK column header");
        header[0].erase(0, 1);
        unordered_set<string> header_seen;
        for (const auto& h : header)
            if (!header_seen.insert(h).second) throw runtime_error("Duplicate VCOR column " + h);
        auto optional_col = [&](const string& name) {
            auto it = find(header.begin(), header.end(), name);
            return it == header.end() ? -1 : static_cast<int>(it - header.begin());
        };
        auto required_col = [&](const string& name) {
            int c = optional_col(name);
            if (c < 0) throw runtime_error("Missing VCOR column " + name);
            return c;
        };
        const int ca = required_col("ID_A"), cb = required_col("ID_B");
        const int cr = required_col("PHASED_R2");
        const int cpa = optional_col("POS_A"), cpb = optional_col("POS_B");
        const int cca = optional_col("CHROM_A"), ccb = optional_col("CHROM_B");

        bool numeric_fast_path = header == vector<string>{"ID_A", "ID_B", "PHASED_R2"};
        for (int i = 0; numeric_fast_path && i < n; ++i)
            numeric_fast_path = names[i] == "v" + to_string(i);
        int a = -1, b = -1, prevA = n, prevB = n;
        long double r = 0;
        bool have = false;
        uint64_t pairs_read = 0;
        auto advance = [&]() {
            have = false;
            while (rev.next(line)) {
                if (line.empty() || line[0] == '#') continue;
                vector<string> v;
                if (numeric_fast_path) {
                    numeric_row(line, n, a, b, r);
                } else {
                    v = fields(line);
                    if (v.size() != header.size()) throw runtime_error("Malformed VCOR row");
                    auto ia = ids.find(v[ca]), ib = ids.find(v[cb]);
                    if (ia == ids.end() || ib == ids.end()) throw runtime_error("Unknown VCOR variant ID");
                    a = ia->second;
                    b = ib->second;
                }
                if (a >= b || a > prevA || (a == prevA && b >= prevB))
                    throw runtime_error("VCOR must contain unique upper-triangle pairs ordered by PVAR index (ID_A, ID_B)");
                prevA = a;
                prevB = b;
                if ((cpa >= 0 && integer(v[cpa], "VCOR POS_A") != pos[a]) ||
                    (cpb >= 0 && integer(v[cpb], "VCOR POS_B") != pos[b]) ||
                    (cca >= 0 && normalized_chrom(v[cca]) != chrom) ||
                    (ccb >= 0 && normalized_chrom(v[ccb]) != chrom))
                    throw runtime_error("VCOR coordinate does not match PVAR");
                if (!numeric_fast_path) r = real_number(v[cr], "PHASED_R2");
                if (r < 0 || r > 1) throw runtime_error("PHASED_R2 must be in [0, 1]");
                ++pairs_read;
                have = true;
                return;
            }
        };
        advance();
        vector<PartitionState> states;
        states.reserve(distances.size());
        for (long long d : distances) {
            string out = grid_mode ? (filesystem::path(req("--out-dir")) /
                ("D" + to_string(d / 1000) + "kb") / ("chr" + chrom + ".bed")).string() : req("--out");
            string qc = grid_mode ? (filesystem::path(req("--out-dir")) /
                ("D" + to_string(d / 1000) + "kb") / ("chr" + chrom + ".qc.tsv")).string() : req("--qc");
            states.emplace_back(d, out, qc, n);
        }
        vector<long double> weights(n, 0);
        vector<int> observed;
        for (int i = n - 1; i >= 0; --i) {
            for (auto& state : states) {
                while (pos[state.hi] - pos[i] > state.distance) --state.hi;
                state.window_end[i] = state.hi;
            }
            const int max_hi = states.back().hi;
            fill(weights.begin() + i, weights.begin() + max_hi + 1, 0);
            observed.clear();
            while (have && a == i) {
                if (b <= max_hi) {
                    weights[b] = r;
                    observed.push_back(b); // Reverse input => strictly descending indices.
                }
                advance();
            }
            const uint64_t row_expected_max = static_cast<uint64_t>(max_hi - i);
            if (mode == "strict" && observed.size() != row_expected_max)
                throw runtime_error("Incomplete strict LD table for " + names[i] + ": expected " +
                    to_string(row_expected_max) + " pairs within maximum D, found " + to_string(observed.size()) +
                    ". Do not replace missing LD with zero; check --ld-window-r2 -1, NaN and table truncation.");
            for (auto& state : states) {
                const uint64_t row_expected = static_cast<uint64_t>(state.hi - i);
                const uint64_t row_used = mode == "strict" ? row_expected :
                    static_cast<uint64_t>(observed.end() - lower_bound(observed.begin(), observed.end(),
                        state.hi, greater<int>()));
                state.pairs_used += row_used;
                state.pairs_expected += row_expected;
                // score[j] on entry is W(i+1,j). Add the i-th row to get W(i,j).
                state.score[i] = 0;
                // The paper assumes unique coordinates. Our explicit 'group'
                // extension retains co-located records, preserving every pairwise
                // weight, and forbids a physical cut between them. No jitter or
                // variant deletion is performed. Forced groups can have negative Q.
                const bool singleton_end_allowed = i == n - 1 || pos[i] != pos[i + 1];
                state.best[i] = singleton_end_allowed ? state.best[i + 1] :
                    -numeric_limits<long double>::infinity();
                state.block_end[i] = singleton_end_allowed ? i : -1;
                long double acc = 0;
                for (int j = i + 1; j <= state.hi; ++j) {
                    acc += weights[j] - pen;
                    state.score[j] += acc;
                    const long double candidate = state.score[j] + state.best[j + 1];
                    // Deterministic ties prefer the shortest first block. No numerical
                    // tolerance is allowed to hide a strictly better objective.
                    const bool end_allowed = j == n - 1 || pos[j] != pos[j + 1];
                    if (end_allowed && candidate > state.best[i]) {
                        state.best[i] = candidate;
                        state.block_end[i] = j;
                    }
                }
            }
        }
        if (have) throw runtime_error("Unconsumed LD rows");
        // No output is created until all requested D states and input validations
        // have completed. Callers should move outputs atomically after return 0.
        for (const auto& state : states) {
            if (!isfinite(state.best[0])) throw runtime_error("Nonfinite DP objective");
            for (int i = 0; i < n; i = state.block_end[i] + 1) {
                const int j = state.block_end[i];
                if (j + 1 < n && pos[j] == pos[j + 1])
                    throw runtime_error("Optimal partition splits co-located variants at POS " + to_string(pos[j]) +
                        "; this cannot be represented by nonoverlapping BED intervals. Inspect input representation.");
            }
        }
        for (const auto& state : states) {
            vector<pair<int, int>> blocks;
            long long max_span = 0;
            int max_variants = 0, at_limit = 0, near_limit = 0, frontier_blocks = 0;
            for (int i = 0; i < n; i = state.block_end[i] + 1) {
                const int j = state.block_end[i];
                blocks.emplace_back(i, j);
                const long long span = pos[j] - pos[i];
                max_span = max(max_span, span);
                max_variants = max(max_variants, j - i + 1);
                at_limit += span == state.distance;
                near_limit += state.distance > 0 && static_cast<long double>(span) >= 0.95L * state.distance;
                frontier_blocks += j == state.window_end[i] && j < n - 1;
            }
            if (grid_mode) filesystem::create_directories(filesystem::path(state.out_path).parent_path());
            ofstream out(state.out_path);
            if (!out) throw runtime_error("Cannot create BED");
            for (auto [i, j] : blocks) out << chrom << '\t' << pos[i] - 1 << '\t' << pos[j] << '\n';
            out.close();
            if (!out) throw runtime_error("BED write failed");
            ofstream qc(state.qc_path);
            if (!qc) throw runtime_error("Cannot create QC");
            qc << setprecision(numeric_limits<long double>::max_digits10);
            qc << "chrom\tvariants\tpairs\tblocks\tobjective\tpenalty\tmax_distance_bp\tld_mode"
              "\tpairs_within_D\texpected_pairs_within_D\tmissing_pairs_within_D\tpairs_beyond_D"
              "\tmax_block_span_bp\tmax_block_variant_count\tblocks_at_distance_limit"
              "\tblocks_at_95pct_distance\tblocks_at_search_frontier\tobjective_pair_convention"
              "\tcoincident_policy\tduplicate_position_groups\textra_colocated_records\n"
           << chrom << '\t' << n << '\t' << pairs_read << '\t' << blocks.size() << '\t'
           << state.best[0] << '\t' << pen << '\t' << state.distance << '\t' << mode << '\t'
           << state.pairs_used << '\t' << state.pairs_expected << '\t' << state.pairs_expected - state.pairs_used << '\t'
           << pairs_read - state.pairs_used << '\t' << max_span << '\t' << max_variants << '\t'
           << at_limit << '\t' << near_limit << '\t' << frontier_blocks << "\tunordered\t" << coincident << '\t'
           << duplicate_position_groups << '\t' << extra_colocated_records << '\n';
            qc.close();
            if (!qc) throw runtime_error("QC write failed");
        }
        return 0;
    } catch (const exception& e) {
        cerr << "ld_partition: " << e.what() << '\n';
        return 2;
    }
}
