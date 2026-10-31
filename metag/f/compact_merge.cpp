// SPDX-License-Identifier: GPL-3.0-or-later
// Compact, transactional merge for the pinned MetaCache species classifier.
// Candidate admission, rank aliases and LCA thresholds follow MetaCache
// d7646eca4c4dc131262b16d2910923fce3f5d4fc (André Müller and Robin Kobus).
// This program is distributed WITHOUT ANY WARRANTY. See
// https://www.gnu.org/licenses/gpl-3.0.html for the license terms.
//
// Build: g++ -O3 -std=c++17 -Wall -Wextra -pedantic compact_merge.cpp -o compact_merge
//
// State: 64-byte header, then N records of {u64 tax1,u64 tax2,u32 hits1,u32 hits2}.
// Delta: 96-byte header, then changed records {u64 query_id, StateRecord}.
// All integers, including exported taxids and optional header hashes, are LE.
// UPDATE only reads the committed state. APPLY replays absolute replacements;
// the caller must retain its pending receipt and delta until it commits progress.

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

constexpr std::uint64_t FNV_OFFSET = 14695981039346656037ULL;
constexpr std::uint64_t FNV_PRIME = 1099511628211ULL;
constexpr std::uint64_t DIRTY = 1;
constexpr std::uint64_t DELTA_READY = 1, DELTA_PAIRED = 2, DELTA_HEADERS = 4;
constexpr std::uint32_t VERSION = 1;

struct StateHeader {
    char magic[8];
    std::uint32_t version, record_bytes;
    std::uint64_t fragments, applied_parts, state_id, flags, transaction, reserved;
};
struct Record {
    std::uint64_t tax1, tax2;
    std::uint32_t hits1, hits2;
};
struct DeltaHeader {
    char magic[8];
    std::uint32_t version, record_bytes;
    std::uint64_t fragments, state_id, base_parts, transaction, changed_rows;
    std::uint64_t rows, physical_reads, payload_hash, flags, reserved;
};
struct DeltaRecord { std::uint64_t qid; Record value; };
static_assert(sizeof(StateHeader) == 64, "state header layout");
static_assert(sizeof(Record) == 24, "state record layout");
static_assert(sizeof(DeltaHeader) == 96, "delta header layout");
static_assert(sizeof(DeltaRecord) == 32, "delta record layout");

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }
std::string syserr(const std::string& context) { return context + ": " + std::strerror(errno); }

std::string_view trim(std::string_view text) {
    constexpr std::string_view spaces = " \t\r\n";
    const auto first = text.find_first_not_of(spaces);
    if (first == std::string_view::npos) return {};
    return text.substr(first, text.find_last_not_of(spaces) - first + 1);
}
std::uint64_t number(std::string_view text, const std::string& what) {
    text = trim(text);
    std::uint64_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        fail("Invalid " + what + ": " + std::string(text.substr(0, 100)));
    return value;
}
std::uint64_t hash_bytes(const void* ptr, std::size_t size, std::uint64_t hash = FNV_OFFSET) {
    const auto* bytes = static_cast<const unsigned char*>(ptr);
    for (std::size_t i = 0; i < size; ++i) { hash ^= bytes[i]; hash *= FNV_PRIME; }
    return hash;
}
std::uint64_t random_id() {
    std::random_device random;
    std::uint64_t value = (std::uint64_t(random()) << 32) ^ std::uint64_t(random());
    return value ? value : 1;
}
std::uint64_t byte_size(std::uint64_t count, std::uint64_t width, std::uint64_t header = 0) {
    const auto limit = std::min<std::uint64_t>(std::numeric_limits<std::size_t>::max(),
                                             std::numeric_limits<off_t>::max());
    if (count > (limit - header) / width) fail("Requested file is too large for this platform");
    return header + count * width;
}
void require_little_endian() {
    const std::uint32_t one = 1;
    if (*reinterpret_cast<const unsigned char*>(&one) != 1)
        fail("This binary format requires a little-endian host");
}

class File {
public:
    int fd = -1;
    explicit File(const std::string& path, int flags, mode_t mode = 0600) {
        fd = ::open(path.c_str(), flags | O_CLOEXEC, mode);
        if (fd < 0) fail(syserr("Cannot open " + path));
    }
    ~File() { if (fd >= 0) ::close(fd); }
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    std::uint64_t size() const {
        struct stat st {};
        if (::fstat(fd, &st) != 0) fail(syserr("fstat"));
        if (!S_ISREG(st.st_mode) || st.st_size < 0) fail("State/delta/hash input must be a regular file");
        return static_cast<std::uint64_t>(st.st_size);
    }
    void lock(bool write) const {
        if (::flock(fd, (write ? LOCK_EX : LOCK_SH) | LOCK_NB) != 0)
            fail(syserr("State is already in use"));
    }
    void sync() const { if (::fsync(fd) != 0) fail(syserr("fsync")); }
};
void write_all(int fd, const void* ptr, std::size_t size) {
    const char* bytes = static_cast<const char*>(ptr);
    while (size) {
        const ssize_t wrote = ::write(fd, bytes, size);
        if (wrote < 0) { if (errno == EINTR) continue; fail(syserr("write")); }
        if (wrote == 0) fail("write made no progress");
        bytes += wrote; size -= static_cast<std::size_t>(wrote);
    }
}
void pwrite_all(int fd, const void* ptr, std::size_t size, off_t position) {
    const char* bytes = static_cast<const char*>(ptr);
    while (size) {
        const ssize_t wrote = ::pwrite(fd, bytes, size, position);
        if (wrote < 0) { if (errno == EINTR) continue; fail(syserr("pwrite")); }
        if (wrote == 0) fail("pwrite made no progress");
        bytes += wrote; size -= static_cast<std::size_t>(wrote); position += wrote;
    }
}

class Mapping {
public:
    File file;
    void* data = MAP_FAILED;
    std::size_t bytes = 0;
    Mapping(const std::string& path, bool write, bool lock = false)
        : file(path, write ? O_RDWR : O_RDONLY) {
        if (lock) file.lock(write);
        bytes = static_cast<std::size_t>(file.size());
        if (!bytes) fail("Cannot mmap an empty file: " + path);
        data = ::mmap(nullptr, bytes, PROT_READ | (write ? PROT_WRITE : 0), MAP_SHARED, file.fd, 0);
        if (data == MAP_FAILED) fail(syserr("mmap " + path));
    }
    ~Mapping() { if (data != MAP_FAILED) ::munmap(data, bytes); }
    void sync(std::size_t length = 0) {
        if (::msync(data, length ? length : bytes, MS_SYNC) != 0) fail(syserr("msync"));
        file.sync();
    }
};

void validate_state_header(const Mapping& map, bool allow_dirty = false) {
    if (map.bytes < sizeof(StateHeader)) fail("Truncated state header");
    const auto& h = *static_cast<const StateHeader*>(map.data);
    if (std::memcmp(h.magic, "MCTOP2V1", 8) || h.version != VERSION ||
        h.record_bytes != sizeof(Record) || !h.fragments || !h.state_id || h.reserved)
        fail("Unsupported or corrupt compact state header");
    if (h.flags & ~DIRTY) fail("Unknown state flags");
    if (h.flags && !allow_dirty) fail("State has an unfinished apply; replay its pending delta before continuing");
    if (map.bytes != byte_size(h.fragments, sizeof(Record), sizeof(StateHeader)))
        fail("State size does not match its fragment count");
}
bool valid_record(const Record& value) {
    if ((value.tax1 == 0) != (value.hits1 == 0) || (value.tax2 == 0) != (value.hits2 == 0)) return false;
    if (value.tax1 > std::uint64_t(std::numeric_limits<std::int64_t>::max()) ||
        value.tax2 > std::uint64_t(std::numeric_limits<std::int64_t>::max())) return false;
    if (value.hits1 < value.hits2 || (value.tax2 && value.tax1 == value.tax2)) return false;
    return true;
}
bool same_record(const Record& a, const Record& b) {
    return a.tax1 == b.tax1 && a.tax2 == b.tax2 && a.hits1 == b.hits1 && a.hits2 == b.hits2;
}

// Exact K=2 specialization of upstream's ordered, distinct-taxon insertion.
void insert_candidate(Record& value, std::uint64_t taxid, std::uint32_t hits) {
    if (value.tax2 && value.hits2 >= hits) return;
    if (taxid == value.tax1) {
        if (hits > value.hits1) value.hits1 = hits;
    } else if (taxid == value.tax2) {
        if (hits > value.hits2) {
            value.hits2 = hits;
            if (value.hits2 > value.hits1) {
                std::swap(value.tax1, value.tax2);
                std::swap(value.hits1, value.hits2);
            }
        }
    } else if (!value.tax1 || hits > value.hits1) {
        value.tax2 = value.tax1; value.hits2 = value.hits1;
        value.tax1 = taxid; value.hits1 = hits;
    } else {
        value.tax2 = taxid; value.hits2 = hits;
    }
}

class Bitmap {
    std::vector<std::uint64_t> words;
public:
    explicit Bitmap(std::uint64_t n) : words(static_cast<std::size_t>((n + 63) / 64), 0) {}
    bool mark(std::uint64_t zero_based) {
        const auto bit = std::uint64_t(1) << (zero_based % 64);
        auto& word = words[static_cast<std::size_t>(zero_based / 64)];
        if (word & bit) return false;
        word |= bit; return true;
    }
};

using Args = std::map<std::string, std::string>;
Args arguments(int argc, char** argv) {
    Args args;
    for (int i = 2; i < argc; i += 2) {
        const std::string key = argv[i];
        if (key.rfind("--", 0) != 0 || i + 1 >= argc) fail("Options require --name VALUE");
        if (!args.emplace(key, argv[i + 1]).second) fail("Duplicate option: " + key);
    }
    return args;
}
std::string required(const Args& args, const std::string& key) {
    const auto found = args.find(key);
    if (found == args.end() || found->second.empty()) fail("Missing option " + key);
    return found->second;
}
std::string option(const Args& args, const std::string& key, const std::string& fallback) {
    const auto found = args.find(key);
    return found == args.end() ? fallback : found->second;
}
void allowed(const Args& args, std::initializer_list<std::string> keys) {
    for (const auto& entry : args)
        if (std::find(keys.begin(), keys.end(), entry.first) == keys.end()) fail("Unknown option " + entry.first);
}
void report(const Args& args, const std::string& json) {
    const auto found = args.find("--summary");
    if (found != args.end()) {
        File output(found->second, O_WRONLY | O_CREAT | O_EXCL);
        write_all(output.fd, json.data(), json.size()); write_all(output.fd, "\n", 1); output.sync();
    }
    std::cout << json << '\n';
    if (!std::cout) fail("Cannot write JSON summary");
}

std::string state_summary(const char* command, const StateHeader& h) {
    std::ostringstream out;
    out << "{\"command\":\"" << command << "\",\"schema\":1,\"fragments\":" << h.fragments
        << ",\"applied_parts\":" << h.applied_parts << ",\"state_id\":\"" << h.state_id
        << "\",\"transaction\":\"" << h.transaction << "\",\"dirty\":" << (h.flags ? "true" : "false")
        << ",\"state_bytes\":" << byte_size(h.fragments, sizeof(Record), sizeof(StateHeader)) << "}";
    return out.str();
}
void initialize(const Args& args) {
    allowed(args, {"--state", "--fragments", "--summary"});
    const auto n = number(required(args, "--fragments"), "fragment count");
    if (!n) fail("Cannot initialize a zero-fragment sample");
    File file(required(args, "--state"), O_RDWR | O_CREAT | O_EXCL);
    file.lock(true);
    const auto bytes = byte_size(n, sizeof(Record), sizeof(StateHeader));
    if (::ftruncate(file.fd, static_cast<off_t>(bytes)) != 0) fail(syserr("Initialize state size"));
    StateHeader header {};
    std::memcpy(header.magic, "MCTOP2V1", 8);
    header.version = VERSION; header.record_bytes = sizeof(Record);
    header.fragments = n; header.state_id = random_id();
    pwrite_all(file.fd, &header, sizeof(header), 0); file.sync();
    report(args, state_summary("init", header));
}
void info(const Args& args) {
    allowed(args, {"--state", "--summary"});
    Mapping map(required(args, "--state"), false, true);
    validate_state_header(map, true);
    report(args, state_summary("info", *static_cast<const StateHeader*>(map.data)));
}

void split_columns(std::string_view text, std::vector<std::string_view>& parts) {
    parts.clear();
    std::size_t begin = 0;
    for (;;) {
        const auto end = text.find('|', begin);
        parts.push_back(trim(text.substr(begin, end == std::string_view::npos ? end : end - begin)));
        if (end == std::string_view::npos) return;
        begin = end + 1;
    }
}

void update(const Args& args) {
    allowed(args, {"--state", "--input", "--delta", "--paired", "--headers", "--summary"});
    Mapping state(required(args, "--state"), false, true);
    validate_state_header(state);
    const auto& sh = *static_cast<const StateHeader*>(state.data);
    if (sh.applied_parts == std::numeric_limits<std::uint64_t>::max()) fail("Partition counter overflow");
    const auto* records = reinterpret_cast<const Record*>(static_cast<const char*>(state.data) + sizeof(StateHeader));
    const auto paired_text = option(args, "--paired", "0");
    if (paired_text != "0" && paired_text != "1") fail("--paired must be 0 or 1");
    const bool paired = paired_text == "1";
    const auto expected_reads = byte_size(sh.fragments, paired ? 2 : 1);
    std::unique_ptr<Mapping> headers;
    if (args.count("--headers")) {
        headers = std::make_unique<Mapping>(required(args, "--headers"), false);
        if (headers->bytes != byte_size(sh.fragments, sizeof(std::uint64_t)))
            fail("Header-hash file does not contain exactly N uint64 values");
    }
    const auto* header_hashes = headers ? static_cast<const std::uint64_t*>(headers->data) : nullptr;
    File delta(required(args, "--delta"), O_RDWR | O_CREAT | O_EXCL);
    DeltaHeader dh {};
    std::memcpy(dh.magic, "MCDELTV1", 8);
    dh.version = VERSION; dh.record_bytes = sizeof(DeltaRecord);
    dh.fragments = sh.fragments; dh.state_id = sh.state_id; dh.base_parts = sh.applied_parts;
    dh.transaction = random_id(); dh.payload_hash = FNV_OFFSET;
    write_all(delta.fd, &dh, sizeof(dh));

    const std::string input_name = option(args, "--input", "-");
    std::ifstream input_file;
    std::vector<char> input_buffer(4 * 1024 * 1024);
    if (input_name != "-") {
        input_file.rdbuf()->pubsetbuf(input_buffer.data(), input_buffer.size());
        input_file.open(input_name);
        if (!input_file) fail("Cannot open native candidate input " + input_name);
    }
    std::istream& input = input_name == "-" ? std::cin : input_file;
    Bitmap seen(sh.fragments);
    std::vector<std::string_view> columns;
    std::size_t column_count = 0, qid_col = 0, header_col = 0, hits_col = 0;
    bool have_layout = false, have_summary = false;
    std::uint64_t reported_reads = 0, candidate_count = 0;
    std::string line;
    std::vector<DeltaRecord> output_buffer;
    output_buffer.reserve(32768);
    auto flush_delta = [&]() {
        if (output_buffer.empty()) return;
        const auto bytes = output_buffer.size() * sizeof(DeltaRecord);
        dh.payload_hash = hash_bytes(output_buffer.data(), bytes, dh.payload_hash);
        write_all(delta.fd, output_buffer.data(), bytes); output_buffer.clear();
    };
    while (std::getline(input, line)) {
        const std::string_view text = trim(line);
        if (text.empty()) continue;
        if (text.front() == '#') {
            constexpr std::string_view layout_tag = "# TABLE_LAYOUT:";
            constexpr std::string_view summary_tag = "# queries:";
            if (text.substr(0, layout_tag.size()) == layout_tag) {
                if (have_layout || dh.rows || have_summary) fail("Repeated or misplaced native table layout");
                split_columns(text.substr(layout_tag.size()), columns);
                column_count = columns.size();
                auto locate = [&](std::string_view name) {
                    const auto first = std::find(columns.begin(), columns.end(), name);
                    if (first == columns.end() || std::find(first + 1, columns.end(), name) != columns.end())
                        fail("Native layout requires exactly one " + std::string(name) + " column");
                    return static_cast<std::size_t>(first - columns.begin());
                };
                qid_col = locate("query_id"); header_col = locate("query_header"); hits_col = locate("top_hits");
                if (std::find(columns.begin(), columns.end(), "reference_hits") != columns.end() ||
                    std::find(columns.begin(), columns.end(), "all_hits") != columns.end())
                    fail("Compact mode must omit -refhits and -allhits");
                have_layout = true;
            } else if (text.substr(0, summary_tag.size()) == summary_tag) {
                if (!have_layout || have_summary) fail("Repeated or misplaced native query summary");
                reported_reads = number(text.substr(summary_tag.size()), "native physical-read count");
                have_summary = true;
            }
            continue;
        }
        if (!have_layout || have_summary) fail("Native data row before layout or after query summary");
        split_columns(text, columns);
        if (columns.size() != column_count) fail("Native row column count differs from table layout");
        const auto qid = number(columns[qid_col], "query ID");
        if (qid < 1 || qid > sh.fragments) fail("Native query ID out of range: " + std::to_string(qid));
        if (!seen.mark(qid - 1)) fail("Duplicate native query ID: " + std::to_string(qid));
        if (columns[header_col].empty() || columns[header_col].find_first_of(" \t\r\n") != std::string_view::npos)
            fail("Native query header is empty or contains unsupported whitespace");
        if (header_hashes && hash_bytes(columns[header_col].data(), columns[header_col].size()) != header_hashes[qid - 1])
            fail("Native/FASTQ header mismatch at query ID " + std::to_string(qid));
        const Record before = records[qid - 1];
        if (!valid_record(before)) fail("Corrupt candidate state at query ID " + std::to_string(qid));
        Record after = before;
        std::string_view hits = columns[hits_col];
        std::array<std::uint64_t, 2> local_taxa {};
        unsigned local_count = 0;
        while (!hits.empty()) {
            if (local_count == 2) fail("Native row has more than two candidates; use -maxcand 2");
            const auto comma = hits.find(',');
            const auto candidate = trim(hits.substr(0, comma));
            const auto colon = candidate.find(':');
            if (colon == std::string_view::npos) fail("Malformed native top_hits entry");
            const auto taxid = number(candidate.substr(0, colon), "candidate taxid");
            const auto score = number(candidate.substr(colon + 1), "candidate score");
            if (!taxid || taxid > std::uint64_t(std::numeric_limits<std::int64_t>::max()) ||
                !score || score > std::numeric_limits<std::uint32_t>::max())
                fail("Candidate taxid/score outside supported range");
            if (local_count && taxid == local_taxa[0]) fail("Duplicate taxon in native top_hits");
            local_taxa[local_count++] = taxid;
            insert_candidate(after, taxid, static_cast<std::uint32_t>(score));
            ++candidate_count;
            if (comma == std::string_view::npos) break;
            hits = trim(hits.substr(comma + 1));
            if (hits.empty()) fail("Trailing comma in native top_hits");
        }
        ++dh.rows;
        if (!same_record(before, after)) {
            output_buffer.push_back(DeltaRecord{qid, after}); ++dh.changed_rows;
            if (output_buffer.size() == output_buffer.capacity()) flush_delta();
        }
        if (dh.rows % 5000000 == 0)
            std::cerr << "Compact merge: " << dh.rows << '/' << sh.fragments << " fragments; "
                      << dh.changed_rows << " changed states\n";
    }
    if (!input.eof() || input.bad()) fail("Error reading native candidate stream");
    if (!have_layout || !have_summary || dh.rows != sh.fragments || reported_reads != expected_reads)
        fail("Incomplete native query: rows=" + std::to_string(dh.rows) + ", expected=" +
             std::to_string(sh.fragments) + ", physical reads=" + std::to_string(reported_reads));
    flush_delta();
    dh.physical_reads = reported_reads;
    dh.flags = DELTA_READY | (paired ? DELTA_PAIRED : 0) | (headers ? DELTA_HEADERS : 0);
    pwrite_all(delta.fd, &dh, sizeof(dh), 0); delta.sync();
    std::ostringstream result;
    result << "{\"command\":\"update\",\"schema\":1,\"fragments\":" << sh.fragments
           << ",\"rows\":" << dh.rows << ",\"physical_reads\":" << reported_reads
           << ",\"candidate_count\":" << candidate_count << ",\"changed_rows\":" << dh.changed_rows
           << ",\"base_parts\":" << dh.base_parts << ",\"transaction\":\"" << dh.transaction
           << "\",\"state_id\":\"" << dh.state_id << "\",\"headers_checked\":" << (headers ? "true" : "false")
           << ",\"delta_bytes\":" << byte_size(dh.changed_rows, sizeof(DeltaRecord), sizeof(DeltaHeader)) << "}";
    report(args, result.str());
}

void apply(const Args& args) {
    allowed(args, {"--state", "--delta", "--summary"});
    Mapping delta(required(args, "--delta"), false);
    if (delta.bytes < sizeof(DeltaHeader)) fail("Truncated delta header");
    const auto& dh = *static_cast<const DeltaHeader*>(delta.data);
    if (std::memcmp(dh.magic, "MCDELTV1", 8) || dh.version != VERSION ||
        dh.record_bytes != sizeof(DeltaRecord) || !dh.fragments || !dh.state_id || !dh.transaction || dh.reserved ||
        !(dh.flags & DELTA_READY) || (dh.flags & ~(DELTA_READY | DELTA_PAIRED | DELTA_HEADERS)) ||
        dh.rows != dh.fragments || dh.changed_rows > dh.fragments ||
        dh.physical_reads != byte_size(dh.fragments, (dh.flags & DELTA_PAIRED) ? 2 : 1))
        fail("Incomplete or corrupt compact delta header");
    if (delta.bytes != byte_size(dh.changed_rows, sizeof(DeltaRecord), sizeof(DeltaHeader)))
        fail("Delta size does not match its changed-row count");
    const auto* changes = reinterpret_cast<const DeltaRecord*>(static_cast<const char*>(delta.data) + sizeof(DeltaHeader));
    if (hash_bytes(changes, static_cast<std::size_t>(dh.changed_rows) * sizeof(DeltaRecord)) != dh.payload_hash)
        fail("Delta payload checksum mismatch");
    Bitmap seen(dh.fragments);
    for (std::uint64_t i = 0; i < dh.changed_rows; ++i) {
        const auto& entry = changes[i];
        if (!entry.qid || entry.qid > dh.fragments || !seen.mark(entry.qid - 1) || !valid_record(entry.value))
            fail("Invalid or duplicated delta record");
    }
    Mapping state(required(args, "--state"), true, true);
    validate_state_header(state, true);
    auto& sh = *static_cast<StateHeader*>(state.data);
    if (sh.state_id != dh.state_id || sh.fragments != dh.fragments)
        fail("Delta belongs to another sample state");
    if (dh.base_parts == std::numeric_limits<std::uint64_t>::max()) fail("Invalid delta partition counter");
    if (!sh.flags && sh.applied_parts == dh.base_parts + 1 && sh.transaction == dh.transaction) {
        report(args, state_summary("apply_already_complete", sh)); return;
    }
    // A process can stop between the two scalar stores in the final header
    // update. DIRTY + base+1 + the same transaction still needs safe replay.
    const bool matching_dirty = sh.flags && sh.transaction == dh.transaction &&
        (sh.applied_parts == dh.base_parts || sh.applied_parts == dh.base_parts + 1);
    if ((sh.flags && !matching_dirty) || (!sh.flags && sh.applied_parts != dh.base_parts))
        fail("Delta is out of sequence or does not match the unfinished transaction");
    // Writing the transaction first makes interruption before DIRTY harmless:
    // it leaves a clean base state with a new transaction, accepted on replay.
    sh.transaction = dh.transaction; sh.flags = DIRTY; state.sync(sizeof(StateHeader));
    auto* records = reinterpret_cast<Record*>(static_cast<char*>(state.data) + sizeof(StateHeader));
    for (std::uint64_t i = 0; i < dh.changed_rows; ++i) records[changes[i].qid - 1] = changes[i].value;
    if (dh.changed_rows) state.sync();
    sh.applied_parts = dh.base_parts + 1; sh.flags = 0; state.sync(sizeof(StateHeader));
    report(args, state_summary("apply", sh));
}

enum Rank : unsigned { SEQUENCE, FORM, VARIETY, SUBSPECIES, SPECIES, SUBGENUS, GENUS,
    SUBTRIBE, TRIBE, SUBFAMILY, FAMILY, SUBORDER, ORDER, SUBCLASS, CLASS, SUBPHYLUM,
    PHYLUM, SUBKINGDOM, KINGDOM, DOMAIN, ROOT, NONE };
unsigned rank_of(std::string_view text) {
    static const std::unordered_map<std::string, unsigned> names {
        {"sequence",SEQUENCE},{"genome",SEQUENCE},{"form",FORM},{"forma",FORM},
        {"variety",VARIETY},{"varietas",VARIETY},{"subspecies",SUBSPECIES},{"species",SPECIES},
        {"species group",SUBGENUS},{"species subgroup",SUBGENUS},{"subgenus",SUBGENUS},
        {"genus",GENUS},{"subtribe",SUBTRIBE},{"tribe",TRIBE},{"subfamily",SUBFAMILY},
        {"family",FAMILY},{"superfamily",SUBORDER},{"parvorder",SUBORDER},
        {"infraorder",SUBORDER},{"suborder",SUBORDER},{"order",ORDER},
        {"superorder",SUBCLASS},{"infraclass",SUBCLASS},{"subclass",SUBCLASS},{"class",CLASS},
        {"superclass",SUBPHYLUM},{"subphylum",SUBPHYLUM},{"phylum",PHYLUM},{"division",PHYLUM},
        {"superphylum",SUBKINGDOM},{"subkingdom",SUBKINGDOM},{"kingdom",KINGDOM},
        {"subdomain",KINGDOM},{"superkingdom",DOMAIN},{"domain",DOMAIN},{"root",ROOT}
    };
    std::string normalized(text);
    for (auto& ch : normalized) if (ch >= 'A' && ch <= 'Z') ch = char(ch - 'A' + 'a');
    const auto it = names.find(normalized);
    return it == names.end() ? NONE : it->second;
}
struct Node { std::uint64_t parent; unsigned rank; };
class Taxonomy {
    std::unordered_map<std::uint64_t, Node> nodes;
    using Lineage = std::array<std::uint64_t, ROOT + 1>;
    std::unordered_map<std::uint64_t, Lineage> lineages;
public:
    explicit Taxonomy(const std::string& directory) {
        std::unordered_map<std::uint64_t, std::uint64_t> mergers;
        std::string line;
        std::vector<std::string_view> fields;
        std::ifstream merged(directory + "/merged.dmp");
        while (std::getline(merged, line)) {
            if (trim(line).empty()) continue;
            split_columns(line, fields);
            if (fields.size() < 2) fail("Malformed merged.dmp row");
            const auto old_id = number(fields[0], "merged old taxid");
            const auto new_id = number(fields[1], "merged new taxid");
            if (!old_id || !new_id || old_id == new_id) fail("Invalid taxid merger");
            mergers.emplace(old_id, new_id);
            nodes.emplace(old_id, Node{new_id, NONE});
        }
        if (merged.bad()) fail("Cannot read merged.dmp");
        std::ifstream source(directory + "/nodes.dmp");
        if (!source) fail("Cannot read taxonomy nodes.dmp");
        while (std::getline(source, line)) {
            if (trim(line).empty()) continue;
            split_columns(line, fields);
            if (fields.size() < 3) fail("Malformed nodes.dmp row");
            auto id = number(fields[0], "taxonomy taxid");
            auto parent = number(fields[1], "taxonomy parent");
            if (!id) fail("Zero taxonomy taxid");
            auto changed = mergers.find(id); if (changed != mergers.end()) id = changed->second;
            changed = mergers.find(parent); if (changed != mergers.end()) parent = changed->second;
            nodes.emplace(id, Node{parent, rank_of(fields[2])});
        }
        if (source.bad() || !nodes.count(1)) fail("Incomplete taxonomy or missing root");
        nodes.at(1).rank = ROOT;
    }
    const Node& node(std::uint64_t id) const {
        const auto found = nodes.find(id);
        if (found == nodes.end()) fail("Candidate taxid missing from frozen taxonomy: " + std::to_string(id));
        return found->second;
    }
    const Lineage& lineage(std::uint64_t id) {
        const auto cached = lineages.find(id);
        if (cached != lineages.end()) return cached->second;
        Lineage value {};
        std::unordered_set<std::uint64_t> visited;
        auto current = id;
        while (current) {
            if (!visited.insert(current).second || visited.size() > 4096) fail("Cycle in taxonomy lineage");
            const auto& n = node(current);
            // Match upstream: a farther ancestor overwrites a repeated rank.
            if (n.rank != NONE) value[n.rank] = current;
            if (n.parent == current) break;
            current = n.parent;
        }
        return lineages.emplace(id, value).first->second;
    }
    std::uint64_t classify(const Record& value, std::uint32_t hitmin, float fraction) {
        if (!valid_record(value)) fail("Corrupt candidate state during export");
        if (value.tax1) (void)node(value.tax1);
        if (value.tax2) (void)node(value.tax2);
        if (!value.tax1 || value.hits1 < hitmin) return 0;
        const auto rank = node(value.tax1).rank;
        if (rank < SPECIES) fail("Candidate below species rank; compact mode requires -lowest species");
        if (rank > DOMAIN) return 0;
        const float threshold = value.hits1 > hitmin ? float(value.hits1 - hitmin) * fraction : 0.0f;
        if (!value.tax2 || !(float(value.hits2) > threshold)) return value.tax1;
        if (node(value.tax2).rank < SPECIES) fail("Competing candidate below species rank");
        // unordered_map rehash does not invalidate references to its elements.
        const auto& first = lineage(value.tax1);
        const auto& second = lineage(value.tax2);
        for (unsigned r = rank; r <= DOMAIN; ++r)
            if (first[r] && first[r] == second[r]) return first[r];
        return 0; // Root/unranked-only LCA is above the permitted Domain rank.
    }
    unsigned rank(std::uint64_t id) const { return node(id).rank; }
};

void export_taxids(const Args& args) {
    allowed(args, {"--state", "--taxonomy", "--output", "--hitmin", "--hitdiff", "--summary"});
    Mapping state(required(args, "--state"), false, true);
    validate_state_header(state);
    const auto& sh = *static_cast<const StateHeader*>(state.data);
    if (!sh.applied_parts) fail("Cannot export before at least one complete partition is applied");
    const auto hitmin = number(option(args, "--hitmin", "5"), "minimum hits");
    if (!hitmin || hitmin > std::numeric_limits<std::uint16_t>::max()) fail("--hitmin must be 1..65535");
    const std::string difference = option(args, "--hitdiff", "100");
    std::size_t consumed = 0;
    float fraction = std::stof(difference, &consumed);
    if (consumed != difference.size() || !std::isfinite(fraction) || fraction < 0 || fraction > 100)
        fail("--hitdiff must be finite and between 0 and 100");
    if (fraction > 1) fraction *= 0.01f; // Same conversion as upstream options.cpp.
    Taxonomy taxonomy(required(args, "--taxonomy"));
    File output(required(args, "--output"), O_WRONLY | O_CREAT | O_EXCL);
    const auto* records = reinterpret_cast<const Record*>(static_cast<const char*>(state.data) + sizeof(StateHeader));
    std::vector<std::uint64_t> buffer;
    buffer.reserve(131072);
    std::uint64_t species = 0, higher = 0, unclassified = 0;
    for (std::uint64_t i = 0; i < sh.fragments; ++i) {
        const auto taxid = taxonomy.classify(records[i], static_cast<std::uint32_t>(hitmin), fraction);
        if (!taxid) ++unclassified;
        else if (taxonomy.rank(taxid) == SPECIES) ++species;
        else ++higher;
        buffer.push_back(taxid);
        if (buffer.size() == buffer.capacity()) {
            write_all(output.fd, buffer.data(), buffer.size() * sizeof(std::uint64_t)); buffer.clear();
        }
    }
    if (!buffer.empty()) write_all(output.fd, buffer.data(), buffer.size() * sizeof(std::uint64_t));
    output.sync();
    std::ostringstream result;
    result << "{\"command\":\"export\",\"schema\":1,\"fragments\":" << sh.fragments
           << ",\"applied_parts\":" << sh.applied_parts << ",\"assigned_species\":" << species
           << ",\"assigned_higher_rank\":" << higher << ",\"unclassified\":" << unclassified
           << ",\"output_bytes\":" << byte_size(sh.fragments, sizeof(std::uint64_t)) << "}";
    report(args, result.str());
}

void usage() {
    std::cout << "compact_merge 1 (MetaCache-compatible species top2)\n"
              << "  init   --state FILE --fragments N\n"
              << "  info   --state FILE\n"
              << "  update --state FILE --input FILE_OR_- --delta FILE --paired 0|1 [--headers FILE]\n"
              << "  apply  --state FILE --delta FILE\n"
              << "  export --state FILE --taxonomy DIR --output FILE [--hitmin 5 --hitdiff 100]\n"
              << "All commands may add --summary NEW_JSON_FILE. Existing output files are never overwritten.\n"
              << "Header hashes: N little-endian uint64 FNV1a64(R1 native header token, without @).\n"
              << "Update requires every qid1..N exactly once and the native physical-read summary.\n"
              << "Apply is idempotent only for the same transaction; retain its delta until progress commits.\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::ios::sync_with_stdio(false); std::cin.tie(nullptr);
        require_little_endian();
        if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "--version") {
            usage(); return argc < 2 ? 2 : 0;
        }
        const std::string command = argv[1];
        const auto args = arguments(argc, argv);
        if (command == "init") initialize(args);
        else if (command == "info") info(args);
        else if (command == "update") update(args);
        else if (command == "apply") apply(args);
        else if (command == "export") export_taxids(args);
        else fail("Unknown command: " + command);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR compact_merge: " << error.what() << '\n';
        return 1;
    }
}
