// AlphaGenome Atlas: bounded-memory TSV partitioning and allele-aware VCF merge.
// Build: c++ -O3 -std=c++17 -DNDEBUG alphagenome_merge.cpp -lz -o alphagenome_merge
// Ordinary gzip inputs are decompressed once per successful partition phase.
// All statistics below describe PREDICTION COVERAGE, never population genotypes.
#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;
using U64 = std::uint64_t;
constexpr std::size_t BUFFER_SIZE = 1U << 20;
constexpr std::size_t MAX_LINE_SIZE = 1U << 20;
constexpr const char* PARTITION_SCHEMA = "alphagenome-tsv-partition-1";
constexpr const char* MERGE_SCHEMA = "alphagenome-vcf-merge-1";
const std::array<std::string, 3> KINDS = {"avi", "splicing", "features"};
const std::array<std::vector<std::string>, 3> SCORE_COLUMNS = {{
    {"raw_score", "PHRED"},
    {"alphagenome_splicing"},
    {"MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",
     "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",
     "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",
     "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY",
     "IS_INSERTION", "IS_DELETION"}
}};
const std::array<std::size_t, 3> FIELD_OFFSETS = {0, 2, 3};
const std::array<U64, 25> CONTIG_LENGTHS = {
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569
};

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }

std::string json_quote(std::string_view value) {
    std::string out = "\"";
    static constexpr char HEX[] = "0123456789abcdef";
    for (unsigned char c : value) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 32) { out += "\\u00"; out += HEX[c >> 4]; out += HEX[c & 15]; }
                else out += static_cast<char>(c);
        }
    }
    out += '"';
    return out;
}

U64 parse_uint(std::string_view text, const std::string& name) {
    U64 result = 0;
    auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);
    if (text.empty() || parsed.ec != std::errc() || parsed.ptr != text.data() + text.size())
        fail("Invalid unsigned integer for " + name + ": " + std::string(text));
    return result;
}

// A strict parser for the flat scalar JSON receipts written by this executable.
// No external JSON dependency is needed on the analysis host.
std::map<std::string, std::string> read_receipt(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) fail("Cannot read receipt: " + path.string());
    std::string data((std::istreambuf_iterator<char>(input)), {});
    if (data.size() > 65536) fail("Oversized receipt: " + path.string());
    std::size_t at = 0;
    auto whitespace = [&]() { while (at < data.size() && std::isspace(static_cast<unsigned char>(data[at]))) ++at; };
    auto expect = [&](char c) {
        whitespace();
        if (at >= data.size() || data[at++] != c) fail("Malformed receipt: " + path.string());
    };
    auto string_value = [&]() {
        expect('"');
        std::string value;
        while (at < data.size()) {
            unsigned char c = data[at++];
            if (c == '"') return value;
            if (c < 32) fail("Invalid control character in receipt: " + path.string());
            if (c != '\\') { value += static_cast<char>(c); continue; }
            if (at == data.size()) fail("Incomplete JSON escape in " + path.string());
            switch (data[at++]) {
                case '"': value += '"'; break;
                case '\\': value += '\\'; break;
                case '/': value += '/'; break;
                case 'n': value += '\n'; break;
                case 'r': value += '\r'; break;
                case 't': value += '\t'; break;
                case 'b': value += '\b'; break;
                case 'f': value += '\f'; break;
                case 'u': {
                    if (at + 4 > data.size() || data[at] != '0' || data[at + 1] != '0')
                        fail("Unsupported JSON escape in " + path.string());
                    auto nibble = [&](char h) -> unsigned {
                        if (h >= '0' && h <= '9') return h - '0';
                        if (h >= 'a' && h <= 'f') return h - 'a' + 10;
                        if (h >= 'A' && h <= 'F') return h - 'A' + 10;
                        fail("Invalid JSON escape in " + path.string());
                    };
                    value += static_cast<char>((nibble(data[at + 2]) << 4) | nibble(data[at + 3]));
                    at += 4; break;
                }
                default: fail("Invalid JSON escape in " + path.string());
            }
        }
        fail("Unterminated JSON string in " + path.string());
    };
    std::map<std::string, std::string> values;
    expect('{'); whitespace();
    if (at < data.size() && data[at] == '}') { ++at; }
    else {
        for (;;) {
            const auto key = string_value(); expect(':'); whitespace();
            std::string value;
            if (at < data.size() && data[at] == '"') value = string_value();
            else {
                const auto begin = at;
                while (at < data.size() && data[at] >= '0' && data[at] <= '9') ++at;
                if (at == begin) fail("Invalid scalar in receipt: " + path.string());
                value = data.substr(begin, at - begin);
            }
            if (!values.emplace(key, value).second) fail("Duplicate receipt key: " + key);
            whitespace();
            if (at < data.size() && data[at] == '}') { ++at; break; }
            expect(',');
        }
    }
    whitespace();
    if (at != data.size()) fail("Trailing JSON data in " + path.string());
    return values;
}

void atomic_text(const fs::path& path, const std::string& content) {
    if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
    fs::path temporary = path; temporary += ".partial";
    {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        if (!out) fail("Cannot write " + temporary.string());
        out.write(content.data(), static_cast<std::streamsize>(content.size()));
        out.flush();
        if (!out) fail("Failed writing " + temporary.string());
        out.close();
        if (!out) fail("Failed closing " + temporary.string());
    }
    fs::rename(temporary, path);
}

std::pair<std::string, int> chromosome(std::string_view text) {
    if (text.substr(0, 3) == "chr") text.remove_prefix(3);
    if (text == "X") return {"chrX", 22};
    if (text == "Y") return {"chrY", 23};
    if (text == "M" || text == "MT") return {"chrM", 24};
    if (text.empty() || text.size() > 2 || text[0] == '0')
        fail("Unsupported chromosome: " + std::string(text));
    const auto n = parse_uint(text, "chromosome");
    if (n < 1 || n > 22) fail("Unsupported chromosome: " + std::string(text));
    return {"chr" + std::to_string(n), static_cast<int>(n - 1)};
}

int base_index(char base) {
    switch (base) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; }
    return -1;
}

class GzipReader {
    gzFile file_ = nullptr;
    std::vector<char> buffer_ = std::vector<char>(BUFFER_SIZE);
    std::size_t start_ = 0, end_ = 0;
    bool eof_ = false;
    std::string path_;
public:
    explicit GzipReader(const fs::path& path) : path_(path.string()) {
        std::ifstream magic(path, std::ios::binary);
        unsigned char prefix[2] = {0, 0};
        magic.read(reinterpret_cast<char*>(prefix), 2);
        if (magic.gcount() != 2 || prefix[0] != 0x1f || prefix[1] != 0x8b)
            fail("Expected gzip input: " + path_);
        file_ = gzopen(path_.c_str(), "rb");
        if (!file_) fail("Cannot open gzip input: " + path_);
        if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0) {
            gzclose(file_); file_ = nullptr; fail("Cannot allocate gzip buffer for " + path_);
        }
    }
    GzipReader(const GzipReader&) = delete;
    GzipReader& operator=(const GzipReader&) = delete;
    ~GzipReader() { if (file_) gzclose(file_); }
    bool next(std::string& line) {
        line.clear();
        for (;;) {
            if (start_ < end_) {
                const auto* found = static_cast<const char*>(std::memchr(buffer_.data() + start_, '\n', end_ - start_));
                const auto stop = found ? static_cast<std::size_t>(found - buffer_.data()) : end_;
                line.append(buffer_.data() + start_, stop - start_);
                start_ = stop + (found ? 1 : 0);
                if (line.size() > MAX_LINE_SIZE) fail("Oversized TSV line in " + path_);
                if (found) { if (!line.empty() && line.back() == '\r') line.pop_back(); return true; }
            }
            if (eof_) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return !line.empty();
            }
            const int count = gzread(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));
            int code = Z_OK;
            const char* error = gzerror(file_, &code);
            if (count < 0 || (code != Z_OK && code != Z_STREAM_END))
                fail("gzip read/CRC/truncation error in " + path_ + ": " + (error ? error : "unknown"));
            start_ = 0;
            end_ = count > 0 ? static_cast<std::size_t>(count) : 0;
            if (count == 0) {
                if (!gzeof(file_)) fail("Unexpected gzip EOF in " + path_);
                const int close_code = gzclose(file_); file_ = nullptr;
                if (close_code != Z_OK) fail("gzip close/truncation error in " + path_);
                eof_ = true;
            }
        }
    }
};

class CanonicalSink {
    gzFile file_ = nullptr;
    std::string buffer_;
    std::string path_;
    bool closed_ = false;
public:
    U64 uncompressed_bytes = 0;
    uLong crc = crc32(0L, Z_NULL, 0);
    explicit CanonicalSink(const std::optional<fs::path>& path) {
        buffer_.reserve(BUFFER_SIZE + 4096);
        if (path) {
            path_ = path->string(); file_ = gzopen(path_.c_str(), "wb1");
            if (!file_) fail("Cannot write staging gzip: " + path_);
            if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0) {
                gzclose(file_); file_ = nullptr; fail("Cannot allocate gzip buffer for " + path_);
            }
        }
    }
    ~CanonicalSink() { if (file_) gzclose(file_); }
    void flush() {
        if (buffer_.empty()) return;
        crc = crc32(crc, reinterpret_cast<const Bytef*>(buffer_.data()), static_cast<uInt>(buffer_.size()));
        uncompressed_bytes += buffer_.size();
        if (file_) {
            const auto written = gzwrite(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));
            if (written != static_cast<int>(buffer_.size())) fail("Failed writing staging gzip: " + path_);
        }
        buffer_.clear();
    }
    void append(std::string_view value) { buffer_.append(value.data(), value.size()); if (buffer_.size() >= BUFFER_SIZE) flush(); }
    void finish() {
        if (closed_) return;
        flush();
        if (file_) {
            const int result = gzclose(file_); file_ = nullptr;
            if (result != Z_OK) fail("Failed closing staging gzip: " + path_);
        }
        closed_ = true;
    }
};

void split_fields(std::string_view line, std::vector<std::string_view>& fields) {
    fields.clear();
    std::size_t at = 0;
    while (at < line.size()) {
        while (at < line.size() && (line[at] == ' ' || line[at] == '\t')) ++at;
        if (at == line.size()) break;
        const auto begin = at;
        while (at < line.size() && line[at] != ' ' && line[at] != '\t') ++at;
        fields.emplace_back(line.data() + begin, at - begin);
    }
}

std::string header_for(int kind) {
    std::string line = "#CHROM\tPOS\tREF\tALT";
    for (const auto& column : SCORE_COLUMNS[kind]) { line += '\t'; line += column; }
    line += '\n'; return line;
}

void validate_header(std::string_view line, int kind, const std::string& context) {
    std::vector<std::string_view> fields;
    split_fields(line, fields);
    if (fields.size() != SCORE_COLUMNS[kind].size() + 4 || fields[0] != "#CHROM" ||
        fields[1] != "POS" || fields[2] != "REF" || fields[3] != "ALT")
        fail("Unexpected " + KINDS[kind] + " header in " + context + "; expected " + header_for(kind));
    for (std::size_t j = 0; j < SCORE_COLUMNS[kind].size(); ++j)
        if (fields[j + 4] != SCORE_COLUMNS[kind][j])
            fail("Unexpected " + KINDS[kind] + " score column in " + context + ": " + std::string(fields[j + 4]));
}

struct Numeric { bool missing = false; bool zero = false; };

Numeric validate_number(std::string_view token, const std::string& context) {
    if (token == "." || token == "NA" || token == "NaN" || token == "nan") return {true, false};
    std::size_t at = 0, digits = 0;
    bool nonzero = false;
    if (at < token.size() && (token[at] == '+' || token[at] == '-')) ++at;
    while (at < token.size() && token[at] >= '0' && token[at] <= '9') {
        nonzero |= token[at] != '0'; ++at; ++digits;
    }
    if (at < token.size() && token[at] == '.') {
        ++at;
        while (at < token.size() && token[at] >= '0' && token[at] <= '9') {
            nonzero |= token[at] != '0'; ++at; ++digits;
        }
    }
    if (digits == 0) fail("Invalid numeric token " + json_quote(token) + " in " + context);
    if (at < token.size() && (token[at] == 'e' || token[at] == 'E')) {
        ++at;
        if (at < token.size() && (token[at] == '+' || token[at] == '-')) ++at;
        const auto begin = at;
        while (at < token.size() && token[at] >= '0' && token[at] <= '9') ++at;
        if (at == begin) fail("Invalid exponent in " + context + ": " + std::string(token));
    }
    if (at != token.size()) fail("Invalid numeric token " + json_quote(token) + " in " + context);
    // GCC 8's libstdc++ implements integer from_chars only. Keep the fast
    // floating parser where available; the portable branch preserves decimal
    // text and underflow semantics and never round-trips output through double.
#if defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE >= 11 && !defined(ALPHAGENOME_PORTABLE_FLOAT)
    auto number = token;
    if (!number.empty() && number.front() == '+') number.remove_prefix(1);
    double value = 0;
    const auto result = std::from_chars(number.data(), number.data() + number.size(), value, std::chars_format::general);
    if (result.ec == std::errc::result_out_of_range) {
        // Very small finite decimal text may underflow binary64. Preserve its text
        // and do not count it as numeric zero; only overflow is an error.
        std::string copy(token);
        char* end = nullptr;
        value = std::strtod(copy.c_str(), &end);
        if (end != copy.c_str() + copy.size() || !std::isfinite(value))
            fail("Non-finite/overflow numeric token " + json_quote(token) + " in " + context);
    } else if (result.ec != std::errc() || result.ptr != number.data() + number.size() || !std::isfinite(value))
        fail("Non-finite numeric token " + json_quote(token) + " in " + context);
#else
    std::string copy(token);
    char* end = nullptr;
    const double value = std::strtod(copy.c_str(), &end);
    if (end != copy.c_str() + copy.size() || !std::isfinite(value))
        fail("Non-finite/overflow numeric token " + json_quote(token) + " in " + context);
#endif
    return {false, !nonzero};
}

struct ParsedRow {
    std::string chrom;
    int chrom_index = 0;
    U64 pos = 0;
    char ref = 'N', alt = 'N';
    int alt_index = -1;
};

ParsedRow parse_row(const std::vector<std::string_view>& fields, int kind, const std::string& context,
                    std::vector<Numeric>* numbers = nullptr) {
    if (fields.size() != SCORE_COLUMNS[kind].size() + 4)
        fail("Wrong field count in " + context + ": expected " + std::to_string(SCORE_COLUMNS[kind].size() + 4) +
             ", got " + std::to_string(fields.size()));
    ParsedRow row;
    auto chr = chromosome(fields[0]); row.chrom = std::move(chr.first); row.chrom_index = chr.second;
    row.pos = parse_uint(fields[1], "POS in " + context);
    if (!row.pos || row.pos > CONTIG_LENGTHS[row.chrom_index])
        fail("Position outside GRCh38 primary chromosome in " + context);
    if (fields[2].size() != 1 || fields[3].size() != 1 || base_index(fields[2][0]) < 0 ||
        base_index(fields[3][0]) < 0 || fields[2][0] == fields[3][0])
        fail("Expected canonical A/C/G/T single-base REF != ALT in " + context);
    row.ref = fields[2][0]; row.alt = fields[3][0]; row.alt_index = base_index(row.alt);
    if (numbers) { numbers->clear(); numbers->reserve(SCORE_COLUMNS[kind].size()); }
    for (std::size_t j = 4; j < fields.size(); ++j) {
        auto numeric = validate_number(fields[j], context);
        if (numbers) numbers->push_back(numeric);
    }
    return row;
}

struct PositionOrder {
    U64 pos = 0;
    char ref = 'N';
    unsigned alt_mask = 0;
    void accept(const ParsedRow& row, const std::string& context) {
        if (row.pos < pos) fail("Decreasing position in " + context);
        if (row.pos != pos) { pos = row.pos; ref = row.ref; alt_mask = 0; }
        else if (row.ref != ref) fail("Inconsistent REF at one position in " + context);
        const auto bit = 1U << row.alt_index;
        if (alt_mask & bit) fail("Duplicate CHROM/POS/REF/ALT in " + context);
        alt_mask |= bit;
    }
};

struct ChromReceipt {
    std::string chrom, kind, token;
    U64 rows = 0, first_pos = 0, last_pos = 0, gzip_bytes = 0, uncompressed_bytes = 0, crc = 0;
    std::string json() const {
        std::ostringstream out;
        out << "{\n  \"schema\": " << json_quote(PARTITION_SCHEMA)
            << ",\n  \"chrom\": " << json_quote(chrom) << ",\n  \"kind\": " << json_quote(kind)
            << ",\n  \"token\": " << json_quote(token) << ",\n  \"rows\": " << rows
            << ",\n  \"first_pos\": " << first_pos << ",\n  \"last_pos\": " << last_pos
            << ",\n  \"gzip_bytes\": " << gzip_bytes << ",\n  \"uncompressed_bytes\": " << uncompressed_bytes
            << ",\n  \"uncompressed_crc32\": " << crc << "\n}\n";
        return out.str();
    }
    static ChromReceipt load(const fs::path& path) {
        auto values = read_receipt(path);
        auto get = [&](const std::string& key) -> std::string {
            auto it = values.find(key);
            if (it == values.end()) fail("Missing receipt field " + key + " in " + path.string());
            return it->second;
        };
        if (get("schema") != PARTITION_SCHEMA) fail("Unsupported receipt schema in " + path.string());
        ChromReceipt receipt;
        receipt.chrom = get("chrom"); receipt.kind = get("kind"); receipt.token = get("token");
        receipt.rows = parse_uint(get("rows"), "receipt rows");
        receipt.first_pos = parse_uint(get("first_pos"), "receipt first_pos");
        receipt.last_pos = parse_uint(get("last_pos"), "receipt last_pos");
        receipt.gzip_bytes = parse_uint(get("gzip_bytes"), "receipt gzip_bytes");
        receipt.uncompressed_bytes = parse_uint(get("uncompressed_bytes"), "receipt uncompressed_bytes");
        receipt.crc = parse_uint(get("uncompressed_crc32"), "receipt CRC");
        return receipt;
    }
};

bool staging_quick_check(const fs::path& path, const ChromReceipt& receipt) {
    if (!fs::is_regular_file(path) || fs::file_size(path) != receipt.gzip_bytes || receipt.gzip_bytes < 18) return false;
    std::ifstream input(path, std::ios::binary);
    unsigned char bytes[8] = {};
    input.seekg(-8, std::ios::end); input.read(reinterpret_cast<char*>(bytes), 8);
    if (!input) return false;
    auto little = [&](int offset) -> U64 {
        return U64(bytes[offset]) | (U64(bytes[offset + 1]) << 8) | (U64(bytes[offset + 2]) << 16) | (U64(bytes[offset + 3]) << 24);
    };
    return little(0) == receipt.crc && little(4) == (receipt.uncompressed_bytes & 0xffffffffU);
}

struct Arguments {
    std::string mode;
    std::map<std::string, std::string> values;
    explicit Arguments(int argc, char** argv) {
        if (argc < 2) fail("Specify partition or merge; use --help for usage");
        mode = argv[1];
        if (mode == "--help" || mode == "-h") return;
        if (mode != "partition" && mode != "merge") fail("Unknown mode: " + mode);
        const std::set<std::string> allowed = mode == "partition"
            ? std::set<std::string>{"--kind", "--input", "--out-dir", "--token", "--progress-every", "--file-prefix"}
            : std::set<std::string>{"--chrom", "--avi", "--splicing", "--features", "--stats", "--token", "--progress-every"};
        for (int i = 2; i < argc; i += 2) {
            if (i + 1 == argc) fail("Missing value for " + std::string(argv[i]));
            if (!allowed.count(argv[i])) fail("Unknown option: " + std::string(argv[i]));
            if (!values.emplace(argv[i], argv[i + 1]).second) fail("Duplicate option: " + std::string(argv[i]));
        }
    }
    std::string required(const std::string& name) const {
        auto found = values.find(name);
        if (found == values.end() || found->second.empty()) fail("Required option: " + name);
        return found->second;
    }
    U64 progress() const {
        const auto found = values.find("--progress-every");
        return found == values.end() ? 10000000 : parse_uint(found->second, "--progress-every");
    }
};

int partition(const Arguments& args) {
    const auto kind_name = args.required("--kind");
    auto found = std::find(KINDS.begin(), KINDS.end(), kind_name);
    if (found == KINDS.end()) fail("--kind must be avi, splicing, or features");
    const int kind = static_cast<int>(found - KINDS.begin());
    const fs::path input = args.required("--input"), out_dir = args.required("--out-dir");
    const auto prefix = args.values.count("--file-prefix") ? args.values.at("--file-prefix") : std::string{};
    if (!prefix.empty() && prefix != kind_name + ".") fail("--file-prefix must match kind plus '.'");
    const auto token = args.required("--token");
    const U64 progress = args.progress();
    fs::create_directories(out_dir);
    GzipReader reader(input);
    std::string line;
    if (!reader.next(line)) fail("Empty input without header: " + input.string());
    validate_header(line, kind, input.string());
    const auto header = header_for(kind);
    std::vector<std::string_view> fields; fields.reserve(22);
    std::string normalized; normalized.reserve(2048);
    std::set<std::string> visited;
    std::vector<ChromReceipt> completed;
    std::unique_ptr<CanonicalSink> sink;
    std::optional<ChromReceipt> previous;
    ChromReceipt receipt;
    PositionOrder order;
    fs::path destination, temporary, done_path;
    U64 source_rows = 0, line_number = 1, reused_chromosomes = 0;
    auto finish_chromosome = [&]() {
        if (!sink) return;
        sink->finish();
        receipt.uncompressed_bytes = sink->uncompressed_bytes;
        receipt.crc = sink->crc;
        if (previous) {
            if (receipt.rows != previous->rows || receipt.first_pos != previous->first_pos ||
                receipt.last_pos != previous->last_pos || receipt.uncompressed_bytes != previous->uncompressed_bytes ||
                receipt.crc != previous->crc)
                fail("Source content differs from completed chromosome despite matching token: " + receipt.chrom);
            receipt.gzip_bytes = previous->gzip_bytes; ++reused_chromosomes;
            std::cerr << "partition " << kind_name << ' ' << receipt.chrom << " reused; rows=" << receipt.rows << '\n';
        } else {
            receipt.gzip_bytes = fs::file_size(temporary);
            fs::rename(temporary, destination);
            atomic_text(done_path, receipt.json());
            std::cerr << "partition " << kind_name << ' ' << receipt.chrom << " completed; rows=" << receipt.rows
                      << " gzip_bytes=" << receipt.gzip_bytes << '\n';
        }
        completed.push_back(receipt); sink.reset(); previous.reset();
    };
    while (reader.next(line)) {
        ++line_number;
        const std::string context = input.string() + ":" + std::to_string(line_number);
        split_fields(line, fields);
        const auto row = parse_row(fields, kind, context);
        if (!sink || row.chrom != receipt.chrom) {
            finish_chromosome();
            if (!visited.insert(row.chrom).second) fail("Chromosome block revisited in " + context + ": " + row.chrom);
            receipt = ChromReceipt{}; receipt.chrom = row.chrom; receipt.kind = kind_name; receipt.token = token;
            order = PositionOrder{};
            destination = out_dir / (prefix + row.chrom + ".tsv.gz"); temporary = destination; temporary += ".partial";
            done_path = out_dir / (prefix + row.chrom + ".done.json");
            if (fs::exists(done_path)) {
                const auto existing = ChromReceipt::load(done_path);
                if (existing.token != token || existing.kind != kind_name || existing.chrom != row.chrom)
                    fail("Incompatible staging receipt; use a distinct/clean staging directory: " + done_path.string());
                if (staging_quick_check(destination, existing)) previous = existing;
                else std::cerr << "partition " << kind_name << ' ' << row.chrom << " regenerating missing/damaged staging output\n";
            }
            sink = std::make_unique<CanonicalSink>(previous ? std::nullopt : std::optional<fs::path>(temporary));
            sink->append(header);
        }
        order.accept(row, context);
        if (receipt.rows++ == 0) receipt.first_pos = row.pos;
        receipt.last_pos = row.pos;
        normalized.clear(); normalized += row.chrom;
        for (std::size_t j = 1; j < fields.size(); ++j) { normalized += '\t'; normalized.append(fields[j]); }
        normalized += '\n'; sink->append(normalized);
        ++source_rows;
        if (progress && source_rows % progress == 0)
            std::cerr << "partition " << kind_name << " rows=" << source_rows << " at=" << row.chrom << ':' << row.pos << '\n';
    }
    finish_chromosome();
    // Old source.done is never sufficient to reuse outputs: the orchestrator also
    // checks source fingerprints and every chromosome receipt. An interrupted
    // ordinary-gzip source can require rereading its prefix on the next attempt.
    std::ostringstream summary;
    summary << "{\n  \"schema\": " << json_quote(PARTITION_SCHEMA) << ",\n  \"kind\": " << json_quote(kind_name)
            << ",\n  \"token\": " << json_quote(token) << ",\n  \"source_rows\": " << source_rows
            << ",\n  \"reused_chromosomes\": " << reused_chromosomes << ",\n  \"chromosomes\": [";
    for (std::size_t j = 0; j < completed.size(); ++j) {
        if (j) summary << ',';
        summary << '\n' << completed[j].json();
    }
    summary << "\n  ]\n}\n";
    atomic_text(out_dir / (prefix + "source.done.json"), summary.str());
    std::cerr << "partition " << kind_name << " finished; source_rows=" << source_rows << " chromosomes=" << completed.size() << '\n';
    return 0;
}

struct ScoreRow {
    ParsedRow locus;
    std::vector<std::string> values;
    std::vector<Numeric> numbers;
};

struct Group {
    U64 pos = 0;
    char ref = 'N';
    std::array<std::optional<ScoreRow>, 4> alleles;
};

class SourceStream {
    int kind_;
    std::string expected_chrom_, path_, line_;
    GzipReader reader_;
    PositionOrder order_;
    std::vector<std::string_view> fields_;
    U64 line_number_ = 1;
public:
    U64 rows = 0;
    std::optional<ScoreRow> current;
    SourceStream(int kind, const std::string& path, const std::string& chrom)
        : kind_(kind), expected_chrom_(chrom), path_(path), reader_(path) {
        fields_.reserve(22);
        if (!reader_.next(line_)) fail("Empty input without header: " + path_);
        validate_header(line_, kind_, path_); advance();
    }
    void advance() {
        if (!reader_.next(line_)) { current.reset(); return; }
        ++line_number_;
        const std::string context = path_ + ":" + std::to_string(line_number_);
        split_fields(line_, fields_);
        ScoreRow row;
        row.locus = parse_row(fields_, kind_, context, &row.numbers);
        if (row.locus.chrom != expected_chrom_) fail("Unexpected chromosome in " + context + ": " + row.locus.chrom);
        order_.accept(row.locus, context);
        row.values.reserve(SCORE_COLUMNS[kind_].size());
        for (std::size_t j = 4; j < fields_.size(); ++j) row.values.emplace_back(fields_[j]);
        current = std::move(row); ++rows;
    }
    Group take_group() {
        if (!current) fail("Internal error: group requested after EOF");
        Group group; group.pos = current->locus.pos; group.ref = current->locus.ref;
        do {
            const int alt = current->locus.alt_index;
            group.alleles[alt] = std::move(*current);
            advance();
        } while (current && current->locus.pos == group.pos);
        return group;
    }
};

struct FieldStats {
    std::string name;
    int source = 0;
    U64 numeric_present = 0, numeric_missing = 0, numeric_zero = 0, absent_source_alleles = 0;
    U64 sites_with_three_alt_scores = 0;
};

struct SourceStats {
    U64 covered_alleles = 0, covered_sites = 0, sites_with_three_alt_rows = 0;
    U64 sites_with_three_alt_complete_scores = 0;
    U64 numeric_present = 0, numeric_missing = 0, numeric_zero = 0;
};

std::array<FieldStats, 21> make_field_stats() {
    std::array<FieldStats, 21> fields;
    fields[0].name = "AG_AVI_RAW"; fields[0].source = 0;
    fields[1].name = "AG_AVI_PHRED"; fields[1].source = 0;
    fields[2].name = "AG_SPLICE"; fields[2].source = 1;
    for (std::size_t j = 0; j < SCORE_COLUMNS[2].size(); ++j) {
        fields[j + 3].name = "AG_SHAP_" + SCORE_COLUMNS[2][j]; fields[j + 3].source = 2;
    }
    return fields;
}

void vcf_header(const std::string& chrom, const std::string& token, const std::array<FieldStats, 21>& fields) {
    const auto chr = chromosome(chrom);
    std::cout << "##fileformat=VCFv4.3\n##source=allele3_alphagenome_merge_1\n##reference=GRCh38\n"
              << "##AG_schema=" << MERGE_SCHEMA << "\n##AG_recipe_token=" << json_quote(token) << '\n'
              << "##AG_allele_semantics=\"ALT contains published predicted substitutions; presence is not evidence of population observation\"\n"
              << "##AG_missing_semantics=\"Dot means unavailable; AG_SOURCE_MASK distinguishes absent rows from explicit missing scores\"\n"
              << "##AG_feature_semantics=\"AG_SHAP fields are published signed feature contributions, not raw molecular-effect measurements\"\n"
              << "##contig=<ID=" << chrom << ",length=" << CONTIG_LENGTHS[chr.second] << ",assembly=GRCh38>\n";
    for (std::size_t j = 0; j < fields.size(); ++j) {
        std::string description;
        if (j == 0) description = "Published AlphaGenome variant impact raw score; not a genotype statistic";
        else if (j == 1) description = "Published AlphaGenome variant impact PHRED score; not genotype quality";
        else if (j == 2) description = "Published combined AlphaGenome splicing score; distinct from SHAP MERGED_SPLICING";
        else description = "Published signed AVI feature-importance (SHAP) value for " + SCORE_COLUMNS[2][j - 3] + "; not the raw feature score";
        std::cout << "##INFO=<ID=" << fields[j].name << ",Number=A,Type=Float,Description=" << json_quote(description) << ">\n";
    }
    std::cout << "##INFO=<ID=AG_SOURCE_MASK,Number=A,Type=Integer,Description=\"Published row presence: AVI=1, splicing=2, feature importance=4; bitwise union\">\n"
              << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}

int merge(const Arguments& args) {
    const auto chrom = chromosome(args.required("--chrom")).first;
    const auto token = args.required("--token");
    const fs::path stats_path = args.required("--stats");
    const U64 progress = args.progress();
    std::array<std::unique_ptr<SourceStream>, 3> inputs;
    bool supplied = false;
    for (int kind = 0; kind < 3; ++kind) {
        const auto path = args.required("--" + KINDS[kind]);
        if (path != "NONE") { inputs[kind] = std::make_unique<SourceStream>(kind, path, chrom); supplied = true; }
    }
    if (!supplied) fail("At least one source file must be supplied; all three were NONE");
    auto fields = make_field_stats();
    std::array<SourceStats, 3> source_stats;
    std::array<U64, 4> sites_by_nalt{};
    std::array<U64, 8> source_mask_alleles{};
    U64 sites = 0, union_alt_rows = 0, first_pos = 0, last_pos = 0, next_progress = progress;
    vcf_header(chrom, token, fields);
    std::string output; output.reserve(4096);
    constexpr char BASES[] = "ACGT";
    for (;;) {
        U64 pos = std::numeric_limits<U64>::max();
        for (const auto& source : inputs) if (source && source->current) pos = std::min(pos, source->current->locus.pos);
        if (pos == std::numeric_limits<U64>::max()) break;
        std::array<std::optional<Group>, 3> groups;
        char ref = 'N';
        std::array<unsigned, 4> masks{};
        for (int kind = 0; kind < 3; ++kind) {
            if (!inputs[kind] || !inputs[kind]->current || inputs[kind]->current->locus.pos != pos) continue;
            groups[kind] = inputs[kind]->take_group();
            if (ref == 'N') ref = groups[kind]->ref;
            else if (ref != groups[kind]->ref)
                fail("Cross-source REF mismatch at " + chrom + ":" + std::to_string(pos));
            U64 allele_count = 0;
            bool complete = true;
            for (int alt = 0; alt < 4; ++alt) {
                if (!groups[kind]->alleles[alt]) continue;
                ++allele_count; masks[alt] |= 1U << kind;
                const auto& row = *groups[kind]->alleles[alt];
                for (std::size_t j = 0; j < row.numbers.size(); ++j) {
                    auto& field = fields[FIELD_OFFSETS[kind] + j];
                    if (row.numbers[j].missing) { ++field.numeric_missing; ++source_stats[kind].numeric_missing; complete = false; }
                    else {
                        ++field.numeric_present; ++source_stats[kind].numeric_present;
                        if (row.numbers[j].zero) { ++field.numeric_zero; ++source_stats[kind].numeric_zero; }
                    }
                }
            }
            source_stats[kind].covered_alleles += allele_count;
            ++source_stats[kind].covered_sites;
            if (allele_count == 3) {
                ++source_stats[kind].sites_with_three_alt_rows;
                if (complete) ++source_stats[kind].sites_with_three_alt_complete_scores;
                for (std::size_t j = 0; j < SCORE_COLUMNS[kind].size(); ++j) {
                    bool field_complete = true;
                    for (int alt = 0; alt < 4; ++alt)
                        if (groups[kind]->alleles[alt] && groups[kind]->alleles[alt]->numbers[j].missing) field_complete = false;
                    if (field_complete) ++fields[FIELD_OFFSETS[kind] + j].sites_with_three_alt_scores;
                }
            }
        }
        unsigned nalt = 0;
        for (unsigned mask : masks) if (mask) { ++nalt; ++source_mask_alleles[mask]; }
        if (nalt < 1 || nalt > 3) fail("Invalid ALT union at " + chrom + ":" + std::to_string(pos));
        if (++sites == 1) first_pos = pos;
        last_pos = pos;
        ++sites_by_nalt[nalt]; union_alt_rows += nalt;
        output.clear(); output += chrom; output += '\t'; output += std::to_string(pos);
        output += "\t.\t"; output += ref; output += '\t';
        bool first = true;
        for (int alt = 0; alt < 4; ++alt) if (masks[alt]) {
            if (!first) output += ',';
            output += BASES[alt]; first = false;
        }
        output += "\t.\t.\t";
        for (std::size_t field_index = 0; field_index < fields.size(); ++field_index) {
            auto& field = fields[field_index];
            if (field_index) output += ';';
            output += field.name; output += '='; first = true;
            const auto offset = field_index - FIELD_OFFSETS[field.source];
            for (int alt = 0; alt < 4; ++alt) if (masks[alt]) {
                if (!first) output += ',';
                first = false;
                if (!groups[field.source] || !groups[field.source]->alleles[alt]) {
                    output += '.'; ++field.absent_source_alleles;
                } else {
                    const auto& row = *groups[field.source]->alleles[alt];
                    if (row.numbers[offset].missing) output += '.';
                    else output += row.values[offset];
                }
            }
        }
        output += ";AG_SOURCE_MASK="; first = true;
        for (unsigned mask : masks) if (mask) {
            if (!first) output += ',';
            output += static_cast<char>('0' + mask); first = false;
        }
        output += '\n';
        std::cout.write(output.data(), static_cast<std::streamsize>(output.size()));
        if (!std::cout) fail("VCF stdout write failed (downstream compressor may have failed)");
        if (progress && union_alt_rows >= next_progress) {
            std::cerr << "merge " << chrom << " sites=" << sites << " union_alt_rows=" << union_alt_rows << " pos=" << pos << '\n';
            next_progress = union_alt_rows + progress;
        }
    }
    std::cout.flush();
    if (!std::cout) fail("VCF stdout flush failed");
    std::ostringstream stats;
    stats << "{\n  \"schema\": " << json_quote(MERGE_SCHEMA) << ",\n  \"chrom\": " << json_quote(chrom)
          << ",\n  \"token\": " << json_quote(token) << ",\n  \"assembly\": \"GRCh38\",\n  \"source_rows\": {";
    for (int kind = 0; kind < 3; ++kind) { if (kind) stats << ','; stats << json_quote(KINDS[kind]) << ':' << (inputs[kind] ? inputs[kind]->rows : 0); }
    stats << "},\n  \"union_alt_rows\": " << union_alt_rows << ",\n  \"sites\": " << sites
          << ",\n  \"first_pos\": " << first_pos << ",\n  \"last_pos\": " << last_pos
          << ",\n  \"sites_by_nalt\": {\"1\":" << sites_by_nalt[1] << ",\"2\":" << sites_by_nalt[2] << ",\"3\":" << sites_by_nalt[3] << '}'
          << ",\n  \"prediction_base_coverage\": {\"exactly_2_bases\":" << sites_by_nalt[1] << ",\"exactly_3_bases\":" << sites_by_nalt[2]
          << ",\"exactly_4_bases\":" << sites_by_nalt[3] << ",\"at_least_3_bases\":" << sites_by_nalt[2] + sites_by_nalt[3] << '}'
          << ",\n  \"source_mask_alleles\": {";
    for (unsigned mask = 1; mask <= 7; ++mask) { if (mask > 1) stats << ','; stats << json_quote(std::to_string(mask)) << ':' << source_mask_alleles[mask]; }
    stats << "},\n  \"per_source\": {";
    for (int kind = 0; kind < 3; ++kind) {
        if (kind) stats << ',';
        const auto& s = source_stats[kind];
        stats << '\n' << json_quote(KINDS[kind]) << ":{\"covered_alleles\":" << s.covered_alleles << ",\"covered_sites\":" << s.covered_sites
              << ",\"sites_with_three_alt_rows\":" << s.sites_with_three_alt_rows
              << ",\"sites_with_three_alt_complete_scores\":" << s.sites_with_three_alt_complete_scores
              << ",\"numeric_present\":" << s.numeric_present << ",\"numeric_missing\":" << s.numeric_missing << ",\"numeric_zero\":" << s.numeric_zero << '}';
    }
    stats << "\n  },\n  \"per_field\": {";
    for (std::size_t j = 0; j < fields.size(); ++j) {
        if (j) stats << ',';
        const auto& f = fields[j];
        stats << '\n' << json_quote(f.name) << ":{\"source\":" << json_quote(KINDS[f.source])
              << ",\"numeric_present\":" << f.numeric_present << ",\"numeric_missing\":" << f.numeric_missing
              << ",\"numeric_zero\":" << f.numeric_zero << ",\"absent_source_alleles\":" << f.absent_source_alleles
              << ",\"sites_with_three_alt_scores\":" << f.sites_with_three_alt_scores << '}';
    }
    stats << "\n  },\n  \"observed_tri_allelic_snvs\": null,\n  \"observed_quad_allelic_snvs\": null,\n"
          << "  \"observed_allele_statistics_reason\": \"Unavailable: Atlas predicts possible substitutions; population AC/AN/AF or genotypes are required to identify observed tri-/quad-allelic SNVs.\",\n"
          << "  \"prediction_coverage_definition\": \"Number of unique published ALT bases plus the shared reference base; no unlisted ALT is invented.\"\n}\n";
    atomic_text(stats_path, stats.str());
    std::cerr << "merge " << chrom << " finished; sites=" << sites << " union_alt_rows=" << union_alt_rows
              << " predicted_ref_plus_2alt=" << sites_by_nalt[2] << " predicted_ref_plus_3alt=" << sites_by_nalt[3] << '\n';
    return 0;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    try {
        const Arguments args(argc, argv);
        if (args.mode == "--help" || args.mode == "-h") {
            std::cout << "AlphaGenome Atlas prediction coverage converter (C++17 + zlib)\n"
                      << "partition --kind avi|splicing|features --input FILE.gz --out-dir DIR --token TOKEN [--progress-every N]\n"
                      << "merge --chrom CHR --avi FILE.gz|NONE --splicing FILE.gz|NONE --features FILE.gz|NONE --stats FILE.json --token TOKEN [--progress-every N]\n"
                      << "Merge writes VCFv4.3 to stdout; progress/errors go to stderr. Inputs must contain sorted positions in contiguous chromosome blocks.\n"
                      << "Missing tokens: '.', 'NA', 'NaN', 'nan'. Zero and signed finite score text are retained.\n"
                      << "All allele summaries describe predicted substitutions, not observed population alleles.\n";
            return 0;
        }
        return args.mode == "partition" ? partition(args) : merge(args);
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 2;
    }
}
