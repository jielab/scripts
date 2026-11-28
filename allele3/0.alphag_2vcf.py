#!/usr/bin/env python3
"""Merge three AlphaGenome SNV releases into indexed per-chromosome VCFs.

Entry point using the shared scripts/0f submission engine.
Scientific sources and shared runtime snapshots are checked on reuse.

The C++ worker reads ordinary gzip sources once into chromosome gzip parts,
then performs a bounded-memory outer join. Python manages jobs, compression,
indexes, provenance, restart and allele summaries. Inputs are never deleted.
The VCF is a prediction catalogue, not observed population variant calls.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "0f"))
from workflow_paths import configured_path
from workflow_mode import resolve_mode


def entry_defaults(a):
    """Standalone replacement for the resource shell entry; no data is touched."""
    if a.submit == 'auto':
        a.submit = resolve_mode()
    if a.jobs is None:
        a.jobs = int(os.environ.get('jobs_N', '8' if a.submit == 'hpc' else '1'))
    atlas = Path(os.environ.get('ALLELE_ATLAS_DIR') or
                 configured_path('DIR_ANNOT', 'alphagenome/38', default_mode=a.submit)).expanduser().resolve()
    a.input_dir = a.input_dir or str(atlas / 'raw')
    a.out_dir = a.out_dir or str(atlas)
    if a.worker and not a.python:
        a.python = sys.executable
    if not a.python:
        a.python = os.environ.get('ALLELE_PYTHON')
        if not a.python and a.submit == 'hpc':
            a.python = str(Path(configured_path('DIR_WORK', default_mode=a.submit)) /
                           '.conda/envs/allele3-evolution/bin/python')
        elif not a.python:
            conda = Path(os.environ.get('CONDA_PREFIX', '')) / 'bin/python'
            a.python = str(conda) if os.environ.get('CONDA_PREFIX') and os.access(conda, os.X_OK) else sys.executable


# Native scientific sources are embedded. Submission helpers are snapshotted
# from the shared implementation; no private scheduler is embedded here.
EMBEDDED_RESOURCES = {'allele3/f/alphagenome_merge.cpp': '// AlphaGenome Atlas: bounded-memory TSV partitioning and allele-aware VCF merge.\n// Build: c++ -O3 -std=c++17 -DNDEBUG alphagenome_merge.cpp -lz -o alphagenome_merge\n// Ordinary gzip inputs are decompressed once per successful partition phase.\n// All statistics below describe PREDICTION COVERAGE, never population genotypes.\n#include <algorithm>\n#include <array>\n#include <cerrno>\n#include <charconv>\n#include <chrono>\n#include <cmath>\n#include <cstdint>\n#include <cstdlib>\n#include <cstring>\n#include <filesystem>\n#include <fstream>\n#include <iomanip>\n#include <iostream>\n#include <limits>\n#include <map>\n#include <memory>\n#include <optional>\n#include <set>\n#include <sstream>\n#include <stdexcept>\n#include <string>\n#include <string_view>\n#include <vector>\n#include <zlib.h>\n\nnamespace fs = std::filesystem;\nusing U64 = std::uint64_t;\nconstexpr std::size_t BUFFER_SIZE = 1U << 20;\nconstexpr std::size_t MAX_LINE_SIZE = 1U << 20;\nconstexpr const char* PARTITION_SCHEMA = "alphagenome-tsv-partition-1";\nconstexpr const char* MERGE_SCHEMA = "alphagenome-vcf-merge-1";\nconst std::array<std::string, 3> KINDS = {"avi", "splicing", "features"};\nconst std::array<std::vector<std::string>, 3> SCORE_COLUMNS = {{\n    {"raw_score", "PHRED"},\n    {"alphagenome_splicing"},\n    {"MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",\n     "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",\n     "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",\n     "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY",\n     "IS_INSERTION", "IS_DELETION"}\n}};\nconst std::array<std::size_t, 3> FIELD_OFFSETS = {0, 2, 3};\nconst std::array<U64, 25> CONTIG_LENGTHS = {\n    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,\n    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,\n    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,\n    58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569\n};\n\n[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }\n\nstd::string json_quote(std::string_view value) {\n    std::string out = "\\"";\n    static constexpr char HEX[] = "0123456789abcdef";\n    for (unsigned char c : value) {\n        switch (c) {\n            case \'"\': out += "\\\\\\""; break;\n            case \'\\\\\': out += "\\\\\\\\"; break;\n            case \'\\n\': out += "\\\\n"; break;\n            case \'\\r\': out += "\\\\r"; break;\n            case \'\\t\': out += "\\\\t"; break;\n            default:\n                if (c < 32) { out += "\\\\u00"; out += HEX[c >> 4]; out += HEX[c & 15]; }\n                else out += static_cast<char>(c);\n        }\n    }\n    out += \'"\';\n    return out;\n}\n\nU64 parse_uint(std::string_view text, const std::string& name) {\n    U64 result = 0;\n    auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);\n    if (text.empty() || parsed.ec != std::errc() || parsed.ptr != text.data() + text.size())\n        fail("Invalid unsigned integer for " + name + ": " + std::string(text));\n    return result;\n}\n\n// A strict parser for the flat scalar JSON receipts written by this executable.\n// No external JSON dependency is needed on the analysis host.\nstd::map<std::string, std::string> read_receipt(const fs::path& path) {\n    std::ifstream input(path, std::ios::binary);\n    if (!input) fail("Cannot read receipt: " + path.string());\n    std::string data((std::istreambuf_iterator<char>(input)), {});\n    if (data.size() > 65536) fail("Oversized receipt: " + path.string());\n    std::size_t at = 0;\n    auto whitespace = [&]() { while (at < data.size() && std::isspace(static_cast<unsigned char>(data[at]))) ++at; };\n    auto expect = [&](char c) {\n        whitespace();\n        if (at >= data.size() || data[at++] != c) fail("Malformed receipt: " + path.string());\n    };\n    auto string_value = [&]() {\n        expect(\'"\');\n        std::string value;\n        while (at < data.size()) {\n            unsigned char c = data[at++];\n            if (c == \'"\') return value;\n            if (c < 32) fail("Invalid control character in receipt: " + path.string());\n            if (c != \'\\\\\') { value += static_cast<char>(c); continue; }\n            if (at == data.size()) fail("Incomplete JSON escape in " + path.string());\n            switch (data[at++]) {\n                case \'"\': value += \'"\'; break;\n                case \'\\\\\': value += \'\\\\\'; break;\n                case \'/\': value += \'/\'; break;\n                case \'n\': value += \'\\n\'; break;\n                case \'r\': value += \'\\r\'; break;\n                case \'t\': value += \'\\t\'; break;\n                case \'b\': value += \'\\b\'; break;\n                case \'f\': value += \'\\f\'; break;\n                case \'u\': {\n                    if (at + 4 > data.size() || data[at] != \'0\' || data[at + 1] != \'0\')\n                        fail("Unsupported JSON escape in " + path.string());\n                    auto nibble = [&](char h) -> unsigned {\n                        if (h >= \'0\' && h <= \'9\') return h - \'0\';\n                        if (h >= \'a\' && h <= \'f\') return h - \'a\' + 10;\n                        if (h >= \'A\' && h <= \'F\') return h - \'A\' + 10;\n                        fail("Invalid JSON escape in " + path.string());\n                    };\n                    value += static_cast<char>((nibble(data[at + 2]) << 4) | nibble(data[at + 3]));\n                    at += 4; break;\n                }\n                default: fail("Invalid JSON escape in " + path.string());\n            }\n        }\n        fail("Unterminated JSON string in " + path.string());\n    };\n    std::map<std::string, std::string> values;\n    expect(\'{\'); whitespace();\n    if (at < data.size() && data[at] == \'}\') { ++at; }\n    else {\n        for (;;) {\n            const auto key = string_value(); expect(\':\'); whitespace();\n            std::string value;\n            if (at < data.size() && data[at] == \'"\') value = string_value();\n            else {\n                const auto begin = at;\n                while (at < data.size() && data[at] >= \'0\' && data[at] <= \'9\') ++at;\n                if (at == begin) fail("Invalid scalar in receipt: " + path.string());\n                value = data.substr(begin, at - begin);\n            }\n            if (!values.emplace(key, value).second) fail("Duplicate receipt key: " + key);\n            whitespace();\n            if (at < data.size() && data[at] == \'}\') { ++at; break; }\n            expect(\',\');\n        }\n    }\n    whitespace();\n    if (at != data.size()) fail("Trailing JSON data in " + path.string());\n    return values;\n}\n\nvoid atomic_text(const fs::path& path, const std::string& content) {\n    if (!path.parent_path().empty()) fs::create_directories(path.parent_path());\n    fs::path temporary = path; temporary += ".partial";\n    {\n        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);\n        if (!out) fail("Cannot write " + temporary.string());\n        out.write(content.data(), static_cast<std::streamsize>(content.size()));\n        out.flush();\n        if (!out) fail("Failed writing " + temporary.string());\n        out.close();\n        if (!out) fail("Failed closing " + temporary.string());\n    }\n    fs::rename(temporary, path);\n}\n\nstd::pair<std::string, int> chromosome(std::string_view text) {\n    if (text.substr(0, 3) == "chr") text.remove_prefix(3);\n    if (text == "X") return {"chrX", 22};\n    if (text == "Y") return {"chrY", 23};\n    if (text == "M" || text == "MT") return {"chrM", 24};\n    if (text.empty() || text.size() > 2 || text[0] == \'0\')\n        fail("Unsupported chromosome: " + std::string(text));\n    const auto n = parse_uint(text, "chromosome");\n    if (n < 1 || n > 22) fail("Unsupported chromosome: " + std::string(text));\n    return {"chr" + std::to_string(n), static_cast<int>(n - 1)};\n}\n\nint base_index(char base) {\n    switch (base) { case \'A\': return 0; case \'C\': return 1; case \'G\': return 2; case \'T\': return 3; }\n    return -1;\n}\n\nclass GzipReader {\n    gzFile file_ = nullptr;\n    std::vector<char> buffer_ = std::vector<char>(BUFFER_SIZE);\n    std::size_t start_ = 0, end_ = 0;\n    bool eof_ = false;\n    std::string path_;\npublic:\n    explicit GzipReader(const fs::path& path) : path_(path.string()) {\n        std::ifstream magic(path, std::ios::binary);\n        unsigned char prefix[2] = {0, 0};\n        magic.read(reinterpret_cast<char*>(prefix), 2);\n        if (magic.gcount() != 2 || prefix[0] != 0x1f || prefix[1] != 0x8b)\n            fail("Expected gzip input: " + path_);\n        file_ = gzopen(path_.c_str(), "rb");\n        if (!file_) fail("Cannot open gzip input: " + path_);\n        if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0) {\n            gzclose(file_); file_ = nullptr; fail("Cannot allocate gzip buffer for " + path_);\n        }\n    }\n    GzipReader(const GzipReader&) = delete;\n    GzipReader& operator=(const GzipReader&) = delete;\n    ~GzipReader() { if (file_) gzclose(file_); }\n    bool next(std::string& line) {\n        line.clear();\n        for (;;) {\n            if (start_ < end_) {\n                const auto* found = static_cast<const char*>(std::memchr(buffer_.data() + start_, \'\\n\', end_ - start_));\n                const auto stop = found ? static_cast<std::size_t>(found - buffer_.data()) : end_;\n                line.append(buffer_.data() + start_, stop - start_);\n                start_ = stop + (found ? 1 : 0);\n                if (line.size() > MAX_LINE_SIZE) fail("Oversized TSV line in " + path_);\n                if (found) { if (!line.empty() && line.back() == \'\\r\') line.pop_back(); return true; }\n            }\n            if (eof_) {\n                if (!line.empty() && line.back() == \'\\r\') line.pop_back();\n                return !line.empty();\n            }\n            const int count = gzread(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));\n            int code = Z_OK;\n            const char* error = gzerror(file_, &code);\n            if (count < 0 || (code != Z_OK && code != Z_STREAM_END))\n                fail("gzip read/CRC/truncation error in " + path_ + ": " + (error ? error : "unknown"));\n            start_ = 0;\n            end_ = count > 0 ? static_cast<std::size_t>(count) : 0;\n            if (count == 0) {\n                if (!gzeof(file_)) fail("Unexpected gzip EOF in " + path_);\n                const int close_code = gzclose(file_); file_ = nullptr;\n                if (close_code != Z_OK) fail("gzip close/truncation error in " + path_);\n                eof_ = true;\n            }\n        }\n    }\n};\n\nclass CanonicalSink {\n    gzFile file_ = nullptr;\n    std::string buffer_;\n    std::string path_;\n    bool closed_ = false;\npublic:\n    U64 uncompressed_bytes = 0;\n    uLong crc = crc32(0L, Z_NULL, 0);\n    explicit CanonicalSink(const std::optional<fs::path>& path) {\n        buffer_.reserve(BUFFER_SIZE + 4096);\n        if (path) {\n            path_ = path->string(); file_ = gzopen(path_.c_str(), "wb1");\n            if (!file_) fail("Cannot write staging gzip: " + path_);\n            if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0) {\n                gzclose(file_); file_ = nullptr; fail("Cannot allocate gzip buffer for " + path_);\n            }\n        }\n    }\n    ~CanonicalSink() { if (file_) gzclose(file_); }\n    void flush() {\n        if (buffer_.empty()) return;\n        crc = crc32(crc, reinterpret_cast<const Bytef*>(buffer_.data()), static_cast<uInt>(buffer_.size()));\n        uncompressed_bytes += buffer_.size();\n        if (file_) {\n            const auto written = gzwrite(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));\n            if (written != static_cast<int>(buffer_.size())) fail("Failed writing staging gzip: " + path_);\n        }\n        buffer_.clear();\n    }\n    void append(std::string_view value) { buffer_.append(value.data(), value.size()); if (buffer_.size() >= BUFFER_SIZE) flush(); }\n    void finish() {\n        if (closed_) return;\n        flush();\n        if (file_) {\n            const int result = gzclose(file_); file_ = nullptr;\n            if (result != Z_OK) fail("Failed closing staging gzip: " + path_);\n        }\n        closed_ = true;\n    }\n};\n\nvoid split_fields(std::string_view line, std::vector<std::string_view>& fields) {\n    fields.clear();\n    std::size_t at = 0;\n    while (at < line.size()) {\n        while (at < line.size() && (line[at] == \' \' || line[at] == \'\\t\')) ++at;\n        if (at == line.size()) break;\n        const auto begin = at;\n        while (at < line.size() && line[at] != \' \' && line[at] != \'\\t\') ++at;\n        fields.emplace_back(line.data() + begin, at - begin);\n    }\n}\n\nstd::string header_for(int kind) {\n    std::string line = "#CHROM\\tPOS\\tREF\\tALT";\n    for (const auto& column : SCORE_COLUMNS[kind]) { line += \'\\t\'; line += column; }\n    line += \'\\n\'; return line;\n}\n\nvoid validate_header(std::string_view line, int kind, const std::string& context) {\n    std::vector<std::string_view> fields;\n    split_fields(line, fields);\n    if (fields.size() != SCORE_COLUMNS[kind].size() + 4 || fields[0] != "#CHROM" ||\n        fields[1] != "POS" || fields[2] != "REF" || fields[3] != "ALT")\n        fail("Unexpected " + KINDS[kind] + " header in " + context + "; expected " + header_for(kind));\n    for (std::size_t j = 0; j < SCORE_COLUMNS[kind].size(); ++j)\n        if (fields[j + 4] != SCORE_COLUMNS[kind][j])\n            fail("Unexpected " + KINDS[kind] + " score column in " + context + ": " + std::string(fields[j + 4]));\n}\n\nstruct Numeric { bool missing = false; bool zero = false; };\n\nNumeric validate_number(std::string_view token, const std::string& context) {\n    if (token == "." || token == "NA" || token == "NaN" || token == "nan") return {true, false};\n    std::size_t at = 0, digits = 0;\n    bool nonzero = false;\n    if (at < token.size() && (token[at] == \'+\' || token[at] == \'-\')) ++at;\n    while (at < token.size() && token[at] >= \'0\' && token[at] <= \'9\') {\n        nonzero |= token[at] != \'0\'; ++at; ++digits;\n    }\n    if (at < token.size() && token[at] == \'.\') {\n        ++at;\n        while (at < token.size() && token[at] >= \'0\' && token[at] <= \'9\') {\n            nonzero |= token[at] != \'0\'; ++at; ++digits;\n        }\n    }\n    if (digits == 0) fail("Invalid numeric token " + json_quote(token) + " in " + context);\n    if (at < token.size() && (token[at] == \'e\' || token[at] == \'E\')) {\n        ++at;\n        if (at < token.size() && (token[at] == \'+\' || token[at] == \'-\')) ++at;\n        const auto begin = at;\n        while (at < token.size() && token[at] >= \'0\' && token[at] <= \'9\') ++at;\n        if (at == begin) fail("Invalid exponent in " + context + ": " + std::string(token));\n    }\n    if (at != token.size()) fail("Invalid numeric token " + json_quote(token) + " in " + context);\n    // GCC 8\'s libstdc++ implements integer from_chars only. Keep the fast\n    // floating parser where available; the portable branch preserves decimal\n    // text and underflow semantics and never round-trips output through double.\n#if defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE >= 11 && !defined(ALPHAGENOME_PORTABLE_FLOAT)\n    auto number = token;\n    if (!number.empty() && number.front() == \'+\') number.remove_prefix(1);\n    double value = 0;\n    const auto result = std::from_chars(number.data(), number.data() + number.size(), value, std::chars_format::general);\n    if (result.ec == std::errc::result_out_of_range) {\n        // Very small finite decimal text may underflow binary64. Preserve its text\n        // and do not count it as numeric zero; only overflow is an error.\n        std::string copy(token);\n        char* end = nullptr;\n        value = std::strtod(copy.c_str(), &end);\n        if (end != copy.c_str() + copy.size() || !std::isfinite(value))\n            fail("Non-finite/overflow numeric token " + json_quote(token) + " in " + context);\n    } else if (result.ec != std::errc() || result.ptr != number.data() + number.size() || !std::isfinite(value))\n        fail("Non-finite numeric token " + json_quote(token) + " in " + context);\n#else\n    std::string copy(token);\n    char* end = nullptr;\n    const double value = std::strtod(copy.c_str(), &end);\n    if (end != copy.c_str() + copy.size() || !std::isfinite(value))\n        fail("Non-finite/overflow numeric token " + json_quote(token) + " in " + context);\n#endif\n    return {false, !nonzero};\n}\n\nstruct ParsedRow {\n    std::string chrom;\n    int chrom_index = 0;\n    U64 pos = 0;\n    char ref = \'N\', alt = \'N\';\n    int alt_index = -1;\n};\n\nParsedRow parse_row(const std::vector<std::string_view>& fields, int kind, const std::string& context,\n                    std::vector<Numeric>* numbers = nullptr) {\n    if (fields.size() != SCORE_COLUMNS[kind].size() + 4)\n        fail("Wrong field count in " + context + ": expected " + std::to_string(SCORE_COLUMNS[kind].size() + 4) +\n             ", got " + std::to_string(fields.size()));\n    ParsedRow row;\n    auto chr = chromosome(fields[0]); row.chrom = std::move(chr.first); row.chrom_index = chr.second;\n    row.pos = parse_uint(fields[1], "POS in " + context);\n    if (!row.pos || row.pos > CONTIG_LENGTHS[row.chrom_index])\n        fail("Position outside GRCh38 primary chromosome in " + context);\n    if (fields[2].size() != 1 || fields[3].size() != 1 || base_index(fields[2][0]) < 0 ||\n        base_index(fields[3][0]) < 0 || fields[2][0] == fields[3][0])\n        fail("Expected canonical A/C/G/T single-base REF != ALT in " + context);\n    row.ref = fields[2][0]; row.alt = fields[3][0]; row.alt_index = base_index(row.alt);\n    if (numbers) { numbers->clear(); numbers->reserve(SCORE_COLUMNS[kind].size()); }\n    for (std::size_t j = 4; j < fields.size(); ++j) {\n        auto numeric = validate_number(fields[j], context);\n        if (numbers) numbers->push_back(numeric);\n    }\n    return row;\n}\n\nstruct PositionOrder {\n    U64 pos = 0;\n    char ref = \'N\';\n    unsigned alt_mask = 0;\n    void accept(const ParsedRow& row, const std::string& context) {\n        if (row.pos < pos) fail("Decreasing position in " + context);\n        if (row.pos != pos) { pos = row.pos; ref = row.ref; alt_mask = 0; }\n        else if (row.ref != ref) fail("Inconsistent REF at one position in " + context);\n        const auto bit = 1U << row.alt_index;\n        if (alt_mask & bit) fail("Duplicate CHROM/POS/REF/ALT in " + context);\n        alt_mask |= bit;\n    }\n};\n\nstruct ChromReceipt {\n    std::string chrom, kind, token;\n    U64 rows = 0, first_pos = 0, last_pos = 0, gzip_bytes = 0, uncompressed_bytes = 0, crc = 0;\n    std::string json() const {\n        std::ostringstream out;\n        out << "{\\n  \\"schema\\": " << json_quote(PARTITION_SCHEMA)\n            << ",\\n  \\"chrom\\": " << json_quote(chrom) << ",\\n  \\"kind\\": " << json_quote(kind)\n            << ",\\n  \\"token\\": " << json_quote(token) << ",\\n  \\"rows\\": " << rows\n            << ",\\n  \\"first_pos\\": " << first_pos << ",\\n  \\"last_pos\\": " << last_pos\n            << ",\\n  \\"gzip_bytes\\": " << gzip_bytes << ",\\n  \\"uncompressed_bytes\\": " << uncompressed_bytes\n            << ",\\n  \\"uncompressed_crc32\\": " << crc << "\\n}\\n";\n        return out.str();\n    }\n    static ChromReceipt load(const fs::path& path) {\n        auto values = read_receipt(path);\n        auto get = [&](const std::string& key) -> std::string {\n            auto it = values.find(key);\n            if (it == values.end()) fail("Missing receipt field " + key + " in " + path.string());\n            return it->second;\n        };\n        if (get("schema") != PARTITION_SCHEMA) fail("Unsupported receipt schema in " + path.string());\n        ChromReceipt receipt;\n        receipt.chrom = get("chrom"); receipt.kind = get("kind"); receipt.token = get("token");\n        receipt.rows = parse_uint(get("rows"), "receipt rows");\n        receipt.first_pos = parse_uint(get("first_pos"), "receipt first_pos");\n        receipt.last_pos = parse_uint(get("last_pos"), "receipt last_pos");\n        receipt.gzip_bytes = parse_uint(get("gzip_bytes"), "receipt gzip_bytes");\n        receipt.uncompressed_bytes = parse_uint(get("uncompressed_bytes"), "receipt uncompressed_bytes");\n        receipt.crc = parse_uint(get("uncompressed_crc32"), "receipt CRC");\n        return receipt;\n    }\n};\n\nbool staging_quick_check(const fs::path& path, const ChromReceipt& receipt) {\n    if (!fs::is_regular_file(path) || fs::file_size(path) != receipt.gzip_bytes || receipt.gzip_bytes < 18) return false;\n    std::ifstream input(path, std::ios::binary);\n    unsigned char bytes[8] = {};\n    input.seekg(-8, std::ios::end); input.read(reinterpret_cast<char*>(bytes), 8);\n    if (!input) return false;\n    auto little = [&](int offset) -> U64 {\n        return U64(bytes[offset]) | (U64(bytes[offset + 1]) << 8) | (U64(bytes[offset + 2]) << 16) | (U64(bytes[offset + 3]) << 24);\n    };\n    return little(0) == receipt.crc && little(4) == (receipt.uncompressed_bytes & 0xffffffffU);\n}\n\nstruct Arguments {\n    std::string mode;\n    std::map<std::string, std::string> values;\n    explicit Arguments(int argc, char** argv) {\n        if (argc < 2) fail("Specify partition or merge; use --help for usage");\n        mode = argv[1];\n        if (mode == "--help" || mode == "-h") return;\n        if (mode != "partition" && mode != "merge") fail("Unknown mode: " + mode);\n        const std::set<std::string> allowed = mode == "partition"\n            ? std::set<std::string>{"--kind", "--input", "--out-dir", "--token", "--progress-every", "--file-prefix"}\n            : std::set<std::string>{"--chrom", "--avi", "--splicing", "--features", "--stats", "--token", "--progress-every"};\n        for (int i = 2; i < argc; i += 2) {\n            if (i + 1 == argc) fail("Missing value for " + std::string(argv[i]));\n            if (!allowed.count(argv[i])) fail("Unknown option: " + std::string(argv[i]));\n            if (!values.emplace(argv[i], argv[i + 1]).second) fail("Duplicate option: " + std::string(argv[i]));\n        }\n    }\n    std::string required(const std::string& name) const {\n        auto found = values.find(name);\n        if (found == values.end() || found->second.empty()) fail("Required option: " + name);\n        return found->second;\n    }\n    U64 progress() const {\n        const auto found = values.find("--progress-every");\n        return found == values.end() ? 10000000 : parse_uint(found->second, "--progress-every");\n    }\n};\n\nint partition(const Arguments& args) {\n    const auto kind_name = args.required("--kind");\n    auto found = std::find(KINDS.begin(), KINDS.end(), kind_name);\n    if (found == KINDS.end()) fail("--kind must be avi, splicing, or features");\n    const int kind = static_cast<int>(found - KINDS.begin());\n    const fs::path input = args.required("--input"), out_dir = args.required("--out-dir");\n    const auto prefix = args.values.count("--file-prefix") ? args.values.at("--file-prefix") : std::string{};\n    if (!prefix.empty() && prefix != kind_name + ".") fail("--file-prefix must match kind plus \'.\'");\n    const auto token = args.required("--token");\n    const U64 progress = args.progress();\n    fs::create_directories(out_dir);\n    GzipReader reader(input);\n    std::string line;\n    if (!reader.next(line)) fail("Empty input without header: " + input.string());\n    validate_header(line, kind, input.string());\n    const auto header = header_for(kind);\n    std::vector<std::string_view> fields; fields.reserve(22);\n    std::string normalized; normalized.reserve(2048);\n    std::set<std::string> visited;\n    std::vector<ChromReceipt> completed;\n    std::unique_ptr<CanonicalSink> sink;\n    std::optional<ChromReceipt> previous;\n    ChromReceipt receipt;\n    PositionOrder order;\n    fs::path destination, temporary, done_path;\n    U64 source_rows = 0, line_number = 1, reused_chromosomes = 0;\n    auto finish_chromosome = [&]() {\n        if (!sink) return;\n        sink->finish();\n        receipt.uncompressed_bytes = sink->uncompressed_bytes;\n        receipt.crc = sink->crc;\n        if (previous) {\n            if (receipt.rows != previous->rows || receipt.first_pos != previous->first_pos ||\n                receipt.last_pos != previous->last_pos || receipt.uncompressed_bytes != previous->uncompressed_bytes ||\n                receipt.crc != previous->crc)\n                fail("Source content differs from completed chromosome despite matching token: " + receipt.chrom);\n            receipt.gzip_bytes = previous->gzip_bytes; ++reused_chromosomes;\n            std::cerr << "partition " << kind_name << \' \' << receipt.chrom << " reused; rows=" << receipt.rows << \'\\n\';\n        } else {\n            receipt.gzip_bytes = fs::file_size(temporary);\n            fs::rename(temporary, destination);\n            atomic_text(done_path, receipt.json());\n            std::cerr << "partition " << kind_name << \' \' << receipt.chrom << " completed; rows=" << receipt.rows\n                      << " gzip_bytes=" << receipt.gzip_bytes << \'\\n\';\n        }\n        completed.push_back(receipt); sink.reset(); previous.reset();\n    };\n    while (reader.next(line)) {\n        ++line_number;\n        const std::string context = input.string() + ":" + std::to_string(line_number);\n        split_fields(line, fields);\n        const auto row = parse_row(fields, kind, context);\n        if (!sink || row.chrom != receipt.chrom) {\n            finish_chromosome();\n            if (!visited.insert(row.chrom).second) fail("Chromosome block revisited in " + context + ": " + row.chrom);\n            receipt = ChromReceipt{}; receipt.chrom = row.chrom; receipt.kind = kind_name; receipt.token = token;\n            order = PositionOrder{};\n            destination = out_dir / (prefix + row.chrom + ".tsv.gz"); temporary = destination; temporary += ".partial";\n            done_path = out_dir / (prefix + row.chrom + ".done.json");\n            if (fs::exists(done_path)) {\n                const auto existing = ChromReceipt::load(done_path);\n                if (existing.token != token || existing.kind != kind_name || existing.chrom != row.chrom)\n                    fail("Incompatible staging receipt; use a distinct/clean staging directory: " + done_path.string());\n                if (staging_quick_check(destination, existing)) previous = existing;\n                else std::cerr << "partition " << kind_name << \' \' << row.chrom << " regenerating missing/damaged staging output\\n";\n            }\n            sink = std::make_unique<CanonicalSink>(previous ? std::nullopt : std::optional<fs::path>(temporary));\n            sink->append(header);\n        }\n        order.accept(row, context);\n        if (receipt.rows++ == 0) receipt.first_pos = row.pos;\n        receipt.last_pos = row.pos;\n        normalized.clear(); normalized += row.chrom;\n        for (std::size_t j = 1; j < fields.size(); ++j) { normalized += \'\\t\'; normalized.append(fields[j]); }\n        normalized += \'\\n\'; sink->append(normalized);\n        ++source_rows;\n        if (progress && source_rows % progress == 0)\n            std::cerr << "partition " << kind_name << " rows=" << source_rows << " at=" << row.chrom << \':\' << row.pos << \'\\n\';\n    }\n    finish_chromosome();\n    // Old source.done is never sufficient to reuse outputs: the orchestrator also\n    // checks source fingerprints and every chromosome receipt. An interrupted\n    // ordinary-gzip source can require rereading its prefix on the next attempt.\n    std::ostringstream summary;\n    summary << "{\\n  \\"schema\\": " << json_quote(PARTITION_SCHEMA) << ",\\n  \\"kind\\": " << json_quote(kind_name)\n            << ",\\n  \\"token\\": " << json_quote(token) << ",\\n  \\"source_rows\\": " << source_rows\n            << ",\\n  \\"reused_chromosomes\\": " << reused_chromosomes << ",\\n  \\"chromosomes\\": [";\n    for (std::size_t j = 0; j < completed.size(); ++j) {\n        if (j) summary << \',\';\n        summary << \'\\n\' << completed[j].json();\n    }\n    summary << "\\n  ]\\n}\\n";\n    atomic_text(out_dir / (prefix + "source.done.json"), summary.str());\n    std::cerr << "partition " << kind_name << " finished; source_rows=" << source_rows << " chromosomes=" << completed.size() << \'\\n\';\n    return 0;\n}\n\nstruct ScoreRow {\n    ParsedRow locus;\n    std::vector<std::string> values;\n    std::vector<Numeric> numbers;\n};\n\nstruct Group {\n    U64 pos = 0;\n    char ref = \'N\';\n    std::array<std::optional<ScoreRow>, 4> alleles;\n};\n\nclass SourceStream {\n    int kind_;\n    std::string expected_chrom_, path_, line_;\n    GzipReader reader_;\n    PositionOrder order_;\n    std::vector<std::string_view> fields_;\n    U64 line_number_ = 1;\npublic:\n    U64 rows = 0;\n    std::optional<ScoreRow> current;\n    SourceStream(int kind, const std::string& path, const std::string& chrom)\n        : kind_(kind), expected_chrom_(chrom), path_(path), reader_(path) {\n        fields_.reserve(22);\n        if (!reader_.next(line_)) fail("Empty input without header: " + path_);\n        validate_header(line_, kind_, path_); advance();\n    }\n    void advance() {\n        if (!reader_.next(line_)) { current.reset(); return; }\n        ++line_number_;\n        const std::string context = path_ + ":" + std::to_string(line_number_);\n        split_fields(line_, fields_);\n        ScoreRow row;\n        row.locus = parse_row(fields_, kind_, context, &row.numbers);\n        if (row.locus.chrom != expected_chrom_) fail("Unexpected chromosome in " + context + ": " + row.locus.chrom);\n        order_.accept(row.locus, context);\n        row.values.reserve(SCORE_COLUMNS[kind_].size());\n        for (std::size_t j = 4; j < fields_.size(); ++j) row.values.emplace_back(fields_[j]);\n        current = std::move(row); ++rows;\n    }\n    Group take_group() {\n        if (!current) fail("Internal error: group requested after EOF");\n        Group group; group.pos = current->locus.pos; group.ref = current->locus.ref;\n        do {\n            const int alt = current->locus.alt_index;\n            group.alleles[alt] = std::move(*current);\n            advance();\n        } while (current && current->locus.pos == group.pos);\n        return group;\n    }\n};\n\nstruct FieldStats {\n    std::string name;\n    int source = 0;\n    U64 numeric_present = 0, numeric_missing = 0, numeric_zero = 0, absent_source_alleles = 0;\n    U64 sites_with_three_alt_scores = 0;\n};\n\nstruct SourceStats {\n    U64 covered_alleles = 0, covered_sites = 0, sites_with_three_alt_rows = 0;\n    U64 sites_with_three_alt_complete_scores = 0;\n    U64 numeric_present = 0, numeric_missing = 0, numeric_zero = 0;\n};\n\nstd::array<FieldStats, 21> make_field_stats() {\n    std::array<FieldStats, 21> fields;\n    fields[0].name = "AG_AVI_RAW"; fields[0].source = 0;\n    fields[1].name = "AG_AVI_PHRED"; fields[1].source = 0;\n    fields[2].name = "AG_SPLICE"; fields[2].source = 1;\n    for (std::size_t j = 0; j < SCORE_COLUMNS[2].size(); ++j) {\n        fields[j + 3].name = "AG_SHAP_" + SCORE_COLUMNS[2][j]; fields[j + 3].source = 2;\n    }\n    return fields;\n}\n\nvoid vcf_header(const std::string& chrom, const std::string& token, const std::array<FieldStats, 21>& fields) {\n    const auto chr = chromosome(chrom);\n    std::cout << "##fileformat=VCFv4.3\\n##source=allele3_alphagenome_merge_1\\n##reference=GRCh38\\n"\n              << "##AG_schema=" << MERGE_SCHEMA << "\\n##AG_recipe_token=" << json_quote(token) << \'\\n\'\n              << "##AG_allele_semantics=\\"ALT contains published predicted substitutions; presence is not evidence of population observation\\"\\n"\n              << "##AG_missing_semantics=\\"Dot means unavailable; AG_SOURCE_MASK distinguishes absent rows from explicit missing scores\\"\\n"\n              << "##AG_feature_semantics=\\"AG_SHAP fields are published signed feature contributions, not raw molecular-effect measurements\\"\\n"\n              << "##contig=<ID=" << chrom << ",length=" << CONTIG_LENGTHS[chr.second] << ",assembly=GRCh38>\\n";\n    for (std::size_t j = 0; j < fields.size(); ++j) {\n        std::string description;\n        if (j == 0) description = "Published AlphaGenome variant impact raw score; not a genotype statistic";\n        else if (j == 1) description = "Published AlphaGenome variant impact PHRED score; not genotype quality";\n        else if (j == 2) description = "Published combined AlphaGenome splicing score; distinct from SHAP MERGED_SPLICING";\n        else description = "Published signed AVI feature-importance (SHAP) value for " + SCORE_COLUMNS[2][j - 3] + "; not the raw feature score";\n        std::cout << "##INFO=<ID=" << fields[j].name << ",Number=A,Type=Float,Description=" << json_quote(description) << ">\\n";\n    }\n    std::cout << "##INFO=<ID=AG_SOURCE_MASK,Number=A,Type=Integer,Description=\\"Published row presence: AVI=1, splicing=2, feature importance=4; bitwise union\\">\\n"\n              << "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n";\n}\n\nint merge(const Arguments& args) {\n    const auto chrom = chromosome(args.required("--chrom")).first;\n    const auto token = args.required("--token");\n    const fs::path stats_path = args.required("--stats");\n    const U64 progress = args.progress();\n    std::array<std::unique_ptr<SourceStream>, 3> inputs;\n    bool supplied = false;\n    for (int kind = 0; kind < 3; ++kind) {\n        const auto path = args.required("--" + KINDS[kind]);\n        if (path != "NONE") { inputs[kind] = std::make_unique<SourceStream>(kind, path, chrom); supplied = true; }\n    }\n    if (!supplied) fail("At least one source file must be supplied; all three were NONE");\n    auto fields = make_field_stats();\n    std::array<SourceStats, 3> source_stats;\n    std::array<U64, 4> sites_by_nalt{};\n    std::array<U64, 8> source_mask_alleles{};\n    U64 sites = 0, union_alt_rows = 0, first_pos = 0, last_pos = 0, next_progress = progress;\n    vcf_header(chrom, token, fields);\n    std::string output; output.reserve(4096);\n    constexpr char BASES[] = "ACGT";\n    for (;;) {\n        U64 pos = std::numeric_limits<U64>::max();\n        for (const auto& source : inputs) if (source && source->current) pos = std::min(pos, source->current->locus.pos);\n        if (pos == std::numeric_limits<U64>::max()) break;\n        std::array<std::optional<Group>, 3> groups;\n        char ref = \'N\';\n        std::array<unsigned, 4> masks{};\n        for (int kind = 0; kind < 3; ++kind) {\n            if (!inputs[kind] || !inputs[kind]->current || inputs[kind]->current->locus.pos != pos) continue;\n            groups[kind] = inputs[kind]->take_group();\n            if (ref == \'N\') ref = groups[kind]->ref;\n            else if (ref != groups[kind]->ref)\n                fail("Cross-source REF mismatch at " + chrom + ":" + std::to_string(pos));\n            U64 allele_count = 0;\n            bool complete = true;\n            for (int alt = 0; alt < 4; ++alt) {\n                if (!groups[kind]->alleles[alt]) continue;\n                ++allele_count; masks[alt] |= 1U << kind;\n                const auto& row = *groups[kind]->alleles[alt];\n                for (std::size_t j = 0; j < row.numbers.size(); ++j) {\n                    auto& field = fields[FIELD_OFFSETS[kind] + j];\n                    if (row.numbers[j].missing) { ++field.numeric_missing; ++source_stats[kind].numeric_missing; complete = false; }\n                    else {\n                        ++field.numeric_present; ++source_stats[kind].numeric_present;\n                        if (row.numbers[j].zero) { ++field.numeric_zero; ++source_stats[kind].numeric_zero; }\n                    }\n                }\n            }\n            source_stats[kind].covered_alleles += allele_count;\n            ++source_stats[kind].covered_sites;\n            if (allele_count == 3) {\n                ++source_stats[kind].sites_with_three_alt_rows;\n                if (complete) ++source_stats[kind].sites_with_three_alt_complete_scores;\n                for (std::size_t j = 0; j < SCORE_COLUMNS[kind].size(); ++j) {\n                    bool field_complete = true;\n                    for (int alt = 0; alt < 4; ++alt)\n                        if (groups[kind]->alleles[alt] && groups[kind]->alleles[alt]->numbers[j].missing) field_complete = false;\n                    if (field_complete) ++fields[FIELD_OFFSETS[kind] + j].sites_with_three_alt_scores;\n                }\n            }\n        }\n        unsigned nalt = 0;\n        for (unsigned mask : masks) if (mask) { ++nalt; ++source_mask_alleles[mask]; }\n        if (nalt < 1 || nalt > 3) fail("Invalid ALT union at " + chrom + ":" + std::to_string(pos));\n        if (++sites == 1) first_pos = pos;\n        last_pos = pos;\n        ++sites_by_nalt[nalt]; union_alt_rows += nalt;\n        output.clear(); output += chrom; output += \'\\t\'; output += std::to_string(pos);\n        output += "\\t.\\t"; output += ref; output += \'\\t\';\n        bool first = true;\n        for (int alt = 0; alt < 4; ++alt) if (masks[alt]) {\n            if (!first) output += \',\';\n            output += BASES[alt]; first = false;\n        }\n        output += "\\t.\\t.\\t";\n        for (std::size_t field_index = 0; field_index < fields.size(); ++field_index) {\n            auto& field = fields[field_index];\n            if (field_index) output += \';\';\n            output += field.name; output += \'=\'; first = true;\n            const auto offset = field_index - FIELD_OFFSETS[field.source];\n            for (int alt = 0; alt < 4; ++alt) if (masks[alt]) {\n                if (!first) output += \',\';\n                first = false;\n                if (!groups[field.source] || !groups[field.source]->alleles[alt]) {\n                    output += \'.\'; ++field.absent_source_alleles;\n                } else {\n                    const auto& row = *groups[field.source]->alleles[alt];\n                    if (row.numbers[offset].missing) output += \'.\';\n                    else output += row.values[offset];\n                }\n            }\n        }\n        output += ";AG_SOURCE_MASK="; first = true;\n        for (unsigned mask : masks) if (mask) {\n            if (!first) output += \',\';\n            output += static_cast<char>(\'0\' + mask); first = false;\n        }\n        output += \'\\n\';\n        std::cout.write(output.data(), static_cast<std::streamsize>(output.size()));\n        if (!std::cout) fail("VCF stdout write failed (downstream compressor may have failed)");\n        if (progress && union_alt_rows >= next_progress) {\n            std::cerr << "merge " << chrom << " sites=" << sites << " union_alt_rows=" << union_alt_rows << " pos=" << pos << \'\\n\';\n            next_progress = union_alt_rows + progress;\n        }\n    }\n    std::cout.flush();\n    if (!std::cout) fail("VCF stdout flush failed");\n    std::ostringstream stats;\n    stats << "{\\n  \\"schema\\": " << json_quote(MERGE_SCHEMA) << ",\\n  \\"chrom\\": " << json_quote(chrom)\n          << ",\\n  \\"token\\": " << json_quote(token) << ",\\n  \\"assembly\\": \\"GRCh38\\",\\n  \\"source_rows\\": {";\n    for (int kind = 0; kind < 3; ++kind) { if (kind) stats << \',\'; stats << json_quote(KINDS[kind]) << \':\' << (inputs[kind] ? inputs[kind]->rows : 0); }\n    stats << "},\\n  \\"union_alt_rows\\": " << union_alt_rows << ",\\n  \\"sites\\": " << sites\n          << ",\\n  \\"first_pos\\": " << first_pos << ",\\n  \\"last_pos\\": " << last_pos\n          << ",\\n  \\"sites_by_nalt\\": {\\"1\\":" << sites_by_nalt[1] << ",\\"2\\":" << sites_by_nalt[2] << ",\\"3\\":" << sites_by_nalt[3] << \'}\'\n          << ",\\n  \\"prediction_base_coverage\\": {\\"exactly_2_bases\\":" << sites_by_nalt[1] << ",\\"exactly_3_bases\\":" << sites_by_nalt[2]\n          << ",\\"exactly_4_bases\\":" << sites_by_nalt[3] << ",\\"at_least_3_bases\\":" << sites_by_nalt[2] + sites_by_nalt[3] << \'}\'\n          << ",\\n  \\"source_mask_alleles\\": {";\n    for (unsigned mask = 1; mask <= 7; ++mask) { if (mask > 1) stats << \',\'; stats << json_quote(std::to_string(mask)) << \':\' << source_mask_alleles[mask]; }\n    stats << "},\\n  \\"per_source\\": {";\n    for (int kind = 0; kind < 3; ++kind) {\n        if (kind) stats << \',\';\n        const auto& s = source_stats[kind];\n        stats << \'\\n\' << json_quote(KINDS[kind]) << ":{\\"covered_alleles\\":" << s.covered_alleles << ",\\"covered_sites\\":" << s.covered_sites\n              << ",\\"sites_with_three_alt_rows\\":" << s.sites_with_three_alt_rows\n              << ",\\"sites_with_three_alt_complete_scores\\":" << s.sites_with_three_alt_complete_scores\n              << ",\\"numeric_present\\":" << s.numeric_present << ",\\"numeric_missing\\":" << s.numeric_missing << ",\\"numeric_zero\\":" << s.numeric_zero << \'}\';\n    }\n    stats << "\\n  },\\n  \\"per_field\\": {";\n    for (std::size_t j = 0; j < fields.size(); ++j) {\n        if (j) stats << \',\';\n        const auto& f = fields[j];\n        stats << \'\\n\' << json_quote(f.name) << ":{\\"source\\":" << json_quote(KINDS[f.source])\n              << ",\\"numeric_present\\":" << f.numeric_present << ",\\"numeric_missing\\":" << f.numeric_missing\n              << ",\\"numeric_zero\\":" << f.numeric_zero << ",\\"absent_source_alleles\\":" << f.absent_source_alleles\n              << ",\\"sites_with_three_alt_scores\\":" << f.sites_with_three_alt_scores << \'}\';\n    }\n    stats << "\\n  },\\n  \\"observed_tri_allelic_snvs\\": null,\\n  \\"observed_quad_allelic_snvs\\": null,\\n"\n          << "  \\"observed_allele_statistics_reason\\": \\"Unavailable: Atlas predicts possible substitutions; population AC/AN/AF or genotypes are required to identify observed tri-/quad-allelic SNVs.\\",\\n"\n          << "  \\"prediction_coverage_definition\\": \\"Number of unique published ALT bases plus the shared reference base; no unlisted ALT is invented.\\"\\n}\\n";\n    atomic_text(stats_path, stats.str());\n    std::cerr << "merge " << chrom << " finished; sites=" << sites << " union_alt_rows=" << union_alt_rows\n              << " predicted_ref_plus_2alt=" << sites_by_nalt[2] << " predicted_ref_plus_3alt=" << sites_by_nalt[3] << \'\\n\';\n    return 0;\n}\n\nint main(int argc, char** argv) {\n    std::ios::sync_with_stdio(false);\n    std::cin.tie(nullptr);\n    try {\n        const Arguments args(argc, argv);\n        if (args.mode == "--help" || args.mode == "-h") {\n            std::cout << "AlphaGenome Atlas prediction coverage converter (C++17 + zlib)\\n"\n                      << "partition --kind avi|splicing|features --input FILE.gz --out-dir DIR --token TOKEN [--progress-every N]\\n"\n                      << "merge --chrom CHR --avi FILE.gz|NONE --splicing FILE.gz|NONE --features FILE.gz|NONE --stats FILE.json --token TOKEN [--progress-every N]\\n"\n                      << "Merge writes VCFv4.3 to stdout; progress/errors go to stderr. Inputs must contain sorted positions in contiguous chromosome blocks.\\n"\n                      << "Missing tokens: \'.\', \'NA\', \'NaN\', \'nan\'. Zero and signed finite score text are retained.\\n"\n                      << "All allele summaries describe predicted substitutions, not observed population alleles.\\n";\n            return 0;\n        }\n        return args.mode == "partition" ? partition(args) : merge(args);\n    } catch (const std::exception& error) {\n        std::cerr << "ERROR: " << error.what() << \'\\n\';\n        return 2;\n    }\n}\n', 'allele3/f/analysis_integrity.py': '"""Content identities and atomic, locked publication for derived analyses.\n\nPermanent inputs are read only. A failed refresh invalidates current evidence,\nwhile preserving the previous files and receipt in an attempt history.\n"""\nfrom __future__ import annotations\nimport copy\nfrom dataclasses import dataclass\nimport datetime as dt\nimport fcntl\nimport functools\nimport hashlib\nimport inspect\nimport json\nimport os\nfrom pathlib import Path\nimport sys\nimport tempfile\nimport uuid\n\nTERMINAL = {\n    \'complete\', \'complete_descriptive_imputed_panel\', \'complete_descriptive_all_allele_panel\',\n    \'complete_descriptive_unverified\', \'complete_adjustment_unavailable\', \'planned\',\n    \'no_usable_exact_high_quality_contrast\', \'no_usable_contrast\', \'not_estimable\',\n    \'not_available\', \'not_run\', \'complete_with_exclusions\',\n}\nEVIDENCE = TERMINAL - {\'planned\', \'not_available\', \'not_run\'}\nSMALL_LIMIT = 64 * 1024 * 1024\n\n\ndef sha256(path):\n    digest = hashlib.sha256()\n    with Path(path).open(\'rb\') as handle:\n        for block in iter(lambda: handle.read(1024 * 1024), b\'\'):\n            digest.update(block)\n    return digest.hexdigest()\n\n\ndef identity(path, full=False):\n    path = Path(path).expanduser().resolve()\n    if not path.is_file():\n        return {\'path\': str(path), \'status\': \'missing\'}\n    s = path.stat()\n    result = {\'path\': str(path), \'size\': s.st_size, \'mtime_ns\': s.st_mtime_ns}\n    if full:\n        result.update(sha256=sha256(path), verification=\'full_sha256_critical_input_only\')\n    else:\n        result.update(verification=\'size_mtime_only\')\n    if path.stat().st_mtime_ns != s.st_mtime_ns or path.stat().st_size != s.st_size:\n        raise ValueError(f\'Input changed while hashing: {path}\')\n    return result\n\n\ndef atomic_json(path, value):\n    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)\n    temporary = path.with_name(path.name + \'.\' + uuid.uuid4().hex + \'.partial\')\n    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False, default=str) + \'\\n\')\n    os.replace(temporary, path)\n\n\ndef code_identity():\n    root = Path(__file__).resolve().parent.parent\n    files = {str(p.relative_to(root)): {\'size\': p.stat().st_size, \'mtime_ns\': p.stat().st_mtime_ns}\n             for p in sorted(root.rglob(\'*\')) if p.is_file() and p.suffix in {\'.py\', \'.sh\', \'.env\', \'.yml\', \'.html\', \'.js\'} and \'__pycache__\' not in p.parts}\n    return {\'files\': files, \'identity_policy\': \'file_size_mtime; no routine source hashing\',\n            \'python\': sys.version, \'executable\': sys.executable, \'git_commit\': \'unavailable_local_checkout_without_git\'}\n\n\ndef inventory(root, manifest):\n    root = Path(root)\n    return {str(p.relative_to(root)): {k: v for k, v in identity(p).items() if k != \'path\'}\n            for p in sorted(root.rglob(\'*\')) if p.is_file() and p != root / manifest}\n\n\ndef validate_manifest(root, name=\'manifest.json\', require_evidence=False):\n    root = Path(root)\n    result = json.loads((root / name).read_text())\n    allowed = EVIDENCE if require_evidence else TERMINAL\n    if result.get(\'status\') not in allowed:\n        raise ValueError(f\'Nonterminal/unusable result: {root}: {result.get("status", "missing")}\')\n    outputs = result.get(\'output_inventory\')\n    if not isinstance(outputs, dict) or not outputs:\n        raise ValueError(f\'Missing verified output inventory: {root / name}\')\n    for rel, expected in outputs.items():\n        path = root / rel\n        if Path(rel).is_absolute() or \'..\' in Path(rel).parts or path.is_symlink():\n            raise ValueError(f\'Unsafe inventory entry: {rel}\')\n        actual = identity(path, full=\'sha256\' in expected)\n        keys = (\'size\', \'sha256\') if \'sha256\' in expected else (\'size\', \'mtime_ns\')\n        if any(actual.get(k) != expected.get(k) for k in keys):\n            raise ValueError(f\'Output identity changed: {path}\')\n    return result\n\n\ndef _relocate(value, source, destination):\n    if isinstance(value, dict):\n        return {_relocate(k, source, destination): _relocate(v, source, destination) for k, v in value.items()}\n    if isinstance(value, list): return [_relocate(v, source, destination) for v in value]\n    if isinstance(value, str) and (value == source or value.startswith(source + \'/\')):\n        return destination + value[len(source):]\n    return value\n\n\n@dataclass(frozen=True)\nclass InputPath:\n    path: str | Path\n    kind: str = "file"\n\n\ndef check_output_paths(output, inputs):\n    """Call before mkdir, lock files, running manifests, staging or rotation."""\n    out = Path(output).expanduser().resolve()\n    if out == Path(out.anchor) or out.is_file():\n        raise ValueError(f"Invalid output directory: {out}")\n    for item in inputs:\n        p = Path(item.path).expanduser().resolve()\n        if out == p or out in p.parents or p in out.parents:\n            raise ValueError(f"Output overlaps protected {item.kind}: {out} <-> {p}")\n    return out\n\n\ndef declared_inputs(raw, *, input_dirs=(), input_files=(), input_prefixes=(), resolver=None):\n    result = []\n    for kind, keys in (("directory", input_dirs), ("file", input_files)):\n        for key in keys:\n            if raw.get(key): result.append(InputPath(raw[key], kind))\n    for key in input_prefixes:\n        if raw.get(key):\n            prefix = Path(raw[key]).expanduser()\n            result += [InputPath(str(prefix) + suffix, "file") for suffix in (".pgen", ".pvar", ".pvar.zst", ".psam")]\n    # Preserve file tracking and protect directory-valued inputs before writes.\n    excluded = {\'out\', \'out_dir\', \'dir_out\', \'analysis_root\', \'dir_analysis_root\'}\n    for key, value in raw.items():\n        if key in excluded or not isinstance(value, (str, Path)) or not value:\n            continue\n        path = Path(value).expanduser()\n        if path.is_file(): result.append(InputPath(path, \'file\'))\n        elif path.is_dir(): result.append(InputPath(path, \'directory\'))\n    # Shards and annotation roots can live outside the scan root after relocation.\n    modern = raw.get(\'modern_dir\') or raw.get(\'modern\')\n    base = Path(modern).expanduser().resolve() if modern else None\n    if base is not None and (base / \'run.json\').is_file():\n        from scan_store import load_run_metadata\n        run = load_run_metadata(base)\n        result.append(InputPath(base / \'run.json\', \'file\'))\n        for entry in run.get(\'shards\', []):\n            shard = Path(entry).expanduser()\n            if not shard.is_absolute(): shard = base / shard\n            result.append(InputPath(shard, \'directory\'))\n            for filename in (\'COMPLETE.json\', \'scan.json\', \'multiallelic_sites.tsv.gz\'):\n                if (shard / filename).is_file(): result.append(InputPath(shard / filename, \'file\'))\n        options = run.get(\'options\', run.get(\'identity\', {}).get(\'options\', {}))\n        for key in (\'architecture_annotation_cache\', \'ld_block_dir\'):\n            if options.get(key): result.append(InputPath(options[key], \'directory\'))\n        for key in (\'architecture_lcr_bed\', \'architecture_segdup_bed\', \'architecture_low_mappability_bed\'):\n            if options.get(key): result.append(InputPath(options[key], \'file\'))\n    if resolver is not None:\n        result.extend(resolver(raw))\n    dedup = {}\n    for item in result:\n        p = Path(item.path).expanduser().resolve()\n        dedup[(p, item.kind)] = InputPath(p, item.kind)\n    return tuple(dedup.values())\n\n\ndef transactional(manifest, out_arg="out", *, input_dirs=(), input_files=(),\n                  input_prefixes=("pgen_prefix",), resolve_inputs=None):\n    """Publish derived outputs after read-only input/output path validation."""\n    def decorate(function):\n        signature = inspect.signature(function)\n        @functools.wraps(function)\n        def wrapped(*args, **kwargs):\n            bound = signature.bind(*args, **kwargs); bound.apply_defaults()\n            namespace = out_arg not in bound.arguments\n            original = bound.arguments["args"] if namespace else None\n            raw = dict(vars(original)) if namespace else dict(bound.arguments)\n            output_value = raw[out_arg]\n            protected = declared_inputs(raw, input_dirs=input_dirs, input_files=input_files,\n                                        input_prefixes=input_prefixes, resolver=resolve_inputs)\n            output = check_output_paths(output_value, protected)\n            # No filesystem mutation has occurred above this line.\n            output.parent.mkdir(parents=True, exist_ok=True)\n            with (output.parent / ("." + output.name + ".lock")).open("a") as lock:\n                try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)\n                except BlockingIOError as exc: raise RuntimeError(f"Output in use: {output}") from exc\n                # Re-evaluate aliases/dynamic source paths while holding the lock.\n                protected = declared_inputs(raw, input_dirs=input_dirs, input_files=input_files,\n                                            input_prefixes=input_prefixes, resolver=resolve_inputs)\n                check_output_paths(output, protected)\n                input_files_now = {str(Path(p.path)): identity(p.path) for p in protected if Path(p.path).is_file()}\n                run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + uuid.uuid4().hex[:8]\n                history = output.parent / ("." + output.name + ".history")\n                history.mkdir(exist_ok=True)\n                staging = Path(tempfile.mkdtemp(prefix=run_id + "-", dir=history))\n                old_manifest = output / manifest\n                if old_manifest.is_file():\n                    (staging / "previous_manifest.receipt").write_bytes(old_manifest.read_bytes())\n                request = {k: str(v) for k, v in raw.items()}\n                started = {"status": "running", "run_id": run_id, "request": request,\n                           "inputs": input_files_now, "protected_inputs": [dict(path=str(p.path), kind=p.kind) for p in protected],\n                           "code": code_identity(), "validation_policy": "declared input paths before writes; per-run identities"}\n                output.mkdir(parents=True, exist_ok=True)\n                atomic_json(old_manifest, started)\n                try:\n                    if namespace:\n                        cloned = copy.copy(original); setattr(cloned, out_arg, str(staging)); bound.arguments["args"] = cloned\n                    else:\n                        bound.arguments[out_arg] = staging\n                    result = function(*bound.args, **bound.kwargs)\n                    if not isinstance(result, dict): result = json.loads((staging / manifest).read_text())\n                    if result.get("status") not in TERMINAL:\n                        raise ValueError(f"Stage did not finish: {result.get(\'status\', \'missing\')}")\n                    if any(identity(path) != saved for path, saved in input_files_now.items()):\n                        raise ValueError("Scientific inputs changed during analysis")\n                    for path in staging.rglob("*.json"):\n                        try: value = json.loads(path.read_text())\n                        except (ValueError, UnicodeError): continue\n                        atomic_json(path, _relocate(value, str(staging), str(output)))\n                    result = _relocate(result, str(staging), str(output))\n                    result.update(run_id=run_id, request=request, inputs=input_files_now,\n                                  protected_inputs=started["protected_inputs"], code=started["code"],\n                                  output_inventory=inventory(staging, manifest))\n                    atomic_json(staging / manifest, result)\n                    archive = history / (run_id + "-previous")\n                    os.replace(output, archive)\n                    receipt = staging / "previous_manifest.receipt"\n                    if receipt.is_file(): (archive / manifest).write_bytes(receipt.read_bytes())\n                    os.replace(staging, output)\n                    return result\n                except BaseException as exc:\n                    failure = dict(started, status="failed", error=str(exc), diagnostics=str(staging))\n                    atomic_json(output / manifest, failure)\n                    atomic_json(staging / "failure.json", failure)\n                    raise\n        return wrapped\n    return decorate\n\n', 'allele3/f/alphagenome_partition_stdin.cpp': '// AlphaGenome Atlas: bounded-memory TSV partitioning and allele-aware VCF merge.\n// Build: c++ -O3 -std=c++17 -DNDEBUG alphagenome_merge.cpp -lz -o alphagenome_merge\n// Ordinary gzip inputs are decompressed once per successful partition phase.\n// All statistics below describe PREDICTION COVERAGE, never population genotypes.\n#include <algorithm>\n#include <array>\n#include <cerrno>\n#include <charconv>\n#include <chrono>\n#include <cmath>\n#include <cstdint>\n#include <cstdlib>\n#include <cstring>\n#include <filesystem>\n#include <fstream>\n#include <iomanip>\n#include <iostream>\n#include <limits>\n#include <map>\n#include <memory>\n#include <optional>\n#include <set>\n#include <sstream>\n#include <stdexcept>\n#include <string>\n#include <string_view>\n#include <vector>\n#include <zlib.h>\n#include <unistd.h>\n\nnamespace fs = std::filesystem;\nusing U64 = std::uint64_t;\nconstexpr std::size_t BUFFER_SIZE = 1U << 20;\nconstexpr std::size_t MAX_LINE_SIZE = 1U << 20;\nconstexpr const char* PARTITION_SCHEMA = "alphagenome-tsv-partition-1";\nconstexpr const char* MERGE_SCHEMA = "alphagenome-vcf-merge-1";\nconst std::array<std::string, 3> KINDS = {"avi", "splicing", "features"};\nconst std::array<std::vector<std::string>, 3> SCORE_COLUMNS = {{\n    {"raw_score", "PHRED"},\n    {"alphagenome_splicing"},\n    {"MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",\n     "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",\n     "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",\n     "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY",\n     "IS_INSERTION", "IS_DELETION"}\n}};\nconst std::array<std::size_t, 3> FIELD_OFFSETS = {0, 2, 3};\nconst std::array<U64, 25> CONTIG_LENGTHS = {\n    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,\n    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,\n    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,\n    58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569\n};\n\n[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }\n\nstd::string json_quote(std::string_view value) {\n    std::string out = "\\"";\n    static constexpr char HEX[] = "0123456789abcdef";\n    for (unsigned char c : value) {\n        switch (c) {\n            case \'"\': out += "\\\\\\""; break;\n            case \'\\\\\': out += "\\\\\\\\"; break;\n            case \'\\n\': out += "\\\\n"; break;\n            case \'\\r\': out += "\\\\r"; break;\n            case \'\\t\': out += "\\\\t"; break;\n            default:\n                if (c < 32) { out += "\\\\u00"; out += HEX[c >> 4]; out += HEX[c & 15]; }\n                else out += static_cast<char>(c);\n        }\n    }\n    out += \'"\';\n    return out;\n}\n\nU64 parse_uint(std::string_view text, const std::string& name) {\n    U64 result = 0;\n    auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);\n    if (text.empty() || parsed.ec != std::errc() || parsed.ptr != text.data() + text.size())\n        fail("Invalid unsigned integer for " + name + ": " + std::string(text));\n    return result;\n}\n\n// A strict parser for the flat scalar JSON receipts written by this executable.\n// No external JSON dependency is needed on the analysis host.\nstd::map<std::string, std::string> read_receipt(const fs::path& path) {\n    std::ifstream input(path, std::ios::binary);\n    if (!input) fail("Cannot read receipt: " + path.string());\n    std::string data((std::istreambuf_iterator<char>(input)), {});\n    if (data.size() > 65536) fail("Oversized receipt: " + path.string());\n    std::size_t at = 0;\n    auto whitespace = [&]() { while (at < data.size() && std::isspace(static_cast<unsigned char>(data[at]))) ++at; };\n    auto expect = [&](char c) {\n        whitespace();\n        if (at >= data.size() || data[at++] != c) fail("Malformed receipt: " + path.string());\n    };\n    auto string_value = [&]() {\n        expect(\'"\');\n        std::string value;\n        while (at < data.size()) {\n            unsigned char c = data[at++];\n            if (c == \'"\') return value;\n            if (c < 32) fail("Invalid control character in receipt: " + path.string());\n            if (c != \'\\\\\') { value += static_cast<char>(c); continue; }\n            if (at == data.size()) fail("Incomplete JSON escape in " + path.string());\n            switch (data[at++]) {\n                case \'"\': value += \'"\'; break;\n                case \'\\\\\': value += \'\\\\\'; break;\n                case \'/\': value += \'/\'; break;\n                case \'n\': value += \'\\n\'; break;\n                case \'r\': value += \'\\r\'; break;\n                case \'t\': value += \'\\t\'; break;\n                case \'b\': value += \'\\b\'; break;\n                case \'f\': value += \'\\f\'; break;\n                case \'u\': {\n                    if (at + 4 > data.size() || data[at] != \'0\' || data[at + 1] != \'0\')\n                        fail("Unsupported JSON escape in " + path.string());\n                    auto nibble = [&](char h) -> unsigned {\n                        if (h >= \'0\' && h <= \'9\') return h - \'0\';\n                        if (h >= \'a\' && h <= \'f\') return h - \'a\' + 10;\n                        if (h >= \'A\' && h <= \'F\') return h - \'A\' + 10;\n                        fail("Invalid JSON escape in " + path.string());\n                    };\n                    value += static_cast<char>((nibble(data[at + 2]) << 4) | nibble(data[at + 3]));\n                    at += 4; break;\n                }\n                default: fail("Invalid JSON escape in " + path.string());\n            }\n        }\n        fail("Unterminated JSON string in " + path.string());\n    };\n    std::map<std::string, std::string> values;\n    expect(\'{\'); whitespace();\n    if (at < data.size() && data[at] == \'}\') { ++at; }\n    else {\n        for (;;) {\n            const auto key = string_value(); expect(\':\'); whitespace();\n            std::string value;\n            if (at < data.size() && data[at] == \'"\') value = string_value();\n            else {\n                const auto begin = at;\n                while (at < data.size() && data[at] >= \'0\' && data[at] <= \'9\') ++at;\n                if (at == begin) fail("Invalid scalar in receipt: " + path.string());\n                value = data.substr(begin, at - begin);\n            }\n            if (!values.emplace(key, value).second) fail("Duplicate receipt key: " + key);\n            whitespace();\n            if (at < data.size() && data[at] == \'}\') { ++at; break; }\n            expect(\',\');\n        }\n    }\n    whitespace();\n    if (at != data.size()) fail("Trailing JSON data in " + path.string());\n    return values;\n}\n\nvoid atomic_text(const fs::path& path, const std::string& content) {\n    if (!path.parent_path().empty()) fs::create_directories(path.parent_path());\n    fs::path temporary = path; temporary += ".partial";\n    {\n        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);\n        if (!out) fail("Cannot write " + temporary.string());\n        out.write(content.data(), static_cast<std::streamsize>(content.size()));\n        out.flush();\n        if (!out) fail("Failed writing " + temporary.string());\n        out.close();\n        if (!out) fail("Failed closing " + temporary.string());\n    }\n    fs::rename(temporary, path);\n}\n\nstd::pair<std::string, int> chromosome(std::string_view text) {\n    if (text.substr(0, 3) == "chr") text.remove_prefix(3);\n    if (text == "X") return {"chrX", 22};\n    if (text == "Y") return {"chrY", 23};\n    if (text == "M" || text == "MT") return {"chrM", 24};\n    if (text.empty() || text.size() > 2 || text[0] == \'0\')\n        fail("Unsupported chromosome: " + std::string(text));\n    const auto n = parse_uint(text, "chromosome");\n    if (n < 1 || n > 22) fail("Unsupported chromosome: " + std::string(text));\n    return {"chr" + std::to_string(n), static_cast<int>(n - 1)};\n}\n\nint base_index(char base) {\n    switch (base) { case \'A\': return 0; case \'C\': return 1; case \'G\': return 2; case \'T\': return 3; }\n    return -1;\n}\n\nclass GzipReader {\n    gzFile file_ = nullptr;\n    std::vector<char> buffer_ = std::vector<char>(BUFFER_SIZE);\n    std::size_t start_ = 0, end_ = 0;\n    bool eof_ = false;\n    std::string path_;\npublic:\n    explicit GzipReader(const fs::path& path) : path_(path.string()) {\n        if (path == "-") {\n            const int fd = dup(STDIN_FILENO);\n            if (fd < 0) fail("Cannot duplicate indexed input stream");\n            file_ = gzdopen(fd, "rb");\n            if (!file_) { close(fd); fail("Cannot open indexed input stream"); }\n            if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0)\n                fail("Cannot allocate indexed input buffer");\n            return;\n        }\n        std::ifstream magic(path, std::ios::binary);\n        unsigned char prefix[2] = {0, 0};\n        magic.read(reinterpret_cast<char*>(prefix), 2);\n        if (magic.gcount() != 2 || prefix[0] != 0x1f || prefix[1] != 0x8b)\n            fail("Expected gzip input: " + path_);\n        file_ = gzopen(path_.c_str(), "rb");\n        if (!file_) fail("Cannot open gzip input: " + path_);\n        if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0) {\n            gzclose(file_); file_ = nullptr; fail("Cannot allocate gzip buffer for " + path_);\n        }\n    }\n    GzipReader(const GzipReader&) = delete;\n    GzipReader& operator=(const GzipReader&) = delete;\n    ~GzipReader() { if (file_) gzclose(file_); }\n    bool next(std::string& line) {\n        line.clear();\n        for (;;) {\n            if (start_ < end_) {\n                const auto* found = static_cast<const char*>(std::memchr(buffer_.data() + start_, \'\\n\', end_ - start_));\n                const auto stop = found ? static_cast<std::size_t>(found - buffer_.data()) : end_;\n                line.append(buffer_.data() + start_, stop - start_);\n                start_ = stop + (found ? 1 : 0);\n                if (line.size() > MAX_LINE_SIZE) fail("Oversized TSV line in " + path_);\n                if (found) { if (!line.empty() && line.back() == \'\\r\') line.pop_back(); return true; }\n            }\n            if (eof_) {\n                if (!line.empty() && line.back() == \'\\r\') line.pop_back();\n                return !line.empty();\n            }\n            const int count = gzread(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));\n            int code = Z_OK;\n            const char* error = gzerror(file_, &code);\n            if (count < 0 || (code != Z_OK && code != Z_STREAM_END))\n                fail("gzip read/CRC/truncation error in " + path_ + ": " + (error ? error : "unknown"));\n            start_ = 0;\n            end_ = count > 0 ? static_cast<std::size_t>(count) : 0;\n            if (count == 0) {\n                if (!gzeof(file_)) fail("Unexpected gzip EOF in " + path_);\n                const int close_code = gzclose(file_); file_ = nullptr;\n                if (close_code != Z_OK) fail("gzip close/truncation error in " + path_);\n                eof_ = true;\n            }\n        }\n    }\n};\n\nclass CanonicalSink {\n    gzFile file_ = nullptr;\n    std::string buffer_;\n    std::string path_;\n    bool closed_ = false;\npublic:\n    U64 uncompressed_bytes = 0;\n    uLong crc = crc32(0L, Z_NULL, 0);\n    explicit CanonicalSink(const std::optional<fs::path>& path) {\n        buffer_.reserve(BUFFER_SIZE + 4096);\n        if (path) {\n            path_ = path->string(); file_ = gzopen(path_.c_str(), "wb1");\n            if (!file_) fail("Cannot write staging gzip: " + path_);\n            if (gzbuffer(file_, static_cast<unsigned>(BUFFER_SIZE)) != 0) {\n                gzclose(file_); file_ = nullptr; fail("Cannot allocate gzip buffer for " + path_);\n            }\n        }\n    }\n    ~CanonicalSink() { if (file_) gzclose(file_); }\n    void flush() {\n        if (buffer_.empty()) return;\n        crc = crc32(crc, reinterpret_cast<const Bytef*>(buffer_.data()), static_cast<uInt>(buffer_.size()));\n        uncompressed_bytes += buffer_.size();\n        if (file_) {\n            const auto written = gzwrite(file_, buffer_.data(), static_cast<unsigned>(buffer_.size()));\n            if (written != static_cast<int>(buffer_.size())) fail("Failed writing staging gzip: " + path_);\n        }\n        buffer_.clear();\n    }\n    void append(std::string_view value) { buffer_.append(value.data(), value.size()); if (buffer_.size() >= BUFFER_SIZE) flush(); }\n    void finish() {\n        if (closed_) return;\n        flush();\n        if (file_) {\n            const int result = gzclose(file_); file_ = nullptr;\n            if (result != Z_OK) fail("Failed closing staging gzip: " + path_);\n        }\n        closed_ = true;\n    }\n};\n\nvoid split_fields(std::string_view line, std::vector<std::string_view>& fields) {\n    fields.clear();\n    std::size_t at = 0;\n    while (at < line.size()) {\n        while (at < line.size() && (line[at] == \' \' || line[at] == \'\\t\')) ++at;\n        if (at == line.size()) break;\n        const auto begin = at;\n        while (at < line.size() && line[at] != \' \' && line[at] != \'\\t\') ++at;\n        fields.emplace_back(line.data() + begin, at - begin);\n    }\n}\n\nstd::string header_for(int kind) {\n    std::string line = "#CHROM\\tPOS\\tREF\\tALT";\n    for (const auto& column : SCORE_COLUMNS[kind]) { line += \'\\t\'; line += column; }\n    line += \'\\n\'; return line;\n}\n\nvoid validate_header(std::string_view line, int kind, const std::string& context) {\n    std::vector<std::string_view> fields;\n    split_fields(line, fields);\n    if (fields.size() != SCORE_COLUMNS[kind].size() + 4 || fields[0] != "#CHROM" ||\n        fields[1] != "POS" || fields[2] != "REF" || fields[3] != "ALT")\n        fail("Unexpected " + KINDS[kind] + " header in " + context + "; expected " + header_for(kind));\n    for (std::size_t j = 0; j < SCORE_COLUMNS[kind].size(); ++j)\n        if (fields[j + 4] != SCORE_COLUMNS[kind][j])\n            fail("Unexpected " + KINDS[kind] + " score column in " + context + ": " + std::string(fields[j + 4]));\n}\n\nstruct Numeric { bool missing = false; bool zero = false; };\n\nNumeric validate_number(std::string_view token, const std::string& context) {\n    if (token == "." || token == "NA" || token == "NaN" || token == "nan") return {true, false};\n    std::size_t at = 0, digits = 0;\n    bool nonzero = false;\n    if (at < token.size() && (token[at] == \'+\' || token[at] == \'-\')) ++at;\n    while (at < token.size() && token[at] >= \'0\' && token[at] <= \'9\') {\n        nonzero |= token[at] != \'0\'; ++at; ++digits;\n    }\n    if (at < token.size() && token[at] == \'.\') {\n        ++at;\n        while (at < token.size() && token[at] >= \'0\' && token[at] <= \'9\') {\n            nonzero |= token[at] != \'0\'; ++at; ++digits;\n        }\n    }\n    if (digits == 0) fail("Invalid numeric token " + json_quote(token) + " in " + context);\n    if (at < token.size() && (token[at] == \'e\' || token[at] == \'E\')) {\n        ++at;\n        if (at < token.size() && (token[at] == \'+\' || token[at] == \'-\')) ++at;\n        const auto begin = at;\n        while (at < token.size() && token[at] >= \'0\' && token[at] <= \'9\') ++at;\n        if (at == begin) fail("Invalid exponent in " + context + ": " + std::string(token));\n    }\n    if (at != token.size()) fail("Invalid numeric token " + json_quote(token) + " in " + context);\n    // GCC 8\'s libstdc++ implements integer from_chars only. Keep the fast\n    // floating parser where available; the portable branch preserves decimal\n    // text and underflow semantics and never round-trips output through double.\n#if defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE >= 11 && !defined(ALPHAGENOME_PORTABLE_FLOAT)\n    auto number = token;\n    if (!number.empty() && number.front() == \'+\') number.remove_prefix(1);\n    double value = 0;\n    const auto result = std::from_chars(number.data(), number.data() + number.size(), value, std::chars_format::general);\n    if (result.ec == std::errc::result_out_of_range) {\n        // Very small finite decimal text may underflow binary64. Preserve its text\n        // and do not count it as numeric zero; only overflow is an error.\n        std::string copy(token);\n        char* end = nullptr;\n        value = std::strtod(copy.c_str(), &end);\n        if (end != copy.c_str() + copy.size() || !std::isfinite(value))\n            fail("Non-finite/overflow numeric token " + json_quote(token) + " in " + context);\n    } else if (result.ec != std::errc() || result.ptr != number.data() + number.size() || !std::isfinite(value))\n        fail("Non-finite numeric token " + json_quote(token) + " in " + context);\n#else\n    std::string copy(token);\n    char* end = nullptr;\n    const double value = std::strtod(copy.c_str(), &end);\n    if (end != copy.c_str() + copy.size() || !std::isfinite(value))\n        fail("Non-finite/overflow numeric token " + json_quote(token) + " in " + context);\n#endif\n    return {false, !nonzero};\n}\n\nstruct ParsedRow {\n    std::string chrom;\n    int chrom_index = 0;\n    U64 pos = 0;\n    char ref = \'N\', alt = \'N\';\n    int alt_index = -1;\n};\n\nParsedRow parse_row(const std::vector<std::string_view>& fields, int kind, const std::string& context,\n                    std::vector<Numeric>* numbers = nullptr) {\n    if (fields.size() != SCORE_COLUMNS[kind].size() + 4)\n        fail("Wrong field count in " + context + ": expected " + std::to_string(SCORE_COLUMNS[kind].size() + 4) +\n             ", got " + std::to_string(fields.size()));\n    ParsedRow row;\n    auto chr = chromosome(fields[0]); row.chrom = std::move(chr.first); row.chrom_index = chr.second;\n    row.pos = parse_uint(fields[1], "POS in " + context);\n    if (!row.pos || row.pos > CONTIG_LENGTHS[row.chrom_index])\n        fail("Position outside GRCh38 primary chromosome in " + context);\n    if (fields[2].size() != 1 || fields[3].size() != 1 || base_index(fields[2][0]) < 0 ||\n        base_index(fields[3][0]) < 0 || fields[2][0] == fields[3][0])\n        fail("Expected canonical A/C/G/T single-base REF != ALT in " + context);\n    row.ref = fields[2][0]; row.alt = fields[3][0]; row.alt_index = base_index(row.alt);\n    if (numbers) { numbers->clear(); numbers->reserve(SCORE_COLUMNS[kind].size()); }\n    for (std::size_t j = 4; j < fields.size(); ++j) {\n        auto numeric = validate_number(fields[j], context);\n        if (numbers) numbers->push_back(numeric);\n    }\n    return row;\n}\n\nstruct PositionOrder {\n    U64 pos = 0;\n    char ref = \'N\';\n    unsigned alt_mask = 0;\n    void accept(const ParsedRow& row, const std::string& context) {\n        if (row.pos < pos) fail("Decreasing position in " + context);\n        if (row.pos != pos) { pos = row.pos; ref = row.ref; alt_mask = 0; }\n        else if (row.ref != ref) fail("Inconsistent REF at one position in " + context);\n        const auto bit = 1U << row.alt_index;\n        if (alt_mask & bit) fail("Duplicate CHROM/POS/REF/ALT in " + context);\n        alt_mask |= bit;\n    }\n};\n\nstruct ChromReceipt {\n    std::string chrom, kind, token;\n    U64 rows = 0, first_pos = 0, last_pos = 0, gzip_bytes = 0, uncompressed_bytes = 0, crc = 0;\n    std::string json() const {\n        std::ostringstream out;\n        out << "{\\n  \\"schema\\": " << json_quote(PARTITION_SCHEMA)\n            << ",\\n  \\"chrom\\": " << json_quote(chrom) << ",\\n  \\"kind\\": " << json_quote(kind)\n            << ",\\n  \\"token\\": " << json_quote(token) << ",\\n  \\"rows\\": " << rows\n            << ",\\n  \\"first_pos\\": " << first_pos << ",\\n  \\"last_pos\\": " << last_pos\n            << ",\\n  \\"gzip_bytes\\": " << gzip_bytes << ",\\n  \\"uncompressed_bytes\\": " << uncompressed_bytes\n            << ",\\n  \\"uncompressed_crc32\\": " << crc << "\\n}\\n";\n        return out.str();\n    }\n    static ChromReceipt load(const fs::path& path) {\n        auto values = read_receipt(path);\n        auto get = [&](const std::string& key) -> std::string {\n            auto it = values.find(key);\n            if (it == values.end()) fail("Missing receipt field " + key + " in " + path.string());\n            return it->second;\n        };\n        if (get("schema") != PARTITION_SCHEMA) fail("Unsupported receipt schema in " + path.string());\n        ChromReceipt receipt;\n        receipt.chrom = get("chrom"); receipt.kind = get("kind"); receipt.token = get("token");\n        receipt.rows = parse_uint(get("rows"), "receipt rows");\n        receipt.first_pos = parse_uint(get("first_pos"), "receipt first_pos");\n        receipt.last_pos = parse_uint(get("last_pos"), "receipt last_pos");\n        receipt.gzip_bytes = parse_uint(get("gzip_bytes"), "receipt gzip_bytes");\n        receipt.uncompressed_bytes = parse_uint(get("uncompressed_bytes"), "receipt uncompressed_bytes");\n        receipt.crc = parse_uint(get("uncompressed_crc32"), "receipt CRC");\n        return receipt;\n    }\n};\n\nbool staging_quick_check(const fs::path& path, const ChromReceipt& receipt) {\n    if (!fs::is_regular_file(path) || fs::file_size(path) != receipt.gzip_bytes || receipt.gzip_bytes < 18) return false;\n    std::ifstream input(path, std::ios::binary);\n    unsigned char bytes[8] = {};\n    input.seekg(-8, std::ios::end); input.read(reinterpret_cast<char*>(bytes), 8);\n    if (!input) return false;\n    auto little = [&](int offset) -> U64 {\n        return U64(bytes[offset]) | (U64(bytes[offset + 1]) << 8) | (U64(bytes[offset + 2]) << 16) | (U64(bytes[offset + 3]) << 24);\n    };\n    return little(0) == receipt.crc && little(4) == (receipt.uncompressed_bytes & 0xffffffffU);\n}\n\nstruct Arguments {\n    std::string mode;\n    std::map<std::string, std::string> values;\n    explicit Arguments(int argc, char** argv) {\n        if (argc < 2) fail("Specify partition or merge; use --help for usage");\n        mode = argv[1];\n        if (mode == "--help" || mode == "-h") return;\n        if (mode != "partition" && mode != "merge") fail("Unknown mode: " + mode);\n        const std::set<std::string> allowed = mode == "partition"\n            ? std::set<std::string>{"--kind", "--input", "--out-dir", "--token", "--progress-every", "--file-prefix"}\n            : std::set<std::string>{"--chrom", "--avi", "--splicing", "--features", "--stats", "--token", "--progress-every"};\n        for (int i = 2; i < argc; i += 2) {\n            if (i + 1 == argc) fail("Missing value for " + std::string(argv[i]));\n            if (!allowed.count(argv[i])) fail("Unknown option: " + std::string(argv[i]));\n            if (!values.emplace(argv[i], argv[i + 1]).second) fail("Duplicate option: " + std::string(argv[i]));\n        }\n    }\n    std::string required(const std::string& name) const {\n        auto found = values.find(name);\n        if (found == values.end() || found->second.empty()) fail("Required option: " + name);\n        return found->second;\n    }\n    U64 progress() const {\n        const auto found = values.find("--progress-every");\n        return found == values.end() ? 10000000 : parse_uint(found->second, "--progress-every");\n    }\n};\n\nint partition(const Arguments& args) {\n    const auto kind_name = args.required("--kind");\n    auto found = std::find(KINDS.begin(), KINDS.end(), kind_name);\n    if (found == KINDS.end()) fail("--kind must be avi, splicing, or features");\n    const int kind = static_cast<int>(found - KINDS.begin());\n    const fs::path input = args.required("--input"), out_dir = args.required("--out-dir");\n    const auto prefix = args.values.count("--file-prefix") ? args.values.at("--file-prefix") : std::string{};\n    if (!prefix.empty() && prefix != kind_name + ".") fail("--file-prefix must match kind plus \'.\'");\n    const auto token = args.required("--token");\n    const U64 progress = args.progress();\n    fs::create_directories(out_dir);\n    GzipReader reader(input);\n    std::string line;\n    if (!reader.next(line)) fail("Empty input without header: " + input.string());\n    validate_header(line, kind, input.string());\n    const auto header = header_for(kind);\n    std::vector<std::string_view> fields; fields.reserve(22);\n    std::string normalized; normalized.reserve(2048);\n    std::set<std::string> visited;\n    std::vector<ChromReceipt> completed;\n    std::unique_ptr<CanonicalSink> sink;\n    std::optional<ChromReceipt> previous;\n    ChromReceipt receipt;\n    PositionOrder order;\n    fs::path destination, temporary, done_path;\n    U64 source_rows = 0, line_number = 1, reused_chromosomes = 0;\n    auto finish_chromosome = [&]() {\n        if (!sink) return;\n        sink->finish();\n        receipt.uncompressed_bytes = sink->uncompressed_bytes;\n        receipt.crc = sink->crc;\n        if (previous) {\n            if (receipt.rows != previous->rows || receipt.first_pos != previous->first_pos ||\n                receipt.last_pos != previous->last_pos || receipt.uncompressed_bytes != previous->uncompressed_bytes ||\n                receipt.crc != previous->crc)\n                fail("Source content differs from completed chromosome despite matching token: " + receipt.chrom);\n            receipt.gzip_bytes = previous->gzip_bytes; ++reused_chromosomes;\n            std::cerr << "partition " << kind_name << \' \' << receipt.chrom << " reused; rows=" << receipt.rows << \'\\n\';\n        } else {\n            receipt.gzip_bytes = fs::file_size(temporary);\n            fs::rename(temporary, destination);\n            atomic_text(done_path, receipt.json());\n            std::cerr << "partition " << kind_name << \' \' << receipt.chrom << " completed; rows=" << receipt.rows\n                      << " gzip_bytes=" << receipt.gzip_bytes << \'\\n\';\n        }\n        completed.push_back(receipt); sink.reset(); previous.reset();\n    };\n    while (reader.next(line)) {\n        ++line_number;\n        const std::string context = input.string() + ":" + std::to_string(line_number);\n        split_fields(line, fields);\n        const auto row = parse_row(fields, kind, context);\n        if (!sink || row.chrom != receipt.chrom) {\n            finish_chromosome();\n            if (!visited.insert(row.chrom).second) fail("Chromosome block revisited in " + context + ": " + row.chrom);\n            receipt = ChromReceipt{}; receipt.chrom = row.chrom; receipt.kind = kind_name; receipt.token = token;\n            order = PositionOrder{};\n            destination = out_dir / (prefix + row.chrom + ".tsv.gz"); temporary = destination; temporary += ".partial";\n            done_path = out_dir / (prefix + row.chrom + ".done.json");\n            if (fs::exists(done_path)) {\n                const auto existing = ChromReceipt::load(done_path);\n                if (existing.token != token || existing.kind != kind_name || existing.chrom != row.chrom)\n                    fail("Incompatible staging receipt; use a distinct/clean staging directory: " + done_path.string());\n                if (staging_quick_check(destination, existing)) previous = existing;\n                else std::cerr << "partition " << kind_name << \' \' << row.chrom << " regenerating missing/damaged staging output\\n";\n            }\n            sink = std::make_unique<CanonicalSink>(previous ? std::nullopt : std::optional<fs::path>(temporary));\n            sink->append(header);\n        }\n        order.accept(row, context);\n        if (receipt.rows++ == 0) receipt.first_pos = row.pos;\n        receipt.last_pos = row.pos;\n        normalized.clear(); normalized += row.chrom;\n        for (std::size_t j = 1; j < fields.size(); ++j) { normalized += \'\\t\'; normalized.append(fields[j]); }\n        normalized += \'\\n\'; sink->append(normalized);\n        ++source_rows;\n        if (progress && source_rows % progress == 0)\n            std::cerr << "partition " << kind_name << " rows=" << source_rows << " at=" << row.chrom << \':\' << row.pos << \'\\n\';\n    }\n    finish_chromosome();\n    // Old source.done is never sufficient to reuse outputs: the orchestrator also\n    // checks source fingerprints and every chromosome receipt. An interrupted\n    // ordinary-gzip source can require rereading its prefix on the next attempt.\n    std::ostringstream summary;\n    summary << "{\\n  \\"schema\\": " << json_quote(PARTITION_SCHEMA) << ",\\n  \\"kind\\": " << json_quote(kind_name)\n            << ",\\n  \\"token\\": " << json_quote(token) << ",\\n  \\"source_rows\\": " << source_rows\n            << ",\\n  \\"reused_chromosomes\\": " << reused_chromosomes << ",\\n  \\"chromosomes\\": [";\n    for (std::size_t j = 0; j < completed.size(); ++j) {\n        if (j) summary << \',\';\n        summary << \'\\n\' << completed[j].json();\n    }\n    summary << "\\n  ]\\n}\\n";\n    atomic_text(out_dir / (prefix + "source.done.json"), summary.str());\n    std::cerr << "partition " << kind_name << " finished; source_rows=" << source_rows << " chromosomes=" << completed.size() << \'\\n\';\n    return 0;\n}\n\nstruct ScoreRow {\n    ParsedRow locus;\n    std::vector<std::string> values;\n    std::vector<Numeric> numbers;\n};\n\nstruct Group {\n    U64 pos = 0;\n    char ref = \'N\';\n    std::array<std::optional<ScoreRow>, 4> alleles;\n};\n\nclass SourceStream {\n    int kind_;\n    std::string expected_chrom_, path_, line_;\n    GzipReader reader_;\n    PositionOrder order_;\n    std::vector<std::string_view> fields_;\n    U64 line_number_ = 1;\npublic:\n    U64 rows = 0;\n    std::optional<ScoreRow> current;\n    SourceStream(int kind, const std::string& path, const std::string& chrom)\n        : kind_(kind), expected_chrom_(chrom), path_(path), reader_(path) {\n        fields_.reserve(22);\n        if (!reader_.next(line_)) fail("Empty input without header: " + path_);\n        validate_header(line_, kind_, path_); advance();\n    }\n    void advance() {\n        if (!reader_.next(line_)) { current.reset(); return; }\n        ++line_number_;\n        const std::string context = path_ + ":" + std::to_string(line_number_);\n        split_fields(line_, fields_);\n        ScoreRow row;\n        row.locus = parse_row(fields_, kind_, context, &row.numbers);\n        if (row.locus.chrom != expected_chrom_) fail("Unexpected chromosome in " + context + ": " + row.locus.chrom);\n        order_.accept(row.locus, context);\n        row.values.reserve(SCORE_COLUMNS[kind_].size());\n        for (std::size_t j = 4; j < fields_.size(); ++j) row.values.emplace_back(fields_[j]);\n        current = std::move(row); ++rows;\n    }\n    Group take_group() {\n        if (!current) fail("Internal error: group requested after EOF");\n        Group group; group.pos = current->locus.pos; group.ref = current->locus.ref;\n        do {\n            const int alt = current->locus.alt_index;\n            group.alleles[alt] = std::move(*current);\n            advance();\n        } while (current && current->locus.pos == group.pos);\n        return group;\n    }\n};\n\nstruct FieldStats {\n    std::string name;\n    int source = 0;\n    U64 numeric_present = 0, numeric_missing = 0, numeric_zero = 0, absent_source_alleles = 0;\n    U64 sites_with_three_alt_scores = 0;\n};\n\nstruct SourceStats {\n    U64 covered_alleles = 0, covered_sites = 0, sites_with_three_alt_rows = 0;\n    U64 sites_with_three_alt_complete_scores = 0;\n    U64 numeric_present = 0, numeric_missing = 0, numeric_zero = 0;\n};\n\nstd::array<FieldStats, 21> make_field_stats() {\n    std::array<FieldStats, 21> fields;\n    fields[0].name = "AG_AVI_RAW"; fields[0].source = 0;\n    fields[1].name = "AG_AVI_PHRED"; fields[1].source = 0;\n    fields[2].name = "AG_SPLICE"; fields[2].source = 1;\n    for (std::size_t j = 0; j < SCORE_COLUMNS[2].size(); ++j) {\n        fields[j + 3].name = "AG_SHAP_" + SCORE_COLUMNS[2][j]; fields[j + 3].source = 2;\n    }\n    return fields;\n}\n\nvoid vcf_header(const std::string& chrom, const std::string& token, const std::array<FieldStats, 21>& fields) {\n    const auto chr = chromosome(chrom);\n    std::cout << "##fileformat=VCFv4.3\\n##source=allele3_alphagenome_merge_1\\n##reference=GRCh38\\n"\n              << "##AG_schema=" << MERGE_SCHEMA << "\\n##AG_recipe_token=" << json_quote(token) << \'\\n\'\n              << "##AG_allele_semantics=\\"ALT contains published predicted substitutions; presence is not evidence of population observation\\"\\n"\n              << "##AG_missing_semantics=\\"Dot means unavailable; AG_SOURCE_MASK distinguishes absent rows from explicit missing scores\\"\\n"\n              << "##AG_feature_semantics=\\"AG_SHAP fields are published signed feature contributions, not raw molecular-effect measurements\\"\\n"\n              << "##contig=<ID=" << chrom << ",length=" << CONTIG_LENGTHS[chr.second] << ",assembly=GRCh38>\\n";\n    for (std::size_t j = 0; j < fields.size(); ++j) {\n        std::string description;\n        if (j == 0) description = "Published AlphaGenome variant impact raw score; not a genotype statistic";\n        else if (j == 1) description = "Published AlphaGenome variant impact PHRED score; not genotype quality";\n        else if (j == 2) description = "Published combined AlphaGenome splicing score; distinct from SHAP MERGED_SPLICING";\n        else description = "Published signed AVI feature-importance (SHAP) value for " + SCORE_COLUMNS[2][j - 3] + "; not the raw feature score";\n        std::cout << "##INFO=<ID=" << fields[j].name << ",Number=A,Type=Float,Description=" << json_quote(description) << ">\\n";\n    }\n    std::cout << "##INFO=<ID=AG_SOURCE_MASK,Number=A,Type=Integer,Description=\\"Published row presence: AVI=1, splicing=2, feature importance=4; bitwise union\\">\\n"\n              << "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n";\n}\n\nint merge(const Arguments& args) {\n    const auto chrom = chromosome(args.required("--chrom")).first;\n    const auto token = args.required("--token");\n    const fs::path stats_path = args.required("--stats");\n    const U64 progress = args.progress();\n    std::array<std::unique_ptr<SourceStream>, 3> inputs;\n    bool supplied = false;\n    for (int kind = 0; kind < 3; ++kind) {\n        const auto path = args.required("--" + KINDS[kind]);\n        if (path != "NONE") { inputs[kind] = std::make_unique<SourceStream>(kind, path, chrom); supplied = true; }\n    }\n    if (!supplied) fail("At least one source file must be supplied; all three were NONE");\n    auto fields = make_field_stats();\n    std::array<SourceStats, 3> source_stats;\n    std::array<U64, 4> sites_by_nalt{};\n    std::array<U64, 8> source_mask_alleles{};\n    U64 sites = 0, union_alt_rows = 0, first_pos = 0, last_pos = 0, next_progress = progress;\n    vcf_header(chrom, token, fields);\n    std::string output; output.reserve(4096);\n    constexpr char BASES[] = "ACGT";\n    for (;;) {\n        U64 pos = std::numeric_limits<U64>::max();\n        for (const auto& source : inputs) if (source && source->current) pos = std::min(pos, source->current->locus.pos);\n        if (pos == std::numeric_limits<U64>::max()) break;\n        std::array<std::optional<Group>, 3> groups;\n        char ref = \'N\';\n        std::array<unsigned, 4> masks{};\n        for (int kind = 0; kind < 3; ++kind) {\n            if (!inputs[kind] || !inputs[kind]->current || inputs[kind]->current->locus.pos != pos) continue;\n            groups[kind] = inputs[kind]->take_group();\n            if (ref == \'N\') ref = groups[kind]->ref;\n            else if (ref != groups[kind]->ref)\n                fail("Cross-source REF mismatch at " + chrom + ":" + std::to_string(pos));\n            U64 allele_count = 0;\n            bool complete = true;\n            for (int alt = 0; alt < 4; ++alt) {\n                if (!groups[kind]->alleles[alt]) continue;\n                ++allele_count; masks[alt] |= 1U << kind;\n                const auto& row = *groups[kind]->alleles[alt];\n                for (std::size_t j = 0; j < row.numbers.size(); ++j) {\n                    auto& field = fields[FIELD_OFFSETS[kind] + j];\n                    if (row.numbers[j].missing) { ++field.numeric_missing; ++source_stats[kind].numeric_missing; complete = false; }\n                    else {\n                        ++field.numeric_present; ++source_stats[kind].numeric_present;\n                        if (row.numbers[j].zero) { ++field.numeric_zero; ++source_stats[kind].numeric_zero; }\n                    }\n                }\n            }\n            source_stats[kind].covered_alleles += allele_count;\n            ++source_stats[kind].covered_sites;\n            if (allele_count == 3) {\n                ++source_stats[kind].sites_with_three_alt_rows;\n                if (complete) ++source_stats[kind].sites_with_three_alt_complete_scores;\n                for (std::size_t j = 0; j < SCORE_COLUMNS[kind].size(); ++j) {\n                    bool field_complete = true;\n                    for (int alt = 0; alt < 4; ++alt)\n                        if (groups[kind]->alleles[alt] && groups[kind]->alleles[alt]->numbers[j].missing) field_complete = false;\n                    if (field_complete) ++fields[FIELD_OFFSETS[kind] + j].sites_with_three_alt_scores;\n                }\n            }\n        }\n        unsigned nalt = 0;\n        for (unsigned mask : masks) if (mask) { ++nalt; ++source_mask_alleles[mask]; }\n        if (nalt < 1 || nalt > 3) fail("Invalid ALT union at " + chrom + ":" + std::to_string(pos));\n        if (++sites == 1) first_pos = pos;\n        last_pos = pos;\n        ++sites_by_nalt[nalt]; union_alt_rows += nalt;\n        output.clear(); output += chrom; output += \'\\t\'; output += std::to_string(pos);\n        output += "\\t.\\t"; output += ref; output += \'\\t\';\n        bool first = true;\n        for (int alt = 0; alt < 4; ++alt) if (masks[alt]) {\n            if (!first) output += \',\';\n            output += BASES[alt]; first = false;\n        }\n        output += "\\t.\\t.\\t";\n        for (std::size_t field_index = 0; field_index < fields.size(); ++field_index) {\n            auto& field = fields[field_index];\n            if (field_index) output += \';\';\n            output += field.name; output += \'=\'; first = true;\n            const auto offset = field_index - FIELD_OFFSETS[field.source];\n            for (int alt = 0; alt < 4; ++alt) if (masks[alt]) {\n                if (!first) output += \',\';\n                first = false;\n                if (!groups[field.source] || !groups[field.source]->alleles[alt]) {\n                    output += \'.\'; ++field.absent_source_alleles;\n                } else {\n                    const auto& row = *groups[field.source]->alleles[alt];\n                    if (row.numbers[offset].missing) output += \'.\';\n                    else output += row.values[offset];\n                }\n            }\n        }\n        output += ";AG_SOURCE_MASK="; first = true;\n        for (unsigned mask : masks) if (mask) {\n            if (!first) output += \',\';\n            output += static_cast<char>(\'0\' + mask); first = false;\n        }\n        output += \'\\n\';\n        std::cout.write(output.data(), static_cast<std::streamsize>(output.size()));\n        if (!std::cout) fail("VCF stdout write failed (downstream compressor may have failed)");\n        if (progress && union_alt_rows >= next_progress) {\n            std::cerr << "merge " << chrom << " sites=" << sites << " union_alt_rows=" << union_alt_rows << " pos=" << pos << \'\\n\';\n            next_progress = union_alt_rows + progress;\n        }\n    }\n    std::cout.flush();\n    if (!std::cout) fail("VCF stdout flush failed");\n    std::ostringstream stats;\n    stats << "{\\n  \\"schema\\": " << json_quote(MERGE_SCHEMA) << ",\\n  \\"chrom\\": " << json_quote(chrom)\n          << ",\\n  \\"token\\": " << json_quote(token) << ",\\n  \\"assembly\\": \\"GRCh38\\",\\n  \\"source_rows\\": {";\n    for (int kind = 0; kind < 3; ++kind) { if (kind) stats << \',\'; stats << json_quote(KINDS[kind]) << \':\' << (inputs[kind] ? inputs[kind]->rows : 0); }\n    stats << "},\\n  \\"union_alt_rows\\": " << union_alt_rows << ",\\n  \\"sites\\": " << sites\n          << ",\\n  \\"first_pos\\": " << first_pos << ",\\n  \\"last_pos\\": " << last_pos\n          << ",\\n  \\"sites_by_nalt\\": {\\"1\\":" << sites_by_nalt[1] << ",\\"2\\":" << sites_by_nalt[2] << ",\\"3\\":" << sites_by_nalt[3] << \'}\'\n          << ",\\n  \\"prediction_base_coverage\\": {\\"exactly_2_bases\\":" << sites_by_nalt[1] << ",\\"exactly_3_bases\\":" << sites_by_nalt[2]\n          << ",\\"exactly_4_bases\\":" << sites_by_nalt[3] << ",\\"at_least_3_bases\\":" << sites_by_nalt[2] + sites_by_nalt[3] << \'}\'\n          << ",\\n  \\"source_mask_alleles\\": {";\n    for (unsigned mask = 1; mask <= 7; ++mask) { if (mask > 1) stats << \',\'; stats << json_quote(std::to_string(mask)) << \':\' << source_mask_alleles[mask]; }\n    stats << "},\\n  \\"per_source\\": {";\n    for (int kind = 0; kind < 3; ++kind) {\n        if (kind) stats << \',\';\n        const auto& s = source_stats[kind];\n        stats << \'\\n\' << json_quote(KINDS[kind]) << ":{\\"covered_alleles\\":" << s.covered_alleles << ",\\"covered_sites\\":" << s.covered_sites\n              << ",\\"sites_with_three_alt_rows\\":" << s.sites_with_three_alt_rows\n              << ",\\"sites_with_three_alt_complete_scores\\":" << s.sites_with_three_alt_complete_scores\n              << ",\\"numeric_present\\":" << s.numeric_present << ",\\"numeric_missing\\":" << s.numeric_missing << ",\\"numeric_zero\\":" << s.numeric_zero << \'}\';\n    }\n    stats << "\\n  },\\n  \\"per_field\\": {";\n    for (std::size_t j = 0; j < fields.size(); ++j) {\n        if (j) stats << \',\';\n        const auto& f = fields[j];\n        stats << \'\\n\' << json_quote(f.name) << ":{\\"source\\":" << json_quote(KINDS[f.source])\n              << ",\\"numeric_present\\":" << f.numeric_present << ",\\"numeric_missing\\":" << f.numeric_missing\n              << ",\\"numeric_zero\\":" << f.numeric_zero << ",\\"absent_source_alleles\\":" << f.absent_source_alleles\n              << ",\\"sites_with_three_alt_scores\\":" << f.sites_with_three_alt_scores << \'}\';\n    }\n    stats << "\\n  },\\n  \\"observed_tri_allelic_snvs\\": null,\\n  \\"observed_quad_allelic_snvs\\": null,\\n"\n          << "  \\"observed_allele_statistics_reason\\": \\"Unavailable: Atlas predicts possible substitutions; population AC/AN/AF or genotypes are required to identify observed tri-/quad-allelic SNVs.\\",\\n"\n          << "  \\"prediction_coverage_definition\\": \\"Number of unique published ALT bases plus the shared reference base; no unlisted ALT is invented.\\"\\n}\\n";\n    atomic_text(stats_path, stats.str());\n    std::cerr << "merge " << chrom << " finished; sites=" << sites << " union_alt_rows=" << union_alt_rows\n              << " predicted_ref_plus_2alt=" << sites_by_nalt[2] << " predicted_ref_plus_3alt=" << sites_by_nalt[3] << \'\\n\';\n    return 0;\n}\n\nint main(int argc, char** argv) {\n    std::ios::sync_with_stdio(false);\n    std::cin.tie(nullptr);\n    try {\n        const Arguments args(argc, argv);\n        if (args.mode == "--help" || args.mode == "-h") {\n            std::cout << "AlphaGenome Atlas prediction coverage converter (C++17 + zlib)\\n"\n                      << "partition --kind avi|splicing|features --input FILE.gz --out-dir DIR --token TOKEN [--progress-every N]\\n"\n                      << "merge --chrom CHR --avi FILE.gz|NONE --splicing FILE.gz|NONE --features FILE.gz|NONE --stats FILE.json --token TOKEN [--progress-every N]\\n"\n                      << "Merge writes VCFv4.3 to stdout; progress/errors go to stderr. Inputs must contain sorted positions in contiguous chromosome blocks.\\n"\n                      << "Missing tokens: \'.\', \'NA\', \'NaN\', \'nan\'. Zero and signed finite score text are retained.\\n"\n                      << "All allele summaries describe predicted substitutions, not observed population alleles.\\n";\n            return 0;\n        }\n        return args.mode == "partition" ? partition(args) : merge(args);\n    } catch (const std::exception& error) {\n        std::cerr << "ERROR: " << error.what() << \'\\n\';\n        return 2;\n    }\n}\n'}
NATIVE_SOURCE = None

def runtime_resources():
    resources = dict(EMBEDDED_RESOURCES)
    for relative in ('allele3/f/batch_runner.py', '0f/qiming_resources.py',
                     '0f/lsf_submit.py', '0f/submission_io.py', '0f/workflow_lock.py',
                     '0f/workflow_cpu.py', '0f/0os.sh'):
        resources[relative] = (ROOT.parent / relative).read_text(encoding='utf-8')
    return resources


def prepare_runtime(out_dir, write=True):
    global NATIVE_SOURCE
    import tempfile
    resources = runtime_resources()
    bundle_id = digest(resources)
    base = Path(out_dir) / '.split2vcf.runtime' / bundle_id
    NATIVE_SOURCE = base / 'allele3/f/alphagenome_merge.cpp'
    if not write:
        return
    for relative, content in resources.items():
        path = base / relative
        data = content.encode('utf-8')
        if path.exists():
            if path.is_symlink() or path.read_bytes() != data:
                raise ValueError('Embedded runtime changed: ' + str(path))
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, temp = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
        try:
            with os.fdopen(fd, 'wb') as handle:
                handle.write(data)
            os.replace(temp, path)
        finally:
            if os.path.exists(temp):
                os.unlink(temp)
    sys.path.insert(0, str(base / 'allele3/f'))
    NATIVE_SOURCE = base / 'allele3/f/alphagenome_merge.cpp'


SCHEMA = "allele3-alphagenome-conversion-1"
VERSION = "2026-11-07.4-chromosome-dag"
FILES = {
    "avi": "alphagenome_variant_impact_score_snvs.tsv.gz",
    "splicing": "combined_alphagenome_splicing_snvs.tsv.gz",
    "features": "combined_ag_cond_linear_ensemble_20260417_feature_importance_indels_with_am_snvs.tsv.gz",
}
CHROMS = ["chr" + str(n) for n in range(1, 23)] + ["chrX", "chrY", "chrM"]
MAX_JSON = 16 * 1024 * 1024


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def read_json(path):
    path = Path(path)
    if path.stat().st_size > MAX_JSON:
        raise ValueError("Oversized checkpoint: " + str(path))
    with path.open() as handle:
        result = json.load(handle)
    if not isinstance(result, dict):
        raise ValueError("Checkpoint must be a JSON object: " + str(path))
    return result


def atomic_json(path, value):
    from batch_runner import atomic_json as write
    write(path, value)


def identity(path):
    """Cheap restart identity, explicitly not a full SHA of a TB input."""
    path = Path(path).resolve()
    st = path.stat()
    if not path.is_file() or st.st_size == 0:
        raise ValueError("Missing or empty file: " + str(path))
    with path.open("rb") as handle:
        head = handle.read(1024 * 1024)
        handle.seek(max(0, st.st_size - 1024 * 1024))
        tail = handle.read(1024 * 1024)
    if (path.stat().st_size, path.stat().st_mtime_ns) != (st.st_size, st.st_mtime_ns):
        raise ValueError("File changed while fingerprinting: " + str(path))
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                edge_sha256=hashlib.sha256(head + tail).hexdigest(),
                identity_method="size+mtime_ns+sha256(first_and_last_1MiB); not full-file hash")


def chrom_list(value):
    if value.lower() == "all":
        return None
    result = []
    for token in value.split(","):
        token = token.strip().removeprefix("chr")
        if "-" in token:
            lo, hi = map(int, token.split("-", 1))
            if not 1 <= lo <= hi <= 22:
                raise ValueError("Chromosome range must be within 1-22")
            items = [str(n) for n in range(lo, hi + 1)]
        else:
            items = ["M" if token in {"MT", "M"} else token]
        for item in items:
            chrom = "chr" + item
            if chrom not in CHROMS or chrom in result:
                raise ValueError("Invalid or repeated chromosome: " + chrom)
            result.append(chrom)
    if not result:
        raise ValueError("No chromosomes selected")
    return sorted(result, key=CHROMS.index)


class Log:
    def __init__(self, path):
        self.handle = Path(path).open("a", buffering=1)
        self.lock = threading.Lock()

    def __call__(self, message):
        line = f"{now()} {message}"
        with self.lock:
            self.handle.write(line + "\n")
            print(line, file=sys.stderr, flush=True)

    def close(self):
        self.handle.close()


class Processes:
    def __init__(self, log):
        self.log = log
        self.active = set()
        self.threads = []
        self.lock = threading.Lock()
        self.stopped = threading.Event()

    def start(self, argv, **kwargs):
        with self.lock:
            if self.stopped.is_set():
                raise RuntimeError("Conversion cancelled after another worker failed")
            proc = subprocess.Popen([str(x) for x in argv], **kwargs)
            self.active.add(proc)
            return proc

    def stderr(self, proc, label):
        def drain():
            for raw in iter(proc.stderr.readline, b""):
                self.log(f"[{label}] {raw.decode('utf-8', errors='replace').rstrip()}")
            proc.stderr.close()
        thread = threading.Thread(target=drain, daemon=True)
        self.threads.append(thread)
        thread.start()
        return thread

    def finish(self, proc, thread, label):
        code = proc.wait()
        thread.join()
        with self.lock:
            self.active.discard(proc)
        if code:
            raise RuntimeError(f"{label} failed with exit status {code}; see split2vcf.log")

    def run(self, argv, label):
        proc = self.start(argv, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        self.finish(proc, self.stderr(proc, label), label)

    def stop(self):
        self.stopped.set()
        with self.lock:
            active = list(self.active)
        for proc in active:
            if proc.poll() is None:
                proc.terminate()
        for proc in active:
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
        for thread in self.threads:
            thread.join(timeout=5)


@contextmanager
def locked(path, wait=False):
    with Path(path).open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | (0 if wait else fcntl.LOCK_NB))
        except BlockingIOError:
            raise RuntimeError("Another converter owns this output directory") from None
        yield


def parallel(items, jobs, function, processes):
    results = {}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(function, item): item for item in items}
        try:
            for future in as_completed(futures):
                results[futures[future]] = future.result()
        except BaseException:
            processes.stop()
            for future in futures:
                future.cancel()
            raise
    return results


def build_native(a, state, source_sha, processes, log):
    if a.native_binary:
        binary = Path(a.native_binary).expanduser().resolve()
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError("--native-binary is not executable")
        log("Using explicitly supplied native binary " + str(binary))
        return binary
    cxx = shutil.which(a.cxx)
    if not cxx:
        raise ValueError("C++17 compiler missing; install g++ and zlib headers, or use --native-binary")
    tools_dir = state / "tools"
    tools_dir.mkdir(exist_ok=True)
    binary = tools_dir / ("alphagenome_merge-" + source_sha[:16])
    receipt = binary.with_suffix(".json")
    if binary.is_file() and receipt.is_file():
        saved = read_json(receipt)
        if saved.get("source_sha256") == source_sha and saved.get("binary") == identity(binary):
            return binary
        raise ValueError("Native binary checkpoint differs; use a new output directory")
    temporary = binary.with_suffix(".partial")
    log("Compiling C++17 streaming engine")
    version = subprocess.run([cxx, "-dumpversion"], capture_output=True, text=True, timeout=15)
    # GCC 8 exposes C++17 filesystem in a separate support library.
    fs_libs = ["-lstdc++fs"] if version.stdout.strip().split(".")[0] == "8" else []
    processes.run([cxx, "-O3", "-DNDEBUG", "-std=c++17", "-Wall", "-Wextra",
                   NATIVE_SOURCE, "-lz", *fs_libs, "-o", temporary], "compiler")
    os.replace(temporary, binary)
    atomic_json(receipt, dict(source_sha256=source_sha, binary=identity(binary), compiler=cxx))
    return binary


def validate_vcf(vcf, index, chrom, first_pos, last_pos):
    import pysam
    with pysam.TabixFile(str(vcf), index=str(index)) as handle:
        if tuple(handle.contigs) != (chrom,):
            raise ValueError("Unexpected chromosome in index: " + str(index))
        for pos in {first_pos, last_pos}:
            row = next(handle.fetch(chrom, pos - 1, pos), None)
            if row is None or int(row.split("\t", 2)[1]) != pos:
                raise ValueError(f"Index failed exact boundary query: {chrom}:{pos}")


def completed_chrom(out, state, chrom, token, bounds):
    receipt = state / "completed" / (chrom + ".json")
    if not receipt.is_file():
        return None
    saved = read_json(receipt)
    if saved.get("schema") != SCHEMA or saved.get("token") != token:
        raise ValueError("Chromosome receipt belongs to another conversion: " + chrom)
    for expected in saved.get("products", []):
        if identity(expected["path"]) != expected:
            raise ValueError("Completed output changed: " + expected["path"])
    if len(saved.get("products", [])) != 2:
        raise ValueError("Incomplete chromosome receipt: " + chrom)
    validate_vcf(out / (chrom + ".vcf.bgz"), out / (chrom + ".vcf.bgz.tbi"), chrom, *bounds)
    return saved["stats"]


def merge_chrom(chrom, a, out, state, binary, sources, token, processes, log):
    parts = {kind: out / (kind + "." + chrom + ".tsv.gz") for kind in FILES}
    entries = [entry for source in sources.values() for entry in source["chromosomes"] if entry["chrom"] == chrom]
    bounds = min(e["first_pos"] for e in entries), max(e["last_pos"] for e in entries)
    cached = completed_chrom(out, state, chrom, token, bounds)
    if cached is not None:
        log(f"REUSE {chrom}: {cached['sites']:,} indexed prediction sites")
        return cached
    vcf = out / (chrom + ".vcf.bgz")
    index = Path(str(vcf) + ".tbi")
    # The manifest and exclusive lock establish ownership of these filenames.
    # A crash between publishing VCF/index/receipt may leave an uncommitted pair;
    # rebuild only that pair. Completed chromosomes are never overwritten.
    partial = out / (".split2vcf." + chrom + ".partial.vcf.bgz")
    partial_index = Path(str(partial) + ".tbi")
    stats_path = state / "completed" / (chrom + ".stats.partial.json")
    for path in (partial, partial_index, stats_path):
        if path.exists():
            path.unlink()
    cmd = [binary, "merge", "--chrom", chrom, "--stats", stats_path,
           "--token", token, "--progress-every", str(a.progress_every)]
    for kind, path in parts.items():
        present = any(e["chrom"] == chrom for e in sources[kind]["chromosomes"])
        cmd += ["--" + kind, path if present else "NONE"]
    log(f"MERGE {chrom}: outer join, Number=A arrays, BGZF level {a.compression_level}")
    native = processes.start(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    native_err = processes.stderr(native, chrom)
    bgzip = shutil.which("bgzip") if a.compressor == "auto" else None
    try:
        if bgzip:
            extra_threads = max(0, a.cores // min(a.jobs, a.cores) - 2)
            with partial.open("wb") as destination:
                compressor = processes.start([bgzip, "-c", "-l", str(a.compression_level),
                                              "-@", str(extra_threads)],
                                             stdin=native.stdout, stdout=destination, stderr=subprocess.PIPE)
                native.stdout.close()
                compressor_err = processes.stderr(compressor, chrom + "/bgzip")
                processes.finish(native, native_err, chrom + "/merge")
                processes.finish(compressor, compressor_err, chrom + "/bgzip")
                destination.flush()
                os.fsync(destination.fileno())
        else:
            import pysam
            log(f"{chrom}: using single-threaded pysam BGZF fallback (install bgzip for threaded compression)")
            with pysam.BGZFile(str(partial), "wb" + str(a.compression_level)) as destination:
                while block := native.stdout.read(4 * 1024 * 1024):
                    destination.write(block)
            native.stdout.close()
            processes.finish(native, native_err, chrom + "/merge")
        import pysam
        log("INDEX " + chrom)
        pysam.tabix_index(str(partial), preset="vcf", force=False, csi=False)
        validate_vcf(partial, partial_index, chrom, *bounds)
        stats = read_json(stats_path)
        if stats.get("token") != token or stats.get("chrom") != chrom:
            raise ValueError("Native statistics identity mismatch: " + chrom)
        os.replace(partial, vcf)
        os.replace(partial_index, index)
        atomic_json(state / "completed" / (chrom + ".json"),
                    dict(schema=SCHEMA, token=token, completed_at=now(),
                         products=[identity(vcf), identity(index)], stats=stats))
        stats_path.unlink()
        log(f"DONE {chrom}: sites={stats['sites']:,}, predicted ALT rows={stats['union_alt_rows']:,}")
        return stats
    except BaseException:
        if native.poll() is None:
            native.terminate()
        raise


def summaries(out, stats, token, requested, source_ids, log):
    columns = ["chrom", "prediction_sites", "predicted_alt_rows", "exactly_2_possible_bases",
               "exactly_3_possible_bases", "exactly_4_possible_bases", "at_least_3_possible_bases",
               "avi_source_rows", "splicing_source_rows", "features_source_rows"]
    rows = []
    for chrom in sorted(stats, key=CHROMS.index):
        s = stats[chrom]
        cov = s["prediction_base_coverage"]
        rows.append([chrom, s["sites"], s["union_alt_rows"], cov["exactly_2_bases"],
                     cov["exactly_3_bases"], cov["exactly_4_bases"], cov["at_least_3_bases"],
                     *[s["source_rows"][kind] for kind in FILES]])
    total = ["TOTAL"] + [sum(r[i] for r in rows) for i in range(1, len(columns))]
    temporary = out / ".split2vcf.summary.tsv.partial"
    with temporary.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows + [total]:
            handle.write("\t".join(map(str, row)) + "\n")
    os.replace(temporary, out / "split2vcf.summary.tsv")
    summary = dict(schema=SCHEMA, version=VERSION, token=token, status="complete",
                   scope="selected chromosomes present in input union", chromosomes=requested,
                   created_at=now(), totals=dict(zip(columns[1:], total[1:])),
                   observed_tri_allelic_sites=None, observed_quad_allelic_sites=None,
                   observed_allele_counts_reason="Unavailable: prediction tables contain no population GT/AC/AN/AF",
                   source_identity=source_ids, per_chromosome=stats)
    products = [out / (c + suffix) for c in requested for suffix in (".vcf.bgz", ".vcf.bgz.tbi")]
    products.append(out / "split2vcf.summary.tsv")
    summary["output_inventory"] = stat_inventory(out, products)
    atomic_json(out / "split2vcf.summary.json", summary)
    log("ALLELE SUMMARY: " + json.dumps(summary["totals"], sort_keys=True))
    log("3rd allele inclusive (>=3 POSSIBLE bases) = " + str(total[6]))
    log("4th allele (4 POSSIBLE bases) = " + str(total[5]))
    for kind in FILES:
        coverage = {key: sum(s["per_source"][kind][key] for s in stats.values())
                    for key in next(iter(stats.values()))["per_source"][kind]}
        log("SOURCE COVERAGE " + kind + ": " + json.dumps(coverage, sort_keys=True))
    masks = {str(n): sum(s["source_mask_alleles"].get(str(n), 0) for s in stats.values()) for n in range(1, 8)}
    log("SOURCE MASK ALT COUNTS (AVI=1, splicing=2, features=4): " + json.dumps(masks, sort_keys=True))
    log("Observed population tri-/quad-allelic sites and carrier counts: NA (no GT/AC/AN/AF)")
    log("Finished. Compressed staging retained for restart; original downloads unchanged.")
    return summary


def stat_inventory(root, paths):
    return {str(path.relative_to(root)): dict(size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)
            for path in paths}


def checked_partition(directory, kind, token, prefix=""):
    receipt = directory / (kind + ".source.done.json")
    if not receipt.is_file():
        return None
    saved = read_json(receipt)
    if saved.get("schema") != "alphagenome-tsv-partition-1" or saved.get("token") != token or saved.get("kind") != kind:
        raise ValueError("Partition receipt identity mismatch: " + kind)
    seen = set()
    for entry in saved["chromosomes"]:
        chrom = entry["chrom"]
        if chrom not in CHROMS or chrom in seen:
            raise ValueError("Invalid/duplicate chromosome in partition receipt")
        seen.add(chrom)
        part = directory / (prefix + chrom + ".tsv.gz")
        marker = directory / (prefix + chrom + ".done.json")
        if not marker.is_file() or read_json(marker) != entry:
            return None
        if not part.is_file() or part.stat().st_size != entry["gzip_bytes"] or part.stat().st_size < 18:
            return None
        with part.open("rb") as handle:
            if handle.read(2) != b"\x1f\x8b":
                return None
            handle.seek(-8, 2)
            crc, length = struct.unpack("<II", handle.read(8))
        if crc != entry["uncompressed_crc32"] or length != entry["uncompressed_bytes"] % 2**32:
            return None
    if sum(e["rows"] for e in saved["chromosomes"]) != saved["source_rows"]:
        raise ValueError("Partition source/chr row counts disagree")
    return saved


def partition_token(kind, source_ids):
    return digest(dict(source=source_ids[kind], kind=kind,
                       native=hashlib.sha256(NATIVE_SOURCE.read_bytes()).hexdigest()))


def checked_chrom(out, kind, chrom, token):
    marker = out / (kind + '.' + chrom + '.done.json')
    part = out / (kind + '.' + chrom + '.tsv.gz')
    if not marker.is_file():
        return None
    e = read_json(marker)
    if (e.get('schema'), e.get('kind'), e.get('chrom'), e.get('token')) != (
            'alphagenome-tsv-partition-1', kind, chrom, token):
        raise ValueError('Incompatible chromosome receipt: ' + str(marker))
    if not part.is_file() or part.stat().st_size != e['gzip_bytes'] or e['gzip_bytes'] < 18:
        return None
    with part.open('rb') as h:
        if h.read(2) != b'\x1f\x8b':
            return None
        h.seek(-8, 2)
        crc, size = struct.unpack('<II', h.read(8))
    return e if (crc, size) == (e['uncompressed_crc32'], e['uncompressed_bytes'] % 2**32) else None


def dag_context(a):
    import pysam
    out = Path(a.out_dir)
    state = out / '.split2vcf'
    state.mkdir(parents=True, exist_ok=True)
    (state / 'completed').mkdir(exist_ok=True)
    source_ids = {k: identity(getattr(a, k)) for k in FILES}
    native_sha = hashlib.sha256(NATIVE_SOURCE.read_bytes()).hexdigest()
    recipe = dict(schema=SCHEMA, version=VERSION, source_identity=source_ids,
                  native_source_sha256=native_sha, compression_level=a.compression_level,
                  orchestrator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  reference='GRCh38', join='outer CHROM,POS,REF,ALT; one record per POS/REF',
                  explicit_native_binary=identity(a.native_binary) if a.native_binary else None)
    with locked(out / '.split2vcf.lock', wait=True):
        path = state / 'manifest.json'
        if path.exists():
            saved = read_json(path)
            old = saved['recipe']
            if old != recipe:
                raise ValueError('Conversion code/recipe changed; use a new output directory')
            if digest(old) != saved['token']:
                raise ValueError('Invalid original conversion identity')
            token = saved['token']
        else:
            if any(out.glob('chr*.vcf.bgz')):
                raise ValueError('Unowned VCF outputs; use a distinct output directory')
            token = digest(recipe)
            atomic_json(path, dict(token=token, recipe=recipe, created_at=now()))
        orchestrator = state / 'dag.orchestrator.json'
        record = dict(schema='alphagenome-dag-1', conversion_token=token,
                      orchestrator_sha256=recipe['orchestrator_sha256'],
                      partition_method='Per-chromosome Tabix extraction; isolated staging; same-chromosome merge dependencies')
        if orchestrator.exists() and read_json(orchestrator) != record:
            raise ValueError('DAG orchestrator changed; finish or retire its jobs before changing code')
        atomic_json(orchestrator, record)
    contigs = {}
    for k in FILES:
        with pysam.TabixFile(getattr(a, k)) as t:
            contigs[k] = list(t.contigs)
        if not contigs[k] or any(c not in CHROMS for c in contigs[k]):
            raise ValueError('Unsupported/empty chromosome inventory for ' + k)
    return out, state, source_ids, token, contigs


def dag_source(out, kind, source_ids, contigs):
    token = partition_token(kind, source_ids)
    saved = checked_partition(out, kind, token, kind + '.')
    if saved is not None and {e['chrom'] for e in saved['chromosomes']} == set(contigs[kind]):
        return saved
    return None


def finalize_source(out, kind, source_ids, contigs):
    # Only publication of the small source receipt is serialized. Extraction
    # and compression for different chromosomes never share a lock or filename.
    with locked(out / ('.split2vcf.partition.' + kind + '.lock'), wait=True):
        token = partition_token(kind, source_ids)
        rows = [checked_chrom(out, kind, c, token) for c in contigs[kind]]
        if any(e is None for e in rows):
            return
        if dag_source(out, kind, source_ids, contigs) is None:
            atomic_json(out / (kind + '.source.done.json'), dict(
                schema='alphagenome-tsv-partition-1', kind=kind, token=token,
                source_rows=sum(e['rows'] for e in rows), chromosomes=rows))


def dag_partition(a, kind, out, state, source_ids, contigs, processes, log):
    import tempfile
    global NATIVE_SOURCE
    chrom = a.chrom
    if chrom not in contigs[kind]:
        raise ValueError('Partition worker requires --chrom from the source index')
    label = kind + '.' + chrom
    token = partition_token(kind, source_ids)
    with locked(out / ('.split2vcf.partition.' + label + '.lock')):
        if checked_chrom(out, kind, chrom, token) is not None:
            log('SKIP ' + label + ': chromosome part verified')
            finalize_source(out, kind, source_ids, contigs)
            return
        tabix = shutil.which('tabix')
        if not tabix:
            raise ValueError('Indexed partition requires tabix on PATH')
        original = NATIVE_SOURCE
        try:
            NATIVE_SOURCE = original.with_name('alphagenome_partition_stdin.cpp')
            with locked(out / '.split2vcf.compile.lock', wait=True):
                binary = build_native(a, state, hashlib.sha256(NATIVE_SOURCE.read_bytes()).hexdigest(), processes, log)
        finally:
            NATIVE_SOURCE = original
        log('PARTITION ' + label + ': independent indexed chromosome extraction')
        # The native engine also emits a source receipt. Keep every worker in
        # its own directory and publish only after BOTH pipeline processes pass.
        # This prevents a failed Tabix stream from publishing a truncated part.
        staging_root = state / 'partition-staging'
        staging_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=label + '.', dir=staging_root) as temporary:
            staging = Path(temporary)
            extract = processes.start([tabix, '-h', getattr(a, kind), chrom],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            extract_err = processes.stderr(extract, label + '/tabix')
            try:
                native = processes.start([binary, 'partition', '--kind', kind, '--input', '-',
                                          '--out-dir', staging, '--file-prefix', kind + '.', '--token', token,
                                          '--progress-every', str(a.progress_every)],
                                         stdin=extract.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                extract.stdout.close()
                native_err = processes.stderr(native, label + '/partition')
                processes.finish(native, native_err, label + '/partition')
                processes.finish(extract, extract_err, label + '/tabix')
                if identity(getattr(a, kind)) != source_ids[kind]:
                    raise ValueError('Source changed during indexed partition: ' + label)
                saved = checked_partition(staging, kind, token, kind + '.')
                if saved is None or [e['chrom'] for e in saved['chromosomes']] != [chrom]:
                    raise ValueError('Incomplete or unexpected indexed chromosome: ' + label)
                # Remove an invalid old receipt before replacing its data. A
                # valid completed part was already skipped above.
                (out / (label + '.done.json')).unlink(missing_ok=True)
                os.replace(staging / (label + '.tsv.gz'), out / (label + '.tsv.gz'))
                os.replace(staging / (label + '.done.json'), out / (label + '.done.json'))
            except BaseException:
                processes.stop()
                raise
        finalize_source(out, kind, source_ids, contigs)
        log('DONE partition ' + label)


def dag_worker(a):
    if a.submit == 'hpc' and not os.environ.get('LSB_JOBID'):
        raise ValueError('Worker requires an LSF allocation')
    if a.submit == 'local' and shutil.which('bsub') and not os.environ.get('LSB_JOBID'):
        raise ValueError('Local worker refused on an LSF login host')
    os.environ['PATH'] = str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', '')
    out, state, source_ids, token, contigs = dag_context(a)
    log = Log(out / 'split2vcf.log')
    processes = Processes(log)
    try:
        if a.stage == 'partition':
            dag_partition(a, a.kind, out, state, source_ids, contigs, processes, log)
        elif a.stage == 'merge':
            chrom = a.chrom
            with locked(out / ('.split2vcf.merge.' + chrom + '.lock')):
                sources = {}
                for kind in FILES:
                    e = checked_chrom(out, kind, chrom, partition_token(kind, source_ids)) if chrom in contigs[kind] else None
                    if chrom in contigs[kind] and e is None:
                        raise ValueError('Required part not complete: ' + kind + '.' + chrom)
                    sources[kind] = dict(chromosomes=[e] if e else [])
                with locked(out / '.split2vcf.compile.lock', wait=True):
                    binary = build_native(a, state, hashlib.sha256(NATIVE_SOURCE.read_bytes()).hexdigest(), processes, log)
                merge_chrom(chrom, a, out, state, binary, sources, token, processes, log)
        elif a.stage == 'summary':
            with locked(out / '.split2vcf.summary.lock'):
                requested = chrom_list(a.chromosomes) or sorted(set().union(*map(set, contigs.values())), key=CHROMS.index)
                stats = {}
                for c in requested:
                    es = [checked_chrom(out, k, c, partition_token(k, source_ids)) for k in FILES if c in contigs[k]]
                    if not es or any(e is None for e in es):
                        raise ValueError('Missing partition: ' + c)
                    result = completed_chrom(out, state, c, token, (min(e['first_pos'] for e in es), max(e['last_pos'] for e in es)))
                    if result is None:
                        raise ValueError('Missing verified VCF/index: ' + c)
                    stats[c] = result
                summaries(out, stats, token, requested, source_ids, log)
        else:
            raise ValueError('Worker requires an explicit DAG stage')
    finally:
        processes.stop()
        log.close()


def dag_plan(a):
    out, state, source_ids, token, contigs = dag_context(a)
    requested = chrom_list(a.chromosomes) or sorted(set().union(*map(set, contigs.values())), key=CHROMS.index)
    if any(c not in set().union(*map(set, contigs.values())) for c in requested):
        raise ValueError('Requested chromosome absent from all sources')
    py = a.python or sys.executable
    base = [py, str(Path(__file__).resolve()), '--worker', '--submit', a.submit,
            '--input-dir', a.input_dir, '--out-dir', a.out_dir, '--jobs', '1']
    for key in (*FILES, 'cores', 'compression_level', 'compressor', 'progress_every', 'cxx', 'native_binary'):
        if getattr(a, key) is not None:
            base += ['--' + key.replace('_', '-'), str(getattr(a, key))]
    jobs = []
    parts = {kind: {c: checked_chrom(out, kind, c, partition_token(kind, source_ids))
                    for c in contigs[kind]} for kind in FILES}
    if a.phase != 'merge':
        for kind in FILES:
            for c in requested:
                if c not in parts[kind] or parts[kind][c] is not None:
                    continue
                jobs.append(dict(id='partition_' + kind + '.' + c,
                    argv=base+['--stage', 'partition', '--kind', kind, '--chrom', c],
                    cwd=str(out), resources=dict(queue=a.partition_queue or a.queue,
                        cores=a.partition_cores, memory_gib=a.partition_memory_gb),
                    outputs=[str(out/(kind+'.'+c+suffix)) for suffix in ('.tsv.gz', '.done.json')]))
    if a.phase != 'partition':
        merge_jobs = []
        for c in requested:
            deps = []
            entries = []
            for kind in FILES:
                if c not in parts[kind]:
                    continue
                entry = parts[kind][c]
                if entry is None:
                    if a.phase == 'merge':
                        raise ValueError('Merge-only requires complete part: ' + kind + '.' + c)
                    deps.append('partition_' + kind + '.' + c)
                else:
                    entries.append(entry)
            if not deps and completed_chrom(out, state, c, token,
                    (min(e['first_pos'] for e in entries), max(e['last_pos'] for e in entries))) is not None:
                continue
            jobs.append(dict(id=c, argv=base+['--stage', 'merge', '--chrom', c], cwd=str(out), deps=deps,
                outputs=[str(out/(c+s)) for s in ('.vcf.bgz', '.vcf.bgz.tbi')] + [str(state/'completed'/(c+'.json'))]))
            merge_jobs.append(c)
        jobs.append(dict(id='summary', argv=base+['--stage','summary','--chromosomes',','.join(requested)],
                         cwd=str(out), deps=merge_jobs, outputs=[str(out/'split2vcf.summary.json'), str(out/'split2vcf.summary.tsv')]))
    cores = a.job_cores or a.cores
    if cores < a.cores:
        raise ValueError('--job-cores cannot be smaller than --cores')
    plan = dict(name='alphagenome', driver_python=py, max_running=max(1, len(jobs)),
                resources=dict(queue=a.queue, cores=cores, memory_gib=a.job_memory_gb, hours=a.job_hours),
                local_memory_gib=a.job_memory_gb,
                inputs=[dict(path=str(Path(__file__).resolve()),sha256=True),dict(path=str(NATIVE_SOURCE),sha256=True)] +
                       [dict(path=getattr(a,k),sha256=False) for k in FILES] +
                       [dict(path=getattr(a,k)+'.tbi',sha256=True) for k in FILES],
                jobs=jobs)
    key = digest(dict(base=base, phase=a.phase, chromosomes=requested, resources=plan['resources'], partition_queue=a.partition_queue, partition_cores=a.partition_cores, partition_memory_gb=a.partition_memory_gb, token=token, code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))[:16]
    job_dir = Path(a.job_dir or str(out)+'.split2vcf.jobs/dag_'+key)
    saved_plan = job_dir / 'dag.plan.json'
    if saved_plan.is_file():
        plan = read_json(saved_plan)
    elif not a.dry_run:
        job_dir.mkdir(parents=True, exist_ok=True)
        atomic_json(saved_plan, plan)
    return plan, job_dir


def dag_submit(a):
    with locked(Path(a.out_dir) / ".split2vcf.submit.lock"):
        return dag_submit_locked(a)


def dag_submit_locked(a):
    plan, job_dir = dag_plan(a)
    if a.dry_run:
        print(json.dumps(dict(status='plan_only',job_dir=str(job_dir),plan=plan),indent=2))
        return 0
    if not plan['jobs']:
        print(json.dumps(dict(status='already_complete', job_dir=str(job_dir))))
        return 0
    import batch_runner
    original = batch_runner.shell_text
    def named(plan, job, out):
        import re
        text = original(plan, job, out)
        label = 'alphagenome.' + job['id'].replace('partition_', 'split.')
        return re.sub(r'^#BSUB -J .*$', '#BSUB -J "'+label+'"', text, count=1, flags=re.MULTILINE)
    batch_runner.shell_text = named
    try:
        result = batch_runner.run_plan(plan, job_dir, mode=a.submit, write_only=a.write_cmd,
                                      jobs_n=a.jobs, wait=a.wait=='TRUE', retry_failed=a.retry_failed)
    finally:
        batch_runner.shell_text = original
    print(json.dumps(dict(job_dir=str(job_dir),result=result),indent=2,default=str))
    return 0


def parser():
    class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    examples = r"""Example:
  cd /work/sph-huangj/data/annot/alphagenome
  /work/sph-huangj/.conda/envs/allele3-evolution/bin/python alphag2vcf.py \
    --submit hpc --queue v3-64

  --write-cmd     Generate the DAG without submission.
  --retry-failed  Retry failed jobs in the same frozen DAG.
  --phase status Show completed source/chromosome receipts.

DAG: independent missing source/chromosome splits -> same-chromosome merge -> summary.
Completed splits and VCFs skip. Tabix reads only each missing chromosome.
A merge with all three parts complete starts immediately; other merges depend
only on the missing splits for that chromosome. LSF decides concurrency.
Outputs: avi/splicing/features.chr*.tsv.gz and chr*.vcf.bgz(.tbi) in --out-dir.
A partial chromosome is rebuilt from its beginning. Completed parts are retained.
"""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=HelpFormatter, epilog=examples)
    p.add_argument("--input-dir", help="Original BGZF TSV directory; default ALLELE_ATLAS_DIR/raw or 38/raw beside this script")
    p.add_argument("--out-dir", "--output-dir", dest="out_dir", help="VCF/index/log directory; default ALLELE_ATLAS_DIR or 38 beside this script")
    for kind, filename in FILES.items():
        p.add_argument("--" + kind, help="Override input: " + filename)
    p.add_argument("--chromosomes", default="all", help="all, 1-22,X,Y,M, or subset; split only missing requested chromosomes")
    p.add_argument("--phase", choices=("all", "partition", "merge", "status"), default="all")
    p.add_argument("--cores", "--threads", dest="cores", type=int, default=8, help="CPU budget (>=2 for streaming merge+compress)")
    p.add_argument("--jobs", "--jobs-N", type=int, default=None, help="Local DAG concurrency; HPC creates one job per chromosome")
    p.add_argument("--partition-jobs", type=int, default=3, help="Concurrent source partition streams; 3 for fast storage")
    p.add_argument("--compression-level", type=int, choices=range(1, 10), default=1, help="BGZF zlib compression level")
    p.add_argument("--compressor", choices=("auto", "pysam"), default="auto", help="auto uses bgzip when installed, else pysam")
    p.add_argument("--progress-every", type=int, default=10000000, help="Native progress interval in rows")
    p.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    p.add_argument("--native-binary", help="Explicit precompiled engine built from included source")
    p.add_argument("--python", help="Worker Python with pysam; default current Python")
    p.add_argument("--submit", choices=("auto", "local", "hpc"), default=os.environ.get("submit", os.environ.get("SUBMIT", "auto")))
    p.add_argument("--job-dir", default=os.environ.get("ALLELE_JOB_DIR"))
    p.add_argument("--queue", default=os.environ.get("ALLELE_JOB_QUEUE", os.environ.get("HPC_CPU_QUEUE", "v3-64")))
    p.add_argument("--job-cores", type=int, default=int(os.environ.get("ALLELE_JOB_CORES", "24")))
    p.add_argument("--job-memory-gb", type=int, default=int(os.environ.get("ALLELE_JOB_MEMORY_GIB", "32")))
    p.add_argument("--job-hours", type=int, default=int(os.environ.get("ALLELE_JOB_HOURS", "168")))
    p.add_argument("--wait", choices=("TRUE", "FALSE"), default=os.environ.get("SUBMIT_WAIT", "FALSE"))
    p.add_argument("--write-cmd", action="store_true", default=os.environ.get("ALLELE_WRITE_CMD") == "1")
    p.add_argument("--retry-failed", action="store_true", default=os.environ.get("ALLELE_RETRY_FAILED") == "1")
    p.add_argument("--dry-run", action="store_true", help="Print plan; no compile, decompress or conversion")
    p.add_argument("--partition-queue", help="Optional queue override for source jobs")
    p.add_argument("--partition-cores", type=int, default=2, help="CPU slots per chromosome split (scheduler may round up)")
    p.add_argument("--partition-memory-gb", type=int, default=4, help="GiB per chromosome split")
    p.add_argument("--stage", choices=("partition", "merge", "summary"), help=argparse.SUPPRESS)
    p.add_argument("--kind", choices=tuple(FILES), help=argparse.SUPPRESS)
    p.add_argument("--chrom", choices=CHROMS, help=argparse.SUPPRESS)
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--help-full", action="help", help="Show all conversion options")
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    a = parser().parse_args(['--help'] if argv == ['help'] else argv)
    entry_defaults(a)
    for key in ("cores", "jobs", "partition_jobs", "progress_every", "job_memory_gb", "job_hours", "partition_cores", "partition_memory_gb"):
        if getattr(a, key) < 1:
            raise ValueError("--" + key.replace("_", "-") + " must be positive")
    if a.cores < 2:
        raise ValueError("--cores must be >=2 for native merge and compression")
    a.input_dir = str(Path(a.input_dir).expanduser().resolve())
    a.out_dir = str(Path(a.out_dir or a.input_dir).expanduser().resolve())
    if a.python:
        a.python = shutil.which(a.python) or str(Path(a.python).expanduser().resolve())
    if a.native_binary:
        a.native_binary = str(Path(a.native_binary).expanduser().resolve())
    for kind, filename in FILES.items():
        setattr(a, kind, str(Path(getattr(a, kind) or Path(a.input_dir) / filename).expanduser().resolve()))
    if len({getattr(a, kind) for kind in FILES}) != 3:
        raise ValueError("Three inputs must be distinct files")
    if a.phase == "status":
        state = Path(a.out_dir) / ".split2vcf"
        print(json.dumps(dict(output=a.out_dir,
                              sources_completed=[k for k in FILES if (Path(a.out_dir) / (k + ".source.done.json")).is_file()],
                              chromosomes_completed=[c for c in CHROMS if (state / "completed" / (c + ".json")).is_file()],
                              note="Receipt inventory only; conversion rerun validates identities and indexes"), indent=2))
        return 0
    chrom_list(a.chromosomes)
    if not os.access(a.python, os.X_OK):
        raise ValueError('Missing scientific Python: ' + a.python + '; use --python PATH')
    if Path(a.python).resolve() != Path(sys.executable).resolve():
        os.execv(a.python, [a.python, str(Path(__file__).resolve()), *argv])
    prepare_runtime(a.out_dir, write=True)
    # Fail before creating a plan or contacting LSF for missing/empty inputs.
    if not a.dry_run:
        for kind in FILES:
            path = Path(getattr(a, kind))
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError("Missing or empty " + kind + " input: " + str(path) +
                                 "; set --input-dir to the directory containing the original TSV.gz files")
    if a.worker:
        dag_worker(a)
        return 0
    return dag_submit(a)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError, KeyError) as error:
        raise SystemExit("split2vcf: " + str(error))
