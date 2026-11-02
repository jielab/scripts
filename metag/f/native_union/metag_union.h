/*
 * Metag additions to Kraken 2, distributed under Kraken 2's MIT license.
 * The underlying CompactHashTable and classifier remain the pinned upstream
 * implementation.  This adapter combines lookup evidence, not final calls.
 */
#ifndef METAG_UNION_H_
#define METAG_UNION_H_

#include "compact_hash.h"
#include "taxonomy.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <unistd.h>

namespace metag_union {

using kraken2::hkey_t;
using kraken2::hvalue_t;
using kraken2::KeyValueStore;
using kraken2::Taxonomy;

static std::string hashes_filename;

// Strip only this extension's arguments; retain the upstream option parser and
// its interpretation of -S (interleaved pairs), -T, -g, -Q, -R, -K, and so on.
inline void ExtractArguments(int &argc, char **argv) {
  int output = 1;
  bool positional = false;
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (!positional && arg == "--") positional = true;
    if (!positional && (arg == "--metag-hashes" ||
                         arg.compare(0, 15, "--metag-hashes=") == 0)) {
      if (!hashes_filename.empty())
        errx(EX_USAGE, "--metag-hashes may appear only once");
      if (arg == "--metag-hashes") {
        if (++i == argc) errx(EX_USAGE, "--metag-hashes requires a filename");
        hashes_filename = argv[i];
      } else {
        hashes_filename = arg.substr(15);
      }
      if (hashes_filename.empty()) errx(EX_USAGE, "empty --metag-hashes filename");
    } else {
      argv[output++] = argv[i];
    }
  }
  argc = output;
  argv[argc] = nullptr;
}

inline uint64_t FileSize(const std::string &path) {
  struct stat sb;
  if (stat(path.c_str(), &sb) != 0 || !S_ISREG(sb.st_mode) || sb.st_size < 0)
    errx(EX_NOINPUT, "metag union: not a readable regular file: %s", path.c_str());
  return static_cast<uint64_t>(sb.st_size);
}

inline std::string CanonicalPath(const std::string &path) {
  char *resolved = realpath(path.c_str(), nullptr);
  if (!resolved) errx(EX_NOINPUT, "metag union: cannot resolve %s", path.c_str());
  std::string result(resolved);
  free(resolved);
  return result;
}

inline std::string Sibling(const std::string &path, const char *name) {
  auto slash = path.find_last_of('/');
  return path.substr(0, slash + 1) + name;
}

inline void RequireSameFileBytes(const std::string &a, const std::string &b) {
  const auto size = FileSize(a);
  if (FileSize(b) != size)
    errx(EX_DATAERR, "metag union: incompatible metadata lengths: %s and %s", a.c_str(), b.c_str());
  // Inode identity is sufficient for a shared hardlink/symlink target. Published
  // database files must remain immutable; the caller verifies their runtime
  // state before and after query and before committing any sample outputs.
  struct stat sa, sb;
  if (stat(a.c_str(), &sa) || stat(b.c_str(), &sb))
    errx(EX_NOINPUT, "metag union: metadata disappeared during validation");
  if (sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino) return;
  std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
  if (!fa || !fb) errx(EX_NOINPUT, "metag union: cannot open metadata for comparison");
  std::array<char, 65536> ba, bb;
  uint64_t left = size;
  while (left) {
    const auto count = static_cast<std::streamsize>(std::min<uint64_t>(left, ba.size()));
    fa.read(ba.data(), count); fb.read(bb.data(), count);
    if (fa.gcount() != count || fb.gcount() != count ||
        std::memcmp(ba.data(), bb.data(), static_cast<size_t>(count)))
      errx(EX_DATAERR, "metag union: metadata is not byte-identical: %s and %s", a.c_str(), b.c_str());
    left -= static_cast<uint64_t>(count);
  }
}

struct HashHeader {
  uint64_t capacity, used, key_bits, value_bits;
};

inline HashHeader ReadHashHeader(const std::string &path, uint64_t nodes) {
  static_assert(sizeof(size_t) == 8, "metag union requires the pinned 64-bit Kraken format");
  const uint16_t endian_probe = 1;
  if (*reinterpret_cast<const unsigned char *>(&endian_probe) != 1)
    errx(EX_CONFIG, "metag union requires little-endian native index files");
  const auto size = FileSize(path);
  HashHeader h{};
  std::ifstream input(path, std::ios::binary);
  input.read(reinterpret_cast<char *>(&h), sizeof(h));
  if (!input || size < sizeof(h) || !h.capacity || h.used > h.capacity ||
      !h.key_bits || !h.value_bits || h.value_bits >= 32 ||
      (h.key_bits + h.value_bits != 32 && h.key_bits + h.value_bits != 40))
    errx(EX_DATAERR, "metag union: invalid compact hash header: %s", path.c_str());
  const uint64_t cell_bytes = (h.key_bits + h.value_bits) / 8;
  if (h.capacity > (std::numeric_limits<uint64_t>::max() - sizeof(h)) / cell_bytes ||
      size != sizeof(h) + h.capacity * cell_bytes)
    errx(EX_DATAERR, "metag union: truncated or oversized compact hash: %s", path.c_str());
  if (nodes > (uint64_t(1) << h.value_bits))
    errx(EX_DATAERR, "metag union: taxonomy does not fit compact hash value bits: %s", path.c_str());
  return h;
}

inline uint64_t ValidateTaxonomyFile(const std::string &path) {
  const auto size = FileSize(path);
  char magic[8];
  uint64_t sizes[3] = {0, 0, 0};
  std::ifstream input(path, std::ios::binary);
  input.read(magic, sizeof(magic));
  input.read(reinterpret_cast<char *>(sizes), sizeof(sizes));
  if (!input || std::memcmp(magic, "K2TAXDAT", 8) || sizes[0] < 2 ||
      sizes[0] > (std::numeric_limits<uint64_t>::max() - 32) / sizeof(kraken2::TaxonomyNode))
    errx(EX_DATAERR, "metag union: invalid taxonomy header: %s", path.c_str());
  uint64_t expected = 32 + sizes[0] * sizeof(kraken2::TaxonomyNode);
  for (int i = 1; i < 3; ++i) {
    if (sizes[i] > std::numeric_limits<uint64_t>::max() - expected)
      errx(EX_DATAERR, "metag union: taxonomy length overflow: %s", path.c_str());
    expected += sizes[i];
  }
  if (size != expected)
    errx(EX_DATAERR, "metag union: truncated or oversized taxonomy: %s", path.c_str());
  return sizes[0];
}

inline std::vector<std::string> ValidateInputs(const std::string &taxo,
                                              const std::string &options) {
  if (hashes_filename.empty())
    errx(EX_USAGE, "metag_classify requires --metag-hashes LISTFILE");
  if (FileSize(options) != sizeof(kraken2::IndexOptions))
    errx(EX_DATAERR, "metag union: options file has an incompatible native size: %s", options.c_str());
  const uint64_t nodes = ValidateTaxonomyFile(taxo);
  FileSize(hashes_filename);
  std::ifstream input(hashes_filename);
  if (!input) errx(EX_NOINPUT, "metag union: cannot read hash list");
  std::set<std::string> unique;
  std::vector<std::string> hashes;
  std::string line;
  HashHeader first{};
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] != '/' || line.find('\0') != std::string::npos)
      errx(EX_DATAERR, "metag union: hash list requires one absolute path per nonempty line");
    const auto path = CanonicalPath(line);
    if (!unique.insert(path).second)
      errx(EX_DATAERR, "metag union: duplicate hash path: %s", path.c_str());
    // Use the supplied category directory for sibling metadata, even when the
    // immutable hash itself is a symlink to a generation directory.
    RequireSameFileBytes(taxo, Sibling(line, "taxo.k2d"));
    RequireSameFileBytes(options, Sibling(line, "opts.k2d"));
    const auto header = ReadHashHeader(path, nodes);
    if (!hashes.empty() && (header.key_bits != first.key_bits || header.value_bits != first.value_bits))
      errx(EX_DATAERR, "metag union: category hash cell widths differ");
    if (hashes.empty()) first = header;
    hashes.push_back(path);
  }
  if (!input.eof()) errx(EX_IOERR, "metag union: error reading hash list");
  if (hashes.empty()) errx(EX_DATAERR, "metag union: empty hash list");
  return hashes;
}

class UnionHashTable final : public KeyValueStore {
 public:
  UnionHashTable(const std::vector<std::string> &paths, const Taxonomy &taxonomy)
      : taxonomy_(taxonomy), node_count_(taxonomy.node_count()) {
    for (const auto &path : paths) {
      auto h = ReadHashHeader(path, node_count_);
      if (h.key_bits + h.value_bits == 32)
        tables_.emplace_back(new kraken2::CompactHashTable<kraken2::CompactHashCell>(path, true));
      else
        tables_.emplace_back(new kraken2::CompactHashTable<kraken2::CompactHashCell40>(path, true));
    }
  }

  hvalue_t Get(hkey_t key) const override {
    hvalue_t value = 0;
    for (const auto &table : tables_) value = Combine(value, table->Get(key));
    return value;
  }

  void GetBatch(const hkey_t *keys, hvalue_t *out, size_t n) const override {
    // Each table retains upstream's batched lookahead prefetch. Each worker has
    // one scratch buffer reused across pairs: no FASTQ reread or intermediate
    // per-category call files, and no per-minimizer heap allocation.
    if (!n) return;
    tables_.front()->GetBatch(keys, out, n);
    for (size_t i = 0; i < n; ++i) ValidateValue(out[i]);
    static thread_local std::vector<hvalue_t> scratch;
    scratch.resize(n);
    for (size_t table = 1; table < tables_.size(); ++table) {
      tables_[table]->GetBatch(keys, scratch.data(), n);
      for (size_t i = 0; i < n; ++i) out[i] = Combine(out[i], scratch[i]);
    }
  }

 private:
  void ValidateValue(hvalue_t value) const {
    if (value >= node_count_)
      errx(EX_DATAERR, "metag union: hash value exceeds the shared taxonomy");
  }
  hvalue_t Combine(hvalue_t a, hvalue_t b) const {
    ValidateValue(b);
    if (!b || a == b) return a;
    if (!a) return b;
    return static_cast<hvalue_t>(taxonomy_.LowestCommonAncestor(a, b));
  }
  const Taxonomy &taxonomy_;
  uint64_t node_count_;
  std::vector<std::unique_ptr<KeyValueStore>> tables_;
};

}  // namespace metag_union
#endif
