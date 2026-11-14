// A separate, source-pinned Kraken builder with transactional hash snapshots.
// Upstream minimizer/insertion algorithms are included without alteration.
// The private compact_hash.h patch adds checked writes and mutable-load locks.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <memory>
#include <regex>
#include <spawn.h>
#include <stdexcept>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/wait.h>

#define main kraken_original_build_main
#include "build_db.cc"
#undef main

extern char** environ;

namespace metag {
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Fields = std::map<std::string, std::string>;
static volatile sig_atomic_t stop_signal = 0;
static void on_signal(int signal_number) { stop_signal = signal_number; }
static const std::string SCHEMA = "metag-native-checkpoint-v1";

static void fail(const std::string& message) { throw std::runtime_error(message); }
static void require(bool condition, const std::string& message) { if (!condition) fail(message); }
static std::string number(uint64_t n) { return std::to_string(n); }
static std::string seconds(double n) {
  std::ostringstream out; out << std::fixed << std::setprecision(6) << n; return out.str();
}
static double elapsed(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
static bool hex64(const std::string& s) {
  return s.size() == 64 && s.find_first_not_of("0123456789abcdef") == std::string::npos;
}
static uint64_t integer(const std::string& s, const std::string& label) {
  require(!s.empty() && s.find_first_not_of("0123456789") == std::string::npos,
          "invalid nonnegative integer for " + label);
  size_t used = 0;
  uint64_t value = std::stoull(s, &used);
  require(used == s.size(), "invalid integer for " + label); return value;
}
static std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> result; size_t start = 0;
  for (;;) {
    size_t end = s.find(delim, start);
    result.push_back(s.substr(start, end == std::string::npos ? end : end-start));
    if (end == std::string::npos) break;
    start = end + 1;
  }
  return result;
}
static std::string json_string(const std::string& s) {
  std::ostringstream out; out << '"';
  for (unsigned char c : s) {
    if (c == '"' || c == '\\') out << '\\' << c;
    else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else out << c;
  }
  out << '"'; return out.str();
}
static std::string json_fields(const Fields& f) {
  static const std::set<std::string> numeric = {"completed_files", "total_files", "sequences", "bases",
    "elapsed_seconds", "session_seconds", "hash_bytes", "capacity", "used_cells", "key_bits",
    "value_bits", "durable_files", "file_seconds", "checkpoint_seconds", "signal"};
  std::string out = "{\n"; size_t i = 0;
  for (const auto& kv : f) {
    out += "  " + json_string(kv.first) + ": ";
    out += (numeric.count(kv.first) || kv.first == "complete") ? kv.second : json_string(kv.second);
    out += (++i < f.size() ? ",\n" : "\n");
  }
  return out + "}\n";
}
static void write_all(int fd, const void* buf, size_t bytes) {
  const char* p = static_cast<const char*>(buf);
  while (bytes) {
    ssize_t n = write(fd, p, std::min(bytes, size_t(64) << 20));
    if (n < 0 && errno == EINTR) continue;
    require(n > 0, "checked write failed: " + std::string(strerror(errno)));
    p += n; bytes -= size_t(n);
  }
}
static void fsync_dir(const fs::path& path) {
  int fd = open(path.c_str(), O_RDONLY | O_DIRECTORY | O_NOFOLLOW);
  require(fd >= 0, "cannot open directory for fsync: " + path.string());
  int status = fsync(fd), saved = errno; close(fd);
  require(status == 0, "directory fsync failed: " + path.string() + ": " + strerror(saved));
}
static void new_file(const fs::path& path, const void* data, size_t bytes) {
  int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW, 0664);
  require(fd >= 0, "cannot create new file: " + path.string() + ": " + strerror(errno));
  try {
    write_all(fd, data, bytes);
    require(fsync(fd) == 0, "file fsync failed: " + path.string());
  } catch (...) { close(fd); throw; }
  require(close(fd) == 0, "file close failed: " + path.string());
}
static void new_file(const fs::path& path, const std::string& value) {
  new_file(path, value.data(), value.size());
}
static std::string nonce() {
  auto n = std::chrono::duration_cast<std::chrono::nanoseconds>(
    std::chrono::system_clock::now().time_since_epoch()).count();
  return std::to_string(n) + "-" + std::to_string(getpid());
}
static void atomic_file(const fs::path& path, const std::string& value) {
  fs::path temp = path.string() + ".tmp-" + nonce();
  new_file(temp, value);
  require(rename(temp.c_str(), path.c_str()) == 0, "atomic rename failed: " + path.string());
  fsync_dir(path.parent_path());
}
static std::string small_file(const fs::path& path, size_t limit = 1024*1024) {
  struct stat state;
  require(lstat(path.c_str(), &state) == 0 && S_ISREG(state.st_mode) && state.st_size >= 0 &&
          uint64_t(state.st_size) <= limit, "missing, symlink, or oversized metadata: " + path.string());
  std::ifstream input(path, std::ios::binary);
  require(bool(input), "cannot read " + path.string());
  std::ostringstream out; out << input.rdbuf();
  require(!input.bad() && out.str().size() == size_t(state.st_size), "short metadata read: " + path.string());
  return out.str();
}
static std::string sha256(const fs::path& path) {
  // No shell is involved: arbitrary filenames cannot introduce commands.
  // Never fork a parent retaining a TiB-scale mutable table. Even immediate
  // exec after fork can fail strict overcommit or copy enormous page tables.
  // posix_spawn uses libc's safe launch path; request its vfork implementation
  // explicitly on older GNU libc (including CentOS 7). No user code executes
  // in an unsafe hand-written vfork child.
  int p[2]; require(pipe(p) == 0, "cannot create SHA256 pipe");
  posix_spawn_file_actions_t actions;
  posix_spawnattr_t attributes;
  int code = posix_spawn_file_actions_init(&actions);
  if (code != 0) { close(p[0]); close(p[1]); fail("cannot initialize SHA256 spawn actions: " + std::string(strerror(code))); }
  bool attributes_ready = false;
  auto check_setup = [&](int status) {
    if (status == 0) return;
    posix_spawn_file_actions_destroy(&actions);
    if (attributes_ready) posix_spawnattr_destroy(&attributes);
    close(p[0]); close(p[1]);
    fail("cannot configure SHA256 spawn: " + std::string(strerror(status)));
  };
  check_setup(posix_spawnattr_init(&attributes)); attributes_ready = true;
  short spawn_flags = 0;
#ifdef POSIX_SPAWN_USEVFORK
  spawn_flags |= POSIX_SPAWN_USEVFORK;
#endif
  check_setup(posix_spawnattr_setflags(&attributes, spawn_flags));
  check_setup(posix_spawn_file_actions_addclose(&actions, p[0]));
  check_setup(posix_spawn_file_actions_adddup2(&actions, p[1], STDOUT_FILENO));
  if (p[1] != STDOUT_FILENO) check_setup(posix_spawn_file_actions_addclose(&actions, p[1]));
  char* child_args[] = {const_cast<char*>("sha256sum"), const_cast<char*>("--"),
                        const_cast<char*>(path.c_str()), nullptr};
  pid_t child = -1;
  code = posix_spawnp(&child, "sha256sum", &actions, &attributes, child_args, environ);
  posix_spawn_file_actions_destroy(&actions);
  posix_spawnattr_destroy(&attributes);
  if (code != 0) {
    close(p[0]); close(p[1]);
    fail("cannot spawn SHA256 checker without fork: " + std::string(strerror(code)));
  }
  close(p[1]); std::string result; char buf[4096];
  for (;;) {
    ssize_t n = read(p[0], buf, sizeof(buf));
    if (n < 0 && errno == EINTR) continue;
    if (n <= 0) { if (n < 0) result.clear(); break; }
    if (result.size() < 8192) result.append(buf, size_t(n));
  }
  close(p[0]); int status = 0; pid_t done;
  do { done = waitpid(child, &status, 0); } while (done < 0 && errno == EINTR);
  require(done == child && WIFEXITED(status) && WEXITSTATUS(status) == 0,
          "sha256sum failed: " + path.string());
  // GNU sha256sum prefixes an escape marker for backslashes in filenames.
  if (!result.empty() && result.front() == '\\') result.erase(0, 1);
  require(result.size() >= 66 && hex64(result.substr(0, 64)), "invalid SHA256 result");
  return result.substr(0, 64);
}
static std::string field_text(const Fields& fields) {
  std::string result;
  for (const auto& kv : fields) {
    require(kv.first.find_first_of("\t\n\r") == std::string::npos &&
            kv.second.find_first_of("\t\n\r") == std::string::npos, "invalid metadata field");
    result += kv.first + "\t" + kv.second + "\n";
  }
  return result;
}
static Fields parse_fields(const fs::path& path) {
  Fields fields; std::istringstream input(small_file(path)); std::string line;
  while (std::getline(input, line)) {
    auto row = split(line, '\t');
    require(row.size() == 2 && !row[0].empty() && fields.emplace(row[0], row[1]).second,
            "invalid receipt field: " + path.string());
  }
  return fields;
}
static std::string get(const Fields& fields, const std::string& key) {
  auto found = fields.find(key); require(found != fields.end(), "missing receipt field: " + key);
  return found->second;
}
struct Input {
  fs::path path; uint64_t bytes, records, bases; std::string digest; struct stat state;
};
static bool same_state(const struct stat& a, const struct stat& b) {
  return a.st_dev == b.st_dev && a.st_ino == b.st_ino && a.st_size == b.st_size &&
    a.st_mtim.tv_sec == b.st_mtim.tv_sec && a.st_mtim.tv_nsec == b.st_mtim.tv_nsec &&
    a.st_ctim.tv_sec == b.st_ctim.tv_sec && a.st_ctim.tv_nsec == b.st_ctim.tv_nsec;
}
static void unchanged(const Input& file) {
  struct stat now;
  require(stat(file.path.c_str(), &now) == 0 && same_state(file.state, now),
          "input changed during run: " + file.path.string());
}
static std::vector<Input> read_inputs(const fs::path& manifest) {
  std::vector<Input> inputs; std::set<std::string> seen;
  std::istringstream data(small_file(manifest, 256*1024*1024)); std::string line;
  while (std::getline(data, line)) {
    auto row = split(line, '\t');
    require(row.size() == 5, "inputs TSV needs path, bytes, SHA256, records, bases (no header)");
    Input i; i.path = row[0]; i.bytes = integer(row[1], "input bytes"); i.digest = row[2];
    i.records = integer(row[3], "input records"); i.bases = integer(row[4], "input bases");
    require(i.path.is_absolute() && hex64(i.digest) && i.records && i.bases,
            "inputs must be absolute nonempty FASTA paths with SHA256 and positive counts");
    require(row[0].find_first_of("\r\n") == std::string::npos && seen.insert(row[0]).second,
            "duplicate or invalid input path");
    inputs.push_back(i);
  }
  require(!inputs.empty(), "inputs manifest is empty"); return inputs;
}
static void verify_inputs(std::vector<Input>& inputs) {
  for (size_t n = 0; n < inputs.size(); ++n) {
    auto& file = inputs[n];
    require(stat(file.path.c_str(), &file.state) == 0 && S_ISREG(file.state.st_mode) &&
            uint64_t(file.state.st_size) == file.bytes, "input size/type mismatch: " + file.path.string());
    require(sha256(file.path) == file.digest, "input SHA256 mismatch: " + file.path.string());
    unchanged(file);
    std::cerr << "METAG_VERIFY " << n+1 << "/" << inputs.size() << " " << file.path << std::endl;
  }
}
struct Controls {
  fs::path inputs, checkpoint_dir; std::string context;
  uint64_t checkpoint_seconds = 14400, max_seconds = 396000, stop_after_files = 0;
};
static void parse(int argc, char** argv, Controls& extra, Options& opts) {
  std::vector<char*> native; native.push_back(argv[0]);
  std::set<std::string> flags;
  for (int i = 1; i < argc; ++i) {
    std::string flag = argv[i];
    if (flag.rfind("--metag-", 0) != 0) { native.push_back(argv[i]); continue; }
    require(i + 1 < argc && flags.insert(flag).second, "missing or duplicate option: " + flag);
    std::string value = argv[++i];
    if (flag == "--metag-inputs") extra.inputs = value;
    else if (flag == "--metag-context") extra.context = value;
    else if (flag == "--metag-checkpoint-dir") extra.checkpoint_dir = value;
    else if (flag == "--metag-checkpoint-seconds") extra.checkpoint_seconds = integer(value, flag);
    else if (flag == "--metag-max-seconds") extra.max_seconds = integer(value, flag);
    else if (flag == "--metag-stop-after-files") extra.stop_after_files = integer(value, flag);
    else fail("unknown metag option: " + flag);
  }
  require(!extra.inputs.empty() && !extra.checkpoint_dir.empty() && hex64(extra.context),
          "--metag-inputs, --metag-checkpoint-dir and lowercase 64-hex --metag-context are required");
  require(extra.checkpoint_seconds > 0 && extra.max_seconds > 0, "time budgets must be positive");
  opts = Options{};
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK; opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.num_threads = 1; opts.block_size = DEFAULT_BLOCK_SIZE; opts.subblock_size = DEFAULT_SUBBLOCK_SIZE;
  opts.cht_cell_size = 32; opts.deterministic_build = true;
  native.push_back(nullptr); optind = 1;
  ParseCommandLine(int(native.size())-1, native.data(), opts);
  require(optind == int(native.size())-1, "unexpected positional argument");
  require(opts.deterministic_build && !opts.input_is_protein && !opts.maximum_capacity &&
          opts.cht_cell_size == 32, "metag_build forbids -F, -X, -M and non-32-bit cells");
  require(sizeof(size_t) == 8 && opts.capacity <= (std::numeric_limits<uint64_t>::max()-32)/4,
          "64-bit host and valid capacity required");
  extra.checkpoint_dir = fs::absolute(extra.checkpoint_dir).lexically_normal();
}
struct Lock {
  int fd = -1;
  explicit Lock(const fs::path& root) {
    require(!fs::is_symlink(root), "checkpoint root must not be symlink");
    fs::create_directories(root);
    fd = open((root / ".LOCK").c_str(), O_RDWR | O_CREAT | O_NOFOLLOW, 0664);
    require(fd >= 0 && flock(fd, LOCK_EX | LOCK_NB) == 0, "checkpoint directory is locked or inaccessible");
  }
  ~Lock() { if (fd >= 0) { flock(fd, LOCK_UN); close(fd); } }
};
static bool generation_name(const std::string& s) {
  return std::regex_match(s, std::regex("cp-[0-9]+-[0-9]+"));
}
static bool stage_name(const std::string& s) {
  return std::regex_match(s, std::regex("stage-[0-9]+-[0-9]+"));
}
static void own_root(const fs::path& root, const std::string& identity_text) {
  fs::path owner = root / "OWNER.tsv";
  if (fs::exists(owner)) require(small_file(owner) == identity_text, "checkpoint owner/semantic identity mismatch");
  else {
    for (const auto& item : fs::directory_iterator(root))
      require(item.path().filename() == ".LOCK", "refusing unowned nonempty checkpoint directory");
    new_file(owner, identity_text); fsync_dir(root);
  }
}
static void prune(const fs::path& root, const std::string& current) {
  for (const auto& item : fs::directory_iterator(root)) {
    const auto name = item.path().filename().string();
    if ((generation_name(name) && name != current) || stage_name(name)) {
      require(!item.is_symlink() && item.is_directory(), "unexpected owned checkpoint path type");
      fs::remove_all(item.path());
    }
  }
  fsync_dir(root);
}
static uint64_t hash_bytes(size_t capacity) { return 32 + uint64_t(capacity)*4; }
static void hash_header(const fs::path& path, const Options& opts, size_t value_bits, uint64_t used) {
  struct stat state;
  require(lstat(path.c_str(), &state) == 0 && S_ISREG(state.st_mode) &&
          uint64_t(state.st_size) == hash_bytes(opts.capacity), "checkpoint hash length/type mismatch");
  size_t h[4] = {};
  std::ifstream input(path, std::ios::binary); input.read(reinterpret_cast<char*>(h), sizeof(h));
  require(input.gcount() == sizeof(h) && h[0] == opts.capacity && h[1] == used && used <= h[0] &&
          h[2] == 32-value_bits && h[3] == value_bits, "checkpoint hash header mismatch");
}
static void atomic_link(const fs::path& source, const fs::path& output) {
  require(!fs::is_symlink(output), "output must not be a symlink: " + output.string());
  fs::create_directories(output.parent_path());
  fs::path temp = output.string() + ".tmp-" + nonce();
  // Intentionally no silent multi-terabyte copy. Put final DB and checkpoints on
  // the same persistent filesystem; hardlink is atomic and requires no payload copy.
  require(link(source.c_str(), temp.c_str()) == 0,
          "cannot hardlink final hash; final DB and checkpoints must share filesystem: " + std::string(strerror(errno)));
  require(rename(temp.c_str(), output.c_str()) == 0, "cannot publish final hash");
  fsync_dir(output.parent_path());
}
static void progress(const Controls& c, Fields data, size_t durable, double session) {
  data["durable_files"] = number(durable); data["session_seconds"] = seconds(session);
  atomic_file(c.checkpoint_dir / "progress.json", json_fields(data));
  fs::path csv = c.checkpoint_dir / "progress.csv";
  bool empty = !fs::exists(csv) || fs::file_size(csv) == 0;
  int fd = open(csv.c_str(), O_WRONLY | O_APPEND | O_CREAT | O_NOFOLLOW, 0664);
  require(fd >= 0, "cannot append progress CSV");
  std::string line = empty ? "phase,completed_files,total_files,durable_files,sequences,bases,session_seconds,elapsed_seconds\n" : "";
  line += get(data, "phase") + "," + get(data, "completed_files") + "," + get(data, "total_files") + "," +
    number(durable) + "," + get(data, "sequences") + "," + get(data, "bases") + "," +
    seconds(session) + "," + get(data, "elapsed_seconds") + "\n";
  try { write_all(fd, line.data(), line.size()); require(fsync(fd) == 0, "CSV fsync failed"); }
  catch (...) { close(fd); throw; }
  require(close(fd) == 0, "CSV close failed");
}
static std::string checkpoint(const Controls& c, const Options& opts, size_t value_bits,
                              CompactHashTable<CompactHashCell>& table, Fields& fields) {
  auto start = Clock::now();
  struct statvfs free_space;
  require(statvfs(c.checkpoint_dir.c_str(), &free_space) == 0, "cannot inspect checkpoint free space");
  uint64_t available = uint64_t(free_space.f_bavail) * uint64_t(free_space.f_frsize);
  require(available >= hash_bytes(opts.capacity) + (8ULL<<20),
          "not enough free space for a new snapshot while retaining CURRENT; need bytes=" + number(hash_bytes(opts.capacity)));
  std::string suffix = nonce(), name = "cp-" + suffix;
  fs::path stage = c.checkpoint_dir / ("stage-" + suffix);
  require(fs::create_directory(stage), "cannot create checkpoint generation");
  std::cerr << "METAG_CHECKPOINT_BEGIN files=" << fields["completed_files"] << " bytes=" << hash_bytes(opts.capacity) << std::endl;
  table.WriteTable((stage / "hash.k2d").c_str());
  hash_header(stage / "hash.k2d", opts, value_bits, table.size());
  fields["hash_sha256"] = sha256(stage / "hash.k2d");
  fields["hash_bytes"] = number(hash_bytes(opts.capacity)); fields["capacity"] = number(table.capacity());
  fields["used_cells"] = number(table.size()); fields["key_bits"] = number(table.key_bits());
  fields["value_bits"] = number(table.value_bits()); fields["generation"] = name;
  fields["checkpoint_seconds"] = seconds(elapsed(start));
  fields["elapsed_seconds"] = seconds(std::stod(fields["elapsed_seconds"]) + elapsed(start));
  new_file(stage / "receipt.json", json_fields(fields));
  Fields canonical = fields; canonical["receipt_json_sha256"] = sha256(stage / "receipt.json");
  new_file(stage / "receipt.tsv", field_text(canonical));
  std::string receipt_hash = sha256(stage / "receipt.tsv");
  fsync_dir(stage);
  fs::path committed = c.checkpoint_dir / name;
  require(rename(stage.c_str(), committed.c_str()) == 0, "cannot finalize checkpoint generation");
  fsync_dir(c.checkpoint_dir);
  atomic_file(c.checkpoint_dir / "CURRENT", name + "\t" + receipt_hash + "\n");
  // The new hash, checksums, receipt and pointer are now durable. Only now may
  // the previous payload be reclaimed. A failed/short write never reaches here.
  prune(c.checkpoint_dir, name);
  std::cerr << "METAG_CHECKPOINT_COMMIT files=" << fields["completed_files"] << " sha256=" << fields["hash_sha256"]
            << " seconds=" << elapsed(start) << std::endl;
  return name;
}
static Fields resume(const Controls& c, const Options& opts, size_t value_bits,
                     const Fields& expected, const std::vector<Input>& inputs,
                     const std::string& recovery = "") {
  auto pointer = recovery.empty() ? split(small_file(c.checkpoint_dir / "CURRENT"), '\t') :
    std::vector<std::string>{recovery, sha256(c.checkpoint_dir/recovery/"receipt.tsv") + "\n"};
  require(pointer.size() == 2 && generation_name(pointer[0]) && pointer[1].size() == 65 && pointer[1].back() == '\n',
          "invalid CURRENT pointer");
  pointer[1].pop_back(); require(hex64(pointer[1]), "invalid CURRENT receipt digest");
  fs::path dir = c.checkpoint_dir / pointer[0];
  require(fs::is_directory(dir) && !fs::is_symlink(dir), "invalid checkpoint generation directory");
  require(sha256(dir / "receipt.tsv") == pointer[1], "checkpoint receipt SHA256 mismatch");
  Fields fields = parse_fields(dir / "receipt.tsv");
  require(sha256(dir / "receipt.json") == get(fields, "receipt_json_sha256"), "checkpoint JSON SHA256 mismatch");
  for (const auto& kv : expected)
    require(get(fields, kv.first) == kv.second, "checkpoint semantic mismatch: " + kv.first);
  require(get(fields, "generation") == pointer[0], "checkpoint generation mismatch");
  uint64_t done = integer(get(fields, "completed_files"), "completed_files"), records = 0, bases = 0;
  require(done <= inputs.size() && integer(get(fields, "total_files"), "total_files") == inputs.size(), "checkpoint file count mismatch");
  for (size_t i = 0; i < done; ++i) { records += inputs[i].records; bases += inputs[i].bases; }
  require(integer(get(fields, "sequences"), "sequences") == records && integer(get(fields, "bases"), "bases") == bases,
          "checkpoint cumulative record/base counts mismatch");
  require(get(fields, "complete") == (done == inputs.size() ? "true" : "false"), "checkpoint completion flag mismatch");
  require(integer(get(fields, "capacity"), "capacity") == opts.capacity &&
          integer(get(fields, "key_bits"), "key_bits") == 32-value_bits &&
          integer(get(fields, "value_bits"), "value_bits") == value_bits &&
          integer(get(fields, "hash_bytes"), "hash_bytes") == hash_bytes(opts.capacity), "checkpoint table metadata mismatch");
  hash_header(dir / "hash.k2d", opts, value_bits, integer(get(fields, "used_cells"), "used_cells"));
  require(sha256(dir / "hash.k2d") == get(fields, "hash_sha256"), "checkpoint hash SHA256 mismatch; refusing restart from zero");
  double old_elapsed = std::stod(get(fields, "elapsed_seconds"));
  require(std::isfinite(old_elapsed) && old_elapsed >= 0, "invalid checkpoint elapsed time");
  fields.erase("receipt_json_sha256");
  prune(c.checkpoint_dir, pointer[0]);
  return fields;
}
static std::pair<uint64_t, uint64_t> process_file(const Input& file, const Options& opts,
  const std::map<std::string, taxid_t>& ids, CompactHashTable<CompactHashCell>& table, const Taxonomy& taxonomy) {
  unchanged(file); uint64_t records = 0, bases = 0;
  BatchSequenceReader reader(file.path.c_str()); Sequence* seq;
  while (reader.LoadBlock(opts.block_size)) {
    while ((seq = reader.NextSequence()) != nullptr) {
      require(seq->format == SequenceFormat::FORMAT_FASTA, "prepared input is not FASTA: " + file.path.string());
      taxid_t taxid = 0;
      for (const auto& id : ExtractNCBISequenceIDs(seq->header)) {
        auto found = ids.find(id);
        if (found != ids.end())
          taxid = taxonomy.LowestCommonAncestor(taxid, taxonomy.GetInternalID(found->second));
      }
      require(taxid != 0, "FASTA sequence lacks valid global-map taxid: " + seq->header);
      ProcessSequence(opts, seq->seq, taxid, table, taxonomy);
      ++records; bases += seq->seq.size();
    }
  }
  unchanged(file);
  require(records == file.records && bases == file.bases, "FASTA record/base counts differ from manifest: " + file.path.string());
  return {records, bases};
}

int run(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--metag-version") {
    std::cout << SCHEMA << " kraken-source=8c190b1b668825935dbf6dee5f969227dc8269bb\n"; return 0;
  }
  auto started = Clock::now();
  struct sigaction action {}; action.sa_handler = on_signal; sigemptyset(&action.sa_mask);
  for (int sig : {SIGUSR1, SIGUSR2, SIGTERM}) require(sigaction(sig, &action, nullptr) == 0, "cannot install soft-stop signal handler");
  Controls c; Options opts; parse(argc, argv, c, opts);
  omp_set_num_threads(opts.num_threads);
  Lock lock(c.checkpoint_dir);
  Fields identity = {{"schema", SCHEMA}, {"context", c.context},
    {"input_manifest_sha256", sha256(c.inputs)}, {"native_binary_sha256", sha256(fs::read_symlink("/proc/self/exe"))},
    {"mapping_sha256", sha256(opts.ID_to_taxon_map_filename)},
    {"nodes_sha256", sha256(fs::path(opts.ncbi_taxonomy_directory)/"nodes.dmp")},
    {"names_sha256", sha256(fs::path(opts.ncbi_taxonomy_directory)/"names.dmp")},
    {"k", number(opts.k)}, {"l", number(opts.l)}, {"spaced_seed_mask", number(opts.spaced_seed_mask)},
    {"toggle_mask", number(opts.toggle_mask)}, {"capacity", number(opts.capacity)},
    {"threads", number(opts.num_threads)}, {"block_size", number(opts.block_size)},
    {"subblock_size", number(opts.subblock_size)}, {"requested_value_bits", number(opts.requested_bits_for_taxid)}};
  own_root(c.checkpoint_dir, field_text(identity));
  identity["identity_sha256"] = sha256(c.checkpoint_dir / "OWNER.tsv");
  auto inputs = read_inputs(c.inputs);
  verify_inputs(inputs);
  require(sha256(c.inputs) == identity["input_manifest_sha256"], "input manifest changed during verification");
  std::map<std::string, taxid_t> ids;
  ReadIDToTaxonMap(ids, opts.ID_to_taxon_map_filename);
  require(!ids.empty(), "global sequence-to-taxid map is empty");
  fs::create_directories(fs::absolute(opts.taxonomy_filename).parent_path());
  fs::create_directories(fs::absolute(opts.options_filename).parent_path());
  std::string final_taxonomy = opts.taxonomy_filename;
  opts.taxonomy_filename += ".tmp-" + nonce();
  GenerateTaxonomy(opts, ids);
  {
    int fd = open(opts.taxonomy_filename.c_str(), O_RDONLY | O_NOFOLLOW);
    require(fd >= 0 && fsync(fd) == 0, "cannot sync generated taxonomy"); close(fd);
  }
  Taxonomy taxonomy(opts.taxonomy_filename.c_str(), false);
  taxonomy.GenerateExternalToInternalIDMap();
  size_t value_bits = 1;
  while ((uint64_t(1) << value_bits) < taxonomy.node_count()) ++value_bits;
  require(value_bits < 32 && (!opts.requested_bits_for_taxid || value_bits <= opts.requested_bits_for_taxid), "taxonomy needs more value bits");
  value_bits = std::max(value_bits, opts.requested_bits_for_taxid);
  identity["taxonomy_sha256"] = sha256(opts.taxonomy_filename);
  require(rename(opts.taxonomy_filename.c_str(), final_taxonomy.c_str()) == 0, "cannot publish taxonomy");
  opts.taxonomy_filename = final_taxonomy;
  fsync_dir(fs::absolute(final_taxonomy).parent_path());
  IndexOptions index_opts = {}; // Zero padding: independent categories have byte-identical opts.
  index_opts.k = opts.k; index_opts.l = opts.l; index_opts.spaced_seed_mask = opts.spaced_seed_mask;
  index_opts.toggle_mask = opts.toggle_mask; index_opts.dna_db = true;
  index_opts.minimum_acceptable_hash_value = 0; index_opts.revcom_version = CURRENT_REVCOM_VERSION;
  index_opts.db_version = 0; index_opts.db_type = 0;
  std::string index_blob(reinterpret_cast<char*>(&index_opts), sizeof(index_opts));
  atomic_file(fs::absolute(opts.options_filename), index_blob);
  identity["options_sha256"] = sha256(opts.options_filename);
  Fields fields = identity; fields["total_files"] = number(inputs.size());
  size_t done = 0, durable = 0; uint64_t sequences = 0, bases = 0; double prior_elapsed = 0;
  std::unique_ptr<CompactHashTable<CompactHashCell>> table;
  std::string current;
  // A kill can fall between publishing the very first immutable generation and
  // publishing CURRENT. Recover that one fully validated generation; never
  // treat a corrupted/missing CURRENT as permission to start from zero.
  if (!fs::exists(c.checkpoint_dir / "CURRENT")) {
    std::vector<std::string> orphan_generations;
    for (const auto& item : fs::directory_iterator(c.checkpoint_dir))
      if (generation_name(item.path().filename().string())) orphan_generations.push_back(item.path().filename().string());
    if (!orphan_generations.empty()) {
      require(orphan_generations.size() == 1, "CURRENT missing with multiple generations; investigate before restarting");
      const auto name = orphan_generations.front();
      resume(c, opts, value_bits, identity, inputs, name);
      atomic_file(c.checkpoint_dir/"CURRENT", name + "\t" + sha256(c.checkpoint_dir/name/"receipt.tsv") + "\n");
      std::cerr << "METAG_RECOVER_FIRST_CURRENT generation=" << name << std::endl;
    }
  }
  if (fs::exists(c.checkpoint_dir / "CURRENT")) {
    fields = resume(c, opts, value_bits, identity, inputs);
    done = integer(fields["completed_files"], "completed_files"); durable = done;
    sequences = integer(fields["sequences"], "sequences"); bases = integer(fields["bases"], "bases");
    prior_elapsed = std::stod(fields["elapsed_seconds"]); current = fields["generation"];
    if (done == inputs.size()) {
      atomic_link(c.checkpoint_dir/current/"hash.k2d", fs::absolute(opts.hashtable_filename));
      std::cerr << "METAG_COMPLETE_REUSE files=" << done << std::endl; return 0;
    }
    if (stop_signal || elapsed(started) >= c.max_seconds) {
      std::cerr << "METAG_PAUSED durable_files=" << durable << " before hash load; reuse CURRENT on resubmission\n";
      return 75;
    }
    table.reset(new CompactHashTable<CompactHashCell>((c.checkpoint_dir/current/"hash.k2d").string(), false));
    std::cerr << "METAG_RESUME files=" << done << "/" << inputs.size() << " reused_sequences=" << sequences << std::endl;
  } else {
    // Ownership was established before this deletion; CURRENT absence is not a
    // license to silently discard a generation previously published as valid.
    for (const auto& item : fs::directory_iterator(c.checkpoint_dir))
      require(!generation_name(item.path().filename().string()), "CURRENT missing but committed generation exists; investigate before restarting");
    prune(c.checkpoint_dir, "");
    table.reset(new CompactHashTable<CompactHashCell>(opts.capacity, 32-value_bits, value_bits));
    std::cerr << "METAG_NEW capacity=" << opts.capacity << " hash_bytes=" << hash_bytes(opts.capacity) << std::endl;
  }
  const size_t session_start_file = done;
  auto last_checkpoint = Clock::now();
  for (;;) {
    bool stopped = stop_signal || elapsed(started) >= c.max_seconds ||
      (c.stop_after_files && done - session_start_file >= c.stop_after_files);
    if (done < inputs.size() && !stopped) {
      std::cerr << "METAG_FILE_BEGIN " << done+1 << "/" << inputs.size() << " path=" << inputs[done].path << std::endl;
      auto file_started = Clock::now(); auto counts = process_file(inputs[done], opts, ids, *table, taxonomy);
      ++done; sequences += counts.first; bases += counts.second;
      fields["file_seconds"] = seconds(elapsed(file_started));
      std::cerr << "METAG_FILE_DONE " << done << "/" << inputs.size() << " sequences=" << counts.first
                << " bases=" << counts.second << " seconds=" << fields["file_seconds"] << std::endl;
    }
    fields["completed_files"] = number(done); fields["total_files"] = number(inputs.size());
    fields["sequences"] = number(sequences); fields["bases"] = number(bases);
    fields["elapsed_seconds"] = seconds(prior_elapsed + elapsed(started));
    fields["complete"] = done == inputs.size() ? "true" : "false";
    fields["signal"] = number(stop_signal); fields["phase"] = "building";
    stopped = stop_signal || elapsed(started) >= c.max_seconds ||
      (c.stop_after_files && done - session_start_file >= c.stop_after_files);
    bool complete = done == inputs.size();
    if (complete || stopped || elapsed(last_checkpoint) >= c.checkpoint_seconds) {
      if (done != durable || current.empty()) {
        fields["phase"] = complete ? "complete" : (stopped ? "paused" : "checkpoint");
        current = checkpoint(c, opts, value_bits, *table, fields);
        durable = done; last_checkpoint = Clock::now();
      }
    }
    fields["elapsed_seconds"] = seconds(prior_elapsed + elapsed(started));
    fields["phase"] = complete ? "complete" : (stopped ? "paused" : "building");
    progress(c, fields, durable, elapsed(started));
    if (complete) {
      atomic_link(c.checkpoint_dir/current/"hash.k2d", fs::absolute(opts.hashtable_filename));
      std::cerr << "METAG_COMPLETE files=" << done << " sequences=" << sequences << " bases=" << bases << std::endl;
      return 0;
    }
    if (stopped) {
      std::cerr << "METAG_PAUSED durable_files=" << durable << " signal=" << stop_signal << " resubmit identical inputs/options to resume\n";
      return 75;
    }
  }
}
} // namespace metag

int main(int argc, char** argv) {
  try { return metag::run(argc, argv); }
  catch (const std::exception& error) {
    std::cerr << "METAG_ERROR " << error.what() << "\nExisting published checkpoints were not reset.\n";
    return 1;
  }
}
