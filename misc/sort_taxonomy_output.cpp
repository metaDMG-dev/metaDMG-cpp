#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include <zlib.h>

// New functionality: sort taxonomy-aware text output in depth-first order.
// Please report bugs, unexpected ordering, or unsupported output formats.

namespace {

struct Row {
    std::string text;
    int taxid;
    std::size_t input_order;
    std::size_t sort_order;
};

struct Taxonomy {
    std::map<int, std::vector<int> > children;
    std::map<int, int> parent;
    std::map<int, std::size_t> order;
};

static std::vector<std::string> split_tab(const std::string &line) {
    std::vector<std::string> result;
    std::size_t start = 0;
    while (true) {
        const std::size_t end = line.find('\t', start);
        result.push_back(line.substr(start, end == std::string::npos ? end : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return result;
}

static bool parse_int(const std::string &value, int &out) {
    char *end = NULL;
    errno = 0;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if (value.empty() || end == value.c_str() || *end != '\0' || errno == ERANGE ||
        parsed < std::numeric_limits<int>::min() || parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    out = static_cast<int>(parsed);
    return true;
}

static std::string trim(const std::string &value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return std::string();
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

static bool gz_getline(gzFile file, std::string &line) {
    char buffer[8192];
    line.clear();
    while (true) {
        char *result = gzgets(file, buffer, sizeof(buffer));
        if (result == NULL) return !line.empty();
        line += result;
        if (!line.empty() && line[line.size() - 1] == '\n') {
            line.erase(line.size() - 1);
            if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
            return true;
        }
        if (gzeof(file)) return true;
    }
}

static bool read_text_lines(const std::string &path, std::vector<std::string> &lines) {
    std::ifstream probe(path.c_str(), std::ios::binary);
    if (!probe) {
        fprintf(stderr, "Error: cannot open %s\n", path.c_str());
        return false;
    }
    unsigned char magic[2] = {0, 0};
    probe.read(reinterpret_cast<char *>(magic), 2);
    const bool compressed = probe.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
    probe.close();
    if (!compressed) {
        std::ifstream plain(path.c_str());
        std::string line;
        while (std::getline(plain, line)) lines.push_back(line);
        return plain.good() || plain.eof();
    }
    gzFile file = gzopen(path.c_str(), "rb");
    if (file == NULL) {
        fprintf(stderr, "Error: cannot open %s\n", path.c_str());
        return false;
    }
    std::string line;
    while (gz_getline(file, line)) lines.push_back(line);
    const int close_status = gzclose(file);
    if (close_status != Z_OK) {
        fprintf(stderr, "Error: failed reading compressed file %s\n", path.c_str());
        return false;
    }
    return true;
}

class TextReader {
public:
    explicit TextReader(const std::string &path) : gz_(NULL), plain_(), compressed_(false) {
        std::ifstream probe(path.c_str(), std::ios::binary);
        unsigned char magic[2] = {0, 0};
        probe.read(reinterpret_cast<char *>(magic), 2);
        compressed_ = probe.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
        if (compressed_) gz_ = gzopen(path.c_str(), "rb");
        else plain_.open(path.c_str());
    }
    ~TextReader() { if (gz_ != NULL) gzclose(gz_); }
    bool good() const { return compressed_ ? gz_ != NULL : plain_.good(); }
    bool getline(std::string &line) {
        if (compressed_) return gz_getline(gz_, line);
        return static_cast<bool>(std::getline(plain_, line));
    }
private:
    gzFile gz_;
    std::ifstream plain_;
    bool compressed_;
};

class TextWriter {
public:
    explicit TextWriter(const std::string &path) : gz_(NULL), plain_(), compressed_(path.size() >= 3 && path.substr(path.size() - 3) == ".gz") {
        if (compressed_) gz_ = gzopen(path.c_str(), "wb");
        else plain_.open(path.c_str());
    }
    ~TextWriter() { if (gz_ != NULL) gzclose(gz_); }
    bool good() const { return compressed_ ? gz_ != NULL : plain_.good(); }
    bool write(const std::string &line) {
        if (compressed_) return gzputs(gz_, line.c_str()) >= 0 && gzputc(gz_, '\n') >= 0;
        plain_ << line << '\n';
        return plain_.good();
    }
private:
    gzFile gz_;
    std::ofstream plain_;
    bool compressed_;
};

static bool write_chunk(const std::string &path, std::vector<Row> &rows) {
    struct RowOrder {
        bool operator()(const Row &a, const Row &b) const {
            return a.sort_order != b.sort_order ? a.sort_order < b.sort_order : a.input_order < b.input_order;
        }
    };
    std::sort(rows.begin(), rows.end(), RowOrder());
    std::ofstream file(path.c_str());
    if (!file) return false;
    for (std::size_t i = 0; i < rows.size(); ++i)
        file << rows[i].sort_order << '\t' << rows[i].input_order << '\t' << rows[i].text << '\n';
    rows.clear();
    return file.good();
}

struct MergeItem {
    std::size_t sort_order;
    std::size_t input_order;
    std::size_t file_index;
    std::string text;
};

struct MergeItemGreater {
    bool operator()(const MergeItem &a, const MergeItem &b) const {
        if (a.sort_order != b.sort_order) return a.sort_order > b.sort_order;
        return a.input_order > b.input_order;
    }
};

static bool read_chunk_item(std::ifstream &file, std::size_t index, MergeItem &item) {
    std::string line;
    if (!std::getline(file, line)) return false;
    const std::size_t first_tab = line.find('\t');
    const std::size_t second_tab = line.find('\t', first_tab + 1);
    if (first_tab == std::string::npos || second_tab == std::string::npos) return false;
    int sort_int = 0, input_int = 0;
    if (!parse_int(line.substr(0, first_tab), sort_int) || !parse_int(line.substr(first_tab + 1, second_tab - first_tab - 1), input_int)) return false;
    item.sort_order = static_cast<std::size_t>(sort_int);
    item.input_order = static_cast<std::size_t>(input_int);
    item.file_index = index;
    item.text = line.substr(second_tab + 1);
    return true;
}

static bool read_nodes(const std::string &path, Taxonomy &taxonomy) {
    std::vector<std::string> lines;
    if (!read_text_lines(path, lines)) return false;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::vector<std::string> fields = split_tab(lines[i]);
        // NCBI nodes.dmp is pipe-delimited; tolerate surrounding spaces.
        std::vector<std::string> pipe_fields;
        std::string field;
        std::stringstream stream(lines[i]);
        while (std::getline(stream, field, '|')) pipe_fields.push_back(field);
        if (pipe_fields.size() < 2) continue;
        const std::string id_string = trim(pipe_fields[0]);
        const std::string parent_string = trim(pipe_fields[1]);
        int id = 0, parent = 0;
        if (!parse_int(id_string, id) || !parse_int(parent_string, parent)) continue;
        taxonomy.children[parent].push_back(id);
        taxonomy.parent[id] = parent;
        // Keep a node present even when it has no children.
        if (taxonomy.children.find(id) == taxonomy.children.end()) taxonomy.children[id];
    }
    if (taxonomy.children.empty()) {
        fprintf(stderr, "Error: no taxonomy records found in %s\n", path.c_str());
        return false;
    }
    for (std::map<int, std::vector<int> >::iterator it = taxonomy.children.begin();
         it != taxonomy.children.end(); ++it) {
        std::sort(it->second.begin(), it->second.end());
    }
    std::vector<int> roots;
    for (std::map<int, std::vector<int> >::const_iterator it = taxonomy.children.begin();
         it != taxonomy.children.end(); ++it) {
        const std::map<int, int>::const_iterator parent = taxonomy.parent.find(it->first);
        if (parent == taxonomy.parent.end() || parent->second == it->first) {
            roots.push_back(it->first);
            continue;
        }
    }
    std::sort(roots.begin(), roots.end());
    std::size_t next = 0;
    std::vector<int> stack(roots.rbegin(), roots.rend());
    while (!stack.empty()) {
        const int node = stack.back();
        stack.pop_back();
        if (taxonomy.order.find(node) != taxonomy.order.end()) continue;
        taxonomy.order[node] = next++;
        const std::vector<int> &kids = taxonomy.children[node];
        for (std::vector<int>::const_reverse_iterator kid = kids.rbegin(); kid != kids.rend(); ++kid)
            stack.push_back(*kid);
    }
    return true;
}

enum Format { FORMAT_LCA, FORMAT_DFIT, FORMAT_STAT };

static bool detect_format(const std::string &header, Format &format, std::size_t &taxid_column) {
    const std::vector<std::string> fields = split_tab(header);
    if (fields.size() >= 9 && fields[0] == "queryid" && fields[7] == "lca") {
        format = FORMAT_LCA;
        taxid_column = 7;
        return true;
    }
    if (!fields.empty() && fields[0] == "taxid" && fields.size() >= 4 && fields[1] == "A") {
        format = FORMAT_DFIT;
        taxid_column = 0;
        return true;
    }
    if (fields.size() >= 3 && fields[0] == "taxid" && fields[1] == "name" && fields[2] == "rank") {
        format = FORMAT_STAT;
        taxid_column = 0;
        return true;
    }
    return false;
}

static void usage() {
    fprintf(stderr, "Usage: sort_taxonomy_output --input FILE --nodes nodes.dmp[.gz] [--output FILE]\n"
            "       [--format auto|lca|dfit|stat]\n"
            "Sort text taxonomy outputs in depth-first taxonomic order.\n");
}

} // namespace

int main(int argc, char **argv) {
    fprintf(stderr, "[new functionality] Sorting taxonomy output in depth-first order. "
            "Please report bugs or unexpected ordering.\n");
    std::string input, output, nodes, format_name = "auto";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if ((arg == "--input" || arg == "-i") && i + 1 < argc) input = argv[++i];
        else if ((arg == "--output" || arg == "-o") && i + 1 < argc) output = argv[++i];
        else if ((arg == "--nodes" || arg == "-n") && i + 1 < argc) nodes = argv[++i];
        else if (arg == "--format" && i + 1 < argc) format_name = argv[++i];
        else { usage(); return 2; }
    }
    if (input.empty() || nodes.empty()) { usage(); return 2; }
    if (output.empty()) output = input + ".sorted.gz";
    if (input.size() >= 9 && input.substr(input.size() - 9) == ".bdamage.gz") {
        fprintf(stderr, "Error: .bdamage.gz is a binary calculation input and cannot be sorted as text\n");
        return 1;
    }

    TextReader reader(input);
    if (!reader.good()) {
        fprintf(stderr, "Error: cannot open %s\n", input.c_str());
        return 1;
    }
    std::vector<std::string> prefix;
    std::string line;
    std::size_t taxid_column = 0;
    Format format = FORMAT_STAT;
    bool header_found = false;
    while (reader.getline(line)) {
        if (!line.empty() && line[0] != '#' && detect_format(line, format, taxid_column)) {
            prefix.push_back(line);
            header_found = true;
            break;
        }
        prefix.push_back(line);
    }
    if (!header_found && format_name != "auto") {
        if (prefix.empty()) {
            fprintf(stderr, "Error: input is empty\n");
            return 1;
        }
        header_found = true;
    }
    if (!header_found) {
        fprintf(stderr, "Error: could not autodetect a supported output header; use --format only for a custom header\n");
        return 1;
    }
    if (format_name != "auto") {
        if (format_name == "lca") { format = FORMAT_LCA; taxid_column = 7; }
        else if (format_name == "dfit") { format = FORMAT_DFIT; taxid_column = 0; }
        else if (format_name == "stat" || format_name == "aggregate") { format = FORMAT_STAT; taxid_column = 0; }
        else { fprintf(stderr, "Error: unknown format %s\n", format_name.c_str()); return 2; }
    }

    Taxonomy taxonomy;
    if (!read_nodes(nodes, taxonomy)) return 1;
    std::map<int, std::size_t> unknown;
    std::size_t next = taxonomy.order.size();
    const std::size_t chunk_size = 100000;
    std::vector<Row> rows;
    std::vector<std::string> chunk_paths;
    std::size_t input_order = 0;
    const char *tmpdir = getenv("TMPDIR");
    if (tmpdir == NULL || tmpdir[0] == '\0') tmpdir = "/tmp";
    char temp_dir_template[512];
    snprintf(temp_dir_template, sizeof(temp_dir_template), "%s/sort_taxonomy_output.XXXXXX", tmpdir);
    char *temp_dir = mkdtemp(temp_dir_template);
    if (temp_dir == NULL) {
        fprintf(stderr, "Error: cannot create a secure temporary directory in %s\n", tmpdir);
        return 1;
    }
    const std::string temp_prefix = std::string(temp_dir) + "/chunk";
    while (reader.getline(line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> fields = split_tab(line);
        if (taxid_column >= fields.size()) {
            fprintf(stderr, "Error: row %lu has too few columns\n", static_cast<unsigned long>(input_order + 1));
            return 1;
        }
        int taxid = 0;
        if (!parse_int(fields[taxid_column], taxid)) {
            fprintf(stderr, "Error: non-numeric taxid in row %lu: %s\n",
                    static_cast<unsigned long>(input_order + 1), fields[taxid_column].c_str());
            return 1;
        }
        if (taxonomy.order.find(taxid) == taxonomy.order.end() && unknown.find(taxid) == unknown.end()) unknown[taxid] = next++;
        const std::size_t sort_order = taxonomy.order.find(taxid) != taxonomy.order.end() ? taxonomy.order[taxid] : unknown[taxid];
        Row row = {line, taxid, input_order++, sort_order};
        rows.push_back(row);
        if (rows.size() >= chunk_size) {
            char path_buffer[512];
            snprintf(path_buffer, sizeof(path_buffer), "%s.%lu", temp_prefix.c_str(),
                     (unsigned long)chunk_paths.size());
            const std::string path = path_buffer;
            if (!write_chunk(path, rows)) {
                fprintf(stderr, "Error: cannot write temporary sort file %s\n", path.c_str());
                rmdir(temp_dir);
                return 1;
            }
            chunk_paths.push_back(path);
        }
    }
    if (!rows.empty()) {
        char path_buffer[512];
        snprintf(path_buffer, sizeof(path_buffer), "%s.%lu", temp_prefix.c_str(),
                 (unsigned long)chunk_paths.size());
        const std::string path = path_buffer;
        if (!write_chunk(path, rows)) {
            fprintf(stderr, "Error: cannot write temporary sort file %s\n", path.c_str());
            rmdir(temp_dir);
            return 1;
        }
        chunk_paths.push_back(path);
    }

    TextWriter writer(output);
    if (!writer.good()) {
        fprintf(stderr, "Error: cannot create %s\n", output.c_str());
        for (std::size_t i = 0; i < chunk_paths.size(); ++i) std::remove(chunk_paths[i].c_str());
        rmdir(temp_dir);
        return 1;
    }
    for (std::size_t i = 0; i < prefix.size(); ++i) writer.write(prefix[i]);
    std::vector<std::unique_ptr<std::ifstream> > files;
    std::priority_queue<MergeItem, std::vector<MergeItem>, MergeItemGreater> queue;
    for (std::size_t i = 0; i < chunk_paths.size(); ++i) {
        files.push_back(std::unique_ptr<std::ifstream>(new std::ifstream(chunk_paths[i].c_str())));
        MergeItem item;
        if (read_chunk_item(*files.back(), i, item)) queue.push(item);
    }
    bool success = true;
    while (!queue.empty()) {
        const MergeItem item = queue.top();
        queue.pop();
        if (!writer.write(item.text)) { success = false; break; }
        MergeItem next_item;
        if (read_chunk_item(*files[item.file_index], item.file_index, next_item)) queue.push(next_item);
    }
    for (std::size_t i = 0; i < chunk_paths.size(); ++i) std::remove(chunk_paths[i].c_str());
    rmdir(temp_dir);
    return success ? 0 : 1;
}
