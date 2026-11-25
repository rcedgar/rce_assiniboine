#include "naview.hpp"

#include "asserta.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace {

const double kPi = 3.14159265358979323846;
const double kLenCut = 0.5;
const double kRt2Over2 = 0.7071067811865476;

struct Conn {
    int helix;
    int start_base; // facing base stored as the historical "start"
    int end_base;
    int other;
    int real_gap;
    std::vector<int> gap_bases;
};

struct LoopB {
    int parent;
    int closing;
    std::vector<int> children;
    std::vector<Conn> conns;
    NavKind kind;
    int start;
    int end;
};

struct Builder {
    int n;
    bool circular;
    std::vector<int> mate;
    std::vector<NavHelix> helices;
    std::vector<int> helix_at; // per base, -1 if unpaired
    std::vector<char> entered;
    std::vector<LoopB> loops;
};

int advance_base(const Builder& b, int base) {
    if (base == b.n - 1) {
        return 0;
    }
    return base + 1;
}

// Real unpaired bases strictly between end_base and next_start, walking
// forward. A linear walk from n-1 to 0 crosses the phantom and emits no
// base for it.
void gap_bases(const Builder& b, int end_base, int next_start, std::vector<int>* out) {
    out->clear();
    int cur = end_base;
    int guard = 0;
    while (true) {
        if (!b.circular && cur == b.n - 1) {
            cur = 0;
        } else {
            cur = advance_base(b, cur);
        }
        if (cur == next_start) {
            break;
        }
        out->push_back(cur);
        asserta(++guard <= b.n + 2);
    }
}

bool crosses_origin(int end_base, int next_start) {
    return next_start <= end_base;
}

void append_int(std::string* s, int v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%d", v);
    s->append(buf);
}

char kind_letter(NavKind k) {
    switch (k) {
    case NavKind::Exterior:
        return 'E';
    case NavKind::Circle:
        return 'C';
    case NavKind::Hairpin:
        return 'P';
    case NavKind::Bulge:
        return 'B';
    case NavKind::Interior:
        return 'I';
    case NavKind::Multibranch:
        return 'M';
    default:
        asserta(false);
        return '?';
    }
}

int helix_len(const NavHelix& h) {
    return h.end5 - h.start5 + 1;
}

// --- pairs and helices -------------------------------------------------

bool build_pairs(Builder* b, const char* db, const char** err) {
    b->mate.assign(static_cast<size_t>(b->n), -1);
    std::vector<int> stack;
    stack.reserve(static_cast<size_t>(b->n));
    for (int i = 0; i < b->n; ++i) {
        char c = db[i];
        if (c == '.') {
            continue;
        }
        if (c == '(') {
            stack.push_back(i);
            continue;
        }
        if (c == ')') {
            if (stack.empty()) {
                if (err) {
                    *err = "unbalanced brackets";
                }
                return false;
            }
            int j = stack.back();
            stack.pop_back();
            b->mate[static_cast<size_t>(i)] = j;
            b->mate[static_cast<size_t>(j)] = i;
            continue;
        }
        if (err) {
            *err = "structure character is not ( . or )";
        }
        return false;
    }
    if (!stack.empty()) {
        if (err) {
            *err = "unbalanced brackets";
        }
        return false;
    }
    return true;
}

void build_helices(Builder* b) {
    b->helix_at.assign(static_cast<size_t>(b->n), -1);
    int i = 0;
    while (i < b->n) {
        int j = b->mate[static_cast<size_t>(i)];
        if (j > i && b->helix_at[static_cast<size_t>(i)] < 0) {
            int a = i;
            int c = j;
            int start5 = a;
            int end3 = c;
            while (a + 1 < c - 1 && b->mate[static_cast<size_t>(a + 1)] == c - 1) {
                a += 1;
                c -= 1;
            }
            NavHelix h;
            h.start5 = start5;
            h.end5 = a;
            h.start3 = c;
            h.end3 = end3;
            h.outer_loop = -1;
            h.inner_loop = -1;
            int id = static_cast<int>(b->helices.size());
            b->helices.push_back(h);
            for (int k = 0; k < helix_len(h); ++k) {
                int p = h.start5 + k;
                int q = h.end3 - k;
                b->helix_at[static_cast<size_t>(p)] = id;
                b->helix_at[static_cast<size_t>(q)] = id;
            }
            i = a + 1;
        } else {
            i += 1;
        }
    }
}

// --- loop tree ---------------------------------------------------------

int build_loop(Builder* b, int ibase) {
    int id = static_cast<int>(b->loops.size());
    LoopB loop;
    loop.parent = -1;
    loop.closing = -1;
    loop.kind = NavKind::Exterior;
    loop.start = -1;
    loop.end = -1;
    b->loops.push_back(loop);

    int cursor = ibase;
    int guard = 0;
    do {
        if (cursor >= 0 && b->mate[static_cast<size_t>(cursor)] >= 0) {
            int hid = b->helix_at[static_cast<size_t>(cursor)];
            asserta(hid >= 0);
            const NavHelix& h = b->helices[static_cast<size_t>(hid)];
            if (!b->entered[static_cast<size_t>(hid)]) {
                asserta(cursor == h.start5);
                b->entered[static_cast<size_t>(hid)] = 1;
                int child = build_loop(b, h.end5 + 1);
                b->loops[static_cast<size_t>(child)].parent = id;
                b->loops[static_cast<size_t>(child)].closing = hid;
                b->helices[static_cast<size_t>(hid)].outer_loop = id;
                b->helices[static_cast<size_t>(hid)].inner_loop = child;
                b->loops[static_cast<size_t>(id)].children.push_back(hid);

                Conn outer;
                outer.helix = hid;
                outer.start_base = h.start5;
                outer.end_base = h.end3;
                outer.other = child;
                outer.real_gap = 0;
                b->loops[static_cast<size_t>(id)].conns.push_back(outer);

                Conn inner;
                inner.helix = hid;
                inner.start_base = h.start3;
                inner.end_base = h.end5;
                inner.other = id;
                inner.real_gap = 0;
                b->loops[static_cast<size_t>(child)].conns.push_back(inner);
            }
            cursor = b->mate[static_cast<size_t>(cursor)];
        }
        if (!b->circular && cursor == b->n - 1) {
            cursor = -1;
        } else if (cursor < 0) {
            cursor = 0;
        } else {
            cursor = advance_base(*b, cursor);
        }
        asserta(++guard <= b->n + 3);
    } while (cursor != ibase);
    return id;
}

void fill_gaps(Builder* b) {
    for (size_t li = 0; li < b->loops.size(); ++li) {
        LoopB& loop = b->loops[li];
        int d = static_cast<int>(loop.conns.size());
        for (int i = 0; i < d; ++i) {
            int j = i + 1;
            if (j == d) {
                j = 0;
            }
            Conn& c = loop.conns[static_cast<size_t>(i)];
            const Conn& nxt = loop.conns[static_cast<size_t>(j)];
            gap_bases(*b, c.end_base, nxt.start_base, &c.gap_bases);
            c.real_gap = static_cast<int>(c.gap_bases.size());
        }
    }
}

int loop_depth(const Builder& b, int u, std::vector<char>* mark) {
    const LoopB& loop = b.loops[static_cast<size_t>(u)];
    if (loop.conns.size() <= 1) {
        return 0;
    }
    if ((*mark)[static_cast<size_t>(u)]) {
        return -1;
    }
    (*mark)[static_cast<size_t>(u)] = 1;
    int count = 0;
    int ret = 0;
    for (size_t i = 0; i < loop.conns.size(); ++i) {
        int d = loop_depth(b, loop.conns[i].other, mark);
        if (d >= 0) {
            if (count == 0) {
                ret = d;
            } else if (d < ret) {
                ret = d;
            }
            count += 1;
        }
    }
    (*mark)[static_cast<size_t>(u)] = 0;
    return ret + 1;
}

int drawing_root(const Builder& b) {
    int best = 0;
    int best_conn = -1;
    int best_depth = -1;
    std::vector<char> mark(b.loops.size(), 0);
    for (int i = 0; i < static_cast<int>(b.loops.size()); ++i) {
        int nc = static_cast<int>(b.loops[static_cast<size_t>(i)].conns.size());
        mark.assign(b.loops.size(), 0);
        int dep = loop_depth(b, i, &mark);
        if (nc > best_conn || (nc == best_conn && dep > best_depth)) {
            best = i;
            best_conn = nc;
            best_depth = dep;
        }
    }
    return best;
}

void classify(Builder* b) {
    for (size_t li = 0; li < b->loops.size(); ++li) {
        LoopB& loop = b->loops[li];
        int d = static_cast<int>(loop.conns.size());
        int empty = 0;
        for (int i = 0; i < d; ++i) {
            if (loop.conns[static_cast<size_t>(i)].real_gap == 0) {
                empty += 1;
            }
        }
        if (li == 0 && !b->circular) {
            loop.kind = NavKind::Exterior;
        } else if (d <= 1 && !(li == 0 && b->circular && d == 0)) {
            loop.kind = (d == 0) ? NavKind::Exterior : NavKind::Hairpin;
            if (d == 0 && b->circular) {
                loop.kind = NavKind::Circle;
            }
        } else if (d == 2 && empty == 1) {
            loop.kind = NavKind::Bulge;
        } else if (d == 2) {
            loop.kind = NavKind::Interior;
        } else if (d >= 3) {
            loop.kind = NavKind::Multibranch;
        } else {
            loop.kind = b->circular ? NavKind::Circle : NavKind::Exterior;
        }
    }
}

void set_spans(Builder* b) {
    for (size_t li = 0; li < b->loops.size(); ++li) {
        LoopB& loop = b->loops[li];
        if (li == 0 && !b->circular) {
            loop.start = 0;
            loop.end = b->n - 1;
            continue;
        }
        if (li == 0 && b->circular) {
            int wrap = -1;
            int d = static_cast<int>(loop.conns.size());
            for (int i = 0; i < d; ++i) {
                int j = (i + 1) % d;
                if (crosses_origin(loop.conns[static_cast<size_t>(i)].end_base,
                                   loop.conns[static_cast<size_t>(j)].start_base)) {
                    wrap = i;
                    break;
                }
            }
            if (wrap < 0 || loop.conns[static_cast<size_t>(wrap)].gap_bases.empty()) {
                loop.start = -1;
                loop.end = -1;
            } else {
                const std::vector<int>& g = loop.conns[static_cast<size_t>(wrap)].gap_bases;
                int end_base = loop.conns[static_cast<size_t>(wrap)].end_base;
                int next_start =
                    loop.conns[static_cast<size_t>((wrap + 1) % d)].start_base;
                bool has_high = end_base < b->n - 1 &&
                                std::find(g.begin(), g.end(), end_base + 1) != g.end();
                bool has_low = next_start > 0 &&
                               std::find(g.begin(), g.end(), next_start - 1) != g.end();
                if (has_high && has_low) {
                    loop.start = end_base + 1;
                    loop.end = next_start - 1;
                } else if (has_low) {
                    loop.start = 0;
                    loop.end = next_start - 1;
                } else if (has_high) {
                    loop.start = end_base + 1;
                    loop.end = b->n - 1;
                } else {
                    loop.start = g.front();
                    loop.end = g.back();
                }
            }
            continue;
        }
        const NavHelix& h = b->helices[static_cast<size_t>(loop.closing)];
        loop.start = h.end5;
        loop.end = h.start3;
    }
}

// Rotate connection index so the parent helix is at position 0.
// The root has no parent; its list stays in walk order.
int parent_conn_index(const LoopB& loop) {
    if (loop.closing < 0) {
        return -1;
    }
    for (int i = 0; i < static_cast<int>(loop.conns.size()); ++i) {
        if (loop.conns[static_cast<size_t>(i)].helix == loop.closing) {
            return i;
        }
    }
    asserta(false);
    return 0;
}

std::string naview_graph_string(const Builder& b, int loop_id) {
    const LoopB& loop = b.loops[static_cast<size_t>(loop_id)];
    int d = static_cast<int>(loop.conns.size());
    int pivot = parent_conn_index(loop); // -1 for the root
    std::string s;
    s.push_back(kind_letter(loop.kind));
    if (d == 0) {
        append_int(&s, b.n);
        return s;
    }
    std::vector<int> order(static_cast<size_t>(d));
    if (pivot < 0) {
        for (int i = 0; i < d; ++i) {
            order[static_cast<size_t>(i)] = i;
        }
    } else {
        for (int i = 0; i < d; ++i) {
            order[static_cast<size_t>(i)] = (pivot + i) % d;
        }
    }
    if (pivot < 0) {
        // Gap that precedes each helix, starting at the wrap / base-0 gap.
        for (int i = 0; i < d; ++i) {
            if (i) {
                s.push_back('+');
            }
            int prev = order[static_cast<size_t>((i + d - 1) % d)];
            append_int(&s, loop.conns[static_cast<size_t>(prev)].real_gap);
        }
    } else {
        for (int i = 0; i < d; ++i) {
            if (i) {
                s.push_back('+');
            }
            int idx = order[static_cast<size_t>(i)];
            append_int(&s, loop.conns[static_cast<size_t>(idx)].real_gap);
        }
    }
    int child_begin = (pivot < 0) ? 0 : 1;
    if (child_begin < d) {
        s.push_back('(');
        for (int i = child_begin; i < d; ++i) {
            if (i > child_begin) {
                s.push_back(',');
            }
            int idx = order[static_cast<size_t>(i)];
            int hid = loop.conns[static_cast<size_t>(idx)].helix;
            append_int(&s, helix_len(b.helices[static_cast<size_t>(hid)]));
            s.push_back(':');
            int child_loop = b.helices[static_cast<size_t>(hid)].inner_loop;
            if (child_loop == loop_id) {
                child_loop = b.helices[static_cast<size_t>(hid)].outer_loop;
            }
            s += naview_graph_string(b, child_loop);
        }
        s.push_back(')');
    }
    return s;
}

// --- segments ----------------------------------------------------------

struct BaseTag {
    NavKind kind;
    int helix_id;
    int loop_id;
    NavSide side;
    int gap_id;
    bool unpaired;
};

void build_segments(const Builder& b, NavLayout* out) {
    std::vector<BaseTag> tag(static_cast<size_t>(b.n));
    int next_gap = 1;
    for (size_t hi = 0; hi < b.helices.size(); ++hi) {
        const NavHelix& h = b.helices[hi];
        for (int k = h.start5; k <= h.end5; ++k) {
            BaseTag t;
            t.kind = NavKind::Helix5;
            t.helix_id = static_cast<int>(hi);
            t.loop_id = -1;
            t.side = NavSide::Five;
            t.gap_id = 0;
            t.unpaired = false;
            tag[static_cast<size_t>(k)] = t;
        }
        for (int k = h.start3; k <= h.end3; ++k) {
            BaseTag t;
            t.kind = NavKind::Helix3;
            t.helix_id = static_cast<int>(hi);
            t.loop_id = -1;
            t.side = NavSide::Three;
            t.gap_id = 0;
            t.unpaired = false;
            tag[static_cast<size_t>(k)] = t;
        }
    }
    for (size_t li = 0; li < b.loops.size(); ++li) {
        const LoopB& loop = b.loops[li];
        int pivot = parent_conn_index(loop);
        int d = static_cast<int>(loop.conns.size());
        for (int i = 0; i < d; ++i) {
            const Conn& c = loop.conns[static_cast<size_t>(i)];
            if (c.gap_bases.empty()) {
                continue;
            }
            NavSide side = NavSide::None;
            if (loop.kind == NavKind::Bulge || loop.kind == NavKind::Interior) {
                // Gap after the parent connection is the 5' side.
                side = (i == pivot) ? NavSide::Five : NavSide::Three;
            }
            int gid = next_gap++;
            for (size_t k = 0; k < c.gap_bases.size(); ++k) {
                int base = c.gap_bases[k];
                BaseTag t;
                t.kind = loop.kind;
                t.helix_id = -1;
                t.loop_id = static_cast<int>(li);
                t.side = side;
                t.gap_id = gid;
                t.unpaired = true;
                tag[static_cast<size_t>(base)] = t;
            }
        }
        if (d == 0) {
            int gid = next_gap++;
            for (int base = 0; base < b.n; ++base) {
                BaseTag t;
                t.kind = loop.kind;
                t.helix_id = -1;
                t.loop_id = 0;
                t.side = NavSide::None;
                t.gap_id = gid;
                t.unpaired = true;
                tag[static_cast<size_t>(base)] = t;
            }
        }
    }

    auto same = [](const BaseTag& a, const BaseTag& c) {
        return a.kind == c.kind && a.helix_id == c.helix_id && a.loop_id == c.loop_id &&
               a.side == c.side && a.gap_id == c.gap_id;
    };

    out->segment_of.assign(static_cast<size_t>(b.n), -1);
    bool wrap = b.circular && b.n > 1 && tag[0].unpaired && tag[static_cast<size_t>(b.n - 1)].unpaired &&
                same(tag[0], tag[static_cast<size_t>(b.n - 1)]);
    int scan_from = 0;
    int scan_to = b.n; // exclusive
    if (wrap) {
        int a = 0;
        while (a + 1 < b.n && same(tag[static_cast<size_t>(a + 1)], tag[0])) {
            a += 1;
        }
        int c = b.n - 1;
        while (c - 1 > a && same(tag[static_cast<size_t>(c - 1)], tag[0])) {
            c -= 1;
        }
        NavSegment seg;
        seg.kind = tag[0].kind;
        seg.start = c;
        seg.end = a;
        seg.helix_id = -1;
        seg.loop_id = tag[0].loop_id;
        seg.side = tag[0].side;
        int sid = static_cast<int>(out->segments.size());
        out->segments.push_back(seg);
        for (int base = 0; base <= a; ++base) {
            out->segment_of[static_cast<size_t>(base)] = sid;
        }
        for (int base = c; base < b.n; ++base) {
            out->segment_of[static_cast<size_t>(base)] = sid;
        }
        scan_from = a + 1;
        scan_to = c;
    }
    int base = scan_from;
    while (base < scan_to) {
        int end = base;
        while (end + 1 < scan_to && same(tag[static_cast<size_t>(end + 1)], tag[static_cast<size_t>(base)])) {
            end += 1;
        }
        NavSegment seg;
        seg.kind = tag[static_cast<size_t>(base)].kind;
        seg.start = base;
        seg.end = end;
        seg.helix_id = tag[static_cast<size_t>(base)].helix_id;
        seg.loop_id = tag[static_cast<size_t>(base)].loop_id;
        seg.side = tag[static_cast<size_t>(base)].side;
        int sid = static_cast<int>(out->segments.size());
        out->segments.push_back(seg);
        for (int k = base; k <= end; ++k) {
            out->segment_of[static_cast<size_t>(k)] = sid;
        }
        base = end + 1;
    }
}

// --- RAG ----------------------------------------------------------------

struct Port {
    int edge;
    int gap;
    int next;
    int prev;
    int node;
    bool live;
};

struct REdge {
    int a;
    int b;
    int len;
    int pa;
    int pb;
    bool live;
};

struct RNode {
    int head;
    bool live;
    bool is_end;
    int extra; // unpaired bases kept when the node has no ports
};

struct Rag {
    std::vector<RNode> nodes;
    std::vector<REdge> edges;
    std::vector<Port> ports;
};

int rag_degree(const Rag& g, int node) {
    int head = g.nodes[static_cast<size_t>(node)].head;
    if (head < 0) {
        return 0;
    }
    int c = 0;
    int p = head;
    do {
        c += 1;
        p = g.ports[static_cast<size_t>(p)].next;
        asserta(c <= static_cast<int>(g.ports.size()) + 1);
    } while (p != head);
    return c;
}

int rag_unpaired(const Rag& g, int node) {
    int head = g.nodes[static_cast<size_t>(node)].head;
    if (head < 0) {
        return g.nodes[static_cast<size_t>(node)].extra;
    }
    int sum = 0;
    int p = head;
    do {
        sum += g.ports[static_cast<size_t>(p)].gap;
        p = g.ports[static_cast<size_t>(p)].next;
    } while (p != head);
    return sum;
}

int add_port(Rag* g, int node, int edge, int gap) {
    int id = static_cast<int>(g->ports.size());
    Port p;
    p.edge = edge;
    p.gap = gap;
    p.node = node;
    p.live = true;
    p.next = id;
    p.prev = id;
    g->ports.push_back(p);
    RNode& nd = g->nodes[static_cast<size_t>(node)];
    if (nd.head < 0) {
        nd.head = id;
    } else {
        int head = nd.head;
        int tail = g->ports[static_cast<size_t>(head)].prev;
        g->ports[static_cast<size_t>(tail)].next = id;
        g->ports[static_cast<size_t>(id)].prev = tail;
        g->ports[static_cast<size_t>(id)].next = head;
        g->ports[static_cast<size_t>(head)].prev = id;
    }
    return id;
}

void unlink_port(Rag* g, int port) {
    Port& p = g->ports[static_cast<size_t>(port)];
    RNode& nd = g->nodes[static_cast<size_t>(p.node)];
    if (p.next == port) {
        nd.head = -1;
    } else {
        g->ports[static_cast<size_t>(p.prev)].next = p.next;
        g->ports[static_cast<size_t>(p.next)].prev = p.prev;
        if (nd.head == port) {
            nd.head = p.next;
        }
    }
    p.live = false;
}

// Degree-1 node A is removed. Its unpaired bases are added into the merged
// flank on B when add_extra is set (length-1 fusion). Step-2 deletion leaves
// add_extra false.
void absorb_leaf(Rag* g, int bnode, int port_on_b, int extra) {
    Port& p = g->ports[static_cast<size_t>(port_on_b)];
    if (p.next == port_on_b) {
        g->nodes[static_cast<size_t>(bnode)].extra = p.gap + extra;
        g->nodes[static_cast<size_t>(bnode)].head = -1;
        p.live = false;
        return;
    }
    Port& prev = g->ports[static_cast<size_t>(p.prev)];
    prev.gap += p.gap + extra;
    unlink_port(g, port_on_b);
}

void retarget_edges(Rag* g, int from, int to) {
    for (size_t i = 0; i < g->edges.size(); ++i) {
        REdge& e = g->edges[i];
        if (!e.live) {
            continue;
        }
        if (e.a == from) {
            e.a = to;
        }
        if (e.b == from) {
            e.b = to;
        }
    }
}

void fuse_length_one(Rag* g, int ei) {
    REdge& e = g->edges[static_cast<size_t>(ei)];
    int a = e.a;
    int b = e.b;
    int pa = e.pa;
    int pb = e.pb;
    int da = rag_degree(*g, a);
    int db = rag_degree(*g, b);
    e.live = false;
    if (da == 1 && db == 1) {
        g->nodes[static_cast<size_t>(a)].extra =
            rag_unpaired(*g, a) + rag_unpaired(*g, b);
        g->nodes[static_cast<size_t>(a)].head = -1;
        g->nodes[static_cast<size_t>(a)].is_end =
            g->nodes[static_cast<size_t>(a)].is_end || g->nodes[static_cast<size_t>(b)].is_end;
        g->nodes[static_cast<size_t>(b)].live = false;
        g->nodes[static_cast<size_t>(b)].head = -1;
        g->ports[static_cast<size_t>(pa)].live = false;
        g->ports[static_cast<size_t>(pb)].live = false;
        return;
    }
    if (da == 1) {
        int extra = rag_unpaired(*g, a);
        absorb_leaf(g, b, pb, extra);
        g->nodes[static_cast<size_t>(b)].is_end =
            g->nodes[static_cast<size_t>(b)].is_end || g->nodes[static_cast<size_t>(a)].is_end;
        g->nodes[static_cast<size_t>(a)].live = false;
        g->nodes[static_cast<size_t>(a)].head = -1;
        g->ports[static_cast<size_t>(pa)].live = false;
        return;
    }
    if (db == 1) {
        int extra = rag_unpaired(*g, b);
        absorb_leaf(g, a, pa, extra);
        g->nodes[static_cast<size_t>(a)].is_end =
            g->nodes[static_cast<size_t>(a)].is_end || g->nodes[static_cast<size_t>(b)].is_end;
        g->nodes[static_cast<size_t>(b)].live = false;
        g->nodes[static_cast<size_t>(b)].head = -1;
        g->ports[static_cast<size_t>(pb)].live = false;
        return;
    }
    // Splice B's cycle into A at the contracted edge, same backbone direction.
    int a_prev = g->ports[static_cast<size_t>(pa)].prev;
    int a_next = g->ports[static_cast<size_t>(pa)].next;
    int b_prev = g->ports[static_cast<size_t>(pb)].prev;
    int b_next = g->ports[static_cast<size_t>(pb)].next;
    g->ports[static_cast<size_t>(a_prev)].gap += g->ports[static_cast<size_t>(pb)].gap;
    g->ports[static_cast<size_t>(b_prev)].gap += g->ports[static_cast<size_t>(pa)].gap;
    g->ports[static_cast<size_t>(a_prev)].next = b_next;
    g->ports[static_cast<size_t>(b_next)].prev = a_prev;
    g->ports[static_cast<size_t>(b_prev)].next = a_next;
    g->ports[static_cast<size_t>(a_next)].prev = b_prev;
    int p = b_next;
    int guard = 0;
    while (p != a_next) {
        g->ports[static_cast<size_t>(p)].node = a;
        p = g->ports[static_cast<size_t>(p)].next;
        asserta(++guard <= static_cast<int>(g->ports.size()) + 1);
    }
    g->ports[static_cast<size_t>(pa)].live = false;
    g->ports[static_cast<size_t>(pb)].live = false;
    g->nodes[static_cast<size_t>(a)].head = a_next;
    g->nodes[static_cast<size_t>(a)].is_end =
        g->nodes[static_cast<size_t>(a)].is_end || g->nodes[static_cast<size_t>(b)].is_end;
    g->nodes[static_cast<size_t>(b)].live = false;
    g->nodes[static_cast<size_t>(b)].head = -1;
    retarget_edges(g, b, a);
}

void delete_degree_one(Rag* g, int node) {
    int head = g->nodes[static_cast<size_t>(node)].head;
    asserta(head >= 0);
    int ei = g->ports[static_cast<size_t>(head)].edge;
    REdge& e = g->edges[static_cast<size_t>(ei)];
    int other = (e.a == node) ? e.b : e.a;
    int port_on_other = (e.a == node) ? e.pb : e.pa;
    if (g->ports[static_cast<size_t>(e.pa)].node == other) {
        port_on_other = e.pa;
    }
    if (g->ports[static_cast<size_t>(e.pb)].node == other) {
        port_on_other = e.pb;
    }
    absorb_leaf(g, other, port_on_other, 0);
    e.live = false;
    g->nodes[static_cast<size_t>(node)].live = false;
    g->nodes[static_cast<size_t>(node)].head = -1;
    g->ports[static_cast<size_t>(head)].live = false;
}

void delete_degree_two(Rag* g, int node) {
    int head = g->nodes[static_cast<size_t>(node)].head;
    int p2 = g->ports[static_cast<size_t>(head)].next;
    asserta(p2 != head);
    int e1 = g->ports[static_cast<size_t>(head)].edge;
    int e2 = g->ports[static_cast<size_t>(p2)].edge;
    // Copy out of edges before push_back. A reference would dangle when the
    // vector reallocates, and the later live=false stores would hit freed memory.
    int a1 = g->edges[static_cast<size_t>(e1)].a;
    int b1 = g->edges[static_cast<size_t>(e1)].b;
    int a2 = g->edges[static_cast<size_t>(e2)].a;
    int b2 = g->edges[static_cast<size_t>(e2)].b;
    int len1 = g->edges[static_cast<size_t>(e1)].len;
    int len2 = g->edges[static_cast<size_t>(e2)].len;
    int pa1 = g->edges[static_cast<size_t>(e1)].pa;
    int pb1 = g->edges[static_cast<size_t>(e1)].pb;
    int pa2 = g->edges[static_cast<size_t>(e2)].pa;
    int pb2 = g->edges[static_cast<size_t>(e2)].pb;
    int n1 = (a1 == node) ? b1 : a1;
    int n2 = (a2 == node) ? b2 : a2;
    int q1 = (g->ports[static_cast<size_t>(pa1)].node == n1) ? pa1 : pb1;
    int q2 = (g->ports[static_cast<size_t>(pa2)].node == n2) ? pa2 : pb2;
    int neu = static_cast<int>(g->edges.size());
    REdge ne;
    ne.a = n1;
    ne.b = n2;
    ne.len = len1 + len2;
    ne.pa = q1;
    ne.pb = q2;
    ne.live = true;
    g->edges.push_back(ne);
    g->ports[static_cast<size_t>(q1)].edge = neu;
    g->ports[static_cast<size_t>(q2)].edge = neu;
    g->edges[static_cast<size_t>(e1)].live = false;
    g->edges[static_cast<size_t>(e2)].live = false;
    g->ports[static_cast<size_t>(head)].live = false;
    g->ports[static_cast<size_t>(p2)].live = false;
    g->nodes[static_cast<size_t>(node)].live = false;
    g->nodes[static_cast<size_t>(node)].head = -1;
}

bool qualifies(const Rag& g, int node) {
    if (!g.nodes[static_cast<size_t>(node)].live) {
        return false;
    }
    if (g.nodes[static_cast<size_t>(node)].is_end) {
        return false;
    }
    int d = rag_degree(g, node);
    if (d < 1 || d >= 3) {
        return false;
    }
    return rag_unpaired(g, node) <= 1;
}

void contract(Rag* g) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < static_cast<int>(g->edges.size()); ++i) {
            if (g->edges[static_cast<size_t>(i)].live && g->edges[static_cast<size_t>(i)].len == 1) {
                fuse_length_one(g, i);
                changed = true;
                break;
            }
        }
    }
    changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < static_cast<int>(g->nodes.size()); ++i) {
            if (!qualifies(*g, i)) {
                continue;
            }
            int d = rag_degree(*g, i);
            if (d == 1) {
                delete_degree_one(g, i);
            } else if (d == 2) {
                delete_degree_two(g, i);
            } else {
                continue;
            }
            changed = true;
            break;
        }
    }
}

struct AdjNode {
    std::vector<int> nbr;
    std::vector<int> gap_after; // parallel to nbr, backbone order
    bool is_end;
    int unpaired;
    bool live;
};

std::vector<AdjNode> adjacency(const Rag& g) {
    std::vector<AdjNode> a(g.nodes.size());
    for (size_t i = 0; i < g.nodes.size(); ++i) {
        a[i].live = g.nodes[i].live;
        a[i].is_end = g.nodes[i].is_end;
        a[i].unpaired = g.nodes[i].live ? rag_unpaired(g, static_cast<int>(i)) : 0;
        if (!g.nodes[i].live || g.nodes[i].head < 0) {
            continue;
        }
        int p = g.nodes[i].head;
        do {
            const Port& port = g.ports[static_cast<size_t>(p)];
            const REdge& e = g.edges[static_cast<size_t>(port.edge)];
            int other = (e.a == static_cast<int>(i)) ? e.b : e.a;
            a[i].nbr.push_back(other);
            a[i].gap_after.push_back(port.gap);
            p = port.next;
        } while (p != g.nodes[i].head);
    }
    return a;
}

int eccentricity(const std::vector<AdjNode>& a, int src) {
    std::vector<int> dist(a.size(), -1);
    std::vector<int> q;
    q.push_back(src);
    dist[static_cast<size_t>(src)] = 0;
    int qi = 0;
    int far = 0;
    while (qi < static_cast<int>(q.size())) {
        int u = q[static_cast<size_t>(qi++)];
        far = dist[static_cast<size_t>(u)];
        for (size_t k = 0; k < a[static_cast<size_t>(u)].nbr.size(); ++k) {
            int v = a[static_cast<size_t>(u)].nbr[k];
            if (dist[static_cast<size_t>(v)] < 0) {
                dist[static_cast<size_t>(v)] = far + 1;
                q.push_back(v);
            }
        }
    }
    return far;
}

char leaf_mark(const AdjNode& n) {
    if (n.is_end) {
        return '*';
    }
    return '@';
}

struct TNode {
    std::vector<TNode> kids;
    char mark; // set on leaves
    bool leaf;
};

std::string write_side(const std::vector<AdjNode>& a, int node, int parent, TNode* tree) {
    tree->leaf = false;
    tree->mark = 0;
    tree->kids.clear();
    const AdjNode& n = a[static_cast<size_t>(node)];
    std::vector<int> kids;
    for (size_t i = 0; i < n.nbr.size(); ++i) {
        if (n.nbr[i] != parent) {
            kids.push_back(n.nbr[i]);
        }
    }
    if (kids.empty()) {
        tree->leaf = true;
        tree->mark = leaf_mark(n);
        std::string s;
        s.push_back(tree->mark);
        return s;
    }
    struct Item {
        std::string text;
        TNode node;
    };
    std::vector<Item> items;
    items.reserve(kids.size());
    for (size_t i = 0; i < kids.size(); ++i) {
        Item item;
        int child = kids[i];
        bool child_leaf = true;
        for (size_t k = 0; k < a[static_cast<size_t>(child)].nbr.size(); ++k) {
            if (a[static_cast<size_t>(child)].nbr[k] != node) {
                child_leaf = false;
                break;
            }
        }
        if (child_leaf) {
            item.node.leaf = true;
            item.node.mark = leaf_mark(a[static_cast<size_t>(child)]);
            item.text.assign(1, item.node.mark);
        } else {
            item.text = write_side(a, child, node, &item.node);
        }
        items.push_back(std::move(item));
    }
    std::stable_sort(items.begin(), items.end(), [](const Item& x, const Item& y) {
        return x.text < y.text;
    });
    std::string s = "(";
    for (size_t i = 0; i < items.size(); ++i) {
        if (i) {
            s.push_back(',');
        }
        s += items[i].text;
        tree->kids.push_back(std::move(items[i].node));
    }
    s.push_back(')');
    return s;
}

int count_edges(const std::vector<AdjNode>& a) {
    int c = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].live) {
            c += static_cast<int>(a[i].nbr.size());
        }
    }
    return c / 2;
}

std::string rag_tree_string(const std::vector<AdjNode>& a, TNode* rooted, bool* bicenter,
                            int* c0, int* c1) {
    *bicenter = false;
    *c0 = -1;
    *c1 = -1;
    std::vector<int> live;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].live) {
            live.push_back(static_cast<int>(i));
        }
    }
    asserta(!live.empty());
    if (live.size() == 1) {
        rooted->leaf = true;
        rooted->mark = (a[static_cast<size_t>(live[0])].nbr.empty() && !a[static_cast<size_t>(live[0])].is_end &&
                        a[static_cast<size_t>(live[0])].unpaired > 1 && false)
                           ? '@'
                           : '*';
        // A lone vertex is '*' even when it began as a hairpin that lost its helix.
        if (a[static_cast<size_t>(live[0])].nbr.size() == 1) {
            rooted->mark = leaf_mark(a[static_cast<size_t>(live[0])]);
        } else {
            rooted->mark = '*';
        }
        rooted->kids.clear();
        return std::string(1, rooted->mark);
    }
    int best_e = 1000000000;
    std::vector<int> ecc(a.size(), 0);
    for (size_t i = 0; i < live.size(); ++i) {
        int e = eccentricity(a, live[i]);
        ecc[static_cast<size_t>(live[i])] = e;
        if (e < best_e) {
            best_e = e;
        }
    }
    std::vector<int> centers;
    for (size_t i = 0; i < live.size(); ++i) {
        if (ecc[static_cast<size_t>(live[i])] == best_e) {
            centers.push_back(live[i]);
        }
    }
    asserta(centers.size() == 1 || centers.size() == 2);
    if (centers.size() == 1) {
        *c0 = centers[0];
        return write_side(a, centers[0], -1, rooted);
    }
    *bicenter = true;
    *c0 = centers[0];
    *c1 = centers[1];
    TNode left;
    TNode right;
    std::string ls = write_side(a, centers[0], centers[1], &left);
    std::string rs = write_side(a, centers[1], centers[0], &right);
    if (rs < ls) {
        std::swap(ls, rs);
        std::swap(left, right);
        std::swap(*c0, *c1);
    }
    rooted->leaf = false;
    rooted->kids.clear();
    rooted->kids.push_back(std::move(left));
    rooted->kids.push_back(std::move(right));
    return ls + "-" + rs;
}

bool far_is_hairpin(const std::vector<AdjNode>& a, int node) {
    return a[static_cast<size_t>(node)].nbr.size() == 1 && !a[static_cast<size_t>(node)].is_end;
}

std::string dual_string(const std::vector<AdjNode>& a) {
    int edges = count_edges(a);
    if (edges == 0) {
        return "*";
    }
    if (edges == 1) {
        int u = -1;
        int v = -1;
        for (size_t i = 0; i < a.size(); ++i) {
            if (!a[i].live || a[i].nbr.empty()) {
                continue;
            }
            u = static_cast<int>(i);
            v = a[i].nbr[0];
            break;
        }
        int hairs = 0;
        if (far_is_hairpin(a, u)) {
            hairs += 1;
        }
        if (far_is_hairpin(a, v)) {
            hairs += 1;
        }
        if (hairs >= 2) {
            return "@@";
        }
        if (hairs == 1) {
            return "@";
        }
        return "*";
    }
    // Write the cycle at the unique center. Bicenter graphs that are a single
    // edge are handled above. A bicenter with extra branches uses the side
    // whose rendered cycle is lexicographically smaller.
    bool bic = false;
    int c0 = -1;
    int c1 = -1;
    TNode ignore;
    rag_tree_string(a, &ignore, &bic, &c0, &c1);

    auto render_node = [&](int node, int skip, auto&& self) -> std::string {
        const AdjNode& n = a[static_cast<size_t>(node)];
        std::vector<int> idx;
        for (int i = 0; i < static_cast<int>(n.nbr.size()); ++i) {
            if (n.nbr[static_cast<size_t>(i)] == skip) {
                continue;
            }
            idx.push_back(i);
        }
        int m = static_cast<int>(idx.size());
        asserta(m >= 1);
        std::vector<std::string> pieces(static_cast<size_t>(m));
        std::vector<int> gaps(static_cast<size_t>(m));
        for (int i = 0; i < m; ++i) {
            int at = idx[static_cast<size_t>(i)];
            int far = n.nbr[static_cast<size_t>(at)];
            gaps[static_cast<size_t>(i)] = n.gap_after[static_cast<size_t>(at)];
            bool leaf = true;
            for (size_t k = 0; k < a[static_cast<size_t>(far)].nbr.size(); ++k) {
                if (a[static_cast<size_t>(far)].nbr[k] != node) {
                    leaf = false;
                    break;
                }
            }
            if (leaf) {
                pieces[static_cast<size_t>(i)] = far_is_hairpin(a, far) ? "@" : "*";
            } else {
                pieces[static_cast<size_t>(i)] = self(far, node, self);
            }
        }
        auto join_cycle = [&](const std::vector<std::string>& pc, const std::vector<char>& seps,
                              bool closed) {
            std::string best;
            int nrot = static_cast<int>(pc.size());
            for (int rev = 0; rev < 2; ++rev) {
                std::vector<std::string> p = pc;
                std::vector<char> s = seps;
                if (rev) {
                    std::reverse(p.begin(), p.end());
                    if (closed) {
                        // Forward seps[i] sits between pc[i] and pc[i+1].
                        // Reversed, starting at pc[0]: the closer becomes the first sep.
                        std::vector<char> rs(s.size());
                        int ns = static_cast<int>(s.size());
                        rs[0] = s[static_cast<size_t>(ns - 1)];
                        for (int i = 1; i < ns; ++i) {
                            rs[static_cast<size_t>(i)] = s[static_cast<size_t>(ns - 1 - i)];
                        }
                        // After reversing pieces, piece 0 stays, then old last, etc.
                        // Rebuild pieces as p0, p_{n-1}, p_{n-2}, ... p1
                        std::vector<std::string> rp(p.size());
                        rp[0] = pc[0];
                        for (int i = 1; i < nrot; ++i) {
                            rp[static_cast<size_t>(i)] = pc[static_cast<size_t>(nrot - i)];
                        }
                        p.swap(rp);
                        s.swap(rs);
                    } else if (p.size() == 2 && s.size() == 1) {
                        std::swap(p[0], p[1]);
                    }
                }
                for (int rot = 0; rot < nrot; ++rot) {
                    std::string cur = "(";
                    for (int i = 0; i < nrot; ++i) {
                        int at = (rot + i) % nrot;
                        cur += p[static_cast<size_t>(at)];
                        bool emit = closed || i + 1 < nrot;
                        if (emit) {
                            int sep_at = closed ? at : 0;
                            if (!closed) {
                                cur.push_back(s[0]);
                            } else {
                                cur.push_back(s[static_cast<size_t>(sep_at)]);
                            }
                        }
                    }
                    cur.push_back(')');
                    if (best.empty() || cur < best) {
                        best = cur;
                    }
                }
            }
            return best;
        };
        if (m == 2) {
            int pos = (gaps[0] > 0 ? 1 : 0) + (gaps[1] > 0 ? 1 : 0);
            char sep = pos >= 2 ? '=' : (pos == 1 ? '-' : '.');
            std::vector<char> seps(1, sep);
            return join_cycle(pieces, seps, false);
        }
        std::vector<char> seps(static_cast<size_t>(m));
        for (int i = 0; i < m; ++i) {
            seps[static_cast<size_t>(i)] = gaps[static_cast<size_t>(i)] > 0 ? '-' : '.';
        }
        return join_cycle(pieces, seps, true);
    };

    if (!bic) {
        return render_node(c0, -999999, render_node);
    }
    std::string left = render_node(c0, c1, render_node);
    std::string right = render_node(c1, c0, render_node);
    // The central edge's dual is already the one-edge case when there is
    // nothing else. With extra branches, prefer the lexicographically smaller
    // side cycle. Both sides include their other stems only.
    if (right < left) {
        return right;
    }
    return left;
}

// --- coordinates -------------------------------------------------------

struct PConn {
    int helix;
    int start; // station
    int end;
    double angle;
    double xrad;
    double yrad;
    bool extruded;
};

struct Placer {
    const Builder* b;
    int N;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<PConn> conns;
    int anchor_helix;
    int icroot;
    double xc;
    double yc;
    double radius;

    int to_station(int base) const {
        if (base < 0) {
            return 0;
        }
        return b->circular ? base : base + 1;
    }

    int station_steps(int end_st, int start_st) const {
        if (start_st > end_st) {
            return start_st - end_st;
        }
        return start_st + N - end_st;
    }

    bool connected(int i, int j) const {
        if (conns[static_cast<size_t>(i)].extruded) {
            return true;
        }
        return station_steps(conns[static_cast<size_t>(i)].end, conns[static_cast<size_t>(j)].start) == 1;
    }

    void determine_radius() {
        int d = static_cast<int>(conns.size());
        double mindit = 1.0e10;
        double radius_now = kRt2Over2;
        int imin = 0;
        bool extrude = false;
        do {
            mindit = 1.0e10;
            extrude = false;
            double sumn = 0.0;
            double sumd = 0.0;
            for (int i = 0; i < d; ++i) {
                int j = (i + 1) % d;
                double dt = conns[static_cast<size_t>(j)].angle - conns[static_cast<size_t>(i)].angle;
                if (dt <= 0.0) {
                    dt += 2.0 * kPi;
                }
                double ci;
                if (!conns[static_cast<size_t>(i)].extruded) {
                    ci = static_cast<double>(station_steps(conns[static_cast<size_t>(i)].end,
                                                          conns[static_cast<size_t>(j)].start));
                } else if (dt <= kPi / 2.0) {
                    ci = 2.0;
                } else {
                    ci = 1.5;
                }
                sumn += dt * (1.0 / ci + 1.0);
                sumd += dt * dt / ci;
                double dit = dt / ci;
                if (dit < mindit && !conns[static_cast<size_t>(i)].extruded && ci > 1.0) {
                    mindit = dit;
                    imin = i;
                }
            }
            radius_now = (sumd > 0.0) ? (sumn / sumd) : kRt2Over2;
            if (radius_now < kRt2Over2) {
                radius_now = kRt2Over2;
            }
            if (mindit * radius_now < kLenCut) {
                conns[static_cast<size_t>(imin)].extruded = true;
                extrude = true;
            }
        } while (extrude && mindit * radius_now < kLenCut);
        radius = radius_now;
    }

    void place_on_chord(int start_st, int end_st) {
        double x0 = px[static_cast<size_t>(start_st)];
        double y0 = py[static_cast<size_t>(start_st)];
        double x1 = px[static_cast<size_t>(end_st)];
        double y1 = py[static_cast<size_t>(end_st)];
        int m = station_steps(start_st, end_st);
        double dx = x1 - x0;
        double dy = y1 - y0;
        for (int j = 1; j < m; ++j) {
            int st = start_st + j;
            if (st >= N) {
                st -= N;
            }
            double f = static_cast<double>(j) / static_cast<double>(m);
            px[static_cast<size_t>(st)] = x0 + dx * f;
            py[static_cast<size_t>(st)] = y0 + dy * f;
        }
    }

    // Returns false when the arc solver fails and the bases were put on the chord.
    bool place_on_arc(int start_st, int end_st) {
        double x0 = px[static_cast<size_t>(start_st)];
        double y0 = py[static_cast<size_t>(start_st)];
        double x1 = px[static_cast<size_t>(end_st)];
        double y1 = py[static_cast<size_t>(end_st)];
        int m = station_steps(start_st, end_st);
        double dx = x1 - x0;
        double dy = y1 - y0;
        double b = std::sqrt(dx * dx + dy * dy);
        if (!(b > 1e-8) || b >= static_cast<double>(m)) {
            place_on_chord(start_st, end_st);
            return b + 1e-8 >= static_cast<double>(m);
        }
        double hhi = static_cast<double>(m) / kPi;
        double hlow = -hhi - b / (static_cast<double>(m) + 1e-6 - b);
        if (b < 1.0) {
            hlow = 0.0;
        }
        double h = 0.0;
        double theta = 0.0;
        bool ok = false;
        for (int iter = 0; iter < 500; ++iter) {
            h = 0.5 * (hhi + hlow);
            double rho = std::sqrt(h * h + 0.25 * b * b);
            double disc = 1.0 - 0.5 / (rho * rho);
            if (std::fabs(disc) > 1.0) {
                place_on_chord(start_st, end_st);
                return false;
            }
            theta = std::acos(disc);
            double phi = std::acos(h / rho);
            double err = theta * static_cast<double>(m) + 2.0 * phi - 2.0 * kPi;
            if (err > 0.0) {
                hlow = h;
            } else {
                hhi = h;
            }
            if (std::fabs(err) <= 1e-4) {
                ok = true;
                break;
            }
        }
        if (!ok) {
            place_on_chord(start_st, end_st);
            return false;
        }
        dx /= b;
        dy /= b;
        double midx = x0 + dx * b * 0.5;
        double midy = y0 + dy * b * 0.5;
        double cx = midx + h * dy;
        double cy = midy - h * dx;
        double mx = x0 - cx;
        double my = y0 - cy;
        double rho = std::sqrt(mx * mx + my * my);
        double ang = std::atan2(my, mx);
        for (int j = 1; j < m; ++j) {
            int st = start_st + j;
            if (st >= N) {
                st -= N;
            }
            px[static_cast<size_t>(st)] = cx + rho * std::cos(ang + j * theta);
            py[static_cast<size_t>(st)] = cy + rho * std::sin(ang + j * theta);
        }
        return true;
    }

    void construct_extruded(int i) {
        int j = (i + 1) % static_cast<int>(conns.size());
        int start = conns[static_cast<size_t>(i)].end;
        int end = conns[static_cast<size_t>(j)].start;
        int m = station_steps(start, end);
        if (m <= 1) {
            return;
        }
        double astart = conns[static_cast<size_t>(i)].angle;
        double aend = conns[static_cast<size_t>(j)].angle;
        double aend2 = aend;
        if (aend2 < astart) {
            aend2 += 2.0 * kPi;
        }
        double aave = 0.5 * (astart + aend2);
        if (m == 2) {
            place_on_arc(start, end);
            return;
        }
        bool collision = true;
        while (collision && m > 1) {
            collision = false;
            place_on_arc(start, end);
            if (m == 2) {
                break;
            }
            int nstart = start + 1;
            if (nstart >= N) {
                nstart -= N;
            }
            int nend = end - 1;
            if (nend < 0) {
                nend += N;
            }
            double ddx = px[static_cast<size_t>(nstart)] - px[static_cast<size_t>(start)];
            double ddy = py[static_cast<size_t>(nstart)] - py[static_cast<size_t>(start)];
            double a1 = std::atan2(ddy, ddx);
            if (a1 < 0.0) {
                a1 += 2.0 * kPi;
            }
            double dac = a1 - astart;
            if (dac < 0.0) {
                dac += 2.0 * kPi;
            }
            if (dac > kPi) {
                collision = true;
            }
            ddx = px[static_cast<size_t>(nend)] - px[static_cast<size_t>(end)];
            ddy = py[static_cast<size_t>(nend)] - py[static_cast<size_t>(end)];
            double a2 = std::atan2(ddy, ddx);
            if (a2 < 0.0) {
                a2 += 2.0 * kPi;
            }
            dac = aend - a2;
            if (dac < 0.0) {
                dac += 2.0 * kPi;
            }
            if (dac > kPi) {
                collision = true;
            }
            if (collision) {
                double ac = std::min(aave, astart + 0.5);
                px[static_cast<size_t>(nstart)] = px[static_cast<size_t>(start)] + std::cos(ac);
                py[static_cast<size_t>(nstart)] = py[static_cast<size_t>(start)] + std::sin(ac);
                start = nstart;
                ac = std::max(aave, aend2 - 0.5);
                px[static_cast<size_t>(nend)] = px[static_cast<size_t>(end)] + std::cos(ac);
                py[static_cast<size_t>(nend)] = py[static_cast<size_t>(end)] + std::sin(ac);
                end = nend;
                m -= 2;
            }
        }
    }

    int to_base(int station) const {
        if (!b->circular) {
            return station == 0 ? -1 : station - 1;
        }
        return station;
    }

    // Step away from this loop along the helix. The facing pair is already
    // placed; xrad/yrad point from that pair toward the other end.
    void generate_region(int ic) {
        const PConn& cp = conns[static_cast<size_t>(ic)];
        const NavHelix& h = b->helices[static_cast<size_t>(cp.helix)];
        int len = helix_len(h);
        int facing = to_base(cp.start);
        if (facing == h.start5) {
            for (int k = 1; k < len; ++k) {
                int p = to_station(h.start5 + k);
                int q = to_station(h.end3 - k);
                px[static_cast<size_t>(p)] = px[static_cast<size_t>(cp.start)] + k * cp.xrad;
                py[static_cast<size_t>(p)] = py[static_cast<size_t>(cp.start)] + k * cp.yrad;
                px[static_cast<size_t>(q)] = px[static_cast<size_t>(cp.end)] + k * cp.xrad;
                py[static_cast<size_t>(q)] = py[static_cast<size_t>(cp.end)] + k * cp.yrad;
            }
        } else {
            asserta(facing == h.start3);
            for (int k = 1; k < len; ++k) {
                int p = to_station(h.start3 + k);
                int q = to_station(h.start5 + (len - 1 - k));
                px[static_cast<size_t>(p)] = px[static_cast<size_t>(cp.start)] + k * cp.xrad;
                py[static_cast<size_t>(p)] = py[static_cast<size_t>(cp.start)] + k * cp.yrad;
                px[static_cast<size_t>(q)] = px[static_cast<size_t>(cp.end)] + k * cp.xrad;
                py[static_cast<size_t>(q)] = py[static_cast<size_t>(cp.end)] + k * cp.yrad;
            }
        }
    }

    void load_conns(int loop_id) {
        const LoopB& loop = b->loops[static_cast<size_t>(loop_id)];
        conns.clear();
        icroot = -1;
        for (size_t i = 0; i < loop.conns.size(); ++i) {
            const Conn& c = loop.conns[i];
            PConn p;
            p.helix = c.helix;
            p.start = to_station(c.start_base);
            p.end = to_station(c.end_base);
            p.extruded = false;
            double ts = 2.0 * kPi * static_cast<double>(p.start) / static_cast<double>(N);
            double te = 2.0 * kPi * static_cast<double>(p.end) / static_cast<double>(N);
            double xs = -std::sin(ts);
            double ys = std::cos(ts);
            double xe = -std::sin(te);
            double ye = std::cos(te);
            double xn = ye - ys;
            double yn = xs - xe;
            double r = std::sqrt(xn * xn + yn * yn);
            asserta(r > 0.0);
            p.xrad = xn / r;
            p.yrad = yn / r;
            p.angle = std::atan2(p.yrad, p.xrad);
            if (p.angle < 0.0) {
                p.angle += 2.0 * kPi;
            }
            if (c.helix == anchor_helix) {
                icroot = static_cast<int>(conns.size());
            }
            conns.push_back(p);
        }
    }

    int find_middle(int icstart, int icend) const {
        int count = 0;
        int ic = icstart;
        int guard = 0;
        int ret = -1;
        bool done = false;
        while (!done) {
            asserta(++guard < static_cast<int>(conns.size()) * 3 + 2);
            if (icroot >= 0 && ic == icroot && anchor_helix >= 0) {
                ret = ic;
            }
            done = ic == icend;
            ic += 1;
            if (ic >= static_cast<int>(conns.size())) {
                ic = 0;
            }
            count += 1;
        }
        if (ret >= 0) {
            return ret;
        }
        ic = icstart;
        for (int i = 1; i < (count + 1) / 2; ++i) {
            ic += 1;
            if (ic >= static_cast<int>(conns.size())) {
                ic = 0;
            }
        }
        return ic;
    }

    void place_pairs() {
        int d = static_cast<int>(conns.size());
        int icstart = (icroot < 0) ? 0 : icroot;
        int count = 0;
        bool done = false;
        int imax = 0;
        do {
            int j = icstart - 1;
            if (j < 0) {
                j = d - 1;
            }
            if (!connected(j, icstart)) {
                done = true;
            } else {
                icstart = j;
            }
            if (++count > d) {
                double maxang = -1.0;
                for (int ic = 0; ic < d; ++ic) {
                    int jn = (ic + 1) % d;
                    double ac = conns[static_cast<size_t>(jn)].angle - conns[static_cast<size_t>(ic)].angle;
                    if (ac < 0.0) {
                        ac += 2.0 * kPi;
                    }
                    if (ac > maxang) {
                        maxang = ac;
                        imax = ic;
                    }
                }
                icstart = imax + 1;
                if (icstart >= d) {
                    icstart = 0;
                }
                done = true;
            }
        } while (!done);

        bool done_all = false;
        int icstart1 = icstart;
        while (!done_all) {
            int icend = icstart;
            count = 0;
            done = false;
            bool rooted = false;
            while (!done) {
                if (icend == icroot) {
                    rooted = true;
                }
                int j = (icend + 1) % d;
                if (connected(icend, j)) {
                    if (++count >= d) {
                        break;
                    }
                    icend = j;
                } else {
                    done = true;
                }
            }
            int icmiddle = find_middle(icstart, icend);
            int icup = icmiddle;
            int icdown = icmiddle;
            int direction = 0;
            done = false;
            while (!done) {
                int ic = icmiddle;
                if (direction < 0) {
                    ic = icup;
                } else if (direction > 0) {
                    ic = icdown;
                }
                if (ic >= 0) {
                    PConn& cp = conns[static_cast<size_t>(ic)];
                    bool is_anchor = (anchor_helix >= 0 && cp.helix == anchor_helix);
                    if (!is_anchor) {
                        if (direction == 0) {
                            double a0 = cp.angle - std::asin(1.0 / (2.0 * radius));
                            double a1 = cp.angle + std::asin(1.0 / (2.0 * radius));
                            px[static_cast<size_t>(cp.start)] = xc + radius * std::cos(a0);
                            py[static_cast<size_t>(cp.start)] = yc + radius * std::sin(a0);
                            px[static_cast<size_t>(cp.end)] = xc + radius * std::cos(a1);
                            py[static_cast<size_t>(cp.end)] = yc + radius * std::sin(a1);
                        } else if (direction < 0) {
                            int j = (ic + 1) % d;
                            PConn& next = conns[static_cast<size_t>(j)];
                            double ac = 0.5 * (cp.angle + next.angle);
                            if (cp.angle > next.angle) {
                                ac -= kPi;
                            }
                            double cnx = std::cos(ac);
                            double cny = std::sin(ac);
                            double lnx = cny;
                            double lny = -cnx;
                            double da = next.angle - cp.angle;
                            if (da < 0.0) {
                                da += 2.0 * kPi;
                            }
                            double rl = 1.0;
                            if (cp.extruded) {
                                rl = (da <= kPi / 2.0) ? 2.0 : 1.5;
                            }
                            px[static_cast<size_t>(cp.end)] =
                                px[static_cast<size_t>(next.start)] + rl * lnx;
                            py[static_cast<size_t>(cp.end)] =
                                py[static_cast<size_t>(next.start)] + rl * lny;
                            px[static_cast<size_t>(cp.start)] =
                                px[static_cast<size_t>(cp.end)] + cp.yrad;
                            py[static_cast<size_t>(cp.start)] =
                                py[static_cast<size_t>(cp.end)] - cp.xrad;
                        } else {
                            int j = ic - 1;
                            if (j < 0) {
                                j = d - 1;
                            }
                            PConn& prev = conns[static_cast<size_t>(j)];
                            double ac = 0.5 * (prev.angle + cp.angle);
                            if (prev.angle > cp.angle) {
                                ac -= kPi;
                            }
                            double cnx = std::cos(ac);
                            double cny = std::sin(ac);
                            double lnx = -cny;
                            double lny = cnx;
                            double da = cp.angle - prev.angle;
                            if (da < 0.0) {
                                da += 2.0 * kPi;
                            }
                            double rl = 1.0;
                            if (prev.extruded) {
                                rl = (da <= kPi / 2.0) ? 2.0 : 1.5;
                            }
                            px[static_cast<size_t>(cp.start)] =
                                px[static_cast<size_t>(prev.end)] + rl * lnx;
                            py[static_cast<size_t>(cp.start)] =
                                py[static_cast<size_t>(prev.end)] + rl * lny;
                            px[static_cast<size_t>(cp.end)] =
                                px[static_cast<size_t>(cp.start)] - cp.yrad;
                            py[static_cast<size_t>(cp.end)] =
                                py[static_cast<size_t>(cp.start)] + cp.xrad;
                            (void)cnx;
                        }
                    }
                }
                if (direction < 0) {
                    if (icdown == icend) {
                        icdown = -1;
                    } else if (icdown >= 0) {
                        icdown += 1;
                        if (icdown >= d) {
                            icdown = 0;
                        }
                    }
                    direction = 1;
                } else {
                    if (icup == icstart) {
                        icup = -1;
                    } else if (icup >= 0) {
                        icup -= 1;
                        if (icup < 0) {
                            icup = d - 1;
                        }
                    }
                    direction = -1;
                }
                done = icup == -1 && icdown == -1;
            }
            int icnext = (icend + 1) % d;
            bool full_turn = (icstart == icstart1 && icnext == icstart1);
            if (icend != icstart && !full_turn) {
                PConn& a = conns[static_cast<size_t>(icstart)];
                PConn& c = conns[static_cast<size_t>(icend)];
                double dx = px[static_cast<size_t>(c.end)] - px[static_cast<size_t>(a.start)];
                double dy = py[static_cast<size_t>(c.end)] - py[static_cast<size_t>(a.start)];
                double midx = px[static_cast<size_t>(a.start)] + dx * 0.5;
                double midy = py[static_cast<size_t>(a.start)] + dy * 0.5;
                double rr = std::sqrt(dx * dx + dy * dy);
                if (rr > 1e-12) {
                    double mx = dx / rr;
                    double my = dy / rr;
                    double vx = xc - midx;
                    double vy = yc - midy;
                    vx /= rr;
                    vy /= rr;
                    double dot = vx * mx + vy * my;
                    double nrx = dot * mx - vx;
                    double nry = dot * my - vy;
                    double nr = std::sqrt(nrx * nrx + nry * nry);
                    if (nr > 1e-12) {
                        nrx /= nr;
                        nry /= nr;
                        double ac = std::atan2(py[static_cast<size_t>(a.start)] - yc,
                                               px[static_cast<size_t>(a.start)] - xc);
                        if (ac < 0.0) {
                            ac += 2.0 * kPi;
                        }
                        double acn = std::atan2(py[static_cast<size_t>(c.end)] - yc,
                                                px[static_cast<size_t>(c.end)] - xc);
                        if (acn < 0.0) {
                            acn += 2.0 * kPi;
                        }
                        if (acn < ac) {
                            acn += 2.0 * kPi;
                        }
                        int sign = (acn - ac > kPi) ? -1 : 1;
                        double nmidx = xc + sign * radius * nrx;
                        double nmidy = yc + sign * radius * nry;
                        double sx = nmidx - midx;
                        double sy = nmidy - midy;
                        if (rooted) {
                            xc -= sx;
                            yc -= sy;
                        } else {
                            int ic = icstart;
                            while (true) {
                                PConn& cp = conns[static_cast<size_t>(ic)];
                                px[static_cast<size_t>(cp.start)] += sx;
                                py[static_cast<size_t>(cp.start)] += sy;
                                px[static_cast<size_t>(cp.end)] += sx;
                                py[static_cast<size_t>(cp.end)] += sy;
                                if (ic == icend) {
                                    break;
                                }
                                ic += 1;
                                if (ic >= d) {
                                    ic = 0;
                                }
                            }
                        }
                    }
                }
            }
            icstart = icnext;
            done_all = icstart == icstart1;
        }
    }

    bool crossed_and_marked() {
        int d = static_cast<int>(conns.size());
        for (int i = 0; i < d; ++i) {
            int j = (i + 1) % d;
            PConn& cp = conns[static_cast<size_t>(i)];
            PConn& next = conns[static_cast<size_t>(j)];
            double dx = px[static_cast<size_t>(cp.end)] - xc;
            double dy = py[static_cast<size_t>(cp.end)] - yc;
            double ac = std::atan2(dy, dx);
            if (ac < 0.0) {
                ac += 2.0 * kPi;
            }
            dx = px[static_cast<size_t>(next.start)] - xc;
            dy = py[static_cast<size_t>(next.start)] - yc;
            double acn = std::atan2(dy, dx);
            if (acn < 0.0) {
                acn += 2.0 * kPi;
            }
            if (acn < ac) {
                acn += 2.0 * kPi;
            }
            double dan = acn - ac;
            double dcp = next.angle - cp.angle;
            if (dcp <= 0.0) {
                dcp += 2.0 * kPi;
            }
            if (std::fabs(dan - dcp) > kPi && !cp.extruded) {
                int steps = station_steps(cp.end, next.start);
                bool eligible = steps > 1;
                if (!b->circular) {
                    eligible = steps != 1;
                }
                if (eligible) {
                    cp.extruded = true;
                    return true;
                }
            }
        }
        return false;
    }

    void place_unpaired() {
        int d = static_cast<int>(conns.size());
        for (int i = 0; i < d; ++i) {
            int j = (i + 1) % d;
            if (conns[static_cast<size_t>(i)].extruded) {
                construct_extruded(i);
                continue;
            }
            int end = conns[static_cast<size_t>(i)].end;
            int start = conns[static_cast<size_t>(j)].start;
            int m = station_steps(end, start);
            if (m <= 1) {
                continue;
            }
            double dx = px[static_cast<size_t>(end)] - xc;
            double dy = py[static_cast<size_t>(end)] - yc;
            double rc = std::sqrt(dx * dx + dy * dy);
            double ac = std::atan2(dy, dx);
            if (ac < 0.0) {
                ac += 2.0 * kPi;
            }
            dx = px[static_cast<size_t>(start)] - xc;
            dy = py[static_cast<size_t>(start)] - yc;
            double rcn = std::sqrt(dx * dx + dy * dy);
            double acn = std::atan2(dy, dx);
            if (acn < 0.0) {
                acn += 2.0 * kPi;
            }
            if (acn < ac) {
                acn += 2.0 * kPi;
            }
            double dan = acn - ac;
            if (dan < 1e-12) {
                place_on_chord(end, start);
                continue;
            }
            double step = dan / static_cast<double>(m);
            for (int k = 1; k < m; ++k) {
                int st = end + k;
                if (st >= N) {
                    st -= N;
                }
                double ang = ac + k * step;
                double rr = rc + (rcn - rc) * (ang - ac) / dan;
                px[static_cast<size_t>(st)] = xc + rr * std::cos(ang);
                py[static_cast<size_t>(st)] = yc + rr * std::sin(ang);
            }
        }
    }

    // Heap stack. A nested molecule can be thousands of helices deep, and this
    // frame is too large to recurse that far on the default stack.
    struct TravFrame {
        int loop_id;
        int anchor;
        int next_ic;
        int d;
        int root_ic;
        bool placed;
        std::vector<PConn> saved;
    };

    void traverse(int loop_id, int anchor) {
        std::vector<TravFrame> stack;
        stack.push_back(TravFrame{loop_id, anchor, 0, 0, -1, false, {}});
        while (!stack.empty()) {
            TravFrame& fr = stack.back();
            if (!fr.placed) {
                anchor_helix = fr.anchor;
                load_conns(fr.loop_id);
                asserta(!conns.empty());
                bool again = true;
                int spins = 0;
                while (again) {
                    asserta(++spins <= static_cast<int>(conns.size()) + 2);
                    determine_radius();
                    if (fr.anchor < 0) {
                        xc = 0.0;
                        yc = 0.0;
                    } else {
                        const PConn& acp = conns[static_cast<size_t>(icroot)];
                        double xo = 0.5 * (px[static_cast<size_t>(acp.start)] + px[static_cast<size_t>(acp.end)]);
                        double yo = 0.5 * (py[static_cast<size_t>(acp.start)] + py[static_cast<size_t>(acp.end)]);
                        xc = xo - radius * acp.xrad;
                        yc = yo - radius * acp.yrad;
                    }
                    place_pairs();
                    again = crossed_and_marked();
                }
                place_unpaired();
                fr.d = static_cast<int>(conns.size());
                fr.root_ic = icroot;
                fr.saved = conns;
                fr.placed = true;
                fr.next_ic = 0;
            }
            int ic = fr.next_ic;
            while (ic < fr.d && ic == fr.root_ic) {
                ic += 1;
            }
            if (ic >= fr.d) {
                stack.pop_back();
                continue;
            }
            fr.next_ic = ic + 1;
            int parent_loop = fr.loop_id;
            conns = fr.saved;
            generate_region(ic);
            fr.saved = conns;
            int hid = fr.saved[static_cast<size_t>(ic)].helix;
            int child = b->helices[static_cast<size_t>(hid)].inner_loop;
            if (child == parent_loop) {
                child = b->helices[static_cast<size_t>(hid)].outer_loop;
            }
            asserta(static_cast<int>(stack.size()) < static_cast<int>(b->loops.size()));
            stack.push_back(TravFrame{child, hid, 0, 0, -1, false, {}});
        }
    }
};

void place_no_pairs(NavLayout* out) {
    int n = static_cast<int>(out->x.size());
    if (n == 1) {
        out->x[0] = 0.0;
        out->y[0] = 0.0;
        return;
    }
    double r = 1.0 / (2.0 * std::sin(kPi / static_cast<double>(n)));
    for (int i = 0; i < n; ++i) {
        double t = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(n);
        out->x[static_cast<size_t>(i)] = r * std::cos(t);
        out->y[static_cast<size_t>(i)] = r * std::sin(t);
    }
}

void place_coordinates(const Builder& b, int root, NavLayout* out) {
    if (b.helices.empty()) {
        place_no_pairs(out);
        return;
    }
    Placer p;
    p.b = &b;
    p.N = b.circular ? b.n : b.n + 1;
    p.px.assign(static_cast<size_t>(p.N), std::nan(""));
    p.py.assign(static_cast<size_t>(p.N), std::nan(""));
    p.traverse(root, -1);
    for (int i = 0; i < b.n; ++i) {
        int st = p.to_station(i);
        asserta(std::isfinite(p.px[static_cast<size_t>(st)]));
        asserta(std::isfinite(p.py[static_cast<size_t>(st)]));
        out->x[static_cast<size_t>(i)] = p.px[static_cast<size_t>(st)];
        out->y[static_cast<size_t>(i)] = p.py[static_cast<size_t>(st)];
    }
}

// --- graph pictures ----------------------------------------------------

double node_weight(const TNode& n) {
    if (n.kids.empty()) {
        return 1.0;
    }
    double w = 0.0;
    for (size_t i = 0; i < n.kids.size(); ++i) {
        w += node_weight(n.kids[i]);
    }
    return w;
}

void place_wedge(const TNode& n, double x, double y, double a0, double a1, NavPicture* pic) {
    if (n.kids.empty()) {
        return;
    }
    double W = node_weight(n);
    double t = a0;
    double span = a1 - a0;
    for (size_t i = 0; i < n.kids.size(); ++i) {
        double w = span * node_weight(n.kids[i]) / W;
        double bis = t + 0.5 * w;
        double cx = x + std::cos(bis);
        double cy = y + std::sin(bis);
        NavLine line;
        line.x0 = x;
        line.y0 = y;
        line.x1 = cx;
        line.y1 = cy;
        pic->lines.push_back(line);
        NavPoint pt;
        pt.x = cx;
        pt.y = cy;
        pic->vertices.push_back(pt);
        place_wedge(n.kids[i], cx, cy, t, t + w, pic);
        t += w;
    }
}

void layout_tree(const TNode& root, bool bicenter, NavPicture* pic) {
    pic->vertices.clear();
    pic->lines.clear();
    if (bicenter) {
        asserta(root.kids.size() == 2);
        NavPoint a;
        a.x = 0.0;
        a.y = 0.0;
        NavPoint b;
        b.x = 1.0;
        b.y = 0.0;
        pic->vertices.push_back(a);
        pic->vertices.push_back(b);
        NavLine line;
        line.x0 = 0.0;
        line.y0 = 0.0;
        line.x1 = 1.0;
        line.y1 = 0.0;
        pic->lines.push_back(line);
        place_wedge(root.kids[0], 0.0, 0.0, kPi / 2.0, 3.0 * kPi / 2.0, pic);
        place_wedge(root.kids[1], 1.0, 0.0, -kPi / 2.0, kPi / 2.0, pic);
        return;
    }
    NavPoint origin;
    origin.x = 0.0;
    origin.y = 0.0;
    pic->vertices.push_back(origin);
    place_wedge(root, 0.0, 0.0, 0.0, 2.0 * kPi, pic);
}

struct DStem {
    char mark;
    int loops; // 1 or 2
};

struct DPiece {
    bool leaf;
    DStem stem;
    std::vector<DPiece> kids;
    std::vector<char> seps;
    bool closed;
};

bool parse_dual(const std::string& s, size_t* i, DPiece* out) {
    if (*i >= s.size()) {
        return false;
    }
    if (s[*i] == '@' || s[*i] == '*') {
        out->leaf = true;
        out->stem.mark = s[*i];
        out->stem.loops = 1;
        *i += 1;
        if (out->stem.mark == '@' && *i < s.size() && s[*i] == '@' && out->kids.empty()) {
            // '@@' is one stem only when it is the entire string. The caller
            // handles that before parsing. A second '@' here is the next piece.
        }
        return true;
    }
    if (s[*i] != '(') {
        return false;
    }
    *i += 1;
    out->leaf = false;
    out->closed = false;
    while (*i < s.size() && s[*i] != ')') {
        DPiece kid;
        if (!parse_dual(s, i, &kid)) {
            return false;
        }
        out->kids.push_back(std::move(kid));
        if (*i < s.size() && (s[*i] == '-' || s[*i] == '=' || s[*i] == '.')) {
            out->seps.push_back(s[*i]);
            *i += 1;
        } else {
            break;
        }
    }
    if (*i >= s.size() || s[*i] != ')') {
        return false;
    }
    *i += 1;
    out->closed = out->seps.size() == out->kids.size();
    return true;
}

void collect_stems(const DPiece& p, std::vector<DStem>* stems) {
    if (p.leaf) {
        stems->push_back(p.stem);
        return;
    }
    for (size_t i = 0; i < p.kids.size(); ++i) {
        collect_stems(p.kids[i], stems);
    }
}

void add_stem_line(NavPicture* pic, const NavPoint& a, const NavPoint& b) {
    NavLine line;
    line.x0 = a.x;
    line.y0 = a.y;
    line.x1 = b.x;
    line.y1 = b.y;
    pic->lines.push_back(line);
}

void add_parallel(NavPicture* pic, const NavPoint& a, const NavPoint& b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double len = std::sqrt(dx * dx + dy * dy);
    asserta(len > 0.0);
    dx /= len;
    dy /= len;
    double px = -dy;
    double py = dx;
    NavPoint a0;
    a0.x = a.x + 0.2 * px;
    a0.y = a.y + 0.2 * py;
    NavPoint b0;
    b0.x = b.x + 0.2 * px;
    b0.y = b.y + 0.2 * py;
    NavPoint a1;
    a1.x = a.x - 0.2 * px;
    a1.y = a.y - 0.2 * py;
    NavPoint b1;
    b1.x = b.x - 0.2 * px;
    b1.y = b.y - 0.2 * py;
    add_stem_line(pic, a0, b0);
    add_stem_line(pic, a1, b1);
}

void add_triangle(NavPicture* pic, const NavPoint& v, double ux, double uy) {
    double px = -uy;
    double py = ux;
    NavPoint A;
    A.x = v.x + 0.45 * ux;
    A.y = v.y + 0.45 * uy;
    NavPoint B;
    B.x = v.x + 0.25 * ux + 0.25 * px;
    B.y = v.y + 0.25 * uy + 0.25 * py;
    NavPoint C;
    C.x = v.x + 0.25 * ux - 0.25 * px;
    C.y = v.y + 0.25 * uy - 0.25 * py;
    add_stem_line(pic, A, B);
    add_stem_line(pic, B, C);
    add_stem_line(pic, C, A);
}

void draw_seps(const DPiece& piece, const std::vector<NavPoint>& stems, int* base, NavPicture* pic) {
    if (piece.leaf) {
        *base += 1;
        return;
    }
    std::vector<int> first(piece.kids.size());
    std::vector<int> last(piece.kids.size());
    for (size_t i = 0; i < piece.kids.size(); ++i) {
        first[i] = *base;
        draw_seps(piece.kids[i], stems, base, pic);
        last[i] = *base - 1;
    }
    for (size_t i = 0; i < piece.seps.size(); ++i) {
        char sep = piece.seps[i];
        if (sep == '.') {
            continue;
        }
        int ia = last[i];
        int ib = (i + 1 < piece.kids.size()) ? first[i + 1] : first[0];
        if (sep == '=') {
            add_parallel(pic, stems[static_cast<size_t>(ia)], stems[static_cast<size_t>(ib)]);
        } else {
            add_stem_line(pic, stems[static_cast<size_t>(ia)], stems[static_cast<size_t>(ib)]);
        }
    }
}

void layout_dual(const std::string& text, NavPicture* pic) {
    pic->vertices.clear();
    pic->lines.clear();
    if (text == "*" || text.empty()) {
        NavPoint p;
        p.x = 0.0;
        p.y = 0.0;
        pic->vertices.push_back(p);
        return;
    }
    if (text == "@" || text == "@@") {
        NavPoint p;
        p.x = 0.0;
        p.y = 0.0;
        pic->vertices.push_back(p);
        add_triangle(pic, p, 1.0, 0.0);
        if (text == "@@") {
            add_triangle(pic, p, -1.0, 0.0);
        }
        return;
    }
    size_t i = 0;
    DPiece root;
    asserta(parse_dual(text, &i, &root));
    asserta(i == text.size());
    std::vector<DStem> stems;
    collect_stems(root, &stems);
    int m = static_cast<int>(stems.size());
    std::vector<NavPoint> pos(static_cast<size_t>(m));
    if (m == 1) {
        pos[0].x = 0.0;
        pos[0].y = 0.0;
    } else {
        for (int k = 0; k < m; ++k) {
            double t = 2.0 * kPi * static_cast<double>(k) / static_cast<double>(m);
            pos[static_cast<size_t>(k)].x = std::cos(t);
            pos[static_cast<size_t>(k)].y = std::sin(t);
        }
    }
    for (int k = 0; k < m; ++k) {
        pic->vertices.push_back(pos[static_cast<size_t>(k)]);
    }
    int base = 0;
    draw_seps(root, pos, &base, pic);
    for (int k = 0; k < m; ++k) {
        if (stems[static_cast<size_t>(k)].mark != '@') {
            continue;
        }
        double ux;
        double uy;
        if (m == 1) {
            ux = 1.0;
            uy = 0.0;
        } else {
            ux = pos[static_cast<size_t>(k)].x;
            uy = pos[static_cast<size_t>(k)].y;
        }
        add_triangle(pic, pos[static_cast<size_t>(k)], ux, uy);
        if (stems[static_cast<size_t>(k)].loops >= 2) {
            add_triangle(pic, pos[static_cast<size_t>(k)], -ux, -uy);
        }
    }
}

TNode naview_tree(const Builder& b, int loop_id) {
    const LoopB& loop = b.loops[static_cast<size_t>(loop_id)];
    TNode node;
    node.leaf = loop.children.empty();
    node.mark = 0;
    int pivot = parent_conn_index(loop);
    int d = static_cast<int>(loop.conns.size());
    int begin = 0;
    std::vector<int> order(static_cast<size_t>(d));
    if (pivot < 0) {
        for (int i = 0; i < d; ++i) {
            order[static_cast<size_t>(i)] = i;
        }
    } else {
        for (int i = 0; i < d; ++i) {
            order[static_cast<size_t>(i)] = (pivot + i) % d;
        }
        begin = 1;
    }
    for (int i = begin; i < d; ++i) {
        int hid = loop.conns[static_cast<size_t>(order[static_cast<size_t>(i)])].helix;
        int child = b.helices[static_cast<size_t>(hid)].inner_loop;
        if (child == loop_id) {
            child = b.helices[static_cast<size_t>(hid)].outer_loop;
        }
        node.kids.push_back(naview_tree(b, child));
    }
    node.leaf = node.kids.empty();
    return node;
}

void publish(const Builder& b, NavLayout* out) {
    out->mate = b.mate;
    out->helices = b.helices;
    out->loops.clear();
    for (size_t i = 0; i < b.loops.size(); ++i) {
        const LoopB& src = b.loops[i];
        NavLoop dst;
        dst.kind = src.kind;
        dst.start = src.start;
        dst.end = src.end;
        dst.parent = src.parent;
        dst.closing_helix = src.closing;
        dst.child_helices = src.children;
        out->loops.push_back(dst);
    }
    out->x.assign(static_cast<size_t>(b.n), 0.0);
    out->y.assign(static_cast<size_t>(b.n), 0.0);
}

}  // namespace

bool naview_layout(const char* db, int n, bool circular, NavLayout* out, const char** err) {
    if (err) {
        *err = 0;
    }
    if (!out || !db || n < 1) {
        if (err) {
            *err = "empty structure";
        }
        return false;
    }
    Builder b;
    b.n = n;
    b.circular = circular;
    if (!build_pairs(&b, db, err)) {
        return false;
    }
    build_helices(&b);
    b.entered.assign(b.helices.size(), 0);
    int exterior = -1;
    if (b.helices.empty()) {
        LoopB loop;
        loop.parent = -1;
        loop.closing = -1;
        loop.kind = circular ? NavKind::Circle : NavKind::Exterior;
        loop.start = 0;
        loop.end = n - 1;
        b.loops.push_back(loop);
        exterior = 0;
    } else {
        exterior = build_loop(&b, circular ? 0 : -1);
        asserta(exterior == 0);
        fill_gaps(&b);
        classify(&b);
        set_spans(&b);
    }
    publish(b, out);
    build_segments(b, out);
    out->naview_graph = naview_graph_string(b, 0);

    Rag rag;
    rag.nodes.resize(b.loops.size());
    for (size_t i = 0; i < b.loops.size(); ++i) {
        rag.nodes[i].head = -1;
        rag.nodes[i].live = true;
        rag.nodes[i].is_end = (i == 0 && !circular);
        rag.nodes[i].extra = 0;
    }
    rag.edges.resize(b.helices.size());
    for (size_t hi = 0; hi < b.helices.size(); ++hi) {
        rag.edges[hi].a = b.helices[hi].outer_loop;
        rag.edges[hi].b = b.helices[hi].inner_loop;
        rag.edges[hi].len = helix_len(b.helices[hi]);
        rag.edges[hi].pa = -1;
        rag.edges[hi].pb = -1;
        rag.edges[hi].live = true;
    }
    for (size_t li = 0; li < b.loops.size(); ++li) {
        const LoopB& loop = b.loops[li];
        for (size_t ci = 0; ci < loop.conns.size(); ++ci) {
            const Conn& c = loop.conns[ci];
            int port = add_port(&rag, static_cast<int>(li), c.helix, c.real_gap);
            REdge& e = rag.edges[static_cast<size_t>(c.helix)];
            if (e.pa < 0) {
                e.pa = port;
            } else {
                e.pb = port;
            }
        }
        if (loop.conns.empty()) {
            rag.nodes[li].extra = n;
        }
    }
    contract(&rag);
    std::vector<AdjNode> adj = adjacency(rag);
    TNode rag_root;
    bool bicenter = false;
    int c0 = -1;
    int c1 = -1;
    out->rag_tree = rag_tree_string(adj, &rag_root, &bicenter, &c0, &c1);
    out->rag_dual = dual_string(adj);

    if (b.helices.empty()) {
        place_no_pairs(out);
    } else {
        int root = drawing_root(b);
        place_coordinates(b, root, out);
    }

    TNode nv = naview_tree(b, 0);
    layout_tree(nv, false, &out->naview_picture);
    layout_tree(rag_root, bicenter, &out->rag_picture);
    layout_dual(out->rag_dual, &out->dual_picture);
    return true;
}
