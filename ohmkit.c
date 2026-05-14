/*
 * ohmkit — 电阻工具箱
 *
 * 功能：
 *   串并联等效电阻    (分数精确运算)
 *   Δ-Y 变换          (桥式电路求解)
 *   色环电阻          (编码 / 解码)
 *
 * 内部全部使用分数（long long 分子/分母 + GCD 约分），零浮点误差。
 *
 * 用法：
 *   ./ohmkit                                   交互 REPL
 *   ./ohmkit "4+((9||27)+(1||27))||27"         表达式计算
 *   ./ohmkit -v "(10+20)||30"                   详细步骤
 *   ./ohmkit delta 10 20 30                    Δ→Y 变换
 *   ./ohmkit wye 10 20 30                      Y→Δ 变换
 *   ./ohmkit bridge 1 2 3 4 5                  桥式电路
 *   ./ohmkit color red violet yellow            色环 → 阻值
 *   ./ohmkit findcolor 4700                     阻值 → 色环
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_EXPR_LEN 1024

/* ═══════════════════════════════════════════════════════════════
   第一部分：分数引擎
   ═══════════════════════════════════════════════════════════════ */

typedef struct {
    long long num;
    long long den;
} Frac;

static long long gcd(long long a, long long b) {
    a = llabs(a); b = llabs(b);
    while (b) { long long t = b; b = a % b; a = t; }
    return a;
}

static Frac frac_make(long long num, long long den) {
    if (den == 0)      { num = 0; den = 1; }
    if (den < 0)       { num = -num; den = -den; }
    long long g = gcd(num, den);
    return (Frac){ num / g, den / g };
}

static Frac frac_add(Frac a, Frac b) {
    return frac_make(a.num * b.den + b.num * a.den, a.den * b.den);
}

static Frac frac_mul(Frac a, Frac b) {
    return frac_make(a.num * b.num, a.den * b.den);
}

static Frac frac_div(Frac a, Frac b) {
    return frac_make(a.num * b.den, a.den * b.num);
}

/* 并联：a || b = ab / (a+b) */
static Frac frac_par(Frac a, Frac b) {
    if (a.num == 0 || b.num == 0) return frac_make(0, 1);
    return frac_make(a.num * b.num, a.num * b.den + b.num * a.den);
}

static int    frac_is_zero(Frac a)      { return a.num == 0; }
static double frac_to_double(Frac a)    { return (double)a.num / (double)a.den; }

static void frac_print(Frac a) {
    if (a.den == 1) printf("%lld", a.num);
    else            printf("%lld/%lld", a.num, a.den);
}

static Frac frac_from_double(double dv) {
    long long int_part = (long long)dv;
    double frac_part = dv - (double)int_part;
    if (frac_part < 0) { frac_part = -frac_part; int_part = -int_part; }
    if (fabs(frac_part) < 1e-12) return frac_make(int_part, 1);

    long long den = 1000000000LL;
    long long num = (long long)(frac_part * den + 0.5);
    Frac f = frac_make(num, den);
    return frac_add(f, frac_make(int_part, 1));
}

/* ═══════════════════════════════════════════════════════════════
   第二部分：Δ-Y 变换 & 桥式电路
   ═══════════════════════════════════════════════════════════════ */

/*
 *        R1                         Ra
 *       / \                        / \
 *      /   \                      /   \
 *   R3 ----- R2      ←→        Rc --- Rb
 *    (Delta)                    (Wye / 星形)
 *
 *  桥式电路（Wheatstone Bridge）:
 *      ┌──R1──┬──R3──┐
 *    A─┤      R5      ├─B
 *      └──R2──┴──R4──┘
 *   左三角 (R1, R2, R5) Δ→Y 后即化为纯串并联。
 */

typedef struct { Frac a, b, c; } DeltaWye;

/* Δ → Y：三边电阻 → 星形三臂 */
static DeltaWye delta_to_wye(Frac R1, Frac R2, Frac R3) {
    Frac sum = frac_add(frac_add(R1, R2), R3);
    if (frac_is_zero(sum)) return (DeltaWye){{0,1},{0,1},{0,1}};
    return (DeltaWye){
        frac_div(frac_mul(R1, R2), sum),   // Ra = R1·R2 / (R1+R2+R3)
        frac_div(frac_mul(R2, R3), sum),   // Rb = R2·R3 / (R1+R2+R3)
        frac_div(frac_mul(R3, R1), sum),   // Rc = R3·R1 / (R1+R2+R3)
    };
}

/* Y → Δ：星形三臂 → 三边电阻 */
static DeltaWye wye_to_delta(Frac Ra, Frac Rb, Frac Rc) {
    Frac num = frac_add(frac_add(frac_mul(Ra, Rb), frac_mul(Rb, Rc)),
                        frac_mul(Rc, Ra));
    if (frac_is_zero(Rc) || frac_is_zero(Ra) || frac_is_zero(Rb))
        return (DeltaWye){{0,1},{0,1},{0,1}};
    return (DeltaWye){
        frac_div(num, Rc),     // R1 = (Ra·Rb+Rb·Rc+Rc·Ra) / Rc
        frac_div(num, Ra),     // R2 = ... / Ra
        frac_div(num, Rb),     // R3 = ... / Rb
    };
}

/* 桥式电路求解：5 个电阻 R1~R5，取左三角 (R1,R2,R5) 做 Δ→Y */
static Frac bridge_solve(Frac R1, Frac R2, Frac R3, Frac R4, Frac R5) {
    /* 左三角 (R1,R2,R5) Δ→Y：Ra 接 A, Rb 接 D, Rc 接 C */
    DeltaWye dw = delta_to_wye(R1, R2, R5);

    /* 变换后：A—Ra—O, O—Rc—C—R3—B, O—Rb—D—R4—B */
    Frac path_c = frac_add(dw.c, R3);   // Rc + R3
    Frac path_d = frac_add(dw.b, R4);   // Rb + R4
    Frac mid    = frac_par(path_c, path_d);

    return frac_add(dw.a, mid);          // Ra + (path_c || path_d)
}

/* ═══════════════════════════════════════════════════════════════
   第三部分：色环电阻
   ═══════════════════════════════════════════════════════════════ */

static const char *digit_names[] = {
    "黑", "棕", "红", "橙", "黄", "绿", "蓝", "紫", "灰", "白"
};
static const char *digit_en[] = {
    "black", "brown", "red", "orange", "yellow",
    "green", "blue", "violet", "grey", "white"
};

/* 乘数：black=1, brown=10, ..., white=1G, gold=0.1, silver=0.01 */
static const double mult_table[] = {
    1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,   // 0-9
    0.1, 0.01                                            // 10=gold, 11=silver
};
static const char *mult_names[] = {
    "黑","棕","红","橙","黄","绿","蓝","紫","灰","白","金","银"
};

/* 误差：棕=1%, 红=2%, 绿=0.5%, 蓝=0.25%, 紫=0.1%, 金=5%, 银=10% */
static const double tol_table[]   = { 1, 2, 0.5, 0.25, 0.1, 5, 10 };
static const int   tol_index[]    = { 1, 2, 4, 5, 6, 10, 11 };
/* tol_index 对应 digit_en / digit_names 中的索引；金=10, 银=11 */

/* 识别颜色名（支持中英文、全称和部分缩写） */
static int match_color(const char *s) {
    if (!s || !*s) return -1;
    char buf[32]; int i = 0;
    while (*s && i < 31) buf[i++] = (char)tolower((unsigned char)*s++);
    buf[i] = '\0';

    /* 英文全称 / 缩写 */
    if (strncmp(buf, "black",  5)==0 || strcmp(buf, "bk")==0) return 0;
    if (strncmp(buf, "brown",  5)==0 || strcmp(buf, "br")==0) return 1;
    if (strncmp(buf, "red",    3)==0 || strcmp(buf, "rd")==0) return 2;
    if (strncmp(buf, "orange", 6)==0 || strcmp(buf, "or")==0) return 3;
    if (strncmp(buf, "yellow", 6)==0 || strcmp(buf, "ye")==0) return 4;
    if (strncmp(buf, "green",  5)==0 || strcmp(buf, "gr")==0) return 5;
    if (strncmp(buf, "blue",   4)==0 || strcmp(buf, "bl")==0) return 6;
    if (strncmp(buf, "violet", 6)==0 || strcmp(buf, "vi")==0) return 7;
    if (strncmp(buf, "grey",   4)==0 || strcmp(buf, "gy")==0) return 8;
    if (strncmp(buf, "gray",   4)==0)                           return 8;
    if (strncmp(buf, "white",  5)==0 || strcmp(buf, "wh")==0) return 9;
    if (strncmp(buf, "gold",   4)==0 || strcmp(buf, "gd")==0) return 10;
    if (strncmp(buf, "silver", 6)==0 || strcmp(buf, "sv")==0) return 11;

    /* 中文 */
    if (strcmp(buf, "黑")==0) return 0;
    if (strcmp(buf, "棕")==0) return 1;
    if (strcmp(buf, "红")==0) return 2;
    if (strcmp(buf, "橙")==0) return 3;
    if (strcmp(buf, "黄")==0) return 4;
    if (strcmp(buf, "绿")==0) return 5;
    if (strcmp(buf, "蓝")==0) return 6;
    if (strcmp(buf, "紫")==0) return 7;
    if (strcmp(buf, "灰")==0) return 8;
    if (strcmp(buf, "白")==0) return 9;
    if (strcmp(buf, "金")==0) return 10;
    if (strcmp(buf, "银")==0) return 11;

    return -1;
}

/* 色环 → 阻值 */
static void color_decode(int n_bands, const char **bands) {
    if (n_bands < 3 || n_bands > 5) {
        fprintf(stderr, "需要 3~5 个色环（3/4 环或 5 环）\n"); return;
    }

    int digits[5], idx;
    for (int i = 0; i < n_bands; i++) {
        idx = match_color(bands[i]);
        if (idx < 0) { fprintf(stderr, "未知颜色: %s\n", bands[i]); return; }
        digits[i] = idx;
    }

    double value, tol = 20.0;
    if (n_bands == 5) {
        /* 5 环：D1 D2 D3 MUL TOL */
        value = (digits[0]*100.0 + digits[1]*10.0 + digits[2]) * mult_table[digits[3]];
        for (int i = 0; i < 7; i++)
            if (tol_index[i] == digits[4]) { tol = tol_table[i]; }
    } else {
        /* 3/4 环：D1 D2 [MUL] [TOL] */
        if (n_bands == 3) {
            value = (digits[0]*10.0 + digits[1]) * mult_table[digits[2]];
        } else {
            value = (digits[0]*10.0 + digits[1]) * mult_table[digits[2]];
            for (int i = 0; i < 7; i++)
                if (tol_index[i] == digits[3]) { tol = tol_table[i]; }
        }
    }

    /* 转分数输出 */
    Frac fv = frac_from_double(value);

    printf("色环: ");
    for (int i = 0; i < n_bands; i++) {
        if (i > 0) printf(" ");
        if (digits[i] < 10)
            printf("%s", digit_names[digits[i]]);
        else
            printf("%s", mult_names[digits[i]]);
    }
    printf("\n阻值: "); frac_print(fv);

    if (value >= 1e9)      printf(" (%.4g GΩ)", value / 1e9);
    else if (value >= 1e6) printf(" (%.4g MΩ)", value / 1e6);
    else if (value >= 1e3) printf(" (%.4g kΩ)", value / 1e3);
    else                   printf(" (%.4g Ω)", value);

    printf("  误差: ±%g%%\n", tol);
}

/* 阻值 → 最近标准色环 */
static void color_encode(double target) {
    /* E24 标准值数列 */
    static const int e24[] = {
        10,11,12,13,15,16,18,20,22,24,27,30,33,36,39,43,47,51,56,62,68,75,82,91
    };

    /* 找量级和最接近的 E24 值 */
    double absv = fabs(target);
    if (absv < 1e-12) { printf("0 Ω = 黑 黑 黑\n"); return; }

    int exp10 = 0;
    double mant = absv;
    while (mant >= 100.0) { mant /= 10.0; exp10++; }
    while (mant < 1.0)    { mant *= 10.0; exp10--; }

    int best_i = 0; double best_diff = 1e30;
    for (int i = 0; i < 24; i++) {
        double diff = fabs(mant - e24[i]);
        if (diff < best_diff) { best_diff = diff; best_i = i; }
    }

    int val = e24[best_i];
    /* val 是两位数，如 47；实际值 = val × 10^exp10 */
    int d1 = val / 10;       // 第一位
    int d2 = val % 10;       // 第二位
    int mul = exp10;          // 乘数指数

    /* 输出 4 环（金 = ±5%） */
    printf("%s %s %s %s  →  ",
           digit_en[d1], digit_en[d2], mult_names[mul], "gold");
    printf("%s %s %s %s  =  %g Ω ±5%%\n",
           digit_names[d1], digit_names[d2], mult_names[mul], "金",
           val * pow(10, mul));
}

/* ═══════════════════════════════════════════════════════════════
   第四部分：表达式解析器（同 v3）
   ═══════════════════════════════════════════════════════════════ */

static Frac parse_expr(const char **s, int verbose, int depth);
static void skip_spaces(const char **s);
static void indent(int d) { while (d--) printf("  "); }

static Frac parse_atom(const char **s, int verbose, int depth) {
    skip_spaces(s);

    if (**s == '(') {
        (*s)++;
        Frac r = parse_expr(s, verbose, depth);
        skip_spaces(s);
        if (**s != ')') { fprintf(stderr, "语法错误：缺少 ')'\n"); exit(1); }
        (*s)++;
        return r;
    }

    if (isdigit((unsigned char)**s) || **s == '.') {
        char *end;
        double dv = strtod(*s, &end);
        if (end == *s) { fprintf(stderr, "语法错误：无效数字 '%c'\n", **s); exit(1); }
        *s = end;
        return frac_from_double(dv);
    }

    fprintf(stderr, "语法错误：意外字符 '%c'\n", **s); exit(1);
}

static Frac parse_term(const char **s, int verbose, int depth) {
    Frac result = parse_atom(s, verbose, depth);
    for (;;) {
        skip_spaces(s);
        if (**s == '|' && *(*s+1) == '|') {
            *s += 2;
            Frac rhs = parse_atom(s, verbose, depth);
            Frac nr = frac_par(result, rhs);
            if (verbose) {
                indent(depth); frac_print(result); printf(" || "); frac_print(rhs);
                printf(" = "); frac_print(nr);
                printf(" (≈ %.4f Ω)\n", frac_to_double(nr));
            }
            result = nr; continue;
        }
        break;
    }
    return result;
}

static Frac parse_expr(const char **s, int verbose, int depth) {
    Frac result = parse_term(s, verbose, depth);
    for (;;) {
        skip_spaces(s);
        if (**s == '+') {
            (*s)++;
            Frac rhs = parse_term(s, verbose, depth);
            Frac nr = frac_add(result, rhs);
            if (verbose) {
                indent(depth); frac_print(result); printf(" + "); frac_print(rhs);
                printf(" = "); frac_print(nr);
                printf(" (≈ %.4f Ω)\n", frac_to_double(nr));
            }
            result = nr; continue;
        }
        break;
    }
    return result;
}

static void skip_spaces(const char **s) {
    while (**s && isspace((unsigned char)**s)) (*s)++;
}

/* ═══════════════════════════════════════════════════════════════
   第五部分：交互模式（电路树构建）
   ═══════════════════════════════════════════════════════════════ */

typedef enum { NODE_SERIES, NODE_PARALLEL, NODE_LEAF } NodeType;

typedef struct Node {
    NodeType   type;
    Frac       value;
    int        child_count;
    struct Node *children[32];
} Node;

static Node *alloc_node(NodeType t) { Node *n = calloc(1, sizeof(Node)); n->type = t; return n; }
static void free_tree(Node *n) {
    if (!n) return;
    for (int i = 0; i < n->child_count; i++) free_tree(n->children[i]);
    free(n);
}

static Frac eval_tree_raw(Node *n) {
    if (!n) return frac_make(0,1);
    if (n->type == NODE_LEAF) return n->value;
    Frac acc = n->type == NODE_SERIES ? frac_make(0,1) : n->children[0]->value;
    for (int i = (n->type == NODE_SERIES ? 0 : 1); i < n->child_count; i++) {
        Frac cv = eval_tree_raw(n->children[i]);
        acc = (n->type == NODE_SERIES) ? frac_add(acc, cv) : frac_par(acc, cv);
    }
    return acc;
}

static Frac eval_tree(Node *n, int verbose, int depth) {
    if (!n) return frac_make(0,1);
    if (n->type == NODE_LEAF) {
        if (verbose) { indent(depth); printf("R = "); frac_print(n->value); printf("\n"); }
        return n->value;
    }
    Frac acc = frac_make(0,1); int has_zero = 0;
    for (int i = 0; i < n->child_count; i++) {
        Frac cv = eval_tree(n->children[i], verbose, depth+1);
        if (n->type == NODE_SERIES) acc = frac_add(acc, cv);
        else {
            if (frac_is_zero(cv)) has_zero = 1;
            else if (i == 0) acc = cv;
            else acc = frac_par(acc, cv);
        }
    }
    if (n->type == NODE_PARALLEL && has_zero) acc = frac_make(0,1);
    if (verbose) {
        indent(depth); printf("%s(", n->type==NODE_SERIES?"串联":"并联");
        for (int i = 0; i < n->child_count; i++) {
            if (i) printf(" + ");
            frac_print(eval_tree_raw(n->children[i]));
        }
        printf(") = "); frac_print(acc); printf("\n");
    }
    return acc;
}

static Node *build_subtree(void) {
    printf("\n  子电路类型？\n  1.单个电阻  2.串联组  3.并联组\n  选择: ");
    fflush(stdout);
    char buf[128];
    if (!fgets(buf, sizeof(buf), stdin)) exit(0);
    int ch = atoi(buf);
    if (ch == 1) {
        Node *n = alloc_node(NODE_LEAF);
        printf("  电阻值 (Ω): "); fflush(stdout);
        if (!fgets(buf, sizeof(buf), stdin)) exit(0);
        n->value = frac_from_double(atof(buf)); return n;
    }
    if (ch == 2 || ch == 3) {
        Node *n = alloc_node(ch == 2 ? NODE_SERIES : NODE_PARALLEL);
        printf("  包含几个子电路？"); fflush(stdout);
        if (!fgets(buf, sizeof(buf), stdin)) exit(0);
        int cnt = atoi(buf); if (cnt < 1) cnt = 1; if (cnt > 32) cnt = 32;
        for (int i = 0; i < cnt; i++) {
            printf("\n  [子电路 %d/%d]", i+1, cnt);
            n->children[n->child_count++] = build_subtree();
        }
        return n;
    }
    return alloc_node(NODE_LEAF);
}

static Node *interactive_build(void) {
    printf("=== 交互式电路构建 ===\n顶层: 1.电阻 2.串联 3.并联\n选择: ");
    fflush(stdout);
    char buf[128];
    if (!fgets(buf, sizeof(buf), stdin)) exit(0);
    int ch = atoi(buf);
    if (ch == 1) {
        Node *n = alloc_node(NODE_LEAF);
        printf("电阻值 (Ω): "); fflush(stdout);
        if (!fgets(buf, sizeof(buf), stdin)) exit(0);
        n->value = frac_from_double(atof(buf)); return n;
    }
    if (ch == 2 || ch == 3) {
        Node *n = alloc_node(ch == 2 ? NODE_SERIES : NODE_PARALLEL);
        printf("包含几个子电路？"); fflush(stdout);
        if (!fgets(buf, sizeof(buf), stdin)) exit(0);
        int cnt = atoi(buf); if (cnt < 1) cnt = 1; if (cnt > 32) cnt = 32;
        for (int i = 0; i < cnt; i++) {
            printf("\n[顶层子电路 %d/%d]", i+1, cnt);
            n->children[n->child_count++] = build_subtree();
        }
        return n;
    }
    return NULL;
}

/* ═══════════════════════════════════════════════════════════════
   第六部分：帮助 & 主入口
   ═══════════════════════════════════════════════════════════════ */

static void usage(const char *prog) {
    printf("ohmkit — 电阻工具箱 (精确分数运算)\n\n");
    printf("用法:\n");
    printf("  %s [表达式]                 串并联计算\n", prog);
    printf("  %s -v [表达式]              显示详细步骤\n", prog);
    printf("  %s delta  R1 R2 R3          Δ → Y 变换\n", prog);
    printf("  %s wye    Ra Rb Rc          Y → Δ 变换\n", prog);
    printf("  %s bridge R1 R2 R3 R4 R5    桥式电路等效电阻\n", prog);
    printf("  %s color  色1 色2 色3 ...   色环 → 阻值\n", prog);
    printf("  %s findcolor 4700           阻值 → 色环\n", prog);
    printf("\n进入交互 REPL 后，可用命令:\n");
    printf("  delta / wye / bridge / color / findcolor / build / quit\n");
}

int main(int argc, char **argv) {
    /* ── 命令行子命令分发 ── */
    if (argc >= 2) {
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
            usage(argv[0]); return 0;
        }

        if (strcmp(argv[1], "delta") == 0 && argc == 5) {
            Frac r1 = frac_from_double(atof(argv[2]));
            Frac r2 = frac_from_double(atof(argv[3]));
            Frac r3 = frac_from_double(atof(argv[4]));
            DeltaWye dw = delta_to_wye(r1, r2, r3);
            printf("Δ (%s, %s, %s) → Y\n",
                   argv[2], argv[3], argv[4]);
            printf("  Ra = "); frac_print(dw.a); printf(" (≈ %.4f Ω)\n", frac_to_double(dw.a));
            printf("  Rb = "); frac_print(dw.b); printf(" (≈ %.4f Ω)\n", frac_to_double(dw.b));
            printf("  Rc = "); frac_print(dw.c); printf(" (≈ %.4f Ω)\n", frac_to_double(dw.c));
            return 0;
        }

        if (strcmp(argv[1], "wye") == 0 && argc == 5) {
            Frac ra = frac_from_double(atof(argv[2]));
            Frac rb = frac_from_double(atof(argv[3]));
            Frac rc = frac_from_double(atof(argv[4]));
            DeltaWye dw = wye_to_delta(ra, rb, rc);
            printf("Y (%s, %s, %s) → Δ\n",
                   argv[2], argv[3], argv[4]);
            printf("  R1 = "); frac_print(dw.a); printf(" (≈ %.4f Ω)\n", frac_to_double(dw.a));
            printf("  R2 = "); frac_print(dw.b); printf(" (≈ %.4f Ω)\n", frac_to_double(dw.b));
            printf("  R3 = "); frac_print(dw.c); printf(" (≈ %.4f Ω)\n", frac_to_double(dw.c));
            return 0;
        }

        if (strcmp(argv[1], "bridge") == 0 && argc == 7) {
            Frac r1 = frac_from_double(atof(argv[2]));
            Frac r2 = frac_from_double(atof(argv[3]));
            Frac r3 = frac_from_double(atof(argv[4]));
            Frac r4 = frac_from_double(atof(argv[5]));
            Frac r5 = frac_from_double(atof(argv[6]));
            printf("桥式电路 (R1=%s, R2=%s, R3=%s, R4=%s, R5=%s)\n",
                   argv[2],argv[3],argv[4],argv[5],argv[6]);
            printf("  将左三角 (R1, R2, R5) Δ → Y ...\n");
            DeltaWye dw = delta_to_wye(r1, r2, r5);
            printf("  Ra = "); frac_print(dw.a); printf("  Rb = "); frac_print(dw.b);
            printf("  Rc = "); frac_print(dw.c); printf("\n");
            Frac result = bridge_solve(r1, r2, r3, r4, r5);
            printf("  等效电阻 = "); frac_print(result);
            printf(" (≈ %.6f Ω)\n", frac_to_double(result));
            return 0;
        }

        if (strcmp(argv[1], "color") == 0 && argc >= 5) {
            color_decode(argc - 2, (const char **)(argv + 2));
            return 0;
        }

        if (strcmp(argv[1], "findcolor") == 0 && argc == 3) {
            color_encode(atof(argv[2]));
            return 0;
        }

        /* 否则当作表达式：支持 -v */
        int verbose = (strcmp(argv[1], "-v") == 0);
        const char *expr = verbose ? argv[2] : argv[1];
        if (!expr) { usage(argv[0]); return 0; }

        const char *p = expr;
        Frac result = parse_expr(&p, verbose, 0);
        skip_spaces(&p);
        if (*p) fprintf(stderr, "警告：尾部未解析: '%s'\n", p);
        printf("\n等效电阻 = "); frac_print(result);
        printf(" (≈ %.6f Ω)\n", frac_to_double(result));
        return 0;
    }

    /* ── 交互 REPL ── */
    printf("ohmkit — 电阻工具箱 (精确分数运算)\n");
    printf("表达式 | delta | wye | bridge | color | findcolor | build | quit\n\n");

    char line[MAX_EXPR_LEN];
    for (;;) {
        printf(">>> ");
        fflush(stdout);
        if (!fgets(line, sizeof(line), stdin)) break;
        line[strcspn(line, "\r\n")] = '\0';
        if (!*line) continue;
        if (strcmp(line, "quit") == 0 || strcmp(line, "exit") == 0) break;
        if (strcmp(line, "help") == 0) { usage(argv[0]); continue; }

        /* 子命令 */
        char cmd[32] = {0}; double v[5] = {0};
        int n = sscanf(line, "%31s %lf %lf %lf %lf %lf", cmd, &v[0],&v[1],&v[2],&v[3],&v[4]);

        if (strcmp(cmd, "delta") == 0 && n >= 4) {
            Frac r1=frac_from_double(v[0]), r2=frac_from_double(v[1]), r3=frac_from_double(v[2]);
            DeltaWye dw = delta_to_wye(r1, r2, r3);
            printf("Δ→Y:  Ra="); frac_print(dw.a);
            printf("  Rb="); frac_print(dw.b);
            printf("  Rc="); frac_print(dw.c);
            printf("\n     (≈ %.4f, %.4f, %.4f Ω)\n",
                   frac_to_double(dw.a),frac_to_double(dw.b),frac_to_double(dw.c));
            continue;
        }
        if (strcmp(cmd, "wye") == 0 && n >= 4) {
            Frac ra=frac_from_double(v[0]), rb=frac_from_double(v[1]), rc=frac_from_double(v[2]);
            DeltaWye dw = wye_to_delta(ra, rb, rc);
            printf("Y→Δ:  R1="); frac_print(dw.a);
            printf("  R2="); frac_print(dw.b);
            printf("  R3="); frac_print(dw.c);
            printf("\n     (≈ %.4f, %.4f, %.4f Ω)\n",
                   frac_to_double(dw.a),frac_to_double(dw.b),frac_to_double(dw.c));
            continue;
        }
        if (strcmp(cmd, "bridge") == 0 && n >= 6) {
            Frac r1=frac_from_double(v[0]), r2=frac_from_double(v[1]), r3=frac_from_double(v[2]);
            Frac r4=frac_from_double(v[3]), r5=frac_from_double(v[4]);
            Frac result = bridge_solve(r1, r2, r3, r4, r5);
            DeltaWye dw = delta_to_wye(r1, r2, r5);
            printf("左三角 Δ→Y: Ra(接A)="); frac_print(dw.a); printf(" Rb="); frac_print(dw.b);
            printf(" Rc="); frac_print(dw.c); printf("\n");
            printf("等效电阻 = "); frac_print(result);
            printf(" (≈ %.6f Ω)\n", frac_to_double(result));
            continue;
        }
        if (strcmp(cmd, "color") == 0) {
            /* 手动解析颜色名 */
            char *save; const char *bands[5]; int nb = 0;
            char *tok = strtok_r(line, " ", &save);
            tok = strtok_r(NULL, " ", &save);  // 跳过 "color"
            while (tok && nb < 5) { bands[nb++] = tok; tok = strtok_r(NULL, " ", &save); }
            color_decode(nb, bands);
            continue;
        }
        if (strcmp(cmd, "findcolor") == 0 && n >= 2) {
            color_encode(v[0]);
            continue;
        }
        if (strcmp(line, "build") == 0) {
            Node *root = interactive_build();
            if (root) {
                printf("\n=== 计算过程 ===\n");
                Frac r = eval_tree(root, 1, 0);
                printf("\n等效电阻 = "); frac_print(r);
                printf(" (≈ %.6f Ω)\n", frac_to_double(r));
                free_tree(root);
            }
            printf("\n"); continue;
        }

        /* 默认：表达式解析 */
        const char *p = line;
        printf("\n");
        Frac result = parse_expr(&p, 0, 0);
        skip_spaces(&p);
        if (*p) fprintf(stderr, "警告：未解析: '%s'\n", p);
        printf("等效电阻 = "); frac_print(result);
        printf(" (≈ %.6f Ω)\n\n", frac_to_double(result));
    }
    printf("再见。\n");
    return 0;
}
