/**
 * @brief 単一始点最短路問題をDijkstra法で解く
 */


// ********************************************************************************
// Include files
// ********************************************************************************
#include <iostream>
#include <limits>
#include <vector>
#include <queue>


// ********************************************************************************
// Type synonyms
// ********************************************************************************
using weight_t = int32_t;  // 辺(u, v)への重みwを表す型
using index_t  = int32_t;  // 頂点vの添字を表す型


// ********************************************************************************
// Constants
// ********************************************************************************
static constexpr weight_t INF = std::numeric_limits<weight_t>::max() / 3;  // 辺が存在しない場合に使用される値
static constexpr index_t  NIL = std::numeric_limits<index_t>::min()  / 3;  // 先行点が存在しない場合に使用される値


// ********************************************************************************
// Structs
// ********************************************************************************
/**< @brief グラフ用ノード(頂点) */
struct vertex {
    weight_t d;     // 始点sからの距離
    index_t pi;     // 先行頂点(の添字)
};

/**
 * @brief グラフ用エッジ(辺)
 * @note  G=(V, E)を重み関数wを持つ重み付きグラフとすると、
 *        辺(u,v)∈Eの重みはw(u, v)と表される。
 */
struct edge {
    index_t   src;  // 辺の始点u
    index_t   dst;  // 辺の終点v
    weight_t  w;    // 辺(u, v)の重み(コスト)
};

/**< @brief 最短路推定値をもつ頂点 (∈V-SでVはグラフの頂点集合, Sは始点sからの最短路重みが最終的に決定された頂点集合) */
struct state {
    index_t  u; // G.Vに属する頂点u
    weight_t d; // 始点sからの距離d

    // <演算子オーバーロード
    bool operator < (const state& s) const { 
        return d > s.d;     // min優先度付きキューのためのREVERSE
    }

    state(index_t u, weight_t d) : u(u), d(d) {}
};


// ********************************************************************************
// Type synonyms (part2)
// ********************************************************************************
using edges_t       = std::vector<edge>;    // グラフG=(V, E)の辺集合E
using vertices_t    = std::vector<vertex>;  // グラフG=(V, E)の頂点集合V
using graph_t       = std::vector<edges_t>; // グラフGの隣接リスト表現


// ********************************************************************************
// Functions
// ********************************************************************************
/**
 * @brief           最短路推定値dと先行点piを初期化する
 * @param[out]      vertices_t& V      グラフの頂点集合V
 * @param[in]       index_t     s      始点s
 */
void init(vertices_t& V, index_t s)
{
    for (auto&& v : V) {
        v.d     = INF;
        v.pi    = NIL;
    }

    V[s].d = 0;
}

/**
 * @brief           uを経由することで、vへの既知の最短路を改善できるか判定し、
 *                  できるならば、v.dとv.piを更新する。
 *                  また緩和(更新)された場合、優先度付きキューQに挿入しておく。
 *
 * @param[inout]    vertices_t&                 V   グラフの頂点集合V
 * @param[in]       const edge&                 e   辺(u, v)
 * @param[out]      std::priority_queue<state>& Q   min優先度付きキュー
 */
void relax(vertices_t& V, const edge& e, std::priority_queue<state>& Q)
{
    index_t u = e.src;
    index_t v = e.dst;

    // 既知の最短路を改善できるか判定する。
    if (V[u].d + e.w < V[v].d) {   // 改善できるならば、

        // 緩和を行う
        V[v].d      = V[u].d + e.w;
        V[v].pi     = u;

        Q.emplace(v, V[v].d);
    }
}

/**
 * @brief           単一始点最短路問題をDijkstra法で解く。
 * @param[in]       const graph_t&  G  非負の重み付き有向グラフG
 * @param[in]       index_t         s  始点s
 * @return          始点sからの最短路重みが最終的に決定された頂点の集合S
 *
 * @details         Dijkstraのアルゴリズムは、始点sからの最短路重みが最終的に決定された頂点の集合Sを管理する。
 *                  アルゴリズムは繰り返し、最小の最短路推定値を持つ頂点u∈V-Sを選択し、uをSに追加し、
 *                  uから出るすべての辺を緩和する。
 */
vertices_t dijkstra(const graph_t& G, index_t s)
{
    vertices_t S(G.size());
    std::priority_queue<state> Q;

    init(S, s);             // すべての頂点のd値とpi値を初期化する。
    Q.emplace(s, S[s].d);   // このループの最初の実行ではu=sである。
    while (!Q.empty()) {

        // Qから最小の最短路推定値を持つ頂点u(∈V-S)を取得する
        state p = Q.top();
        Q.pop();
        index_t  u = p.u;
        weight_t d = p.d;

        // 頂点uの隣接リストに関してループをまわす
        for (const auto& e : G[u]) {    // 頂点uからでる辺(u, v)をそれぞれ緩和し、
            relax(S, e, Q);             // uを経由することでvへの最短路が改善できる場合にはv.dとv.piを更新する
        }
    }

    return S;
}


// ********************************************************************************
// Entry point
// ********************************************************************************

int main()
{
    graph_t G;
    int n, k, u, v, c;
    
    std::cin >> n;
    G.resize(n);
    for (int i = 0; i < n; i++) {
        std::cin >> u >> k;
        G[u].resize(k);
        for (int j = 0; j < k; j++) {
            std::cin >> v >> c;
            G[u][j].src = u;
            G[u][j].dst = v;
            G[u][j].w   = c;
        }
    }

    vertices_t S = dijkstra(G, 0);
    for (int i = 0; i < n; i++) {
        std::cout << i << " " << S[i].d << std::endl; 
    }
    return 0;
}
