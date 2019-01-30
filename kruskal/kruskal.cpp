/**
 * @brief 最小全域木問題(minimum-spanning tree problem)における
 *        Kruskalのアルゴリズム(Kruskal's algorithms)の実装を行う
 * @date  2016/02/16 ~ 2016/03/12
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "kruskal.hpp"
#include "./disjoint_sets/disjoint_sets.hpp"
#include <iostream>



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の定義
//****************************************

/**
 * @brief  Kruskalのアルゴリズム
 *
 * @note   Kruskalのアルゴリズムでは、集合Aは与えられたグラフの頂点集合を頂点集合とする森である
 *         Aに加える安全な辺は、常に2つの異なる連結成分を連結するグラフの最小重み辺である
 *
 * @note   Kruskalのアルゴリズムでは、成長させる森に付け加える安全な辺は、森に属する2つの木を連結するすべての辺の中で、重みが最小の辺(u, v)である
 *         (u, v)が連結する2つの木をC1およびC2とする. (u, v)はC1を別の木と連結する軽い辺だから、(u, v)はこの森に対して安全な辺である
 *         各ステップで重みが可能な限り小さい辺を森に加えているから、Kruskalのアルゴリズムは貪欲アルゴリズムである         
 *
 * @note   Kruskalのアルゴリズムの総実行時間はΟ(ElgV)である
 *
 * @param  const graph& G グラフG
 * @return 辺集合Aとその重み(最小全域木の重み)
 */
std::pair<edges_t, weight_t> kruskal(const graph_t& G)
{
    const index_t n = G.size();
    disjoint_sets ds(static_cast<std::size_t>(n));  // 互いな素な集合族のためのデータ構造を準備
    struct cmp { bool operator()(const edge& e, const edge& f) { return e.w < f.w; } };
    edges_t E;  // グラフGから集合G.Eを取り出す
    for (auto&& es : G) { for (auto&& e : es) { E.push_back(e); } }

    weight_t w = 0; edges_t A;             // Aを空集合に初期化し、
    for (index_t v = 0; v < n; ++v) {      // 各頂点がそれぞれ1つの木である|V|本の木を生成する
        ds.make_set(v);
    }
    std::sort(E.begin(), E.end(), cmp());  // 重みwの非減少順でG.Eの辺をソートする
    for (auto&& e : E) { // 辺を重みの小さいものから順に検討する
        // このループでは、各辺(u, v)について、端点uとvが同じ木に属するかどうかを調べる
        if (ds.find_set(e.src) != ds.find_set(e.dst)) {  // 両端点が同じ木に属さないならば、
            A.push_back(e); w += e.w;                    // 辺(u, v)をAに加え、
            ds.merge(e.src, e.dst);                      // 2つの木の頂点集合をマージする
        }
    }
    return std::make_pair(A, w);
}



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



//****************************************
// エントリポイント
//****************************************

int main(void)
{
    using namespace std;
    using namespace graph;
    
    const int n = 9;

    graph_t G; G.resize(n);

    int adjmtx[n][n] = {
        // 0,  1,  2,  3,  4,  5,  6,  7,  8
        { -1,  4, -1, -1, -1, -1, -1,  8, -1, },  // 0
        {  4, -1,  8, -1, -1, -1, -1, 11, -1, },  // 1
        { -1,  8, -1,  7, -1,  4, -1, -1,  2, },  // 2
        { -1, -1,  7, -1,  9, 14, -1, -1, -1, },  // 3
        { -1, -1, -1,  9, -1, 10, -1, -1, -1, },  // 4
        { -1, -1,  4, 14, 10, -1,  2, -1, -1, },  // 5
        { -1, -1, -1, -1, -1,  2, -1,  1,  6, },  // 6
        {  8, 11, -1, -1, -1, -1,  1, -1,  7, },  // 7
        { -1, -1,  2, -1, -1, -1,  6,  7, -1, },  // 8
    };

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = adjmtx[i][j];
            if (k != -1) G[i].emplace_back(i, j, k);
        }
    }

    auto mst = kruskal(G);
    cout << mst.second << endl;
    
    return 0;
}

