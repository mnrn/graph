/**
 * @brief  有向非巡回グラフにおける単一始点最短路問題を扱う
 * @date   2016/02/18
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iostream>
#include "../topological_sort/tsort.hpp"
#include "dsp.hpp"



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の定義
//***************************************

/**
 * @brief  重み付き有向非巡回グラフG = (V, E)の辺を頂点のトポロジカルソート順に緩和することで、
 *         単一始点からすべての頂点に至るすべての最短路を全体でΘ(V + E)時間で計算できる
 *         負辺があっても負閉路が存在しないから、最短路は有向非巡回グラフ上では常に明確に定義される
 *
 * @note   このアルゴリズムはまず有効非巡回グラフに対してトポロジカルソートを行い、頂点を線形に順序づける
 *         頂点uから頂点vへの道があれば、トポロジカルソートではuはvに先行する
 *         このアルゴリズムは全頂点をトポロジカルソート順に一度だけ走査する. 各頂点を検討するとき、その頂点から出るすべての辺を緩和する
 *
 * @note   このアルゴリズムの総実行時間はΘ(V + E)である
 *
 * @param  const graph_t& G 重み付き有向非巡回グラフG
 * @param  index_t        s 始点s
 */
std::pair<indices_t, array_t>  dsp(const graph_t& G, index_t s)
{
    index_t n = G.size();
    indices_t pi(n); array_t d(n);

    // Θ(V)の手続きによって最短路推定値と先行点を初期化する
    auto initsinglesource = [](indices_t& pi, array_t& d, index_t s, index_t n) -> void {
        for (index_t i = 0; i < n; i++) { d[i]  = limits::inf; pi[i] = limits::nil;}
        d[s] = 0;
    };
    // 辺(u, v)の緩和(relaxing)はuを経由することでvへの既知の最短路が改善できるか否か判定し、改善できるならばv.dとv.πを更新する
    // 緩和によって最短路推定値v.dが減少し、vの先行点属性v.πが更新されることがある. 以下のコードは、辺(u, v)上の緩和をΟ(1)時間で実行する
    auto relax = [](const edge& e, indices_t& pi, array_t& d) -> void {
        index_t v = e.dst, u = e.src;
        if (d[u] != limits::inf && d[v] > d[u] + e.w) { d[v] = d[u] + e.w; pi[v] = u; }
    };

    
    array_t sorted = tsort(G);      // Gの頂点をトポロジカルソートする
    initsinglesource(pi, d, s, n);  // 最短路推定値と先行点を初期化
    for (auto& u : sorted) {        // for トポロジカルソート順に、各頂点u(頂点ごとに1回ずつ実行)
        for (auto e : G[u]) {       // for 各頂点 v ∈ G.Adj[u](全体として各辺をちょうど1回ずつ)
            relax(e, pi, d);        // 緩和する
        }
    }
    return std::make_pair(pi, d);
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
    const int n = 6;
    graph::graph_t G(n);
    const int m = 10;

    for (int i = 0; i < m; i++) {
        int a, b, c;
        std::cin >> a >> b >> c;
        G[a].emplace_back(a, b, c);
    }
    for (auto& j : G) {
        for (auto& i : j) {
            std::cout << i.src << " " << i.dst << " " << i.w << std::endl;
        }
    }

    auto p = dsp(G, 1);
    for (int i = 0; i < n; i++) {
        std::cout << p.second[i] << std::endl;
    }
    

    return 0;
}
