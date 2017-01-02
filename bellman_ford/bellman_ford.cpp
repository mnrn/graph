/**
 * @brief  最短経路問題(shortest paths problem)におけるBellman-Fordアルゴリズム(Bellman-Ford algorithm)を扱う
 * @date   2016/02/18 ~ 2016/03/12
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iostream>
#include "../graph/relax.hpp"
#include "bellman_ford.hpp"


//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の定義
//****************************************

/**
 * @brief  Bellman-Fordアルゴリズム
 *
 * @note   Bellman-Fordアルゴリズムは負辺の存在を許す一般の単一始点最短路問題を解く
 *         重み関数w: E -> Rを持つ重み付き有向グラフG = (V, E)と始点sが与えられたとき、
 *         始点から到達可能な負閉路が存在するか否かを示すブール値をBellman-Fordアルゴリズムは返す
 *         アルゴリズムは、このような負閉路が存在すれば解が存在しないと報告し、そうでなければ最短路とその重みを生成する
 *
 * @note   Bellman-Fordアルゴリズムは、辺を次々に緩和することで、始点sから各頂点v ∈ Vへの最短路重みの推定値v.dを
 *         実際の最短路重みδ(s, v)に一致するまで徐々に減らす
 *         アルゴリズムが値TRUEを返すのは、グラフの始点から到達可能な負閉路を含まないとき、かつそのときに限る
 *
 * @note   Bellman-FordアルゴリズムはΟ(VE)時間で走る
 *
 * @param  const graph_t& G グラフG
 * @param  index_t        s 始点s
 */
std::pair<bool, vertices_t> bellman_ford(const graph_t& G, index_t s)
{   
    index_t n = G.size();
    vertices_t V(n);
    auto relax_pred = [] (const vertices_t& V, index_t u) -> bool { return V[u].d != limits::inf; };

    
    initialize_single_source(V, s);  // すべての頂点のd値とπ値を初期化する
    // アルゴリズムはグラフのすべての辺を|V| - 1回走査する
    for (index_t i = 0; i < n - 1; i++) {
        for (auto&& es : G) { for (auto&& e : es) {
                relax(V, e, relax_pred);  // グラフの各辺をそれぞれ1回緩和する
            }
        }
    }
    // Gがsから到達可能な負閉路を含まなければ、終了時に、すべての辺(u, v)に対して、
    //   u.d = δ(s, v)
    //      <= δ(s, u) + w(u, v) (∵ 三角不等式)
    //       = u.d + w(u, v)
    // だから、BELLMAN-FORDは値FALSEを返すことはなく、TRUEを返す
    for (auto&& es : G) { for (auto&& e : es) {               // 負の重みを持つ閉路の有無を判定する
            index_t v = e.dst, u = e.src;
            if (V[u].d != limits::inf && V[v].d > V[u].d + e.w) {  // Gが始点sから到達可能な負閉路を含むとき、
                return std::make_pair(false, V);       // FALSEを返す
            }
        }
    }
    // Gがsから到達可能な負閉路を含まなければ、値TRUEを返し、すべての頂点v ∈ Vに対してδ(s, v)が成り立ち、
    // 先行点部分グラフGπはsを根とする最短路木である
    return std::make_pair(true, V);
}



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



