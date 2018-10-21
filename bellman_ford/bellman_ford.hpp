/**
 * @brief  最短経路問題(shortest paths problem)におけるBellman-Fordアルゴリズム(Bellman-Ford algorithm)を扱う
 * @date   2016/02/18 ~ 2016/03/12
 */



//****************************************
// インクルードガード
//****************************************

#ifndef BELLMAN_FORD_HPP
#define BELLMAN_FORD_HPP



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "../graph/graph.hpp"



//****************************************
// 関数の宣言
//****************************************

GRAPH_BEGIN

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
 * @param  const graph_t& G グラフG
 * @param  index_t        s 始点s
 */
std::pair<bool, vertices_t> bellman_ford(const graph_t& G, index_t s);



GRAPH_END



#endif  // end of BELLMANFORD_HPP
