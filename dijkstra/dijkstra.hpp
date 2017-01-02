/**
 * @brief  最短経路問題(shortest paths problem)におけるDijkstraのアルゴリズム(Dijkstra's algorithm)を扱う
 * @date   2016/02/20 ~ 2016/03/12
 */



//****************************************
// インクルードガード
//****************************************

#ifndef __DIJKSTRA_HPP__
#define __DIJKSTRA_HPP__



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "../Graph/graph.hpp"



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の宣言
//****************************************

/**
 * @brief  すべての辺重みが非負であるという仮定の下で、Dijkstra(ダイクストラ)のアルゴリズム(Dijkstra's algorithm)は
 *         重み付き有向グラフG = (V, E)上の単一始点最短路問題を解く. ここでは各辺(u, v) ∈ Eについてw(u, v) >= 0を仮定する
 *
 * @note   Dijkstraのアルゴリズムは、始点sからの最短路重みが最終的に決定された頂点の集合Sを管理する
 *         アルゴリズムは、繰り返し、最小の最短路推定値を持つ頂点u ∈ V - Sを選択し、uをSに追加し、
 *         uから出るすべての辺を緩和する. ここではd値をキーとする頂点のmin優先度付きキューQを用いる
 *
 * @note   優先度付きキューの優先度更新を行わないため、優先度付きキューが空になるまでに行われる挿入の数はΟ(E)であるが、
 *         EXTRACT-MIN呼び出し時に、最短路の更新が行われないならば、無視をすることで、全体としての実行時間をΟ(ElgV)としている
 *
 * @param  const graph_t& G    非負の重み付き有向グラフG
 * @param  index_t        s    始点s
 * @return 始点sからの最短路重みが最終的に決定された頂点の集合S
 */
vertices_t dijkstra(const graph_t& G, index_t s);



/**
 * @brief  すべての辺重みが非負であるという仮定の下で、Dijkstra(ダイクストラ)のアルゴリズム(Dijkstra's algorithm)は
 *         重み付き有向グラフG = (V, E)上の単一始点最短路問題を解く. ここでは各辺(u, v) ∈ Eについてw(u, v) >= 0を仮定する
 *
 * @note   頂点v ∈ Vに対し、それぞれΟ(V)時間の操作を行うので、全体でΟ(V^2)時間を要す
 *
 * @param  const matrix_t& W    非負の重み付き有向グラフW
 * @param  index_t         s    始点s
 * @return 始点sからの最短路重みが最終的に決定された頂点の集合S
 */
vertices_t dijkstra(const matrix_t& W, index_t s);



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of __DIJKSTRA_HPP__
