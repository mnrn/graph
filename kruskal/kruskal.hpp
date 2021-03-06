/**
 * @brief 最小全域木問題(minimum-spanning tree problem)における
 *        Kruskalのアルゴリズム(Kruskal's algorithms)の宣言を行う
 * @date  2016/02/16 ~ 2016/03/12
 */



//****************************************
// インクルードガード
//****************************************

#ifndef KRUSKAL_HPP
#define KRUSKAL_HPP



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "../graph/graph.hpp"



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の宣言
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
std::pair<edges_t, weight_t> kruskal(const graph_t& G);



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of KRUSKAL_HPP
