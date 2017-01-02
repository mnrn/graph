/**
 * @brief トポロジカルソート
 * @date  2016/02/14 ~ 2016/03/12
 */



//****************************************
// インクルードガード
//****************************************

#ifndef __TSORT_HPP__
#define __TSORT_HPP__



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
 * @brief  トポロジカルソートを行います
 *
 * @note   有向非巡回グラフG = (V, E)のトポロジカルソート(topological sort)は頂点集合上の線形順序で、
 *         Gが辺(u, v)を含むならば、この線形順序でuがvより先に現れるものである. (グラフに巡回路があればこのような線形順序は存在しない)
 *         グラフのトポロジカルソートは、すべての有向辺が左から右へ向かう、水平線上での頂点の並べ方である
 *
 * @note   つぎの簡単なアルゴリズムは有向非巡回グラフをトポロジカルソートする
 *
 *         TOPOLOGICAL-SORT(G)
 *         1 各頂点の終了時刻v.fを計算するためにDFS(G)を呼び出す
 *         2 ある頂点の探索が終了するたびに、この頂点を連結リストの先頭に挿入する
 *         3 return 頂点の連結リスト
 *
 * @note   深さ優先探索にΘ(V + E)時間かかり、|V|個の頂点のそれぞれを連結リストの先頭に挿入するのにΟ(1)時間かかるので、
 *         トポロジカルソートはΘ(V + E)時間で実行できる
 *
 * @param  const graph_t& G 有向非巡回グラフ
 * @return 既ソートリスト
 */
array_t tsort(const graph_t& G);



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of __TSORT_HPP__










