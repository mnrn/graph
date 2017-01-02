/**
 * @brief 深さ優先探索に関する物置
 * @date  2016/02/12 ~ 2016/03/12
 */



//****************************************
// インクルードガード
//****************************************

#ifndef __DFS_HPP__
#define __DFS_HPP__



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
 * @brief  深さ優先探索を行います
 * 
 * @note   深さ優先探索は、その名が示すように、可能ならば常にそのグラフの"より深い部分を"探索するという戦略に従う
 *         未探索の外向辺が残る頂点の中で、最後に発見した頂点vから出る辺を深さ優先探索は探索する. vの辺をすべて探索し終えると、
 *         vを発見したときに通った辺を「バックトラック(逆戻り)」し、vの直前の頂点を出る(未探索の)辺の探索に戻る
 *         この処理は視点から到達可能なすべての頂点を発見するまで続く. 未発見の頂点が残されていれば、その1つを新たな始点として
 *         探索を繰り返す. アルゴリズムはすべての頂点を発見するまでこのプロセスを繰り返す
 *
 *         幅優先探索と同様、発見済みの頂点uの隣接リストを走査中に頂点vを発見すると、深さ優先探索はvの先行点属性v.πをuに設定し、
 *         この事象を記録する. 先行点部分グラフが木である幅優先探索と違い、深さ優先探索では複数の始点から探索を繰り返すことがあるから、
 *         先行点部分グラフが複数の木から構成されることがある. そこで、深さ優先探索の先行点部分グラフ(predecessor subgraph)を
 *         幅優先探索と少し違って、Gπ = (V, Eπ)と定義する. ここで
 *           Eπ = {(v.π, v) : v ∈ V かつ v.π != NIL }
 *         である. 深さ優先探索の先行点部分グラフは複数の深さ優先木(depth-first tree)から構成される深さ優先森(depth-first forest)
 *         を形成する. Eπに属する辺を木辺(tree edge)と呼ぶ
 *
 *         幅優先探索と同様、深さ優先探索は頂点を探索状態で示す色で彩色する. 初期状態では各頂点は白であり、探索の中で発見されれば(discovered)
 *         灰に変わり、探索が終了すれば(finished)、すなわちこの頂点の隣接リストが完全に黒に変わる. 深さ優先探索では、
 *         各頂点が1つの深さ優先木の中にだけ現れることが保証されるので、これらの木は互いに共通部分を持たない
 *
 *         深さ優先探索は深さ優先森を構成するとともに各頂点に時刻印(timestamp)を押す. 各頂点vは2種類の時刻印を持つ
 *         第1の時刻印v.dはvを最初に発見し、灰に彩色した時刻を記録する. 第2の時刻印v.fはvの隣接リストを調べ終えて黒に彩色した時刻を記録する
 *         これらの時刻印はグラフ構造に関する重要な情報を供給し、深さ優先探索の振る舞いを理解するのに役立つ
 *
 *         以下に示す手続きDFSでは、頂点uを発見した時刻を属性u.dに記録し、uの処理を終了した時刻を属性u.fに記録する
 *         |V|個の各頂点について発見事象と終了事象はそれぞれ1個しか生起しないから、時刻印は1から2|V|の範囲の整数である
 *         任意の頂点uについて
 *           u.d < u.f
 *         が成立する. 頂点uは時刻u.d以前はWHITE, 時刻u.dと時刻u.fの間はGRAY, 時刻u.f以降はBLACKである
 *
 * @param  グラフG(無向でも有向でもよい)
 * @param  深さ優先森
 */
std::pair<vertices_t, array_t> dfs(const graph_t& G);



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of __DFS_HPP__

