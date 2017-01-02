/**
 * @brief 幅優先探索の実装
 * @date  2016/02/11 ~ 2016/03/13
 */



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include <iostream>
#include <iterator>
#include <functional>
#include <utility>
#include <algorithm>
#include <vector>
#include <queue>

#include "bfs.hpp"



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の定義
//****************************************

/**
 * @brief  幅優先探索を行います
 *
 * @note   幅優先探索(breadth-first-search)は最も単純なグラフ探索アルゴリズムの1つである
 *         グラフG=(V, E)を始点(source vertex)sが与えられたとき、幅優先探索はGの辺を組織的に探索して、
 *         sから到達可能なすべての頂点を"発見し"、すべての到達可能な頂点についてsからの距離(変数の最小値)を計算する
 *         さらに、sを根とし、到達可能な頂点をすべて含む「幅優先木」を構成する. sから到達可能な任意の頂点vについて、
 *         幅優先木におけるsからvへの単純道はGにおけるsからvへの「最短路」、すなわち最小数の辺を含む道に対応する
 *         有向グラフと無向グラフのどちらに対してもこのアルゴリズムは正しく動く
 *
 *         探索頂点と未探索頂点の境界を、境界の幅一杯にわたって一様に拡張することが幅優先探索という名前の由来である
 *         すなわち、このアルゴリズムはsから距離k+1にある頂点を発見する前に距離kにあるすべての頂点を発見する
 *
 * @note   手続きBFSはグラフを探索しながら幅優先木を構築する.この木はπ属性に対応する.形式的に言うと、
 *         sを始点とするグラフG=(V, E)に対して、Gの先行点部分グラフ(predecessor subgraph)をGπ=(Vπ, Eπ)として定義する.ただし、
 *           Vπ = { v ∈ V : v.π != NIL } ∪ {s}
 *         かつ
 *           Eπ = { (v.π, v) : v ∈ Vπ - {s} }
 *         である.Vπがsから到達可能な全頂点から構成され、すべてのv∈Vπに対して、sからvに至る唯一の単純道がGπに存在し、
 *         しかもこれがGにおけるsからvに至る最短路になっているとき、先行点部分グラフGπを幅優先木(breadth-first tree)と呼ぶ
 *         幅優先木は連結で|Eπ| = |Vπ| - 1を満たすから、実際に木である. Eπの辺を木辺(tree edge)と呼ぶ
 *
 * @note   BFSの総実行時間はΟ(V+E)である.したがって、幅優先探索はGの隣接リスト表現のサイズの線形時間で走る
 *
 * @param  const graph& G  グラフG
 * @param  index_t s  始点s
 * @return 幅優先木
 */
vertices_t bfs(const graph_t& G, index_t s)
{
    index_t n = G.size();
    vertices_t V(n);

    for (auto&& u : V) {            // すべての頂点uについて、
        u.color = vcolor::white;    // uを白に彩色し、
        u.d     = limits::inf;      // u.dを無限大に設定し、
        u.pi    = limits::nil;      // uの親をNILに設定する
    }
    // 手続き開始と同時に始点sを発見すると考え、
    V[s].color = vcolor::gray;      // 始点sを灰色に彩色する 
    V[s].d     = 0;                 // s.dを0に初期化し、
    V[s].pi    = limits::nil;       // 始点の先行点をNILに設定する

    std::queue<index_t> Q;
    Q.push(s);                      // sだけを含むようにQを初期化する

    // 以下のwhile文に対して、つぎのループ不変式が成立する
    // while文の条件判定を行う時点ではキューQはすべての灰頂点を含む
    while (!Q.empty()) {
        index_t u = Q.front(); Q.pop();
        for (auto&& e : G[u]) {                  // uの隣接リストに
            index_t v = e.dst;                   // 属する各頂点vを考える
            if (V[v].color == vcolor::white) {   // vが白ならvは未発見である
                V[v].color = vcolor::gray;       // vを灰色に彩色し、
                V[v].d     = V[u].d + 1;         // 距離v.dをu.d+1に設定し、
                V[v].pi    = u;                  // uをその親v.piとして記録し、
                Q.push(v);                       // vをキューQの末尾に置く
            }
        }
        V[u].color = vcolor::black;  // uの隣接リストに属するすべての頂点の探索が完了すると、この頂点を黒に彩色する
    }
    // ある頂点を灰に彩色したときには、この頂点をQへ挿入し、ある頂点をQから削除したときには、この頂点を黒に彩色するので、
    // ループ不変式が保存される

    return V;
}


/**
 * @brief BFSが幅優先木を計算した後でこの手続きを用いれば、sからvへの最短路上の頂点を印刷できる
 */
void printpath(const vertices_t& V, index_t s, index_t v)
{
    if (v == s) {
        std::cout << s << " ";
    }
    else if (V[v].pi == limits::nil) {
        std::cout << "s\"から\"vへの道は存在しない\n"; 
    }
    else {
        printpath(V, s, V[v].pi);
        std::cout << v << " ";
    }
}



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



