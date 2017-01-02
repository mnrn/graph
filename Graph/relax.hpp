/**
 * @brief  単一始点最短路問題(single-source shortest paths problem)に関する(特に、緩和についての)物置
 *
 *
 * @note   最短路問題(shortest paths problem)では、辺を実数値重みに写す関数w:E -> Rを持つ、重み付き有向グラフG=(V,E)が与えられる
 *         道p = <v0,v1,...,vk>の重み(weight)は、道を構成する辺の重みの和
 *           w(p) = Σ[i = 1, k]w(v(i-1),vi)
 *         である. uからvへの最短路重み(shortest path weight) δ(u,v)を
 *           δ(u,v) = { min{w(p): u~p~>v } : uからvへの道が存在するとき, ∞ : それ以外のとき } 
 *         と定義する. このとき、頂点uから頂点vへの最短路(shortest path)は重みがw(p) = δ(u,v)の道pである
 *         時間、コスト、罰金、損失、あるいはそのほかの、道に沿って線形に増大し最小化することが望ましい量など、距離以外の距離関数(metric)を
 *         辺重みと見なすこともできる
 *
 * @note   最短路重みだけではなく最短路上の頂点も計算したい場合が多い. グラフG=(V, E)が与えられたとき、各頂点v ∈ Vに対して、
 *         ある別の頂点かNILを値とする先行点(predecessor)v.πを管理する. 本章の最短路アルゴリズムは、頂点vから始まる先行点の連鎖が、
 *         sからvへの最短路に沿って逆戻りするように属性πを設定する. しかし、最短路アルゴリズムを実行中はπ値が最短路を示すとは限らない
 *         幅優先探索と同様、π値から誘導される先行点部分グラフ(predecessor subgraph)Gπ=(Vπ, Eπ)に興味がある
 *         ここでふたたび、頂点集合VπはNILを先行点としないGの頂点集合に始点sを加えたもの、すなわち
 *           Vπ = {v∈V : v.π != NIL} ∪ {s},
 *         有効辺の集合Eπは、Vπの頂点に対するπ値によって誘導される辺集合、すなわち
 *           Eπ = {(v.π, v) ∈ E : v ∈ Vπ - {s}}
 *         である
 *
 * @note   アルゴリズムが生成するπ値に対応するGπは「最短路木」である. 最短路木は簡単に言うと、始点sから到達可能なすべての頂点について
 *         sからの最短路を含む根付き木である. 最短路木は幅優先木に似ているが、辺数ではなく、辺重みによって定義された始点からの最短路を含んでいる
 *         正確に定義する. G=(V,E)を重み関数w:E -> Rを持つ重み付き有向グラフとし、Gは始点s ∈ Vから到達可能な閉路を含まず、最短路が明確に定義されていると
 *         仮定する. sを根とする最短路木(shortest paths tree)は有向部分グラフG'=(V', E')である. ここで、V'⊆VとE'⊆Eは以下の条件を満たす
 *
 *           1. V'がGにおいてsから到達可能な頂点の集合である
 *           2. G'はsを根とする根付き木である
 *           3. すべてのv∈V'に対して、G'におけるsからvへの単純道がGにおけるsからvへの最短路である
 *
 *         最短路は必ずしも一意ではないから、最短路木も一意に決まらないことがある
 *
 * @note   緩和(relaxation)をアルゴリズムは用いる. 各頂点v∈Vに対して、始点sからvへの最短路重みの上界を属性v.dとして管理する
 *         v.dを最短路推定値(shortest path estimate)と呼ぶ. 以下のΘ(V)時間手続きによって最短路推定値と先行点を初期化する
 *
 *           INITIALIZE-SINGLE-SOURCE(G, s)
 *           1 for each vertex v ∈ G.V
 *           2   v.d = ∞
 *           3   v.π = NIL
 *           4 s.d = 0
 *
 *         初期化の後、すべてのv∈Vについてv.π = NIL、すべてのv∈V-{s}についてv.d = ∞である
 *
 *         辺(u, v)の緩和(relaxing)は、uを経由することでvへの既知の最短路を改善できるか否か判定し、改善できるならばv.dとv.πを更新する
 *         緩和によって最短路推定値v.dが減少し、vの先行点属性v.πが更新されることがある. 辺(u, v)上の緩和をΟ(1)時間で実行する疑似コードを以下に示す
 *
 *           RELAX(u, v, w)
 *           1 if v.d > u.d + w(u, v)
 *           2   v.d = u.d + w(u, v)
 *           3   v.π = u
 *
 * @note   最短路と緩和の性質を明示する
 *
 *           三角不等式
 *             任意の辺(u, v) ∈ Eに対してδ(s, v) <= δ(s, u) + w(u, v)が成立する
 *
 *           上界性
 *             すべての頂点v ∈ Vに対して、v.d >= δ(s, v)が常に成立する. v.dが値δ(s, v)に達するとその後は変化しない
 *
 *           無経路性
 *             頂点sからvへの道がなければ、v.d = δ(s, v) = ∞が常に成立する
 *
 *           収束性
 *             ある辺(u, v) ∈ Eに対して、s ~> u -> vがGの最短路であり、辺(u, v)を緩和前のある時点でu.d = δ(s, u)だったならば、
 *             緩和後は常にv.d = δ(s, v)が成立する
 *
 *           経路緩和性
 *             p = <v0, v1, ..., vk>がs = v0からvkへの最短路で、pの辺が(v0, v1), (v1, v2),...,(v(k-1), vk)の順に緩和されたとき、
 *             vk.d = δ(s, vk)が成立する. この性質は他の任意の緩和とは無関係に成立する. たとえこれらの辺の緩和がpの辺の緩和と混在したとしても
 *             この性質は成立する
 *
 *           先行点部分グラフ性
 *             すべてのv ∈ Vに対してv.d = δ(s, v)が成立すれば、先行点部分グラフはsを根とする最短路木である
 *
 *
 * @date   2016/03/13
 */



//****************************************
// インクルードガード
//****************************************

#ifndef __RELAX_HPP__
#define __RELAX_HPP__



//****************************************
// 必要なヘッダファイルのインクルード
//****************************************

#include "graph.hpp"
#include <functional>



//****************************************
// 名前空間の始端
//****************************************

GRAPH_BEGIN



//****************************************
// 関数の宣言
//****************************************

/**
 * @brief  Θ(V)の手続きによって最短路推定値と先行点を初期化する
 * @note   初期化の後、すべてのv ∈ Vについてv.π = NIL、すべてのv ∈ V - {s}についてv.d = ∞である
 */
static inline void initialize_single_source(vertices_t& V, index_t s)
{
    for (auto&& v : V) {
        v.d     = limits::inf;
        v.pi    = limits::nil;
    }
    V[s].d = 0;
}


/**
 * @brief  Θ(V)の手続きによって最短路推定値と先行点および頂点色を初期化する
 * @note   初期化の後、すべてのv ∈ Vについてv.π = NIL、
 *         すべてのv ∈ V - {s}についてv.d = ∞、 v.color = WHITEである　
 */
static inline void initialize_single_source_with_color(vertices_t& V, index_t s)
{
    for (auto&& v : V) {
        v.d     = limits::inf;
        v.pi    = limits::nil;
        v.color = vcolor::white;
    }
    V[s].d = 0;
    V[s].color = vcolor::gray;
}


/**
 * @brief  Θ(V)の手続きによって最短路推定値と先行点および訪問済みフラグを初期化する
 * @note   初期化の後、すべてのv ∈ Vについてv.π = NIL、v.visited = false
 *         すべてのv ∈ V - {s}についてv.d = ∞である　
 */
static inline void initialize_single_source_with_visitor(vertices_t& V, index_t s)
{
    for (auto&& v : V) {
        v.d       = limits::inf;
        v.pi      = limits::nil;
        v.visited = false;
    }
    V[s].d = 0;
}



/**
 * @brief  辺(u, v)の緩和(relaxing)は、uを経由することでvへの既知の最短路を
 *         改善できるか否かを判定し、改善できるならばv.dとv.πを更新する
 *         緩和によって最短路推定値v.dが減少し、vの先行点属性v.πが更新されることがある
 *
 * @note   以下のコードは、辺(u, v)上の緩和をΟ(1)時間で実行する
 *         ただし、手続きpred(V, u)がΟ(1)で実行されることを仮定する
 *
 * @param  vertices_t& V   頂点集合V
 * @param  index_t u       辺(u, v)の始点u (ただし、u ∈ V)
 * @param  index_t v       辺(u, v)の終点v (ただし、v ∈ V)
 * @param  weight_t w      辺(u, v)の重みw
 * @param  Predicate pred  relax可能な前提条件を記述した述語
 */
static void relax(vertices_t& V,
           index_t u, index_t v, weight_t w,
           std::function<bool(const vertices_t&, index_t)> pred)
{
    if (pred(V, u) && V[v].d > V[u].d + w) {
        V[v].d = V[u].d + w;
        V[v].pi = u;
    }
}


// オーバーロードされたrelax関数群

static inline void relax(vertices_t& V, const edge& e,
           std::function<bool(const vertices_t&, index_t)> pred)
{
    relax(V, e.src, e.dst, e.w, pred);
}
static inline void relax(vertices_t& V, const matrix_t& W,
           index_t u, index_t v,
           std::function<bool(const vertices_t&, index_t)> pred)
{
    relax(V, u, v, W[u][v], pred);
}


/**
 * @brief  辺(u, v)を緩和すると同時に、頂点vおよび道s~>vの重みをmin優先度付きキューQに挿入する
 * 
 * @tparam PriorityQueue min優先度付きキューの型
 * @param vertices_t&    V 頂点集合V
 * @param const edge&    e 辺(u, v)
 * @param PriorityQueue& Q min優先度付きキュー
 */
template<class PriorityQueue>
void relax_with_heap(vertices_t& V, const edge& e, PriorityQueue& Q)
{
    index_t u = e.src, v = e.dst;
    if (V[v].color != vcolor::black && V[v].d > V[u].d + e.w) {
        V[v].d     = V[u].d + e.w;
        V[v].pi    = u;
        V[v].color = vcolor::gray;
        Q.emplace(v, V[v].d);
    }
}



//****************************************
// 名前空間の終端
//****************************************

GRAPH_END



#endif  // end of __RELAX_HPP__

